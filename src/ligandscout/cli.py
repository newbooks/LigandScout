"""Typer entry point for LigandScout."""

import logging
import time
from datetime import datetime, timezone
from pathlib import Path

import httpx
import typer

from ligandscout.accessions import normalize_accession
from ligandscout.report import render_summary
from ligandscout.search import run_simple
from ligandscout.storage import save_results

LOG = logging.getLogger(__name__)

app = typer.Typer(
    help="Find protein ligands from a UniProt accession.",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=False,
)


@app.command()
def search(
    ctx: typer.Context,
    accession_id: str | None = typer.Argument(None, help="Swiss-Prot accession ID."),
    simple: bool = typer.Option(
        False, "--simple", help="Search this accession (default mode)."
    ),
    explore: bool = typer.Option(
        False, "--explore", help="Exploration mode (not implemented)."
    ),
    log_level: str = typer.Option("INFO", "--log-level", "-l", help="Logging level."),
) -> None:
    """Search for binding partners of ACCESSION_ID."""
    if accession_id is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()
    if simple and explore:
        raise typer.BadParameter("--simple and --explore cannot be used together.")
    level = log_level.upper()
    if level not in {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}:
        raise typer.BadParameter(
            "Use DEBUG, INFO, WARNING, ERROR, or CRITICAL.", param_hint="--log-level"
        )
    logging.basicConfig(level=level, format="[%(levelname)s]: %(message)s", force=True)
    if explore:
        typer.echo("Exploration mode is not implemented yet.", err=True)
        raise typer.Exit(2)
    try:
        accession = normalize_accession(accession_id)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="ACCESSION_ID") from exc
    # Keep HTTPX request diagnostics out of the default application log.
    logging.getLogger("httpx").setLevel(
        logging.WARNING if level != "DEBUG" else logging.DEBUG
    )
    try:
        with httpx.Client(timeout=30.0, follow_redirects=False) as client:
            outcome = run_simple(
                accession,
                client=client,
                sleep=time.sleep,
                now=lambda: datetime.now(timezone.utc),
            )
        folder = save_results(outcome, root=Path.cwd() / "results")
    except KeyboardInterrupt:
        typer.echo(
            "Interrupted; any output folder may contain incomplete results.", err=True
        )
        raise typer.Exit(130) from None
    except OSError:
        LOG.debug("Output publication failed", exc_info=True)
        typer.echo(
            "Could not save results; the output folder may contain incomplete or mixed results.",
            err=True,
        )
        raise typer.Exit(1) from None
    typer.echo(render_summary(outcome, folder=folder), nl=False)
    errors = [source.result.error for source in outcome.sources if source.result.error]
    for error in errors:
        LOG.error("%s", error.message)
    if errors:
        raise typer.Exit(1)


def main() -> None:
    """Run the installed command."""
    app()
