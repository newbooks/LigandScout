# Repository Ground Rules

## Branches and approval

- Never edit files or create commits on `main`. Switch to a working branch before making changes.
- On branches other than `main`, proceed with task-related code changes, documentation updates, and tests without requesting additional approval. This standing authorization remains subject to tool permissions. It does not authorize destructive actions, merges into `main`, publishing, or external communications.
- Preserve unrelated user changes.

## Implementation specifications

- Whenever an algorithm changes, update the relevant implementation specification as part of the same change without waiting for a separate request. Keep the documented inputs, outputs, steps, assumptions, and limitations consistent with the code.
- If no implementation specification exists, create one under `docs/specs/` for the affected functionality.

## Tests

- Every code change must include new or updated pytest tests for the affected functionality. Cover meaningful edge cases and add regression tests for bug fixes.
- Run the relevant pytest tests whenever you add or modify code. Resolve failures caused by the change before considering the work complete.
- Report which tests ran and their results. If tests cannot run, explain the blocker and do not claim they passed.
- Documentation-only changes do not require pytest tests.
