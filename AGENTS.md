# Repository Ground Rules

## Branches and approval

- Never edit files or create commits on `main`. Switch to a working branch before making changes.
- On branches other than `main`, proceed with task-related edits, implementation, documentation, and tests without requesting additional approval. This standing authorization does not override tool permissions or authorize destructive actions, merging into `main`, publishing, or external communications.
- Preserve unrelated user changes.

## Implementation specifications

- Whenever an algorithm changes, automatically update the relevant implementation specification in the same change. Keep the documented inputs, outputs, steps, assumptions, and limitations consistent with the code.
- If no implementation specification exists, create one under `docs/specs/` for the affected functionality.

## Tests

- Every code implementation or behavior change must include new or updated pytest tests covering the affected behavior. Include meaningful edge cases and regression coverage for bug fixes.
- Run the relevant pytest tests whenever new code is developed or existing code changes. Resolve failures caused by the change before considering the work complete.
- Report which tests ran and their results. If tests cannot run, explain the blocker and do not claim they passed.
- Documentation-only changes do not require pytest tests.
