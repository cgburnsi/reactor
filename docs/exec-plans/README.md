# Exec Plans

Use this folder for implementation-ready plans that drive substantial changes.

## Folder layout
- `docs/exec-plans/TEMPLATE.md`: reusable template
- `docs/exec-plans/active/`: plans currently in progress
- `docs/exec-plans/completed/`: finished plans

## Naming convention
- `YYYY-MM-DD-<short-kebab-title>.md`
- Example: `2026-02-16-input-parser-v1.md`

## How to use
1. Copy `TEMPLATE.md` into `active/` with a dated name.
2. Complete Sections 1-4 before writing code.
3. Implement by vertical slices (Section 5).
4. Keep Sections 7-9 updated during execution.
5. When done, complete Sections 10-11 and move the file to `completed/`.

## Collaboration loop (agent + user)
At minimum, run these checkpoints in order:
1. Intent checkpoint: confirm success criteria, non-goals, and examples.
2. Design checkpoint: compare options and choose one approach.
3. Slice review: after each slice, verify behavior before continuing.

## Policy reminders
- Respect dependency limits in `AGENTS.md`.
- Keep changes small and deterministic.
- Record any tradeoff in the Decision Log.
