# AGENTS.md - How to work in this repo (Python)

## Purpose
Operational rules for agents working in this repository. Keep this file concise. Put detailed plans in `docs/exec-plans/` and roadmap detail in `docs/roadmap.md`.

## Project Intent
Build a small, correct, modular Python library that implements mesh capabilities for a downstream Finite Volume Method (FVM) solver.

## Source of Truth (read first)
- `README.md` - installation, usage, and project overview
- `ARCHITECTURE.md` - stable technical blueprint (contracts, boundaries, pipeline)
- `docs/roadmap.md` - milestones, sequencing, and strategy direction
- `PLANS.md` - planning/execution workflow for substantial changes
- `docs/exec-plans/README.md` - ExecPlan usage in this repo
- `docs/exec-plans/TEMPLATE.md` - required ExecPlan template

If any source file is missing, create a minimal stub rather than inventing undocumented architecture.

## Execution Rules
- Work milestone-by-milestone per `docs/roadmap.md` unless explicitly approved otherwise.
- Use ExecPlans for algorithmic or cross-cutting work.
- Keep changes small, deterministic, and focused.
- No unrelated refactors in the same change.

## Definition of Done
A change is done only if:
1. It runs without errors for intended use.
2. Tests pass (unit + relevant regression checks).
3. New behavior has tests (or an explicit note why not).
4. Docs are updated when behavior/interfaces/formats change.
5. Behavior is deterministic (within stated tolerances).
6. Dependency rules are respected.

## Numerical Tolerance Policy
Use:
- `abs(err) <= max(abs_tol, rel_tol * max(1.0, abs(reference)))`

## Approval Policy
Loose mode: small and medium scoped edits are allowed autonomously.
Ask before:
- deleting files
- major architecture rewrites
- dependency or tooling policy changes

## Dependency Rules (hard)
Allowed dependencies:
- Python standard library
- NumPy
- Matplotlib

No other packages may be added (runtime or development).

### Dependency Gate
If a task appears to require a new dependency:
1. Stop immediately.
2. Record:
   - blocked feature
   - why stdlib + NumPy + Matplotlib are insufficient
   - two alternatives that keep dependency set unchanged
   - if still justified: proposed dependency, scope, and rationale
3. Do not add the dependency without explicit approval.

## Testing and I/O Guidance
- Use stdlib `unittest` for tests.
- Prefer stdlib `json`, `csv`, `pathlib`, and `logging` for I/O/instrumentation.

## Commit Message Policy
Use Conventional Commits:
- `type(scope): subject`
Preferred types: `docs`, `feat`, `fix`, `test`, `refactor`, `chore`.
