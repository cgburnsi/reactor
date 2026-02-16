# ExecPlan: Input Schema, Parser, and Validator (v1)

## 1) Purpose
- Problem: early mesh development is blocked without a stable, validated way to read user domain and boundary inputs.
- Outcome: a deterministic parser + validator that produces a normalized internal model or actionable errors.
- Why now: this unlocks all downstream meshing work and reduces rework risk in later milestones.

## 2) Scope
- In scope:
  - define v1 input schema for domain geometry, boundary tags, and meshing controls
  - implement parser into normalized Python data structures
  - implement deterministic validation and error reporting
  - add unit tests for happy path and key failure modes
  - add a runnable example that loads input and prints normalized model
- Out of scope:
  - full triangular meshing implementation
  - mesh-quality optimization and smoothing
  - production export/reporting formats beyond parser diagnostics

### Future-compatibility note
- The v1 schema and internal model should be cell-type extensible (`tri` now, `quad` later) without implementing quad meshing in this plan.
- Keep Milestone 1 implementation triangle-first; do not add Milestone 4 logic here.

## 3) Alignment Checkpoint (Must Confirm Before Coding)
- Success criteria:
  - given valid input, parser returns normalized data structure with deterministic ordering/content
  - given invalid input, validator reports clear error messages tied to offending fields
  - tests cover representative valid/invalid cases and pass with stdlib `unittest`
- Non-goals:
  - do not add external dependencies
  - do not implement meshing algorithms in this plan
- Example scenario:
  - Input: simple polygon domain with one boundary tag set and target element size control
  - Expected output/behavior: normalized model printed to stdout; repeated runs produce identical output
- Open questions:
  - [ ] final on-disk input format for v1 (`json` vs `csv` + conventions)
  - [ ] required vs optional meshing-control fields in v1 schema
- Assumptions to confirm:
  - [ ] v1 input can use stdlib-supported formats only
  - [ ] deterministic ordering can be enforced without changing user-authored field order semantics

## 4) Design Checkpoint
- Option 1: strict dataclass-style internal objects with explicit normalization pass
  - Pros: clear ownership and validation phases; easier to test deterministically
  - Cons: more initial boilerplate
- Option 2: dictionary-first parsing with inline validation
  - Pros: faster initial coding
  - Cons: weaker boundaries and harder long-term evolution
- Selected approach: Option 1, to optimize for correctness and maintainability before algorithmic expansion.

## 5) Implementation Plan (Vertical Slices)
1. Slice 1: schema + parser scaffold
   - Files: `src/` parser module(s), `tests/` parser tests, `examples/` parser demo
   - Done when: valid minimal input loads into normalized internal structure with deterministic output
2. Slice 2: validation and error taxonomy
   - Files: `src/` validation module(s), `tests/` validation tests
   - Done when: invalid inputs fail deterministically with actionable field-level errors
3. Slice 3: integration and docs sync
   - Files: `README.md` or `docs/` usage notes, final test cases
   - Done when: documented example runs cleanly and all tests pass

## 6) Validation Plan
- Commands:
  - `python -m unittest`
  - `python <example_parser_script>.py`
- Expected results:
  - unit tests pass
  - example script prints normalized model for valid input
  - invalid test fixtures produce deterministic, readable errors
- Numerical tolerance (if applicable):
  - not applicable for parser/validator logic in this plan

## 7) Progress Log (Update During Execution)
- [ ] Intent checkpoint complete (success criteria, non-goals, examples agreed)
- [ ] Design checkpoint complete (approach selected)
- [ ] Slice 1 implemented + validated
- [ ] Slice 2 implemented + validated
- [ ] Slice 3 implemented + validated
- [ ] Documentation updated

### Status Notes
- 2026-02-16 00:00 - Initial plan created from milestone decomposition.

## 8) Surprises and Discoveries
- None yet.

## 9) Decision Log
- 2026-02-16 - Decision: start with explicit internal model + normalization pass.
  - Context: parser choices affect all downstream meshing code.
  - Chosen: structured model-first approach.
  - Tradeoff: slightly slower startup in exchange for clearer correctness guarantees.

## 10) Outcomes and Retrospective
- Delivered:
  - pending
- Validation evidence:
  - pending
- Remaining risks:
  - pending
- Follow-ups:
  - pending

## 11) Handoff Notes (for a fresh agent)
- Current state: planning complete; implementation has not started.
- Next action: complete Alignment Checkpoint open questions, then implement Slice 1.
- Known constraints: stdlib + NumPy + Matplotlib only; tests via stdlib `unittest`.
