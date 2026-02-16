# ExecPlan: <Short Title>

## 1) Purpose
- Problem: <What problem this plan solves>
- Outcome: <What exists when this plan is done>
- Why now: <Why this is the next best increment>

## 2) Scope
- In scope:
  - <Item 1>
  - <Item 2>
- Out of scope:
  - <Item A>
  - <Item B>

## 3) Alignment Checkpoint (Must Confirm Before Coding)
- Success criteria:
  - <Observable behavior 1>
  - <Observable behavior 2>
- Non-goals:
  - <What we explicitly will not do>
- Example scenario:
  - Input: <example>
  - Expected output/behavior: <example>
- Open questions:
  - [ ] <Question 1>
  - [ ] <Question 2>
- Assumptions to confirm:
  - [ ] <Assumption 1>
  - [ ] <Assumption 2>

## 4) Design Checkpoint
- Option 1: <Approach summary>
  - Pros: <...>
  - Cons: <...>
- Option 2: <Approach summary>
  - Pros: <...>
  - Cons: <...>
- Selected approach: <Option + reason>

## 5) Implementation Plan (Vertical Slices)
1. Slice 1: <small end-to-end piece>
   - Files: `<path1>`, `<path2>`
   - Done when: <behavior + test proof>
2. Slice 2: <next piece>
   - Files: `<path3>`
   - Done when: <behavior + test proof>
3. Slice 3: <next piece>
   - Files: `<path4>`
   - Done when: <behavior + test proof>

## 6) Validation Plan
- Commands:
  - `python -m unittest <test_target>`
  - `<additional command>`
- Expected results:
  - <What output/metrics indicate success>
- Numerical tolerance (if applicable):
  - `abs(err) <= max(abs_tol, rel_tol * max(1.0, abs(reference)))`

## 7) Progress Log (Update During Execution)
- [ ] Intent checkpoint complete (success criteria, non-goals, examples agreed)
- [ ] Design checkpoint complete (approach selected)
- [ ] Slice 1 implemented + validated
- [ ] Slice 2 implemented + validated
- [ ] Slice 3 implemented + validated
- [ ] Documentation updated

### Status Notes
- YYYY-MM-DD HH:MM - <note>

## 8) Surprises and Discoveries
- <What was learned that changed implementation details>

## 9) Decision Log
- YYYY-MM-DD - Decision: <decision>
  - Context: <why decision was needed>
  - Chosen: <what was chosen>
  - Tradeoff: <what was accepted>

## 10) Outcomes and Retrospective
- Delivered:
  - <what shipped>
- Validation evidence:
  - <command/results summary>
- Remaining risks:
  - <known limitations>
- Follow-ups:
  - <next plan candidates>

## 11) Handoff Notes (for a fresh agent)
- Current state: <what is complete right now>
- Next action: <single best next step>
- Known constraints: <dependency/policy constraints>
