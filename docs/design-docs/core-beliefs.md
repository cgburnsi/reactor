# Core Beliefs

This document captures stable engineering priorities used to make tradeoff decisions.

## Priority Order
1. Correctness and validity of mesh/domain outputs
2. Determinism and reproducibility
3. Clear interfaces and modular boundaries
4. Extensibility for additional meshing strategies
5. Performance optimization

## Working Principles
- Prefer simple, testable designs over clever implementations.
- Keep algorithm-specific logic behind shared method interfaces.
- Fail fast on invalid input with actionable errors.
- Avoid broad refactors unrelated to current milestone intent.
- Defer optimization until correctness and validation harness are stable.

## Dependency Principle
This repository remains constrained to:
- Python standard library
- NumPy
- Matplotlib

If additional dependencies appear necessary, follow the Dependency Gate in `AGENTS.md`.

## Determinism Principle
For fixed input and configuration:
- outputs should be reproducible
- ordering should be explicit where iteration order matters
- seeded randomness is required if randomness is introduced

## How to use this doc
When a design tradeoff is unclear, choose the option that best aligns with the priority order above and record the decision in the active ExecPlan.
