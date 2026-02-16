# ARCHITECTURE.md

## Purpose
Define the stable technical blueprint for this repository:
- core data contracts
- module boundaries
- meshing pipeline flow
- validation and determinism requirements
- extension points for multiple meshing methods

This document should remain mostly stable over time. Implementation sequencing and temporary decisions belong in `docs/exec-plans/`.

## Scope and Non-Scope
In scope:
- architecture for 2D unstructured mesh generation suitable for downstream FVM
- triangle-first implementation with future quad/hybrid support
- deterministic behavior and validation requirements

Out of scope:
- detailed algorithm derivations
- task-by-task execution details
- benchmark result logs (store in ExecPlans/docs)

## Design Principles
- Correctness before optimization.
- Deterministic outputs for fixed inputs and parameters.
- Small, explicit module interfaces.
- Strategy-based extensibility (no single hard-coded meshing algorithm).
- Separation of concerns: parsing, generation, validation, quality, export.

## Planned Repository Layout
Target layout (create as milestones progress):
- `src/`
- `src/models/`
- `src/io/`
- `src/meshing/`
- `src/validation/`
- `src/quality/`
- `src/export/`
- `tests/`
- `examples/`
- `docs/`

## Core Contracts

### DomainModel
Represents parsed and validated problem definition.
Expected fields (v1 target):
- geometry (outer boundary + optional holes)
- boundary tags and boundary condition metadata
- meshing controls (global size + optional local controls)
- method selection/config hints

### MeshModel
Represents mesh output independent of algorithm.
Expected fields (v1 target):
- node coordinates
- cell connectivity
- cell type metadata (triangle now; quad-ready later)
- face/edge connectivity
- boundary tag assignments

### MeshingMethod Interface
Common strategy interface for generation methods.
Conceptual contract:
- input: `DomainModel`, method constraints/config
- output: `MeshModel`
- behavior: deterministic for fixed input/config

### MeshValidator Interface
Performs structural and geometric checks.
Examples:
- connectivity consistency
- orientation and inversion checks
- boundary closure and tag consistency

### QualityEvaluator Interface
Computes quality metrics without mutating topology.
Examples:
- triangle skewness/aspect metrics
- quad orthogonality/skew metrics (future)

## Meshing Pipeline (Conceptual)
1. Input parse/load (`io`)
2. Domain validation + normalization (`validation` + `models`)
3. Strategy selection (`meshing` orchestrator)
4. Mesh generation by selected method(s)
5. Mesh validation (`validation`)
6. Quality evaluation (`quality`)
7. Optional quality improvement pass(es) (`meshing`/`quality`)
8. Export/reporting (`export`, `examples`)

## Orchestration Model
A `MeshingOrchestrator` coordinates strategy execution.
Responsibilities:
- select strategy from config
- invoke one or more methods in sequence
- support region-based assignment (future hybrid mode)
- run validator gates after each major step
- return a single validated `MeshModel`

The orchestrator owns workflow policy. Individual strategies own only local generation logic.

## Strategy Roadmap
Near-term:
- `AdvancingFrontTriMethod` (triangle baseline)

Planned:
- `QuadLayerMethod` for near-wall boundary layers
- hybrid workflows combining quad near-wall + tri interior
- explicit tri/quad transition handling and stitching

## Determinism and Reproducibility
Requirements:
- identical input/config must produce identical output ordering and values (within tolerance where applicable)
- iteration order must be defined where unordered containers could introduce drift
- if randomness is introduced in the future, seed must be explicit and logged

## Validation Invariants
Minimum invariants for a valid mesh output:
- all cell indices reference valid nodes
- no duplicate/self-intersecting boundary loops in normalized domain
- no inverted elements in generated mesh
- boundary tags preserved from domain model to mesh model
- topology consistency across adjacency relations

Numerical comparisons follow policy in `AGENTS.md`:
- `abs(err) <= max(abs_tol, rel_tol * max(1.0, abs(reference)))`

## Error Handling
- Parsing and validation errors must be actionable and field-specific.
- Fail fast on invalid input before entering expensive generation stages.
- Validation failures in generation should surface with enough context for debugging.

## Dependencies and Tooling Constraints
Hard constraints (from `AGENTS.md`):
- Python standard library
- NumPy
- Matplotlib

Testing framework:
- stdlib `unittest`

## Relationship to Planning Docs
- `AGENTS.md` defines milestones, policies, and execution priority.
- `ARCHITECTURE.md` defines stable technical design and boundaries.
- `docs/exec-plans/*` define execution steps, decisions, and progress for specific work items.

## Related Design Docs
- `docs/roadmap.md` - milestone outcomes, sequencing, and meshing strategy direction.
- `docs/design-docs/core-beliefs.md` - design priority order and tradeoff policy.
- `docs/design-docs/mesh-data-model.md` - conceptual domain/mesh structures and invariants.
- `docs/design-docs/validation-invariants.md` - required correctness checks before/after generation.

## Change Policy for This File
Update `ARCHITECTURE.md` when:
- core interfaces/contracts change
- module boundaries are redefined
- pipeline stages or ownership responsibilities change

Do not update for routine implementation progress unless architecture actually changed.
