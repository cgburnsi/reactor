# Roadmap

## Purpose
This document tracks capability outcomes, sequencing, and planned execution themes.

## Milestones
1. Input/domain model and validation pipeline
2. Baseline triangular meshing on canonical geometries
3. Quality improvement, boundary fidelity, and export/reporting
4. Quad-capable meshing for boundary-layer-focused grids

### Milestone 1: Input/domain model and validation pipeline
- Goal: deterministically parse user input into a normalized internal domain model with clear validation errors.
- Required:
  - input schema and parser for domain, boundaries, and meshing controls
  - deterministic validation with actionable error messages
  - unit tests covering valid/invalid examples and edge cases
  - one runnable example that prints normalized parsed data
- Proof checks:
  - invalid input fails with deterministic errors
  - valid input produces identical normalized model on repeated runs

### Milestone 2: Baseline triangular meshing on canonical geometries
- Goal: generate valid triangular meshes for simple polygons (including at least one hole case).
- Required:
  - baseline mesh generation algorithm implementation
  - topology/geometry consistency checks (orientation, inversion, connectivity)
  - unit tests and at least one benchmark-style example case
  - deterministic behavior for fixed inputs and parameters
- Proof checks:
  - generated meshes contain no inverted elements on canonical examples
  - boundary loops and adjacency checks pass for generated meshes

### Milestone 3: Quality improvement, boundary fidelity, and export/reporting
- Goal: improve mesh quality while preserving boundary semantics and enabling downstream use.
- Required:
  - quality metrics and at least one quality-improvement pass
  - boundary tag preservation through all transformations
  - export format for nodes/elements/boundary tags
  - plotting/reporting example including quality histogram
- Proof checks:
  - quality metrics improve (or do not regress) after improvement pass
  - exported mesh reloads without repair and retains boundary tags

### Milestone 4: Quad-capable meshing for boundary-layer-focused grids
- Goal: support quadrilateral-dominant meshes where boundary-layer resolution is important.
- Required:
  - quad-mesh generation strategy and controls for boundary-layer spacing
  - tri/quad transition handling where full-quad regions are not possible
  - quality checks specific to quads (skewness, aspect ratio, orthogonality)
  - deterministic tests on canonical boundary-layer geometries
- Proof checks:
  - boundary-layer cases generate valid quad-capable meshes without topological defects
  - boundary tags and near-wall spacing constraints are preserved in output

## Execution Priority
- Work milestone-by-milestone in order; do not start later milestones early unless explicitly approved.
- During Milestone 1, prioritize parser correctness, schema clarity, and deterministic validation.
- During Milestone 2, prioritize geometric/topological correctness over runtime optimization.
- During Milestone 3, prioritize output fidelity and verification/reporting.
- During Milestone 4, prioritize boundary-layer fidelity and robust tri/quad transitions.
- Within each milestone, resolve the largest technical unknown first and implement in vertical slices.
- Keep implementation triangle-first until Milestones 1-3 are complete or explicit approval is given to start Milestone 4 early.

## Meshing Strategy Direction
- Use a pipeline with pluggable meshing strategies, not a single hard-coded algorithm.
- Stable contracts:
  - `DomainModel`
  - `MeshModel`
  - `MeshingMethod`: `generate(domain, constraints) -> MeshModel`
- Use an orchestrator to assign methods by region, run sequence steps, validate between steps, and stitch outputs.
- Hybrid target: quad-focused boundary layers near walls and triangular meshing in interior regions, with explicit tri/quad transitions.

## Planned ExecPlan Sequence (strategy extensibility)
- EP-A: Meshing method interfaces + orchestrator skeleton
- EP-B: Advancing-front triangular meshing strategy
- EP-C: Quad boundary-layer meshing strategy
- EP-D: Tri/quad transition and mesh stitching validation
