# Mesh Data Model (v1 Draft)

Defines the intended structure of domain and mesh representations. This is a design target and may be refined as milestones progress.

## Goals
- Keep model independent from any single meshing algorithm.
- Support triangle-first implementation while remaining quad-ready.
- Enable deterministic validation and downstream FVM use.

## DomainModel (conceptual)
Required fields:
- `geometry`
  - `outer_boundary`: ordered points
  - `holes`: zero or more ordered point loops
- `boundary_tags`
  - mapping from boundary segments/loops to semantic tags
- `meshing_controls`
  - global size target
  - optional local/region controls
- `method_config`
  - selected method identifier and method-specific options

Invariants:
- boundaries are non-self-intersecting
- hole loops are contained within outer boundary
- tag references point to existing boundary entities

## MeshModel (conceptual)
Required fields:
- `nodes`: array/list of coordinates
- `cells`: connectivity entries referencing node indices
- `cell_types`: per-cell type metadata (`tri` now, `quad` planned)
- `faces_or_edges`: adjacency/boundary connectivity representation
- `boundary_tags`: tags attached to boundary edges/faces

Invariants:
- all indices are in bounds
- no duplicate invalid connectivity entries
- boundary entities carry valid tags

## Indexing and Ordering
- Use explicit indexing conventions and document zero-based assumptions.
- Preserve deterministic ordering of nodes/cells where feasible.
- Any reordering pass must be deterministic for fixed input/config.

## Extensibility Notes
- Do not hard-code triangle-only assumptions into shared models.
- New cell types should be additive through `cell_types` and method-specific generation logic.
- Transition-region metadata may be added later for hybrid tri/quad workflows.

## Open Questions
- Exact representation for `faces_or_edges` in v1.
- Minimal metadata set required by downstream FVM assembly.
