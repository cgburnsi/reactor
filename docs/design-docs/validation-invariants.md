# Validation Invariants

Defines baseline checks that must hold for normalized domain input and generated mesh output.

## Domain Validation (pre-meshing)
- Input schema is complete for required fields.
- Geometry loops are well-formed and non-self-intersecting.
- Hole loops are valid and placed inside outer boundary.
- Boundary tags reference valid geometry entities.
- Meshing controls are within allowed numeric ranges.

Failure behavior:
- fail fast before meshing starts
- provide actionable, field-specific error messages

## Mesh Validation (post-generation)
- all connectivity indices reference existing nodes
- no inverted elements
- adjacency/topology relations are internally consistent
- boundary entities are preserved and tagged correctly
- no orphan cells or malformed boundary connectivity

## Determinism Checks
For fixed input/config:
- repeated runs produce identical structural outputs
- ordering-dependent structures maintain stable order

## Numerical Tolerance Policy
For floating-point comparisons, use:
- `abs(err) <= max(abs_tol, rel_tol * max(1.0, abs(reference)))`

## Minimum Validation Gates by Milestone
- Milestone 1: schema + domain normalization checks
- Milestone 2: triangle mesh validity and topology consistency
- Milestone 3: quality metric checks and export integrity
- Milestone 4: quad/hybrid transition consistency and boundary-layer fidelity

## Reporting
Validation failures should report:
- failing invariant
- affected entity identifiers/indices
- minimal contextual values needed for debugging
