"""
core/constants.py
-----------------
Global constants and tolerances for SnapMesh
"""

# Geometric tolerance for point coincidence
# Two points are "the same" if distance < GEOM_TOL
GEOM_TOL = 1e-10

# Tolerance for containment/intersection tests
# Used in ray casting, segment intersection, etc.
CONTAINMENT_TOL = 1e-8

# Tolerance for mesh quality checks
# Volumes/areas smaller than this are problematic
VOLUME_TOL = 1e-12

# Squared versions (for performance - avoid sqrt)
GEOM_TOL_SQ = GEOM_TOL ** 2
CONTAINMENT_TOL_SQ = CONTAINMENT_TOL ** 2