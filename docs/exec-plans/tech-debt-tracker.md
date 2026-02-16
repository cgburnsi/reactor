# Tech Debt Tracker

Track deferred decisions and known shortcomings that are intentionally postponed.

## Status key
- `open`: identified, not yet scheduled
- `planned`: linked to an ExecPlan
- `closed`: addressed and verified

## Items
| ID | Status | Area | Summary | Impact | Planned ExecPlan | Notes |
|---|---|---|---|---|---|---|
| TD-001 | open | Input schema | Finalize v1 input format choice (`json` vs alternatives). | Blocks parser ergonomics finalization. | 2026-02-16-input-schema-parser-validator-v1.md | Decide during Alignment Checkpoint. |
| TD-002 | open | Meshing architecture | Finalize method interface shape for orchestrator plugins. | Affects all strategy implementations. | 2026-02-16-ep-a-meshing-interfaces-orchestrator-skeleton.md | Keep interface minimal and stable. |
| TD-003 | open | Hybrid meshing | Define tri/quad transition policy for v1 hybrid output. | Blocks robust boundary-layer + interior merge. | 2026-02-16-ep-d-tri-quad-transition-and-stitching.md | Requires EP-B and EP-C outputs. |
