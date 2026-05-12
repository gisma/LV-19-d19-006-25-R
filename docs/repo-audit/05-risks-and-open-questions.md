# Risks and Open Questions

## Documentation-Folder Audit Update

- QUESTION (`src/README_1st.md`, "Repository Mapping"): Why are `src/S5_decisions` and `src/lib` documented when the actual folders are `src/S5-decision` and `src/r-libs`?
- QUESTION (`src/S5-decision/doc/S5-meta.qmd` vs. `src/S5-decision/05-9_decision_selector.qmd`): Is S5 only final selection after candidate generation, or does it include candidate generation through per-domain scripts and selector settings?
- QUESTION (`src/S5-decision/doc/S5-meta.qmd`): Why is S5 titled "Validation and Selection Space" when `src/README_1st.md` and the selector tutorial identify S6 as validation space?
- QUESTION (`src/S4_signatures/doc/S4-meta.qmd`): Where are the documented symbolic tokens/controlled compression implemented and registered?
- QUESTION (`src/S4_signatures/doc/04-3_signatures_RF-segs_classification.qmd`): Should RF classification be documented as S4L even though its implementation is under `src/S4_signatures`?
- QUESTION (`src/run/meta.md`): Are the documented `run_S1_all.R`, `run_S2_all.R`, and `run_S3_all.R` planned, removed, or intentionally absent?
- QUESTION (`src/S5-decision/05-9_decision_selector.qmd`): Where is the documented S6 validation implemented?
- QUESTION (`src/r-libs/meta.md`): How should the "no side effects" library rule be interpreted for helper functions that download data or write files?
- QUESTION (`src/_core/meta.md`): How should the "no hard-coded absolute paths" core rule be reconciled with fixed local tool locations in executable workflows?

