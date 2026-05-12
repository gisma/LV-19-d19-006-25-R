# Repository Audit Index

## Scope

- FACT (`README.md`): This repository is the technical backbone for the LV-19-d19-006-25 course and provides reproducible building blocks for environmental data workflows.
- FACT (`src/README_1st.md`): The repository is framed as a structured spatial research pipeline for optimized precipitation network design in Burgwald, with conceptual layers S0-S6.
- FACT (`AGENTS.md`): Audit work must remain documentation-only and must write only below `docs/repo-audit/`.
- FACT (`docs/repo-audit/07-doc-folder-audit.md`): A targeted documentation-folder audit was added after the first audit because documentation under `src/**/doc`, archive docs, module `meta.md` files, and the S5 selector tutorial contain primary architectural and workflow statements.

## Audit Documents

- FACT (`docs/repo-audit/01-repo-map.md`): Repository structure, key files, configuration, and entry points.
- FACT (`docs/repo-audit/02-architecture.md`): Layer architecture, component dependencies, control flow, and data flow.
- FACT (`docs/repo-audit/03-workflows.md`): Setup, execution, build/test status, data processing, and export workflows.
- FACT (`docs/repo-audit/04-project-ideas.md`): Project ideas, domain intent, modelling assumptions, and teaching goals.
- FACT (`docs/repo-audit/05-risks-and-open-questions.md`): Contradictions, unclear paths, missing docs, and open questions.
- FACT (`docs/repo-audit/06-onboarding.md`): 30/60/120 minute onboarding route.
- FACT (`docs/repo-audit/07-doc-folder-audit.md`): Complete inventory and assessment of documentation files under `src/`.

## Recommended Reading Order

1. FACT: Start with `00-index.md`.
2. FACT: Read `07-doc-folder-audit.md` next, because it records the primary documentation sources that refine the earlier audit.
3. FACT: Read `01-repo-map.md` for the physical repository structure.
4. FACT: Read `02-architecture.md` and `03-workflows.md` together; the documented S-layer semantics and the runnable scripts do not always align perfectly.
5. FACT: Read `04-project-ideas.md` for conceptual intent and didactic framing.
6. FACT: Finish with `05-risks-and-open-questions.md` and `06-onboarding.md`.

## Primary Source Set Read for the Documentation-Folder Audit

- FACT: The required inventory commands were run:
  - `find src -type d \( -iname "doc" -o -iname "docs" \) | sort`
  - `find src -type f \( -iname "*.md" -o -iname "*.qmd" -o -iname "README*" \) | sort`
- FACT: Documentation folders found and read as source contexts: `src/S1_observation/doc`, `src/S2_features/doc`, `src/S3_structure/doc`, `src/S4_signatures/doc`, `src/S5-decision/doc`.
- FACT: Documentation files read explicitly:
  - `src/README_1st.md`
  - `src/S1_observation/doc/S1-meta.md`
  - `src/S2_features/doc/S2-meta.qmd`
  - `src/S3_structure/archive/meta.md`
  - `src/S3_structure/archive/segmentation-meta.md`
  - `src/S3_structure/doc/S3-meta.qmd`
  - `src/S4_signatures/doc/04-3_signatures_RF-segs_classification.qmd`
  - `src/S4_signatures/doc/S4-meta.qmd`
  - `src/S5-decision/05-9_decision_selector.qmd`
  - `src/S5-decision/doc/S5-meta.qmd`
  - `src/_core/meta.md`
  - `src/r-libs/meta.md`
  - `src/run/meta.md`
  - `src/tools/meta.md`
