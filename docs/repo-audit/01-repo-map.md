# Repo-Map

## Verzeichnisbaum

```text
.
├── README.md
├── AGENTS.md
├── LV-19-d19-006-25-R.Rproj
├── src/
│   ├── README_1st.md
│   ├── _core/
│   ├── r-libs/
│   ├── tools/
│   ├── run/
│   ├── S1_observation/
│   ├── S2_features/
│   ├── S3_structure/
│   ├── S4_signatures/
│   └── S5-decision/
├── data/
│   ├── raw/
│   ├── productive/
│   ├── tmp/
│   └── deprec/
├── outputs/
│   ├── figures/
│   ├── logs/
│   └── tables/
├── metadata/
├── docs/
├── figures/
├── tmp/
└── ki_lit_recherche.*
```

## Wichtigste Ordner

- FACT (`src/README_1st.md`, "Repository Mapping"): `src/_core` ist als Infrastruktur, `src/S1_observation` bis `src/S5_decisions` als Pipeline-Schichten, `src/lib` als Libraries, `src/tools` als Tools und `src/run` als Orchestrierung beschrieben.
- FACT (tatsaechliche Struktur): Der Codeordner heisst `src/r-libs`, nicht `src/lib`; der S5-Codeordner heisst `src/S5-decision`, nicht `src/S5_decisions`.
- FACT (`src/_core/meta.md`): `_core` soll Projekt-Setup, Pfade, Toolchain-Registrierung, globale Optionen und Konstanten bereitstellen.
- FACT (`src/r-libs/meta.md`): `r-libs` ist fuer wiederverwendbare Funktionen gedacht.
- FACT (`src/tools/meta.md`): `tools` enthaelt externe Glue-Skripte und JavaScript-Helfer; konkrete Dateien sind `NDVI.js`, `evi.js`, `kndvi.js`, `lai.js`, `savi.js`.
- FACT (`src/run/meta.md`): `run` soll explizite Einstiegspunkte bereitstellen; tatsaechlich liegt dort nur `meta.md`, keine `run_S*_all.R`-Datei.
- FACT (`data/productive/*`): Produktive Artefakte folgen den Schichten `S1_observation`, `S2_features`, `S3_structure`, `S4_signatures`, `S5-decision` bzw. historische Varianten wie `S4_decision`, `S5_decisions`.
- FACT (`data/raw/*`): Rohdaten enthalten u.a. AOI-Daten, DGM/DOM-Gemeinden, DWD-Stationen, GADM, CLC und Anbieter-Unterordner.

## Zentrale Dateien und Rolle

- FACT (`README.md`, "Repository content and purpose"): Beschreibt das Repo als kursbezogene Referenzimplementation fuer Earth Observation, Terrain/GIS, Statistik/ML und Helper-Funktionen.
- FACT (`src/README_1st.md`): Definiert die Meta-Architektur "Spatial Research Pipeline" fuer optimiertes Niederschlagsnetzdesign im Burgwald.
- FACT (`LV-19-d19-006-25-R.Rproj`): RStudio-Projekt mit UTF-8-Encoding, 2 Leerzeichen fuer Tabs, Sweave und pdfLaTeX.
- FACT (`src/_core/00-folders.R`, `burgwald_folders()`): Definiert den Ordner-Skeleton fuer `data/raw`, `data/productive`, `outputs`, `metadata`, `docs`, `src`, `tmp`.
- FACT (`src/_core/01-setup-burgwald.R`): Laedt Pakete, erstellt Ordner, kompiliert Pfade aus `outputs.tsv`, definiert AOI, CLC-Legende, Temp-Verzeichnisse und `envrmt`.
- FACT (`src/_core/02-helper.R`): Definiert `createFolders_simple()`, `burgwald_outputs()` und `compile_outputs()`.
- FACT (`src/_core/outputs.tsv`): Zentrale Output-Registry; Schluessel werden in `paths[[...]]` auf absolute produktive Pfade abgebildet.
- FACT (`src/r-libs/01-fun-data-retrieval.R`): Enthaelt Datenabruf-Helfer wie `aoi_with_buffer()`, `get_osm_burgwald_by_key()`, `download_if_missing()`, `run_if_missing()`, DWD- und RADOLAN-Funktionen.
- FACT (`src/r-libs/metrics-fun.R`): Enthaelt OTB/MeanShift/PCA/ARI-Helfer fuer S3-Segmentierung sowie spaeter genutzte Decision-Helfer laut S5-Skripten.

## Konfiguration und Automatisierung

- FACT (`.gitignore`): Ignoriert u.a. `.Rproj.user`, `.Rhistory`, `data`, `tmp`, `figures`, `ki_lit_recherche.*`.
- FACT (`.renvignore`): Ignoriert u.a. `data/`, `tmp/`, `outputs/`, `docs/`, `figures/`, `*.tif`, `*.nc`, `*.rds`, `*.zip`, `*.gz`, `*.gpkg`.
- FACT (Inventur): Keine CI/CD-Dateien gefunden; `.github` existiert nicht.
- FACT (Inventur): Keine Build-Systemdateien wie `Makefile`, `justfile`, `Taskfile`, `Dockerfile`, Compose-Dateien oder Paketmanifestdateien gefunden.
- INFERENCE: Ausfuehrung ist skriptbasiert ueber `Rscript` oder interaktiv in RStudio vorgesehen, nicht ueber ein formales Build- oder Workflow-Tool.

