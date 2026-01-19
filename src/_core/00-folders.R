# src/_core/00-folders.R
# Single source: folder skeleton as relative paths.

burgwald_folders <- function() {
  c(
    # raw
    "data/raw/",
    "data/raw/AOI_Burgwald/",
    "data/raw/AOI_Burgwald/dem/",
    "data/raw/AOI_Burgwald/clc/",
    "data/raw/AOI_Burgwald/osm_by_key/",
    "data/raw/providers/dwd/",
    "data/raw/providers/cdse/",
    "data/raw/providers/gdalcubes/",
    "data/raw/providers/osm/",
    "data/raw/providers/dwd/radolan/",
    
    # productive (canonical modelling artefacts)
    "data/productive/",
    "data/productive/S0_problem/",
    "data/productive/S1_observation/",
    "data/productive/S2_features/",
    "data/productive/S3_structure/",
    "data/productive/S4_signatures/",
    "data/productive/S5_decisions/",
    "data/productive/S5_validation/",
    
    # project misc
    "outputs/figures/",
    "outputs/tables/",
    "outputs/logs/",
    "metadata/",
    "docs/",
    "src/",
    "tmp/"
  )
}
