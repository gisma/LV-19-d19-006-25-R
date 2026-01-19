# src/_core/02-helper.R
# Minimal clone of link2GI::createFolders behaviour.

createFolders_simple <- function(root_folder, folders, create_folders = TRUE) {
  folders <- lapply(folders, function(f) file.path(root_folder, f))
  folders <- folders[!duplicated(folders)]
  names(folders) <- basename(unlist(folders))
  
  tmplt <- unlist(folders)
  while (any(duplicated(names(folders)))) {
    tmplt <- dirname(tmplt)
    dpl <- which(duplicated(names(folders), fromLast = FALSE) | duplicated(names(folders), fromLast = TRUE))
    names(folders)[dpl] <- paste(basename(tmplt)[dpl], names(folders)[dpl], sep = "_")
  }
  
  if (create_folders) {
    for (f in folders) if (!dir.exists(f)) dir.create(f, recursive = TRUE)
  }
  folders
}



burgwald_outputs <- function() {
  utils::read.delim(
    file = file.path(here::here(), "src","_core", "outputs.tsv"),
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "#"
  )
}



compile_outputs <- function(root_folder, outputs_df) {
  rel <- ifelse(
    outputs_df$ext == "dir",
    file.path("data/productive", outputs_df$S, outputs_df$topic, outputs_df$key),
    file.path("data/productive", outputs_df$S, outputs_df$topic,
              paste0(outputs_df$key, ".", outputs_df$ext))
  )
  abs <- file.path(root_folder, rel)
  names(abs) <- outputs_df$key
  abs
}