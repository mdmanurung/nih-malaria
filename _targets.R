library(targets)
library(tarchetypes)
library(tidyselect)
library(tidyverse)
library(future)
library(future.batchtools)
library(future.callr)
library(here)
library(data.table)

# options
set_resources <- function(resources) {
  # helper-function to set resources for SLURM.
  tar_resources(
    future = tar_resources_future(
      plan = tweak(
        batchtools_slurm,
        template = "configs/batchtools.slurm.tmpl",
        resources = resources
      )
    )
  )
}
tar_option_set(
  packages = c(
    "here", "tidyverse", "tibble", "qs", "data.table", "future",
    "patchwork", "scattermore", "ggbeeswarm",
    "glue", "fs", "assertthat",
    "SummarizedExperiment", "SingleCellExperiment",
    "DropletUtils", "SingleR", "scuttle", "scDblFinder",
    "Seurat", "SeuratObject", "SeuratWrappers",
    "scater", "scran", "miloR", "batchelor",
    "Matrix",
    "Trex", "scRepertoire"
  ),
  format = "qs", # default storage format
  memory = "transient",
  garbage_collection = TRUE,
  error = "abridge"
)
future::plan(future.batchtools::batchtools_slurm, template = "configs/batchtools.slurm.tmpl")

# parameters
params <- yaml::read_yaml(here::here("configs", "pipeline_parameters.yaml"))
resources <- yaml::read_yaml(here::here("configs", "resources.yaml"))
path_df <- data.table::fread("data/metadata/path_files.tsv")

# functions
tar_source(here::here("R"))

# lists of steps
filter_droplets <- list(
  tar_map(
    unlist = FALSE,
    values = path_df,
    tar_target(
      sce_filt,
      command = step_filter_cells(batch, path_abs_raw, params),
      resources = set_resources(resources$small)
    ),
    tar_target(
      sce_demux,
      command = step_demultiplex(sce_filt, batch),
      resources = set_resources(resources$small)
    ),
    tar_target(
      sce_norm,
      command = step_normalize(sce_demux, sce_filt),
      resources = set_resources(resources$medium),
    ),
    tar_target(
      sce_qc,
      command = step_annotate_qc(sce_norm),
      resources = set_resources(resources$small)
    ),
    tar_target(
      sce_singlets,
      command = step_keep_singlets_goodcells(sce_qc),
    ),
    names = starts_with("batch")
  )
)


merge_objects <- list(
  tar_combine(
    name = sce_merged_bg,
    tar_select_targets(filter_droplets, starts_with("sce_filt")),
    command = step_merge_background(list(!!!.x)),
    deployment = "main"
  ),
  tar_combine(
    name = sce_merged,
    tar_select_targets(filter_droplets, starts_with("sce_singlets_")),
    command = step_merge_sce(list(!!!.x)),
    deployment = "main"
  )
)

integrate <- list(
  tar_target(
    srt_mnn,
    step_integrate_fastmnn(sce_merged),
    resources = set_resources(resources$small)
  )
)

diff_abundance <- list(
  tar_target(
    milo_age,
    step_DA_milo_age(srt_mnn),
    resources = set_resources(resources$medium)
  )
)

repertoire_analysis <- list(
  tar_target(
    srt_trex,
    step_trex(srt_mnn),
    resources = set_resources(resources$large)
  )
)


# pipeline
list(
  filter_droplets,
  merge_objects,
  integrate,
  diff_abundance,
  repertoire_analysis
)

