
step_DA_milo_age <- function(srt){

  n_neighbors <- 30
  dim <- 30

  srt <- srt[, srt$residence=="rural"]

  set.seed(42)
  srt <- Seurat::FindMultiModalNeighbors(
    srt,
    k.nn = n_neighbors,
    reduction.list = list("rnaMNN", "adtMNN"),
    dims.list = list(1:dim, 1:dim)) |>
    Seurat::RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

  DefaultAssay(srt) <- "RNA"
  cite_milo <- Milo(SingleCellExperiment(
    list(counts=srt@assays$RNA$counts),
    colData = srt@meta.data))
  miloR::graph(cite_milo) <- miloR::graph(buildFromAdjacency(srt@graphs$wsnn, k=n_neighbors))
  reducedDim(cite_milo, "WNN.UMAP") <- srt@reductions$wnn.umap@cell.embeddings

  set.seed(42)
  cite_milo <- makeNhoods(cite_milo,
                          k=n_neighbors,
                          refined=TRUE,
                          prop=0.1,
                          refinement_scheme="graph")
  cite_milo <- countCells(cite_milo,
                          meta.data = data.frame(colData(cite_milo)),
                          samples="demux_id")
cite_milo <- buildNhoodGraph(cite_milo)


  message("Running miloR...")
  design.mat <- srt@meta.data[,c("demux_id", "capture_id", "residence", "age_group")] |>
    mutate(age_group = case_when(
      age_group == "11-17" ~ "adolescents",
      age_group == "18-35" ~ "adults",
      TRUE ~ "older"
    )) |>
    mutate(age_group = ordered(age_group, levels=c("adolescents", "adults", "older"))) |>
    dplyr::select(demux_id, capture_id, age_group)
  design.mat <- distinct(design.mat)
  rownames(design.mat) <- design.mat$demux_id
  design.mat <- design.mat[colnames(nhoodCounts(cite_milo)), , drop=FALSE]


  # trend test
  set.seed(42)
  d_res <-
    testNhoods(
      cite_milo,
      design = ~ capture_id + age_group,
      design.df = design.mat,
      model.contrasts = "age_group.L",
      fdr.weighting = "neighbour-distance",
      reduced.dim = "WNN.UMAP",
      norm.method="logMS",
      robust = TRUE)
  d_res <- annotateNhoods(cite_milo, d_res, coldata_col = "singler_monaco_fine")

  
  # one vs rest
  mod <- model.matrix(~ 0 + age_group + capture_id, data=design.mat)
  mod.contrast <- makeContrasts(
    YoungVsRest = c("age_groupadolescents - (age_groupadults + age_groupolder)/2"),
    AdultsVsRest = c("age_groupadults - (age_groupadolescents + age_groupolder)/2"),
    OldVsRest = c("age_groupolder - (age_groupadults + age_groupadolescents)/2"),
    levels=mod
  )

  d_list <- lapply(1:3, function(i){
    d_res <-
      testNhoods(
        cite_milo,
        design = ~ 0 + age_group + capture_id,
        design.df = design.mat,
        model.contrasts = mod.contrast[,i],
        fdr.weighting = "neighbour-distance",
        reduced.dim = "WNN.UMAP",
        norm.method="logMS",
        robust = TRUE)
    d_res <- annotateNhoods(cite_milo, d_res, coldata_col = "singler_monaco_fine")

    d_res
  })
  names(d_list) <- c("Young", "Adults", "Older")


  list(milo=cite_milo, trend=d_res, one_vs_rest=d_list)
}
