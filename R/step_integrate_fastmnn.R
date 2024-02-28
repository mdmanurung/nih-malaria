

step_integrate_fastmnn <- function(sce){

  ncores <- ifelse(parallel::detectCores() >= 6, 6, 1)
  BPP <- BiocParallel::MulticoreParam(ncores)
  future::plan("multicore")


  # remove low quality cells
  sce <- sce[, !sce$discard.qc]

  message("Correcting RNA...")
  gene.vars <- modelGeneVar(
    sce,
    block = sce$capture_id
  )
  hvgs <- getTopHVGs(gene.vars, n=4000)

  set.seed(42)
  gene.mnn <- fastMNN(
    sce,
    batch = sce$capture_id,
    k = 20,
    d = 50,
    get.variance = TRUE,
    subset.row = hvgs,
    BSPARAM = BiocSingular::IrlbaParam(),
    BNPARAM = BiocNeighbors::AnnoyParam(),
    BPPARAM = BPP
  )
  rownames(sce) <- gsub("_", "-", rownames(logcounts(sce)))

  srt <- CreateSeuratObject(counts = counts(sce),
                            meta.data = as.data.frame(colData(sce)))
  srt@assays$RNA$data <- logcounts(sce)

  srt[["rnaMNN"]] <- SeuratObject::CreateDimReducObject(
    embeddings = reducedDim(gene.mnn, "corrected"),
    key="rnaMNN_", assay="RNA")
  srt <- srt |>
    FindNeighbors(reduction="rnaMNN", dims=1:50, return.neighbor = TRUE, annoy.metric="cosine") |>
    RunUMAP(nn.name="RNA.nn", reduction.name="rnaUMAP")



  ## ADT
  message("Correcting ADT")
  set.seed(42)
  sce <- swapAltExp(sce, "ADT")

  adt.mnn <- fastMNN(
    sce,
    batch = sce$capture_id,
    assay.type = "dsb",
    k = 20,
    d = 50,
    get.variance = TRUE,
    subset.row = !grepl("Isotype", rownames(sce)),
    BSPARAM = BiocSingular::IrlbaParam(),
    BNPARAM = BiocNeighbors::AnnoyParam(),
    BPPARAM = BPP
  )

  srt[["ADT"]] <-  SeuratObject::CreateAssay5Object(
    counts = counts(sce),
    data = assay(sce, "dsb")
  )
  srt[["adtMNN"]] <- SeuratObject::CreateDimReducObject(
    embeddings = reducedDim(adt.mnn, "corrected"),
    key="adtMNN_", assay="ADT")

  DefaultAssay(srt) <- "ADT"
  srt <- srt |>
    FindNeighbors(reduction="adtMNN", dims=1:50, return.neighbor = TRUE, annoy.metric="cosine") |>
    RunUMAP(nn.name="ADT.nn", reduction.name="adtUMAP")


  message("Running WNN...")
  srt <- FindMultiModalNeighbors(
    srt,
    reduction.list = list("rnaMNN", "adtMNN"),
    dims.list = list(1:40, 1:40),
    modality.weight.name = "RNA.weight") |>
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnnUMAP", reduction.key = "wnnUMAP_") |>
    FindClusters(graph.name="wsnn", resolution=seq(.1, 2, .1))


  srt
}

