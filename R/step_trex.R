

step_trex <- function(srt){

  future::plan("multicore", workers = 6)

  message("Combine contigs...")
  path_b1 <- "data/pilot/pilot_pbmc1/per_sample_outs/PBMC1_count/vdj_t/filtered_contig_annotations.csv"
  path_b2 <- "data/pilot/pilot_pbmc2/per_sample_outs/PBMC2_count/vdj_t/filtered_contig_annotations.csv"

  contigs.b1 <- data.table::fread(path_b1) |>
    combineTCR(samples="B1", removeNA=TRUE, removeMulti=TRUE, filterMulti=TRUE)
  contigs.b2 <- data.table::fread(path_b2) |>
    combineTCR(samples="B2", removeNA=TRUE, removeMulti=TRUE, filterMulti=TRUE)
  contigs <- rbind(contigs.b1$B1, contigs.b2$B2)

  ix <- contigs$barcode %in% srt$Barcode
  contigs.filter <- contigs[ix, ]

  srt2 <- combineExpression(contigs.filter, srt)
  srt2 <- srt2[, !is.na(srt2$CTaa)]


  ## rerun MNN without TCR genes
  message("Correcting RNA data...")
  DefaultAssay(srt2) <- "RNA"
  srt2[["RNA"]] <- split(srt2[["RNA"]], f=srt2$capture_id)
  srt2 <- srt2 |>
    FindVariableFeatures(nfeatures=3000)
  VariableFeatures(srt2) <- VariableFeatures(srt2)[!grepl("TR[ABDG][VDJ]", VariableFeatures(srt2))]
  srt2 <- srt2 |>
    ScaleData(do.scale=F, do.center=F) |>
    RunPCA()

  srt2 <- IntegrateLayers(
    object = srt2, method = FastMNNIntegration,
    new.reduction = "newMNN",
    verbose = FALSE
  )
  srt2 <- JoinLayers(srt2)
  srt2 <- srt2 |>
    RunUMAP(reduction="newMNN", dims=1:30)

  ## Trex
  message("Trex embedding...")
  Trex_vectorsb <- maTrex(srt2,
                          chains = "TRB",
                          encoder.model = "VAE",
                          encoder.input = "KF")
  srt2[["trex_kf"]] <- CreateDimReducObject(embeddings=as.matrix(Trex_vectorsb))
  srt2 <- Seurat::FindMultiModalNeighbors(
    srt2,
    k.nn = 50,
    reduction.list = list("newMNN", "adtMNN", "trex_kf"),
    dims.list = list(1:30, 1:30, 1:30)
  ) |>
    RunUMAP(nn.name = "weighted.nn",
            reduction.name = "wnn.umap",
            reduction.key = "wnnUMAP_") |>
    FindClusters(graph.name="wsnn", resolution=seq(.1, 1, .1))

  message("Done!")
  srt2

}
