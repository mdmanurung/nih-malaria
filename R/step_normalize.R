step_normalize <- function(sce_cells, sce_filter) {
    library(scuttle)
    library(dsb)

    ncores <- ifelse(parallel::detectCores() >= 4, 4, 1)
    BPP <- BiocParallel::MulticoreParam(ncores)

    sce.cells <- sce_cells$cells
    sce.bg <- sce_filter$background

    ## normalize RNA -----------------------------------------------------------
    message("Normalizing RNA")
    if (SingleCellExperiment::mainExpName(sce.cells) != "RNA") {
        sce.cells <- SingleCellExperiment::swapAltExp(sce.cells, "RNA")
    }

    set.seed(42)
    clust <- scran::quickCluster(sce.cells, BPPARAM = BPP)
    sce.cells <- scuttle::computePooledFactors(sce.cells, clusters = clust, BPPARAM = BPP, min.mean = 0.1)
    sce.cells <- scuttle::logNormCounts(sce.cells)


    ## normalize ADT  ----------------------------------------------------------
    message("Normalizing ADT")

    sce.cells <- SingleCellExperiment::swapAltExp(sce.cells, "ADT")
    if (SingleCellExperiment::mainExpName(sce.bg) != "ADT") {
        sce.bg <- SingleCellExperiment::swapAltExp(sce.bg, "ADT")
    }

    isotype.controls <- rownames(sce.cells)[grepl("[Ii]sotype", rownames(sce.cells))]
    set.seed(42)
    dsb.norm <- dsb::DSBNormalizeProtein(
        cell_protein_matrix = counts(sce.cells),
        empty_drop_matrix = counts(sce.bg),
        denoise.counts = TRUE,
        use.isotype.control = TRUE,
        isotype.control.name.vec = isotype.controls,
        quantile.clipping = TRUE,
        quantile.clip = c(0.001, 0.9995)
    )
    SummarizedExperiment::assay(sce.cells, "dsb") <- dsb.norm


    ## done --------------------------------------------------------------------
    message("Normalization finished")

    sce.cells <- SingleCellExperiment::swapAltExp(sce.cells, "RNA")
    sce.cells
}
