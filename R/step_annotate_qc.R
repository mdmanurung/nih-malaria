step_annotate_qc <- function(sce.cells) {

    ncores <- ifelse(parallel::detectCores() >= 6, 6, 1)
    BPP <- BiocParallel::MulticoreParam(ncores)

    # unpack object
    # sce.cells <- sce_demuxed$cells

    ## uniquify barcodes ------------------------------------------------------
    sce.cells$Barcode <- paste0(sce.cells$capture_id, "_", sce.cells$Barcode)

    ## doublet annotation ------------------------------------------------------

    message("Annotate doublets using scDblFinder...")
    assertthat::assert_that("soup_id" %in% colnames(colData(sce.cells)))
    doublet_ratio <- mean(sce.cells$soup_id == "Doublet")
    sce.cells <- scDblFinder::scDblFinder(sce.cells, dbr = doublet_ratio, BPPARAM = BPP)

    ## qc metric annotation ----------------------------------------------------

    ### GEX  -------------------------------------------------------------------
    message("Annotate QC metrics...")
    is.mito <- grepl("^MT-", rowData(sce.cells)$Symbol)
    is.ribo <- grepl("^RP[SL]", rowData(sce.cells)$Symbol)

    cellqc <- scuttle::perCellQCMetrics(sce.cells, subsets = list(Mito = is.mito, Ribo = is.ribo), BPPARAM = BPP)
    cellqc$qc_complexity <- cellqc$detected / cellqc$sum
    geneqc <- cellqc |>
        as.data.frame() |>
        dplyr::select(
            qc_libSize = sum,
            qc_nGene = detected,
            qc_complexity,
            qc_pctMito = subsets_Mito_percent,
            qc_pctRibo = subsets_Ribo_percent
        ) |>
        dplyr::mutate(Barcode = sce.cells$Barcode, .before = qc_libSize)
    colData(sce.cells) <- merge(colData(sce.cells), geneqc, by = "Barcode", sort = FALSE)


    ### ADT --------------------------------------------------------------------
    sce.cells <- SingleCellExperiment::swapAltExp(sce.cells, "ADT")
    controls <- grep("[Ii]sotype", rownames(sce.cells))
    qc.stats <- DropletUtils::cleanTagCounts(sce.cells, controls = controls)

    sce.cells$qc_zeroAmbient <- qc.stats$zero.ambient
    sce.cells$qc_highIsotypes <- qc.stats$high.controls

    sce.cells <- swapAltExp(sce.cells, "RNA")


    #### Annotate ----------------------------------------------------------------

    message("Running SingleR...")

    monaco <- celldex::MonacoImmuneData()
    singler.fine.monaco <- SingleR::SingleR(
        test = sce.cells, assay.type.test = 1,
        ref = list(MNC = monaco),
        labels = list(monaco$label.fine),
        BPPARAM = BPP
    )
    sce.cells$singler_monaco_fine <- singler.fine.monaco$pruned.labels


    singler.l2 <- qs::qread("data/singlecell_reference/singlerTrained_azimuth_L2.qs")
    pred.l2 <- SingleR::classifySingleR(sce.cells, singler.l2, assay.type=1, BPPARAM=BPP)
    sce.cells$singler_azimuth_l2 <- pred.l2$pruned.labels
    sce.cells$singler_azimuth_l1 <- case_when(
        sce.cells$singler_azimuth_l2 %in% c("B intermediate", "B memory", "B naive", "Plasmablast") ~ "B cells",
        sce.cells$singler_azimuth_l2 %in% c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC") ~ "MonoDCs",
        sce.cells$singler_azimuth_l2 %in% c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "Treg") ~ "CD4+T cells",
        sce.cells$singler_azimuth_l2 %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM") ~ "CD8+T cells",
        sce.cells$singler_azimuth_l2 %in% c("dnT", "gdT", "MAIT") ~ "UnconvT cells",
        sce.cells$singler_azimuth_l2 %in% c("ILC", "NK Proliferating", "NK_CD56bright", "NK") ~ "NKILC",
        sce.cells$singler_azimuth_l2 %in% c("Doublet", "Platelet", "HSPC") ~ "Other"
    )


    sce.cells
}
