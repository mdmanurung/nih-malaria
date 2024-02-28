step_filter_cells <- function(batch, path_raw, PARAMS) {
    assertthat::assert_that(fs::dir_exists(path_raw))

    # batch <- "B1"
    # path_raw <- path_df$path_abs_raw[1]

    ncores <- ifelse(parallel::detectCores() >= 4, 4, 1)
    BPP <- BiocParallel::MulticoreParam(ncores)

    # read data
    sce <- read10xCounts(path_raw, BPPARAM = BPP)
    sce <- sce[, Matrix::colSums(counts(sce)) > 0] # remove barcodes with 0 counts

    # split count matrices
    rowData(sce)$Type[rowData(sce)$ Type == "Gene Expression"] <- "RNA"
    rowData(sce)$Type[rowData(sce)$ Type == "Antibody Capture"] <- "ADT"
    sce <- splitAltExps(sce, rowData(sce)$Type)

    # get library sizes for RNA and ADT
    md <- data.frame(
        Barcode = sce$Barcode,
        rna.size = Matrix::colSums(counts(sce)),
        prot.size = Matrix::colSums(counts(altExp(sce, "ADT"))),
        n.gene = Matrix::colSums(counts(sce) > 0)
    )

    # pre-filtering to reduce computational load
    keep <- md$rna.size > 10 & md$prot.size > 10
    md <- md[keep, ]
    sce <- sce[, keep]

    # call non-empty droplets
    set.seed(42)
    edrops <- DropletUtils::emptyDrops(counts(sce),
        niters = PARAMS$EDROPS$ITER,
        lower = PARAMS$EDROPS$LOWER,
        test.ambient = TRUE,
        BPPARAM = BPP
    )
    is.cell <- edrops$FDR <= PARAMS$EDROPS$FDR # barcodes of cells


    # further split ADT matrix to ADT and HTO matrix
    sce <- swapAltExp(sce, "ADT")
    rowData(sce)$Type[!grepl("ADT", rowData(sce)$ID)] <- "HTO"
    sce <- splitAltExps(sce, rowData(sce)$Type)
    sce <- swapAltExp(sce, "RNA")


    ## get empty droplets as well
    ## useful for denoising ADT with dsb R package
    droplet.barcode <- md |>
        mutate(rna.size = log10(rna.size), prot.size = log10(prot.size)) |>
        dplyr::filter(rna.size <= PARAMS$EDROPS$BG[[batch]]$upper.rna) |>
        dplyr::filter(prot.size <= PARAMS$EDROPS$BG[[batch]]$upper.adt & prot.size >= PARAMS$EDROPS$BG[[batch]]$lower.adt) |>
        pull(Barcode)
    droplet.barcode.ix <- which(droplet.barcode %in% sce$Barcode)

    droplet.plot <- md |>
        ggplot(aes(log10(prot.size), log10(rna.size))) +
        scattermore::geom_scattermore() +
        annotate("rect",
            xmin = PARAMS$EDROPS$BG[[batch]]$lower.adt,
            xmax = PARAMS$EDROPS$BG[[batch]]$upper.adt,
            ymin = 1, ymax = PARAMS$EDROPS$BG[[batch]]$upper.rna,
            alpha = .15, fill = "blue"
        )
    cowplot::ggsave2(glue::glue("results/step_filter/droplet_plots_{batch}.png"),
        droplet.plot,
        width = 2.5, height = 2.5, dpi = 400
    )


    # final SCE object
    sce.cell <- sce[, which(is.cell)]
    sce.background <- sce[, droplet.barcode.ix]
    sce.cell$capture_id <- batch
    sce.background$capture_id <- batch

    # add ambient profile
    sums.adt <- Matrix::rowSums(counts(altExp(sce.background, "ADT")))
    sums.adt <- sums.adt / sum(sums.adt)
    sums.rna <- Matrix::rowSums(counts(sce.background))
    sums.rna <- sums.rna / sum(sums.rna)
    rowData(sce.cell)$ambient <- sums.rna
    rowData(altExp(sce.cell, "ADT"))$ambient <- sums.adt

    # tidy up feature names
    ## RNA
    rna.unique <- uniquifyFeatureNames(rowData(sce.cell)$ID, rowData(sce.cell)$Symbol)
    rownames(sce.cell) <- rowData(sce.cell)$Symbol.unique <- rna.unique

    ## ADT: clean ADT names
    rd_adt <- SummarizedExperiment::rowData(SingleCellExperiment::altExp(sce.cell, "ADT"))
    is.isotype <- grepl("[Ii]sotype", rd_adt$Symbol)
    adt.recoded <- stringr::str_remove(rd_adt$Symbol, "anti-([Hh]uman)") |>
        stringr::str_remove("^Hu\\.") |>
        stringr::str_remove("^HuMs\\.") |>
        stringr::str_remove("^HuMsRt\\.") |>
        stringr::str_replace_all("-", "_") |>
        make.names()
    ix <- stringr::str_which(adt.recoded, "^CD")
    adt.recoded[ix] <- stringr::str_extract(adt.recoded[ix], "CD[0-9]+")


    rowData(altExp(sce.cell, "ADT"))$is.isotype <- is.isotype
    rownames(altExp(sce.cell, "ADT")) <- rowData(altExp(sce.cell, "ADT"))$Symbol.recoded <- adt.recoded
    rownames(altExp(sce.background, "ADT")) <- rowData(altExp(sce.background, "ADT"))$Symbol.recoded <- adt.recoded

    # save barcodes for vireo
    fwrite(data.table(sce.cell$Barcode),
        here("results", "barcodes_filtered", glue("barcodes_{batch}.tsv.gz")),
        col.names = FALSE
    )

    list(cells = sce.cell, background = sce.background, droplet.plot = droplet.plot)
}
