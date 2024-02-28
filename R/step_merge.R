
clean_string_function <- function(original_string) {
  # Check if the string starts with "Isotype"
  if (grepl("^Isotype_", original_string)) {
    # Return the original string if it starts with "Isotype"
    return(original_string)
  } else {
    # Remove everything after (and including) the underscore for other strings
    return(sub("_.+", "", original_string))
  }
}
cleanv <- Vectorize(clean_string_function, USE.NAMES=F)



step_merge_sce <- function(sce_list){

    # tar_load(starts_with("sce_qc"))
    # sce_list <- list(sce_singlet_B1, sce_singlet_B2)

    ncores <- ifelse(parallel::detectCores() >= 6, 6, 1)
    BPP <- BiocParallel::MulticoreParam(ncores)


    # merge GEX
    message("Merging GEX data")
    keep_coldata <- colnames(colData(sce_list[[1]]))
    captures <- sapply(sce_list, function(x) unique(colData(x)$capture_id))

    sce.merged <- scMerge::sce_cbind(sce_list,
                       method = "intersect",
                       cut_off_batch = 0,
                       cut_off_overall = 0,
                       exprs = c("counts", "logcounts"),
                       colData_names = keep_coldata)
    SingleCellExperiment::mainExpName(sce.merged) <- "GEX"
    SummarizedExperiment::rowData(sce.merged) <- SummarizedExperiment::rowData(sce_list[[1]])
    colnames(sce.merged) <- sce.merged$Barcode

    sce_rowdata <- lapply(sce_list, function(x) rowData(x))
    rowdata_rna <- rowData(sce_list[[1]])
    rowdata_rna$ambient <- NULL
    ambient_rna <- do.call(cbind, lapply(sce_rowdata, function(x) x$ambient))
    colnames(ambient_rna) <- paste0("ambient_", captures)
    rowdata_rna <- cbind(rowdata_rna, ambient_rna)
    rowData(sce.merged) <- rowdata_rna

    # merge ADT
    message("Merging ADT data")
    sce_list_adt <- lapply(sce_list, function(x) SingleCellExperiment::swapAltExp(x, "ADT"))
    sce_list_adt <- lapply(sce_list_adt, function(x){
      rownames(x) <- rowData(x)$Symbol
      x
      })

    sce.merged.adt <- scMerge::sce_cbind(sce_list_adt,
                           method = "intersect",
                           cut_off_batch = 0,
                           cut_off_overall = 0,
                           exprs = c("counts", "dsb"),
                           colData_names = keep_coldata)
    SingleCellExperiment::mainExpName(sce.merged.adt) <- "ADT"
    colnames(sce.merged.adt) <- sce.merged$Barcode

    sce_rowdata_adt <- lapply(sce_list_adt, function(x) rowData(x))
    rowdata_adt <- rowData(sce_list_adt[[1]])
    rowdata_adt$ambient <- NULL
    ambient_adt <- do.call(cbind, lapply(sce_rowdata_adt, function(x) x$ambient))
    colnames(ambient_adt) <- paste0("ambient_", captures)
    rowdata_adt <- cbind(rowdata_adt, ambient_adt)
    rowData(sce.merged.adt) <- rowdata_adt


    # store to main SCE
    SingleCellExperiment::altExp(sce.merged, "ADT") <- sce.merged.adt
    newname <- stringr::str_remove(rowdata_adt$Symbol, "Hu\\.|HuMs\\.|HuMsRt\\.") |>
      cleanv()
    newname[1:130] <- paste0(newname[1:130], '_adt')
    rownames(altExp(sce.merged, "ADT")) <- newname


    ## ADD sample annotation
      meta <- openxlsx::read.xlsx(here::here("data", "pilot", "donorIDs_and_corresponding_hashtags.xlsx"))
      meta$age <- naturalsort::naturalfactor(meta$age)
      ix <- match(sce.merged$demux_id, meta$sample.ID)
      sce.merged$trialID <- meta$trial[ix]
      sce.merged$age_group <- meta$age[ix]
      sce.merged$residence <- meta$X6[ix]

    ## multibatchnorm
    message("Running multibatchnorm")
    sce.merged <- batchelor::multiBatchNorm(
        sce.merged,
        batch = sce.merged$capture_id,
        normalize.all = TRUE,
        BPPARAM = BPP,
        norm.args=list(use_altexps=FALSE)
    )


    # tag qc metrics
    message("Tagging cell-level QC metrics")
    qcres.byLineage <- SummarizedExperiment::colData(sce.merged) |>
        as.data.frame() |>
        dplyr::mutate(
            out.libsize.l1 = as.logical(scuttle::isOutlier(qc_libSize, nmads=3, log=T, type="lower", batch=capture_id)),
            out.ngene.l1 = as.logical(scuttle::isOutlier(qc_nGene, nmads=3, log=T, type="lower", batch=capture_id)),
            out.mito.l1 = as.logical(scuttle::isOutlier(qc_pctMito, nmads=3, log=F, type="higher", batch=capture_id)),
            out.complexity.l1 = as.logical(scuttle::isOutlier(qc_complexity, nmads=3, log=F, type="lower", batch=capture_id)),
            .by = singler_azimuth_l1
        ) |>
        dplyr::select(Barcode, tidyselect::starts_with("out"))
    qcres.byLineage$zeroAmbient <- sce.merged$qc_zeroAmbient
    qcres.byLineage$highIsotypes <- sce.merged$qc_highIsotypes
    qcres.byLineage$discard <-
        with(qcres.byLineage, out.libsize.l1 | out.ngene.l1 | out.mito.l1 | out.complexity.l1 | zeroAmbient | highIsotypes)
    S4Vectors::metadata(sce.merged) <- qcres.byLineage


    sce.merged$discard.qc <- qcres.byLineage$discard

    sce.merged

}


step_merge_background <- function(filter_list){

    bg_list <- lapply(filter_list, function(filter) filter$background)

    merged.bg <- scMerge::sce_cbind(bg_list,
                                    exprs="counts",
                                    cut_off_batch = 0,
                                    cut_off_overall = 0,
                                    colData_names = "capture_id")
    SingleCellExperiment::mainExpName(merged.bg) <- "GEX"

    bg_list_adt <- lapply(bg_list, function(x) SingleCellExperiment::swapAltExp(x, "ADT"))
    merged.bg.adt <- scMerge::sce_cbind(bg_list_adt,
                                        exprs = "counts",
                                        cut_off_batch = 0,
                                        cut_off_overall = 0,
                                        colData_names = "capture_id")
    SingleCellExperiment::mainExpName(merged.bg.adt) <- "ADT"

    SingleCellExperiment::altExp(merged.bg, "ADT") <- merged.bg.adt
    merged.bg
}
