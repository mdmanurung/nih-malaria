step_demultiplex <- function(sce_filter, batch) {
    require(deMULTIplex2)
    require(DropletUtils)

    ncores <- ifelse(parallel::detectCores() >= 6, 6, 1)
    BPP <- BiocParallel::MulticoreParam(ncores)
    future::plan("multicore", workers = ncores - 1)

    # get object
    sce.cells <- sce_filter$cells
    sce.bg <- sce_filter$background


    # genotype demultiplexing ------------------------------------------------------

    message("Adding souporcell results...")
    path <- glue::glue("results/souporcell/soup_{batch}/clusters.tsv")
    df <- data.table::fread(path) |>
        mutate(soup_id = case_when(
            status == "singlet" ~ assignment,
            status == "doublet" ~ "Doublet",
            status == "unassigned" ~ "Negative"
        )) |>
        dplyr::select(soup_id,
            soup_posteriorSinglet = singlet_posterior,
            soup_posteriorDoublet = doublet_posterior,
            soup_logProbSinglet = log_prob_singleton,
            soup_logProbDoublet = log_prob_doublet
        )

    colData(sce.cells) <- cbind(colData(sce.cells), df)


    # HTO ---------------------------------------------------------------------

    message("Running demultiplex2")
    hto.mat <- t(counts(altExp(sce.cells, "HTO")))
    demux.res <- demultiplexTags(hto.mat, plot.umap = "none", plot.diagnostics = FALSE)
    sce.cells$demux_id <- demux.res$final_assign


    # DONE ---------------------------------------------------------------------
    
    list(cells = sce.cells)
}
