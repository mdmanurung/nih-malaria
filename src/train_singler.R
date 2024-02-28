
tar_load(srt_integ)
ref.sce <- qs::qread("data/singlecell_reference/pbmc_multimodal_sce.qs")

sumexp <- SummarizedExperiment(list(logcounts=logcounts(ref.sce)),
                               colData = colData(ref.sce),
                               rowData = data.frame(genes=rownames(ref.sce)))
sumexp$Barcode <- colnames(sumexp)
common <- intersect(rownames(sce_norm_B1), rownames(sumexp))


set.seed(42)
cd <- as.data.table(colData(sumexp))
l2.subsets <- cd %>% dplyr::count(celltype.l2) %>% dplyr::filter(n>100) %>% pull(celltype.l2)
ix_l2 <- cd %>%
  dplyr::filter(celltype.l2 %in% l2.subsets) %>%
  slice_sample(n=1000, by=celltype.l2) %>%
  pull(Barcode)
ix2 <- which(sumexp$Barcode %in% ix_l2)
sumexp_subset <- sumexp[common, ix2]

ncores <- ifelse(parallel::detectCores() >= 6, 6, 1)
BPP <- BiocParallel::MulticoreParam(ncores)



set.seed(42); trained.2 <-
  trainSingleR(sumexp_subset,
               labels=ref.sce$celltype.l2[ix2],
               de.method = "wilcox",
               de.n = 10,
               aggr.ref=F,
               BPPARAM=BPP)

qs::qsave(trained.2, "data/singlecell_reference/singlerTrained_azimuth_L2.qs")

trained.2 <- qs::qread("data/singlerTrained_azimuth_L2.qs")

sce <- as.SingleCellExperiment(srt_integ)
pred2 <- SingleR::classifySingleR(sce_norm_B1, trained.2, assay.type=1, BPPARAM=BPP)
qs::qsave(pred2, "singler_l2.qs")
