

step_keep_singlets_goodcells <- function(sce){


  # remove doublets called by HTO and unassigneds by either SNPs or HTO
  # i.e. use HTO as the "gold standard" of cell assignment

  # doublet.snps <- sce_qc_B1$soup_id %in% c("Doublet")
  doublet.htos <- sce$demux_id %in% c("multiplet")
  unassigned.both <- (sce$demux_id == "negative" | sce$soup_id == "Negative")

  # remove doublets according to HTO
  discard <- !(doublet.htos | unassigned.both)
  sce_singlet <- sce[, discard]

  sce_singlet
}
