################################################################################
# CDH1: mutation-driven vs dosage-driven inactivation
# Supports the manuscript claim that CDH1 is not dosage-coupled because it is
# inactivated by (truncating) mutation rather than copy-number loss.
# Contrast: ERBB2, whose expression IS driven by copy number (dosage).
# Prereq: run mbc_ALL_IN_ONE.R first. Uses: mutations (raw MAF), mutation_matrix,
#   rna_matrix_final, cna_matrix_final, samp0, theme_publication, save_pdf, sub_cols
################################################################################
library(ggplot2)
if (!dir.exists("output")) dir.create("output")
samp <- samp0                                     # samples with CNA + RNA

trunc_classes <- c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
                   "Splice_Site","Translation_Start_Site","Nonstop_Mutation")

## ---- per-sample status for a gene -----------------------------------------
gene_status <- function(gene){
  # expression & copy number
  rna <- if(gene %in% rownames(rna_matrix_final)) as.numeric(rna_matrix_final[gene,samp]) else rep(NA,length(samp))
  cna <- if(gene %in% rownames(cna_matrix_final)) as.numeric(cna_matrix_final[gene,samp]) else rep(NA,length(samp))
  # any mutation (from binary matrix, if the gene is in it)
  mut <- rep(NA_real_, length(samp))
  if(gene %in% colnames(mutation_matrix)){
    inM <- samp %in% rownames(mutation_matrix)
    mut[inM] <- mutation_matrix[samp[inM], gene]
  }
  # truncating mutation (from raw MAF)
  maf_g <- mutations[mutations$Hugo_Symbol==gene,]
  trunc_samples <- unique(maf_g$Tumor_Sample_Barcode[maf_g$Variant_Classification %in% trunc_classes])
  trunc <- as.numeric(samp %in% trunc_samples)
  data.frame(SAMPLE_ID=samp, rna=rna, cna=cna, mut=mut, trunc=trunc)
}

## ---- variance explained: mutation vs copy number --------------------------
compare_drivers <- function(gene){
  d <- gene_status(gene); d <- d[complete.cases(d[,c("rna","cna","mut")]),]
  if(nrow(d) < 15){ cat("  ",gene,": too few complete samples\n"); return(invisible()) }
  r2 <- function(f) summary(lm(f, data=d))$r.squared
  cat(sprintf("\n[%s]  n=%d | mutated=%d (truncating=%d)\n", gene, nrow(d), sum(d$mut>0), sum(d$trunc>0)))
  cat(sprintf("  expr ~ copy-number : R2 = %.3f\n", r2(rna~cna)))
  cat(sprintf("  expr ~ mutation    : R2 = %.3f\n", r2(rna~mut)))
  cat(sprintf("  expr ~ CN + mut    : R2 = %.3f\n", r2(rna~cna+mut)))
  # which effect is significant in the joint model?
  cf <- summary(lm(rna~cna+mut, data=d))$coefficients
  cat("  joint model coefficients (p-values):\n")
  print(round(cf[,c("Estimate","Pr(>|t|)")],4))
  # expression: mutated vs WT
  if(sum(d$mut>0)>=3 && sum(d$mut==0)>=3){
    w <- suppressWarnings(wilcox.test(rna ~ (mut>0), data=d))
    cat(sprintf("  expr mutated vs WT: Wilcoxon p = %.4f (median WT %.2f vs mut %.2f)\n",
                w$p.value, median(d$rna[d$mut==0]), median(d$rna[d$mut>0])))
  }
  invisible(d)
}

cat("==== CDH1 (expected: mutation-driven) vs ERBB2 (expected: dosage-driven) ====\n")
dC <- compare_drivers("CDH1")
dE <- compare_drivers("ERBB2")

## ---- CDH1 mutation spectrum (how many are truncating / LoF?) ---------------
cat("\n==== CDH1 mutation spectrum (all biopsies in MAF) ====\n")
cdh1_maf <- mutations[mutations$Hugo_Symbol=="CDH1",]
if(nrow(cdh1_maf)>0){
  tb <- sort(table(cdh1_maf$Variant_Classification), decreasing=TRUE)
  print(tb)
  cat(sprintf("  Truncating/LoF fraction: %.0f%%\n",
              100*mean(cdh1_maf$Variant_Classification %in% trunc_classes)))
}

## ---- FIGURE: CDH1 expression by mutation status ---------------------------
if(!is.null(dC) && nrow(dC)>0){
  dC$grp <- ifelse(dC$trunc>0,"Truncating mut", ifelse(dC$mut>0,"Other mut","Wild-type"))
  dC$grp <- factor(dC$grp, levels=c("Wild-type","Other mut","Truncating mut"))
  dC2 <- dC[!is.na(dC$rna),]
  if(nlevels(droplevels(dC2$grp))>=2){
    pM <- ggplot(dC2, aes(grp, rna, fill=grp)) +
      geom_boxplot(outlier.shape=NA, alpha=.75) + geom_jitter(width=.15, size=1.6, alpha=.6) +
      scale_fill_manual(values=c("Wild-type"="#4DBBD5","Other mut"="#F39B7F","Truncating mut"="#E64B35"), guide="none") +
      labs(title="CDH1 expression by mutation status",
           subtitle="Lower expression in mutated tumors supports mutation-driven inactivation",
           x=NULL, y="CDH1 (log2 expression)") + theme_publication(12)
    save_pdf(pM, "17_cdh1_mutation.pdf", 7, 6)
  }
}

cat("\n==== INTERPRETATION ====\n")
cat(" If for CDH1: R2(mutation) > R2(copy-number) and mutated<WT expression,\n")
cat("   while for ERBB2: R2(copy-number) >> R2(mutation),\n")
cat("   -> supports 'CDH1 inactivation is mutation-driven, not dosage-driven'.\n")
cat(" (Methylation is not in MBCproject; cite literature for the epigenetic component.)\n")
################################################################################
