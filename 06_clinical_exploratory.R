################################################################################
# CLINICAL-PATHOLOGICAL ASSOCIATIONS (no survival data available in MBCproject)
# Links the dosage-sensitive drivers to clinical/pathological phenotypes.
# Prereq: run mbc_ALL_IN_ONE.R first, so these are in memory:
#   clinical_matched, cna_matrix_final, rna_matrix_final, Y_final, dsc, samp0,
#   theme_publication, save_pdf, sub_cols
################################################################################
library(ggplot2)
if (!dir.exists("output")) dir.create("output")

## ---- helpers --------------------------------------------------------------
clean_clin <- function(x){ x <- as.character(x)
  x[toupper(trimws(x)) %in% c("NOT_FOUND_IN_RECORD","ABSTRACTION_PENDING","DON'T_KNOW",
    "DONT_KNOW","UNKNOWN","NA","","NOT AVAILABLE","NOT APPLICABLE")] <- NA; x }
getclin <- function(col, ids){
  if(!col %in% colnames(clinical_matched)) return(rep(NA, length(ids)))
  clean_clin(clinical_matched[[col]][match(ids, clinical_matched$SAMPLE_ID)]) }
grade_num <- function(x){ x <- toupper(clean_clin(x)); g <- rep(NA_real_, length(x))
  g[grepl("III|POOR|GRADE 3|\\b3\\b", x)] <- 3
  g[is.na(g) & grepl("\\bII\\b|MODER|GRADE 2|\\b2\\b", x)] <- 2
  g[is.na(g) & grepl("\\bI\\b|WELL|GRADE 1|\\b1\\b", x)] <- 1
  g }
rho_line <- function(tag, x, y){ ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<10){ cat(sprintf("  %-28s n<10, skipped\n", tag)); return(invisible()) }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman"))
  cat(sprintf("  %-28s rho=%+.2f | p=%.3f | n=%d\n", tag, ct$estimate, ct$p.value, sum(ok))) }

samp <- samp0                                   # samples with both CNA + RNA

## ===========================================================================
## C1. CDH1 dosage sensitivity vs LOBULAR histology (ILC)
##     CDH1 (E-cadherin) loss defines invasive lobular carcinoma -> strong,
##     non-circular, biologically motivated validation of the CDH1 finding.
## ===========================================================================
cat("\n==== C1: CDH1 vs lobular histology ====\n")
ilc_raw <- getclin("CANCER_IDENTIFICATION_ILC", samp)
ILC <- ifelse(toupper(ilc_raw)=="YES","ILC", ifelse(toupper(ilc_raw)=="NO","non-ILC", NA))
if(sum(!is.na(ILC)) < 10){                       # fallback to histology text
  h <- toupper(paste(getclin("DX_HISTOLOGY",samp), getclin("BX_HISTOLOGY",samp)))
  ILC <- ifelse(grepl("LOBULAR",h),"ILC", ifelse(grepl("DUCTAL",h),"non-ILC",NA)) }
cdh1_cna <- if("CDH1" %in% rownames(cna_matrix_final)) as.numeric(cna_matrix_final["CDH1",samp]) else rep(NA,length(samp))
cdh1_rna <- if("CDH1" %in% rownames(rna_matrix_final)) as.numeric(rna_matrix_final["CDH1",samp]) else rep(NA,length(samp))
d1 <- data.frame(ILC=ILC, cna=cdh1_cna, rna=cdh1_rna)[!is.na(ILC),]
cat("  Histology groups:\n"); print(table(d1$ILC))
if(nlevels(factor(d1$ILC))==2 && nrow(d1)>=10){
  loss <- d1$cna < 0
  cat("  CDH1 loss x histology (Fisher):\n")
  print(tryCatch(fisher.test(table(d1$ILC, loss)), error=function(e) "n/a"))
  cat("  CDH1 expression by histology (Wilcoxon):\n")
  print(tryCatch(wilcox.test(rna ~ ILC, data=d1), error=function(e) "n/a"))
  pC1 <- ggplot(d1[!is.na(d1$rna),], aes(ILC, rna, fill=ILC)) +
    geom_boxplot(outlier.shape=NA, alpha=.7) + geom_jitter(width=.15, size=1.5, alpha=.6) +
    scale_fill_manual(values=c("ILC"="#E64B35","non-ILC"="#4DBBD5"), guide="none") +
    labs(title="CDH1 expression by histology", subtitle="E-cadherin loss defines lobular carcinoma",
         x=NULL, y="CDH1 (log2 expression)") + theme_publication(12)
  save_pdf(pC1, "14_cdh1_lobular.pdf", 6, 6)
} else cat("  Too few labeled samples for the ILC test.\n")

## ===========================================================================
## C2. Dosage drivers vs DISEASE AGGRESSIVENESS (time from dx to metastasis)
##     Shorter time-to-metastasis = more aggressive disease.
## ===========================================================================
cat("\n==== C2: dosage drivers vs time-to-metastasis ====\n")
ttm <- suppressWarnings(as.numeric(getclin("CALC_TIME_TO_METS_DX_MONTHS", samp)))
cat(sprintf("  time-to-mets available: %d samples (median %.0f months)\n",
            sum(is.finite(ttm)), median(ttm, na.rm=TRUE)))
cs <- dsc$Coupling[match(samp, dsc$SAMPLE_ID)]
cat("  (negative rho = higher value -> shorter time -> more aggressive)\n")
for(g in c("ERBB2","CCND1","CDH1"))
  if(g %in% rownames(rna_matrix_final)) rho_line(paste0(g," expr vs time-to-mets"),
                                                 as.numeric(rna_matrix_final[g,samp]), ttm)
rho_line("coupling score vs TTM", cs, ttm)

## ===========================================================================
## C3. CCND1 vs PROLIFERATION (tumor grade, Ki67)
##     CCND1 (cyclin D1) drives cell-cycle progression -> expect grade/Ki67 link.
## ===========================================================================
cat("\n==== C3: dosage drivers vs grade / Ki67 ====\n")
grade <- grade_num(getclin("BX_GRADE", samp)); if(sum(!is.na(grade))<10) grade <- grade_num(getclin("DX_GRADE", samp))
ki67  <- suppressWarnings(as.numeric(gsub("[^0-9.]","", getclin("BX_KI67", samp))))
cat(sprintf("  grade available: %d | Ki67 available: %d\n", sum(is.finite(grade)), sum(is.finite(ki67))))
for(g in c("CCND1","ERBB2","CDH1")) if(g %in% rownames(rna_matrix_final)){
  ex <- as.numeric(rna_matrix_final[g,samp])
  rho_line(paste0(g," vs grade"), ex, grade)
  rho_line(paste0(g," vs Ki67"),  ex, ki67)
}

## ===========================================================================
## C4. (Exploratory) dosage drivers vs metastatic burden (# of involved sites)
## ===========================================================================
cat("\n==== C4 (exploratory): dosage vs metastatic burden ====\n")
site_cols_abst <- grep("^METS_SITES_AT_METS_ABST_", colnames(clinical_matched), value=TRUE)
if(length(site_cols_abst) >= 3){
  burden <- rowSums(sapply(site_cols_abst, function(c)
    toupper(clean_clin(clinical_matched[[c]][match(samp, clinical_matched$SAMPLE_ID)]))=="YES"), na.rm=TRUE)
  cat(sprintf("  metastatic-site burden available: %d (median %d sites)\n",
              sum(burden>0), median(burden[burden>0])))
  for(g in c("CCND1","ERBB2","CDH1")) if(g %in% rownames(rna_matrix_final))
    rho_line(paste0(g," vs #mets sites"), as.numeric(rna_matrix_final[g,samp]), as.numeric(burden))
  rho_line("coupling score vs #mets sites", cs, as.numeric(burden))
} else cat("  metastatic-site columns not found.\n")

cat("\n==== DONE. Report: C1 histology test + figure 14, C2/C3 rho tables. ====\n")
## NOTE on interpretation / honesty:
##  - These are associations in a small, biopsy-level cohort; treat as exploratory.
##  - Correct for multiple testing if you report several (e.g., p.adjust over the set).
##  - Do NOT over-claim causation; frame as "consistent with known biology".
################################################################################
