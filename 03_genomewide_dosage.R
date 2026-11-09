################################################################################
# GENOME-WIDE DOSAGE SENSITIVITY (removes the driver-panel selection bias)
# Computes CNA->RNA cis-coupling for ALL genes, locates the drivers in that
# distribution, and tests cross-cohort reproducibility genome-wide.
# Prereq: run mbc_ALL_IN_ONE.R first. Uses the FULL matrices:
#   cna_matrix (all genes), rna_matrix (all genes), clinical_matched, final_samples
#   For TCGA: tcga_rna, tcga_cna, have_tcga
################################################################################
library(ggplot2)
if (!dir.exists("output")) dir.create("output")
clean_num <- function(x) sub("\\..*$","",x)
one_per_patient <- function(pid){ set.seed(42)
  sort(unlist(tapply(seq_along(pid),pid,function(ii) ii[sample(length(ii),1)]),use.names=FALSE)) }

## ---- fast per-gene Spearman between two aligned genes x samples matrices ----
gw_coupling <- function(CNA, RNA){
  g  <- intersect(rownames(CNA), rownames(RNA))
  CNA<- CNA[g,,drop=FALSE]; RNA<- RNA[g,,drop=FALSE]
  s  <- intersect(colnames(CNA), colnames(RNA)); CNA<-CNA[,s,drop=FALSE]; RNA<-RNA[,s,drop=FALSE]
  # drop genes constant in either layer
  keep <- apply(CNA,1,function(x) sd(x,na.rm=TRUE)>0) & apply(RNA,1,function(x) sd(x,na.rm=TRUE)>0)
  CNA<-CNA[keep,,drop=FALSE]; RNA<-RNA[keep,,drop=FALSE]
  CNr <- t(apply(CNA,1,rank)); RNr <- t(apply(RNA,1,rank))     # ranks -> Spearman
  cz <- CNr-rowMeans(CNr); rz <- RNr-rowMeans(RNr)
  rho <- rowSums(cz*rz)/sqrt(rowSums(cz^2)*rowSums(rz^2))
  n <- ncol(CNr); tt <- rho*sqrt((n-2)/(1-rho^2)); p <- 2*pt(-abs(tt),df=n-2)
  data.frame(Gene=rownames(CNA), rho=as.numeric(rho), p=p,
             padj=p.adjust(p,"BH"), n=n, row.names=NULL) }

## ---- build full, patient-level MBCproject matrices ------------------------
# full CNA (numeric), full RNA (log2), clean gene names, dedupe
CN <- cna_matrix; mode(CN)<-"numeric"; rownames(CN)<-clean_num(rownames(CN)); CN<-CN[!duplicated(rownames(CN)),,drop=FALSE]
RN <- rna_matrix;                      rownames(RN)<-clean_num(rownames(RN)); RN<-RN[!duplicated(rownames(RN)),,drop=FALSE]
S  <- intersect(colnames(CN), colnames(RN)); S <- intersect(S, final_samples)
pidS <- clinical_matched$PATIENT_ID[match(S, clinical_matched$SAMPLE_ID)]
Spt <- S[one_per_patient(pidS)]                                   # patient-level
cat("Genome-wide set: patient-level n =", length(Spt),
    "(from", length(S), "biopsies)\n")

gw <- gw_coupling(CN[,Spt,drop=FALSE], RN[,Spt,drop=FALSE])
gw <- gw[order(-gw$rho),]
cat(sprintf("Genes tested genome-wide: %d\n", nrow(gw)))
cat(sprintf("Dosage-sensitive (FDR<0.05 & rho>0): %d (%.1f%%)\n",
            sum(gw$padj<0.05 & gw$rho>0), 100*mean(gw$padj<0.05 & gw$rho>0)))
cat(sprintf("Median genome-wide rho = %.2f\n", median(gw$rho)))
write.csv(gw, "output/genomewide_dosage_MBC.csv", row.names=FALSE)

## ---- where do the drivers fall in the genome-wide distribution? -----------
drivers <- c("ERBB2","CCND1","CDH1","PIK3CA","TP53","GATA3","FOXA1")
pct <- function(g){ if(!(g %in% gw$Gene)) return(NA); mean(gw$rho <= gw$rho[gw$Gene==g])*100 }
cat("\nDriver position in genome-wide dosage distribution:\n")
for(g in drivers) if(g %in% gw$Gene)
  cat(sprintf("  %-7s rho=%+.2f | percentile %4.1f | FDR %s\n", g,
      gw$rho[gw$Gene==g], pct(g), ifelse(gw$padj[gw$Gene==g]<0.05,"<0.05","ns")))

## ---- Figure: genome-wide rho distribution with drivers marked -------------
lab <- gw[gw$Gene %in% c("ERBB2","CCND1","CDH1","GATA3","FOXA1"),]
pG <- ggplot(gw, aes(rho)) +
  geom_histogram(bins=60, fill="grey80", color="grey60", linewidth=.1) +
  geom_vline(data=lab, aes(xintercept=rho, color=Gene), linewidth=.9) +
  ggrepel::geom_text_repel(data=lab, aes(x=rho, y=0, label=Gene, color=Gene),
                           angle=90, hjust=0, size=3.5, fontface="bold",
                           nudge_y=max(table(cut(gw$rho,60)))*0.3, show.legend=FALSE) +
  scale_color_manual(values=c(ERBB2="#E64B35",CCND1="#D73027",CDH1="#B2182B",
                              GATA3="#4DBBD5",FOXA1="#3C5488"), guide="none") +
  labs(title="Genome-wide copy-number dosage sensitivity",
       subtitle=sprintf("CNA-RNA cis-coupling across %d genes (patient-level, n=%d)", nrow(gw), length(Spt)),
       x="Spearman rho (CNA vs own expression)", y="Number of genes") +
  theme_publication(12)
save_pdf(pG, "15_genomewide_dosage.pdf", 9, 6)

## ---- TCGA genome-wide + cross-cohort reproducibility (unbiased) -----------
if(exists("have_tcga") && isTRUE(have_tcga) && exists("tcga_cna") && exists("tcga_rna")){
  cat("\n==== Genome-wide reproducibility in TCGA-BRCA ====\n")
  gc2 <- if("Hugo_Symbol" %in% colnames(tcga_cna)) "Hugo_Symbol" else colnames(tcga_cna)[1]
  m2  <- intersect(c("Hugo_Symbol","Entrez_Gene_Id"), colnames(tcga_cna))
  Ct  <- as.matrix(tcga_cna[,setdiff(colnames(tcga_cna),m2),drop=FALSE]); mode(Ct)<-"numeric"
  rownames(Ct)<-clean_num(tcga_cna[[gc2]]); Ct<-Ct[!duplicated(rownames(Ct)),,drop=FALSE]
  Rt  <- as.matrix(tcga_rna[,-c(1:2)]); rownames(Rt)<-clean_num(make.names(tcga_rna$Hugo_Symbol,unique=TRUE))
  mode(Rt)<-"numeric"; if(max(Rt,na.rm=TRUE)>100) Rt<-log2(Rt+1); Rt<-Rt[!duplicated(rownames(Rt)),,drop=FALSE]
  St  <- intersect(colnames(Ct), colnames(Rt))
  gwT <- gw_coupling(Ct[,St,drop=FALSE], Rt[,St,drop=FALSE])
  cat(sprintf("TCGA genes tested: %d | n samples: %d\n", nrow(gwT), length(St)))

  mrg <- merge(gw[,c("Gene","rho")], gwT[,c("Gene","rho")], by="Gene", suffixes=c("_mbc","_tcga"))
  ct <- suppressWarnings(cor.test(mrg$rho_mbc, mrg$rho_tcga, method="spearman"))
  cat(sprintf("-> GENOME-WIDE cross-cohort dosage reproducibility: Spearman rho = %.2f | p < %.0e | genes = %d\n",
              ct$estimate, max(ct$p.value,1e-300), nrow(mrg)))
  # unbiased reproducibility figure (2D density of per-gene rho, both cohorts)
  pR <- ggplot(mrg, aes(rho_mbc, rho_tcga)) +
    geom_bin2d(bins=70) + scale_fill_gradient(low="grey85", high="#B2182B", name="genes") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey40") +
    labs(title="Genome-wide reproducibility of dosage sensitivity",
         subtitle=sprintf("Per-gene CNA-RNA coupling, %d genes | Spearman rho=%.2f", nrow(mrg), ct$estimate),
         x="CNA-RNA coupling (MBCproject)", y="CNA-RNA coupling (TCGA-BRCA)") +
    theme_publication(12)
  save_pdf(pR, "16_genomewide_reproducibility.pdf", 7, 6.5)
} else cat("\n(TCGA objects not in memory; run section 5 of mbc_ALL_IN_ONE.R for the genome-wide reproducibility.)\n")

cat("\n==== DONE. Report: % dosage-sensitive genome-wide, driver percentiles,\n",
    "     and the genome-wide cross-cohort rho (the unbiased reproducibility). ====\n")
################################################################################
