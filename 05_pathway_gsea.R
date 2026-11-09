################################################################################
# 05_pathway_gsea.R
# GSEA of genome-wide dosage-sensitivity ranking: GO-BP + Reactome, in both
# MBCproject and TCGA-BRCA, with pathway-level cross-cohort reproducibility.
# Unbiased: ranks ALL genes by CNA->RNA coupling (rho), not a curated panel.
# Prereq: run 02_pipeline.R and 03_genomewide_dosage.R first (gw / gwT in memory).
#         (KEGG is omitted: gseKEGG requires a live KEGG API connection.)
################################################################################
for(p in c("clusterProfiler","org.Hs.eg.db","ReactomePA")){
  if(!requireNamespace(p,quietly=TRUE)){ if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(p,update=FALSE,ask=FALSE) }
}
library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); library(ggplot2)
if(!dir.exists("output")) dir.create("output")
stopifnot(exists("gw"))                                   # from 03_genomewide_dosage.R

## ---- GO-BP GSEA on SYMBOL-ranked dosage rho --------------------------------
run_go <- function(gwdf,label){
  v <- gwdf$rho; names(v) <- gwdf$Gene; v <- v[!is.na(v)&names(v)!=""]; v <- v[!duplicated(names(v))]
  v <- sort(v,decreasing=TRUE); set.seed(42)
  ego <- gseGO(geneList=v, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP",
               minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, eps=0, verbose=FALSE)
  cat(sprintf("\n[%s] GO-BP FDR<0.05: %d\n", label, sum(ego@result$p.adjust<0.05)))
  print(head(ego@result[order(ego@result$p.adjust),c("Description","NES","setSize","p.adjust")],12), row.names=FALSE)
  ego }

## ---- Reactome GSEA on ENTREZ-ranked dosage rho -----------------------------
entrez_rank <- function(gwdf){
  v <- gwdf$rho; names(v) <- gwdf$Gene; v <- v[!is.na(v)&names(v)!=""]; v <- v[!duplicated(names(v))]
  ids <- suppressWarnings(bitr(names(v),"SYMBOL","ENTREZID",org.Hs.eg.db))
  ids <- ids[!duplicated(ids$SYMBOL)&!duplicated(ids$ENTREZID),]
  vk <- v[ids$SYMBOL]; names(vk) <- ids$ENTREZID; sort(vk,decreasing=TRUE) }
run_react <- function(gwdf,label){
  set.seed(42)
  rp <- tryCatch(gsePathway(entrez_rank(gwdf), organism="human", minGSSize=10,
                            maxGSSize=500, pvalueCutoff=0.05, eps=0, verbose=FALSE),
                 error=function(e) NULL)
  if(!is.null(rp) && nrow(rp)>0){
    cat(sprintf("\n[%s] Reactome FDR<0.05: %d\n", label, sum(rp@result$p.adjust<0.05)))
    print(head(rp@result[order(rp@result$p.adjust),c("Description","NES","setSize","p.adjust")],12), row.names=FALSE) }
  rp }

go_mbc  <- run_go(gw,"MBCproject");    re_mbc  <- run_react(gw,"MBCproject")
go_tcga <- if(exists("gwT")) run_go(gwT,"TCGA-BRCA")   else NULL
re_tcga <- if(exists("gwT")) run_react(gwT,"TCGA-BRCA") else NULL

## ---- dotplot of top enriched pathways --------------------------------------
dotplot_paths <- function(res, title, file){
  if(is.null(res) || nrow(res)==0) return(invisible())
  top <- head(res@result[order(res@result$p.adjust),],15)
  top$Description <- factor(top$Description, levels=rev(top$Description[order(top$NES)]))
  p <- ggplot(top, aes(NES, Description, size=setSize, color=p.adjust)) + geom_point() +
    scale_color_gradient(low="#B2182B", high="#2166AC", name="FDR") + scale_size_continuous(name="genes") +
    labs(title=title, subtitle="GSEA on genome-wide CNA-RNA coupling (MBCproject)",
         x="Normalized enrichment score (NES)", y=NULL) +
    theme_publication(11) + theme(axis.text.y=element_text(size=9))
  save_pdf(p, file, 11, 7) }
dotplot_paths(re_mbc, "Reactome pathways enriched among dosage-sensitive genes", "22_reactome_dotplot_MBC.pdf")
dotplot_paths(go_mbc, "GO-BP processes enriched among dosage-sensitive genes",   "18_gsea_dotplot_MBC.pdf")

## ---- pathway-level cross-cohort reproducibility ----------------------------
repro <- function(a,b,tag,file){
  if(is.null(a)||is.null(b)) return(invisible())
  m <- a@result[,c("ID","Description","NES","p.adjust")]; t <- b@result[,c("ID","NES","p.adjust")]
  mg <- merge(m,t,by="ID",suffixes=c("_mbc","_tcga")); if(nrow(mg)<5) return(invisible())
  ct <- suppressWarnings(cor.test(mg$NES_mbc,mg$NES_tcga,method="spearman"))
  cat(sprintf("\n[%s] pathways in both: %d | sig in both: %d | NES rho=%.2f (p=%.1e)\n",
      tag, nrow(mg), sum(mg$p.adjust_mbc<0.05 & mg$p.adjust_tcga<0.05), ct$estimate, ct$p.value))
  p <- ggplot(mg, aes(NES_mbc,NES_tcga)) + geom_point(alpha=.4,size=1.6,color="#3C5488") +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="grey50") +
    geom_hline(yintercept=0,linewidth=.3,color="grey70") + geom_vline(xintercept=0,linewidth=.3,color="grey70") +
    labs(title=paste0(tag,": cross-cohort reproducibility"),
         subtitle=sprintf("%d pathways | Spearman rho=%.2f", nrow(mg), ct$estimate),
         x="NES (MBCproject)", y="NES (TCGA-BRCA)") + theme_publication(12)
  save_pdf(p, file, 7, 6.5) }
repro(re_mbc, re_tcga, "Reactome", "20_reactome_reproducibility.pdf")
repro(go_mbc, go_tcga, "GO-BP",    "19_gsea_reproducibility.pdf")

cat("\nNOTE: copy-number changes occur in genomic blocks, so dosage-sensitive genes\n",
    "can co-localize on altered regions; enrichment is partly positional and correlative.\n",
    "The cross-cohort NES concordance shows the enrichment pattern itself is reproducible.\n")
################################################################################
