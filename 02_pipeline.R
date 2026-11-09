################################################################################
# 02_pipeline.R  — MBCproject multi-omics CONSOLIDATED pipeline
# Run top to bottom. Primary axis = molecular subtype. Site = reference only.
#
# SECTIONS
#   0  Libraries, helpers, config
#   1  Load + preprocess MBCproject (brca_mbcproject_2022)
#   2  Molecular subtype (primary) ; figures 1-6 ; site reference + PERMANOVA
#   3  Validation: BER + patient-disjoint confusion matrix
#   4  Cross-omics convergence: integration value, cis-coupling, coupling score
#   5  External validation (TCGA-BRCA): signature + coupling reproducibility
#   6  Reviewer re-analysis (B1-B5)
#
# OBJECT NAMES (single convention, no X_sub/Y_sub anymore):
#   X          list(Mutations,RNA,CNA)  samples x features (cleaned)
#   Y_final    factor, molecular subtype, named by SAMPLE_ID
#   final_diablo, design, kx
#   mutation_matrix_final / rna_matrix_final / cna_matrix_final  (genes x samples)
#   clinical_matched, cis, cmp
#
# TCGA: set TCGA_DIR below. If the folder is absent, sections 5-6(B4) self-skip.
################################################################################

## ============================ 0. SETUP =====================================
if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("mixOmics","permute","vegan","maftools","ComplexHeatmap","tidyverse",
          "ggplot2","ggrepel","caret","circlize","patchwork","RColorBrewer")
bioc <- c("mixOmics","maftools","ComplexHeatmap")
for (p in pkgs){ if(!require(p,character.only=TRUE,quietly=TRUE)){
  if(p %in% bioc) BiocManager::install(p,update=FALSE,ask=FALSE) else install.packages(p,dependencies=TRUE)
  library(p,character.only=TRUE)}}
set.seed(42)

TCGA_DIR <- "tcga_brca/"            # <-- set to your TCGA-BRCA export folder (or leave; 5-6 self-skip)
ncomp <- 2; pred_dist <- "centroids.dist"
if (!dir.exists("output")) dir.create("output")

sub_cols <- c("HER2_Positive"="#E64B35","HR+/HER2-"="#4DBBD5","Triple_Negative"="#00A087")
theme_publication <- function(bs=13){ theme_classic(bs) + theme(
  plot.title=element_text(size=bs+4,face="bold",hjust=.5,margin=margin(b=6)),
  plot.subtitle=element_text(size=bs,hjust=.5,color="grey40",margin=margin(b=10)),
  axis.title=element_text(size=bs+1,face="bold"), axis.text=element_text(size=bs-1,color="black"),
  axis.line=element_line(color="black",linewidth=.6),
  legend.title=element_text(size=bs,face="bold"), legend.text=element_text(size=bs-2),
  panel.grid.major.y=element_line(color="grey92",linewidth=.4),
  panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
  strip.text=element_text(size=bs,face="bold",color="white"),
  strip.background=element_rect(fill="#3C5488",color=NA), plot.margin=margin(15,20,15,15)) }
save_pdf <- function(pl,fn,w,h){ tryCatch(ggsave(file.path("output",fn),pl,device=cairo_pdf,width=w,height=h,units="in"),
  error=function(e) ggsave(file.path("output",fn),pl,device="pdf",width=w,height=h,units="in")); cat("  Saved:",fn,"\n") }

clean_num <- function(x) sub("\\..*$","",x)
one_per_patient <- function(pid){ set.seed(42)
  sort(unlist(tapply(seq_along(pid),pid,function(ii) ii[sample(length(ii),1)]),use.names=FALSE)) }
# CV-safe block cleaner: binary cols need >= min_ev events in BOTH directions
clean_block <- function(M,min_ev=4){
  M <- M[,apply(M,2,var)>0,drop=FALSE]
  isb <- apply(M,2,function(x) length(unique(x))<=3)
  if(any(isb)){ ones<-colSums(M!=0); zeros<-nrow(M)-ones
    M <- M[,!(isb & (ones<min_ev | zeros<min_ev)),drop=FALSE] }
  nzv <- caret::nearZeroVar(M,freqCut=95/5,uniqueCut=15)
  if(length(nzv)>0) M <- M[,-nzv,drop=FALSE]
  M }
cap <- function(k,M) pmin(k,ncol(M))

## ==================== 1. LOAD + PREPROCESS MBCproject =======================
cat("\n=== 1. LOAD + PREPROCESS ===\n")
req <- c("data_clinical_patient.txt","data_clinical_sample.txt",
         "data_mutations.txt","data_cna.txt","data_mrna_seq_v2_rsem.txt")
for(f in req) if(!file.exists(f)) stop(paste("File not found:",f))
na_s <- c("NA","","[Not Available]","[Not Applicable]","[Unknown]")
clinical_patient <- read.delim("data_clinical_patient.txt",skip=4,stringsAsFactors=FALSE,na.strings=na_s)
clinical_sample  <- read.delim("data_clinical_sample.txt", skip=4,stringsAsFactors=FALSE,na.strings=na_s)
mutations <- read.delim("data_mutations.txt",stringsAsFactors=FALSE,comment.char="#",quote="",fill=TRUE)
cna_data  <- read.delim("data_cna.txt",stringsAsFactors=FALSE,check.names=FALSE)
rna_data  <- read.delim("data_mrna_seq_v2_rsem.txt",stringsAsFactors=FALSE,check.names=FALSE)
if(any(grepl("^TCGA",colnames(rna_data)))) stop("RNA file has TCGA ids -- wrong cohort.")

clinical_full <- clinical_sample %>% left_join(clinical_patient,by="PATIENT_ID",suffix=c("_sample","_patient"))
important_genes <- c("TP53","PIK3CA","ESR1","ERBB2","GATA3","CDH1","MAP3K1","PTEN","AKT1","NF1",
  "RB1","BRCA1","BRCA2","ATM","CHEK2","PALB2","ARID1A","KMT2C","NCOR1","TBX3","MLL3","MAP2K4",
  "CBFB","RUNX1","SF3B1","FOXA1","ZFH4","CASP8","CDKN1B","CCND1")

mutation_matrix <- mutations %>% filter(Hugo_Symbol %in% important_genes) %>%
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode) %>% distinct() %>% mutate(mutated=1) %>%
  pivot_wider(names_from=Hugo_Symbol,values_from=mutated,values_fill=0)
ids <- mutation_matrix$Tumor_Sample_Barcode
mutation_matrix <- as.matrix(mutation_matrix[,-1]); rownames(mutation_matrix) <- ids

gcol <- if("Hugo_Symbol" %in% colnames(cna_data)) "Hugo_Symbol" else colnames(cna_data)[1]
meta <- intersect(c("Hugo_Symbol","Entrez_Gene_Id"),colnames(cna_data))
cna_matrix <- as.matrix(cna_data[,setdiff(colnames(cna_data),meta),drop=FALSE])
rownames(cna_matrix) <- cna_data[[gcol]]
cna_matrix_filtered <- cna_matrix[rownames(cna_matrix) %in% important_genes,,drop=FALSE]
mode(cna_matrix_filtered) <- "numeric"

rna_matrix <- as.matrix(rna_data[,-c(1:2)]); rownames(rna_matrix) <- make.names(rna_data$Hugo_Symbol,unique=TRUE)
mode(rna_matrix) <- "numeric"
if(max(rna_matrix,na.rm=TRUE)>100) rna_matrix <- log2(rna_matrix+1)
rna_matrix <- rna_matrix[rowSums(is.na(rna_matrix))/ncol(rna_matrix)<0.2,]
gv <- apply(rna_matrix,1,var,na.rm=TRUE); rna_matrix <- rna_matrix[!is.na(gv),]; gv <- gv[!is.na(gv)]
rna_matrix_filtered <- rna_matrix[order(gv,decreasing=TRUE)[seq_len(min(5000,length(gv)))],]

common_all <- Reduce(intersect,list(rownames(mutation_matrix),colnames(rna_matrix_filtered),colnames(cna_matrix_filtered)))
use_cna <- length(common_all)>=30
final_samples <- if(use_cna) common_all else intersect(rownames(mutation_matrix),colnames(rna_matrix_filtered))
if(length(final_samples)<30) stop("Too few overlapping samples; check identifier consistency.")
mut_m <- mutation_matrix[final_samples,,drop=FALSE]
rna_m <- rna_matrix_filtered[,final_samples,drop=FALSE]
cna_m <- cna_matrix_filtered[,final_samples,drop=FALSE]
clinical_matched <- clinical_full %>% filter(SAMPLE_ID %in% final_samples) %>% arrange(match(SAMPLE_ID,final_samples))
cat("  Integrable samples:",length(final_samples),"\n")

## ==================== 2. SUBTYPE (PRIMARY) + FIGURES ========================
cat("\n=== 2. MOLECULAR SUBTYPE (primary) ===\n")
clinical_matched <- clinical_matched %>% mutate(
  ER=case_when(toupper(BX_ER)%in%c("POSITIVE","POS","+","1")~"Positive",
               toupper(BX_ER)%in%c("NEGATIVE","NEG","-","0")~"Negative",TRUE~NA_character_),
  PR=case_when(toupper(BX_PR)%in%c("POSITIVE","POS","+","1")~"Positive",
               toupper(BX_PR)%in%c("NEGATIVE","NEG","-","0")~"Negative",TRUE~NA_character_),
  HER2=case_when(toupper(BX_HER2OVERALL)%in%c("POSITIVE","POS","+","1","AMPLIFIED")~"Positive",
                 toupper(BX_HER2OVERALL)%in%c("NEGATIVE","NEG","-","0","NOT AMPLIFIED")~"Negative",TRUE~NA_character_),
  Molecular_Subtype=case_when(HER2=="Positive"~"HER2_Positive",
    (ER=="Positive"|PR=="Positive")&HER2=="Negative"~"HR+/HER2-",
    ER=="Negative"&PR=="Negative"&HER2=="Negative"~"Triple_Negative",TRUE~NA_character_))

Ysub <- clinical_matched$Molecular_Subtype; names(Ysub) <- clinical_matched$SAMPLE_ID
Ysub <- Ysub[!is.na(Ysub)]
ok <- names(table(Ysub))[table(Ysub)>=3]; Ysub <- Ysub[Ysub %in% ok]
subtype_samples <- names(Ysub)
Y_final <- droplevels(factor(Ysub,levels=c("HER2_Positive","HR+/HER2-","Triple_Negative")))
names(Y_final) <- subtype_samples
site_all <- clinical_matched$BX_LOCATION[match(subtype_samples,clinical_matched$SAMPLE_ID)]; names(site_all)<-subtype_samples
mutation_matrix_final <- mut_m[subtype_samples,,drop=FALSE]
rna_matrix_final      <- rna_m[,subtype_samples,drop=FALSE]
cna_matrix_final      <- cna_m[,subtype_samples,drop=FALSE]
cat("  Subtype set n =",length(Y_final),"\n"); print(table(Y_final))

# Figure 1: mutation frequency
mf <- sort(colSums(mutation_matrix_final)/nrow(mutation_matrix_final)*100,decreasing=TRUE)
mf <- mf[!is.na(mf)&!is.na(names(mf))&mf>0]; tn <- min(20,length(mf)); gs <- names(mf)[seq_len(tn)]
mdf <- data.frame(Gene=factor(gs,levels=rev(gs)),Frequency=as.numeric(mf[seq_len(tn)]),row.names=NULL)
save_pdf(ggplot(mdf,aes(Gene,Frequency,fill=Frequency))+
  geom_col(color="black",linewidth=.25,width=.75)+
  geom_text(aes(label=sprintf("%.1f%%",Frequency)),hjust=-.15,size=3.3,fontface="bold")+
  scale_fill_gradientn(colors=c("#4575B4","#FEE090","#D73027"),name="Freq. (%)")+
  scale_y_continuous(expand=expansion(mult=c(0,.18)),labels=function(x)paste0(x,"%"))+coord_flip()+
  labs(title="Somatic Mutation Frequency (subtype set)",subtitle=paste0("Top ",tn," genes | n=",nrow(mutation_matrix_final)),x=NULL,y="Mutation Frequency (%)")+
  theme_publication()+theme(axis.text.y=element_text(face="italic",size=11)),"01_mutation_frequency.pdf",10,8)

# Figure 2: subtype distribution
tb <- sort(table(Y_final),decreasing=TRUE); sdf <- data.frame(Subtype=factor(names(tb),levels=names(tb)),n=as.integer(tb)); sdf$Pct<-sdf$n/sum(sdf$n)*100
save_pdf(ggplot(sdf,aes(Subtype,Pct,fill=Subtype))+geom_col(color="black",linewidth=.3,width=.65)+
  geom_text(aes(label=paste0("n=",n,"\n(",round(Pct,1),"%)")),vjust=-.4,size=4,fontface="bold")+
  scale_fill_manual(values=sub_cols,guide="none")+scale_y_continuous(expand=expansion(mult=c(0,.18)),labels=function(x)paste0(x,"%"))+
  labs(title="Molecular Subtype Distribution",subtitle=paste0("Metastatic Breast Cancer | n=",sum(sdf$n)),x=NULL,y="Percentage (%)")+theme_publication(),"02_subtype_distribution.pdf",8,6)

# Figure 3: RNA PCA by subtype
rna_pca_in <- rna_matrix_final[apply(rna_matrix_final,1,var)>0,,drop=FALSE]
pca <- prcomp(t(rna_pca_in),scale.=TRUE); pv <- summary(pca)$importance[2,1:2]*100
pdf3 <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2],Group=Y_final)
ct3 <- pdf3 %>% group_by(Group) %>% summarise(PC1=mean(PC1),PC2=mean(PC2),.groups="drop")
save_pdf(ggplot(pdf3,aes(PC1,PC2,fill=Group,color=Group))+
  stat_ellipse(geom="polygon",alpha=.12,linewidth=1,level=.95)+
  geom_point(size=2.8,shape=21,color="black",stroke=.4,alpha=.85)+
  geom_point(data=ct3,size=6,shape=23,color="black",stroke=1.8)+
  geom_text_repel(data=ct3,aes(label=Group),color="black",fontface="bold",size=4,show.legend=FALSE)+
  scale_fill_manual(values=sub_cols,name="Molecular Subtype")+scale_color_manual(values=sub_cols,name="Molecular Subtype")+
  labs(title="RNA PCA by Molecular Subtype",subtitle=sprintf("%d genes | %d samples | 95%% ellipses",nrow(rna_matrix_final),ncol(rna_matrix_final)),
       x=sprintf("PC1 (%.1f%%)",pv[1]),y=sprintf("PC2 (%.1f%%)",pv[2]))+
  theme_publication()+theme(panel.grid.major=element_line(color="grey92",linewidth=.4)),"03_rna_pca_subtype.pdf",10,8)
cat(sprintf("  RNA PCA variance: PC1 %.1f%% | PC2 %.1f%%\n",pv[1],pv[2]))

# DIABLO (subtype)
X <- list(Mutations=as.matrix(mutation_matrix_final),RNA=t(rna_matrix_final),CNA=t(cna_matrix_final))
X <- lapply(X,clean_block); cat("  Clean block sizes:\n"); print(sapply(X,ncol))
design <- matrix(0.1,length(X),length(X),dimnames=list(names(X),names(X))); diag(design)<-0
kx <- list(Mutations=cap(c(8,8),X$Mutations),RNA=cap(c(25,50),X$RNA),CNA=cap(c(8,8),X$CNA))
final_diablo <- block.splsda(X,Y_final,ncomp=ncomp,keepX=kx,design=design)
cat("  Subtype DIABLO fitted.\n")

# Figure 4: DIABLO sample plot
sdfp <- do.call(rbind,lapply(seq_along(X),function(i) data.frame(PC1=final_diablo$variates[[i]][,1],PC2=final_diablo$variates[[i]][,2],Group=Y_final,Block=names(X)[i])))
save_pdf(ggplot(sdfp,aes(PC1,PC2,fill=Group,color=Group))+stat_ellipse(geom="polygon",alpha=.12,linewidth=1,level=.95)+
  geom_point(size=2.5,shape=21,color="black",stroke=.4,alpha=.85)+facet_wrap(~Block,nrow=1)+
  scale_fill_manual(values=sub_cols,name="Molecular Subtype")+scale_color_manual(values=sub_cols,name="Molecular Subtype")+
  labs(title="DIABLO Integration by Molecular Subtype",subtitle=paste0(length(Y_final)," samples | ",nlevels(Y_final)," subtypes"),x="Component 1",y="Component 2")+theme_publication(12),"04_diablo_subtype.pdf",14,5)

# Figure 5: loadings
lp <- lapply(names(X),function(bn){ ld<-final_diablo$loadings[[bn]][,1]; t2<-min(15,length(ld)); ld<-sort(ld); ld<-c(head(ld,ceiling(t2/2)),tail(ld,floor(t2/2)))
  d<-data.frame(Feature=factor(names(ld),levels=names(ld)),Loading=as.numeric(ld),Direction=ifelse(ld>0,"Positive","Negative"))
  ggplot(d,aes(Feature,Loading,fill=Direction))+geom_col(color="black",linewidth=.2,width=.75)+geom_hline(yintercept=0,linewidth=.4)+
    scale_fill_manual(values=c(Positive="#E64B35",Negative="#4DBBD5"),name=NULL)+coord_flip()+labs(title=paste("Comp1 Loadings --",bn),x=NULL,y="Loading")+
    theme_publication(11)+theme(legend.position="top",axis.text.y=element_text(size=9,face=ifelse(bn=="RNA","italic","plain"))) })
save_pdf(wrap_plots(lp,ncol=1),"05_diablo_loadings.pdf",9,5*length(lp)+1)

# Figure 6: circos + clean inter-block correlations
draw_circos <- function(cut){ sz<-max(14,min(22,sum(sapply(kx,function(x)x[1]))*0.18))
  cairo_pdf("output/06_circos_correlation.pdf",width=sz,height=sz); par(mar=c(1,1,3,1),bg="white",cex=0.9)
  circosPlot(final_diablo,cutoff=cut,comp=c(1,2),color.blocks=c("#E64B35","#4DBBD5","#00A087")[seq_along(X)],
    color.cor=c("#B2182B","#2166AC"),size.variables=0.7,size.labels=1.1,showIntraLinks=FALSE)
  mtext("Cross-Omics Feature Correlation -- DIABLO (subtype)",3,0.5,cex=1.3,font=2)
  mtext(sprintf("Comp 1-2 | |r| >= %.1f | Red: positive, Blue: negative",cut),3,-0.5,cex=0.9,col="grey40"); dev.off() }
tryCatch({draw_circos(0.5); cat("  Saved: 06_circos_correlation.pdf\n")},
  error=function(e) tryCatch({draw_circos(0.3); cat("  Saved: 06_circos_correlation.pdf (0.3)\n")},
  error=function(e2){try(dev.off(),silent=TRUE); cat("  circos skipped\n")}))
vt <- final_diablo$variates; cmb <- function(a,b) round(cor(vt[[a]][,1],vt[[b]][,1]),2)
cat(sprintf("  Inter-block r (clean): RNA-CNA %.2f | RNA-Mut %.2f | Mut-CNA %.2f\n",cmb("RNA","CNA"),cmb("RNA","Mutations"),cmb("Mutations","CNA")))

# --- SITE REFERENCE (secondary) ---
cat("\n--- Site reference ---\n")
site3 <- ifelse(site_all %in% c("BREAST","AXILLARY LYMPH NODE","BONE"),site_all,"Other")
pdfS <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2],Site=factor(site3,levels=c("BREAST","AXILLARY LYMPH NODE","BONE","Other")))
site_pal <- c("BREAST"="#3C5488","AXILLARY LYMPH NODE"="#F39B7F","BONE"="#8491B4","Other"="grey70")
save_pdf(ggplot(pdfS,aes(PC1,PC2,color=Site))+geom_point(size=2.6,alpha=.85)+scale_color_manual(values=site_pal,name="Biopsy site")+
  labs(title="RNA PCA by Biopsy Site (reference)",subtitle="No clear site structure; supports pooling",
       x=sprintf("PC1 (%.1f%%)",pv[1]),y=sprintf("PC2 (%.1f%%)",pv[2]))+theme_publication()+
  theme(panel.grid.major=element_line(color="grey92",linewidth=.4)),"S1_rna_pca_site_reference.pdf",9,7)
site_full <- clinical_matched$BX_LOCATION[match(final_samples,clinical_matched$SAMPLE_ID)]
okS <- names(table(site_full))[table(site_full)>=3]; selS <- site_full %in% okS
Rs <- t(rna_m[,final_samples[selS],drop=FALSE]); Rs <- Rs[,apply(Rs,2,var)>0,drop=FALSE]
pidS <- clinical_matched$PATIENT_ID[match(final_samples[selS],clinical_matched$SAMPLE_ID)]
Yst <- droplevels(factor(site_full[selS]))
permS <- if(length(pidS)>length(unique(pidS))) how(nperm=999,blocks=factor(pidS)) else 999
cat("  Site PERMANOVA (biopsy-level, n=",length(Yst),"):\n",sep=""); print(adonis2(dist(Rs)~Yst,permutations=permS))

## ==================== 3. VALIDATION (BER + confusion) =======================
cat("\n=== 3. VALIDATION ===\n")
pid <- clinical_matched$PATIENT_ID[match(rownames(X[[1]]),clinical_matched$SAMPLE_ID)]
# CV-safe re-clean: stricter event threshold so sparse binary features do not
# become constant inside a CV fold (otherwise perf() errors on zero variance).
X <- lapply(X, clean_block, min_ev=5)
Y_final <- droplevels(Y_final[rownames(X[[1]])])
kx <- list(Mutations=cap(c(8,8),X$Mutations), RNA=cap(c(25,50),X$RNA), CNA=cap(c(8,8),X$CNA))
final_diablo <- block.splsda(X, Y_final, ncomp=ncomp, keepX=kx, design=design)
pid <- clinical_matched$PATIENT_ID[match(rownames(X[[1]]),clinical_matched$SAMPLE_ID)]
pf <- perf(final_diablo,validation="Mfold",folds=5,nrepeat=50,dist=pred_dist,near.zero.var=TRUE,BPPARAM=BiocParallel::SerialParam())
wv <- pf$WeightedVote.error.rate[[pred_dist]]
cat(sprintf("  Subtype BER (comp %d) = %.3f (chance 0.67)\n",ncomp,wv[grep("BER",rownames(wv)),ncomp]))

cv_confusion <- function(Xl,Y,pid,keepX,k=5,reps=20){
  truth<-pred<-character(0)
  for(r in seq_len(reps)){ up<-unique(pid); fo<-sample(rep_len(seq_len(k),length(up))); names(fo)<-up; folds<-fo[pid]
    for(f in seq_len(k)){ te<-which(folds==f); tr<-which(folds!=f); if(length(unique(Y[tr]))<2) next
      Xtr<-lapply(Xl,function(m) clean_block(m[tr,,drop=FALSE],min_ev=3))
      Xte<-lapply(names(Xtr),function(b) Xl[[b]][te,colnames(Xtr[[b]]),drop=FALSE]); names(Xte)<-names(Xtr)
      if(any(sapply(Xtr,ncol)<2)) next
      kk<-list(Mutations=cap(keepX$Mutations,Xtr$Mutations),RNA=cap(keepX$RNA,Xtr$RNA),CNA=cap(keepX$CNA,Xtr$CNA))
      fit<-tryCatch(block.splsda(Xtr,Y[tr],ncomp=ncomp,keepX=kk,design=design,near.zero.var=TRUE),error=function(e)NULL)
      if(is.null(fit)) next
      pr<-tryCatch(predict(fit,newdata=Xte)$WeightedVote[[pred_dist]][,ncomp],error=function(e)NULL); if(is.null(pr)) next
      truth<-c(truth,as.character(Y[te])); pred<-c(pred,as.character(pr)) } }
  lv<-levels(Y); cm<-table(Truth=factor(truth,lv),Predicted=factor(pred,lv)); rec<-diag(cm)/rowSums(cm)
  list(cm=cm,overall=sum(diag(cm))/sum(cm),balanced=mean(rec,na.rm=TRUE),recall=rec) }
cc <- cv_confusion(X,Y_final,pid,kx)
cat("  Patient-disjoint confusion matrix:\n"); print(cc$cm)
cat(sprintf("  Overall acc %.3f | Balanced acc %.3f (chance 0.33)\n",cc$overall,cc$balanced))
cat("  Per-class recall:\n"); print(round(cc$recall,3))

## ==================== 4. CROSS-OMICS CONVERGENCE ============================
cat("\n=== 4. CONVERGENCE ===\n")
# 4a integration value (single BER, informational; formal test is in section 6-B5)
m_rna  <- mixOmics::splsda(X$RNA,Y_final,ncomp=ncomp,keepX=kx$RNA)
pf_rna <- perf(m_rna,validation="Mfold",folds=5,nrepeat=50,dist=pred_dist,BPPARAM=BiocParallel::SerialParam())
ber_rna <- pf_rna$error.rate$BER[ncomp,pred_dist]; ber_int <- wv[grep("BER",rownames(wv)),ncomp]
cat(sprintf("  Integrated BER %.3f | RNA-only BER %.3f | diff %+.3f\n",ber_int,ber_rna,ber_rna-ber_int))

# 4b cis-coupling (produces `cis`)
samp0 <- intersect(colnames(cna_matrix_final),colnames(rna_matrix_final))
genes0 <- intersect(rownames(cna_matrix_final),rownames(rna_matrix_final))
cis <- do.call(rbind,lapply(genes0,function(g){ cn<-as.numeric(cna_matrix_final[g,samp0]); ex<-as.numeric(rna_matrix_final[g,samp0]); ok<-is.finite(cn)&is.finite(ex)
  if(sum(ok)<10||sd(cn[ok])==0||sd(ex[ok])==0) return(data.frame(Gene=g,rho=NA,p=NA,n=sum(ok)))
  ct<-suppressWarnings(cor.test(cn[ok],ex[ok],method="spearman")); data.frame(Gene=g,rho=unname(ct$estimate),p=ct$p.value,n=sum(ok)) }))
cis <- cis[!is.na(cis$rho),]; cis$padj <- p.adjust(cis$p,"BH"); cis <- cis[order(-cis$rho),]
cat("  cis-coupling (biopsy-level, n=",length(samp0),"):\n",sep=""); print(cis,row.names=FALSE)
write.csv(cis,"output/cis_coupling_CNA_RNA.csv",row.names=FALSE)
cis_plot <- cis; cis_plot$sig <- ifelse(cis_plot$padj<0.05,"FDR<0.05","ns"); cis_plot$Gene <- factor(cis_plot$Gene,levels=cis_plot$Gene[order(cis_plot$rho)])
save_pdf(ggplot(cis_plot,aes(rho,Gene,fill=sig))+geom_col(color="black",linewidth=.2)+geom_vline(xintercept=0,linewidth=.4)+
  scale_fill_manual(values=c("FDR<0.05"="#D73027","ns"="grey75"),name=NULL)+
  labs(title="CNA-RNA cis-coupling (dosage sensitivity)",subtitle=paste0("Spearman rho | n=",length(samp0)),x="Spearman rho (CNA vs RNA)",y=NULL)+theme_publication(12),"10_cis_coupling.pdf",8,8)

# 4c top gene by subtype (ERBB2 = positive control)
top_gene <- as.character(cis$Gene[order(-cis$rho)][1])
dfC <- data.frame(Subtype=Y_final[samp0],Expr=as.numeric(rna_matrix_final[top_gene,samp0])); dfC <- dfC[!is.na(dfC$Subtype),]
save_pdf(ggplot(dfC,aes(Subtype,Expr,fill=Subtype))+geom_boxplot(outlier.shape=NA,alpha=.7)+geom_jitter(width=.15,size=1.5,alpha=.6)+
  scale_fill_manual(values=sub_cols,guide="none")+labs(title=paste0(top_gene," expression by molecular subtype"),
  subtitle="Positive control: copy-number-driven expression aligns with subtype",x=NULL,y=paste0(top_gene," (log2 expression)"))+theme_publication(12),"11_convergence_topgene.pdf",7,6)
cat(sprintf("  Top gene %s ~ subtype: ",top_gene)); print(kruskal.test(Expr~Subtype,data=dfC))

# 4d coupling score (all FDR<0.05 positive genes)
dosage_sig <- as.character(cis$Gene[cis$padj<0.05 & cis$rho>0])
cat("  Dosage-sensitive signature:",paste(dosage_sig,collapse=", "),"\n")
zc <- t(scale(t(cna_matrix_final[dosage_sig,samp0,drop=FALSE]))); zr <- t(scale(t(rna_matrix_final[dosage_sig,samp0,drop=FALSE])))
dsc <- data.frame(SAMPLE_ID=samp0,Coupling=colMeans(zc*zr,na.rm=TRUE),Subtype=Y_final[samp0]); dsc <- dsc[!is.na(dsc$Subtype),]
save_pdf(ggplot(dsc,aes(Coupling,fill=Subtype))+geom_histogram(bins=25,color="black",linewidth=.2,alpha=.85)+
  scale_fill_manual(values=sub_cols,name="Subtype")+labs(title="Sample-level cross-omics coupling score",
  subtitle=paste0("Dosage-sensitive signature (",length(dosage_sig)," genes) | n=",nrow(dsc)),x="Coupling score (mean CNA-RNA z-product)",y="Samples")+theme_publication(12),"12_coupling_score.pdf",8,6)
cat("  Coupling score ~ subtype: "); print(kruskal.test(Coupling~Subtype,data=dsc))

## ==================== 5. EXTERNAL VALIDATION (TCGA) =========================
cat("\n=== 5. EXTERNAL VALIDATION (TCGA-BRCA) ===\n")
have_tcga <- dir.exists(TCGA_DIR) && file.exists(paste0(TCGA_DIR,"data_mrna_seq_v2_rsem.txt"))
if(!have_tcga){ cat("  SKIPPED: set TCGA_DIR to a valid TCGA-BRCA export folder.\n") } else {
  sig_genes <- clean_num(names(final_diablo$loadings$RNA[,1])[final_diablo$loadings$RNA[,1]!=0])
  tcga_rna  <- read.delim(paste0(TCGA_DIR,"data_mrna_seq_v2_rsem.txt"),check.names=FALSE,stringsAsFactors=FALSE)
  tcga_clin <- read.delim(paste0(TCGA_DIR,"data_clinical_sample.txt"),skip=4,stringsAsFactors=FALSE,na.strings=na_s)
  tcga_pat  <- read.delim(paste0(TCGA_DIR,"data_clinical_patient.txt"),skip=4,stringsAsFactors=FALSE,na.strings=na_s)
  pam <- toupper(tcga_pat$SUBTYPE[match(tcga_clin$PATIENT_ID,tcga_pat$PATIENT_ID)])
  tcga_clin$Subtype <- dplyr::case_when(grepl("HER2",pam)~"HER2_Positive",grepl("LUMA|LUMB",pam)~"HR+/HER2-",grepl("BASAL",pam)~"Triple_Negative",TRUE~NA_character_)
  g <- clean_num(make.names(tcga_rna$Hugo_Symbol,unique=TRUE)); M<-as.matrix(tcga_rna[,-c(1:2)]); rownames(M)<-g; mode(M)<-"numeric"
  if(max(M,na.rm=TRUE)>100) M<-log2(M+1); M <- M[rownames(M) %in% sig_genes,,drop=FALSE]
  common <- intersect(colnames(M),tcga_clin$SAMPLE_ID); M<-M[,common,drop=FALSE]
  Ytc <- factor(tcga_clin$Subtype[match(common,tcga_clin$SAMPLE_ID)]); keep<-!is.na(Ytc); M<-M[,keep,drop=FALSE]; Ytc<-droplevels(Ytc[keep])
  okc <- names(table(Ytc))[table(Ytc)>=10]; keep<-Ytc %in% okc; M<-M[,keep,drop=FALSE]; Ytc<-droplevels(Ytc[keep])
  Xt <- t(M); Xt <- Xt[,apply(Xt,2,var)>0,drop=FALSE]
  cat("  TCGA n =",nrow(Xt),"| signature genes found:",ncol(Xt),"/",length(sig_genes),"\n"); print(table(Ytc))
  pm_t <- adonis2(dist(Xt)~Ytc,permutations=999)
  cat(sprintf("  TCGA PERMANOVA (signature): R2 = %.3f | p = %.3f\n",pm_t$R2[1],pm_t$`Pr(>F)`[1]))
  kxt <- min(50,ncol(Xt)); sp<-mixOmics::splsda(Xt,Ytc,ncomp=2,keepX=c(kxt,kxt))
  pf_t<-perf(sp,validation="Mfold",folds=5,nrepeat=50,dist="max.dist",BPPARAM=BiocParallel::SerialParam())
  cat(sprintf("  TCGA BER (comp2) = %.3f (chance %.2f)\n",pf_t$error.rate$BER[2,"max.dist"],1-1/nlevels(Ytc)))

  # coupling reproducibility across cohorts -> cmp (+ Figure 13)
  tcga_cna <- read.delim(paste0(TCGA_DIR,"data_cna.txt"),check.names=FALSE,stringsAsFactors=FALSE)
  gc2 <- if("Hugo_Symbol" %in% colnames(tcga_cna)) "Hugo_Symbol" else colnames(tcga_cna)[1]
  m2 <- intersect(c("Hugo_Symbol","Entrez_Gene_Id"),colnames(tcga_cna))
  Ct <- as.matrix(tcga_cna[,setdiff(colnames(tcga_cna),m2),drop=FALSE]); rownames(Ct)<-tcga_cna[[gc2]]; mode(Ct)<-"numeric"
  Rt <- as.matrix(tcga_rna[,-c(1:2)]); rownames(Rt)<-clean_num(make.names(tcga_rna$Hugo_Symbol,unique=TRUE)); mode(Rt)<-"numeric"
  if(max(Rt,na.rm=TRUE)>100) Rt<-log2(Rt+1)
  g9 <- as.character(cis$Gene); comT <- intersect(colnames(Ct),colnames(Rt))
  tcga_cis <- do.call(rbind,lapply(g9,function(gg){ if(!(gg%in%rownames(Ct))||!(gg%in%rownames(Rt))) return(data.frame(Gene=gg,rho_tcga=NA))
    cn<-as.numeric(Ct[gg,comT]); ex<-as.numeric(Rt[gg,comT]); ok<-is.finite(cn)&is.finite(ex)
    if(sum(ok)<10||sd(cn[ok])==0||sd(ex[ok])==0) return(data.frame(Gene=gg,rho_tcga=NA))
    data.frame(Gene=gg,rho_tcga=cor(cn[ok],ex[ok],method="spearman")) }))
  cmp <- merge(data.frame(Gene=g9,rho_mbc=as.numeric(cis$rho)),tcga_cis,by="Gene"); cmp<-cmp[complete.cases(cmp),]
  ctc <- cor.test(cmp$rho_mbc,cmp$rho_tcga,method="spearman")
  cat(sprintf("  Cross-cohort coupling: Spearman rho = %.2f | p = %.3f (n=%d genes)\n",ctc$estimate,ctc$p.value,nrow(cmp)))
  save_pdf(ggplot(cmp,aes(rho_mbc,rho_tcga,label=Gene))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60")+
    geom_point(size=3,color="#3C5488")+ggrepel::geom_text_repel(size=3.5,fontface="bold")+
    labs(title="Cross-cohort reproducibility of dosage sensitivity",subtitle=sprintf("Driver panel | Spearman rho=%.2f",ctc$estimate),
         x="CNA-RNA coupling (MBCproject)",y="CNA-RNA coupling (TCGA)")+theme_publication(12)+xlim(-.3,.8)+ylim(-.3,.8),"13_coupling_reproducibility.pdf",7,7)
}

## ==================== 6. REVIEWER RE-ANALYSIS (B1-B5) =======================
cat("\n=== 6. REVIEWER RE-ANALYSIS ===\n")
samp <- rownames(X[[1]]); pidA <- clinical_matched$PATIENT_ID[match(samp,clinical_matched$SAMPLE_ID)]
sig_genes <- clean_num(names(final_diablo$loadings$RNA[,1])[final_diablo$loadings$RNA[,1]!=0])

## B1: discovery PERMANOVA on the SAME signature genes (patient-level)
cat("\n[B1] discovery PERMANOVA on signature genes (patient-level)\n")
rn <- clean_num(rownames(rna_matrix_final)); Msig <- rna_matrix_final[rn %in% sig_genes,samp,drop=FALSE]
k1 <- one_per_patient(pidA); Xsig <- t(Msig[,k1,drop=FALSE]); Xsig<-Xsig[,apply(Xsig,2,var)>0,drop=FALSE]
Yk1 <- droplevels(factor(Y_final[samp[k1]]))
pm_sig <- adonis2(dist(Xsig)~Yk1,permutations=999)
cat(sprintf("-> Discovery R2 (signature genes) = %.3f | p = %.3f   [compare to TCGA 25-gene R2]\n",pm_sig$R2[1],pm_sig$`Pr(>F)`[1]))

## B2: cis-coupling at PATIENT level
cat("\n[B2] cis-coupling, one biopsy per patient\n")
samp_pt <- samp0[one_per_patient(clinical_matched$PATIENT_ID[match(samp0,clinical_matched$SAMPLE_ID)])]
cis_pt <- do.call(rbind,lapply(genes0,function(g){ cn<-as.numeric(cna_matrix_final[g,samp_pt]); ex<-as.numeric(rna_matrix_final[g,samp_pt]); ok<-is.finite(cn)&is.finite(ex)
  if(sum(ok)<10||sd(cn[ok])==0||sd(ex[ok])==0) return(data.frame(Gene=g,rho=NA,p=NA,n=sum(ok)))
  ct<-suppressWarnings(cor.test(cn[ok],ex[ok],method="spearman")); data.frame(Gene=g,rho=unname(ct$estimate),p=ct$p.value,n=sum(ok)) }))
cis_pt <- cis_pt[!is.na(cis_pt$rho),]; cis_pt$padj<-p.adjust(cis_pt$p,"BH"); cis_pt<-cis_pt[order(-cis_pt$rho),]
cat("  patient-level n =",length(samp_pt),"\n"); print(cis_pt,row.names=FALSE)

## B3: ERBB2-free coupling score ~ subtype (removes HER2 circularity)
cat("\n[B3] coupling score WITHOUT ERBB2\n")
dosage_noE <- setdiff(dosage_sig,"ERBB2"); cat("  genes:",paste(dosage_noE,collapse=", "),"\n")
if(length(dosage_noE)>=1){
  zc2 <- t(scale(t(cna_matrix_final[dosage_noE,samp0,drop=FALSE]))); zr2 <- t(scale(t(rna_matrix_final[dosage_noE,samp0,drop=FALSE])))
  d3 <- data.frame(SAMPLE_ID=samp0,Score=colMeans(rbind(zc2*zr2),na.rm=TRUE),Subtype=Y_final[samp0]); d3<-d3[!is.na(d3$Subtype),]
  pid3 <- clinical_matched$PATIENT_ID[match(d3$SAMPLE_ID,clinical_matched$SAMPLE_ID)]; d3p<-d3[one_per_patient(pid3),]
  cat(sprintf("  patient-level n = %d | Kruskal-Wallis p = %.4g\n",nrow(d3p),kruskal.test(Score~Subtype,data=d3p)$p.value))
  for(g in dosage_noE){ z<-as.numeric(scale(as.numeric(cna_matrix_final[g,d3$SAMPLE_ID])))*as.numeric(scale(as.numeric(rna_matrix_final[g,d3$SAMPLE_ID])))
    cat("   ",g,"per-sample coupling by subtype:"); print(round(tapply(z,d3$Subtype,mean,na.rm=TRUE),2)) }
} else cat("  (only ERBB2 was dosage-sensitive; report honestly that the axis is HER2-driven)\n")

## B4: frozen-model transfer MBC -> TCGA
cat("\n[B4] frozen-model transfer MBCproject -> TCGA\n")
if(!have_tcga){ cat("  SKIPPED (no TCGA).\n") } else {
  rnM <- clean_num(rownames(rna_matrix_final)); MM <- rna_matrix_final[rnM %in% sig_genes,,drop=FALSE]; rownames(MM)<-clean_num(rownames(MM))
  Xtr <- scale(t(MM)); Ytr <- droplevels(factor(Y_final[rownames(Xtr)]))
  shared <- intersect(colnames(Xtr),colnames(Xt)); shared <- shared[apply(Xtr[,shared,drop=FALSE],2,function(z) all(is.finite(z))&sd(z)>0)]
  cat("  shared signature genes:",length(shared),"\n")
  if(length(shared)>=3){
    fit <- mixOmics::splsda(Xtr[,shared,drop=FALSE],Ytr,ncomp=ncomp,keepX=c(length(shared),length(shared)))
    prT <- predict(fit,newdata=Xt[,shared,drop=FALSE])$class$max.dist[,ncomp]
    cmT <- table(Predicted=factor(prT,levels(Ytc)),True=Ytc); print(cmT)
    recT <- diag(cmT[levels(Ytc),levels(Ytc)])/colSums(cmT)
    cat(sprintf("  Frozen-model balanced acc = %.3f (chance %.2f)\n",mean(recT,na.rm=TRUE),1/nlevels(Ytc)))
  } else cat("  too few shared genes.\n")
}

## B5: integration vs RNA-only, tested across repeats (paired Wilcoxon)
cat("\n[B5] integration vs RNA-only (paired test)\n")
Yc <- droplevels(Y_final[rownames(X[[1]])]); pidC <- clinical_matched$PATIENT_ID[match(rownames(X[[1]]),clinical_matched$SAMPLE_ID)]
BERf <- function(truth,pred,lv){ cm<-table(factor(truth,lv),factor(pred,lv)); mean(1-diag(cm)/rowSums(cm),na.rm=TRUE) }
nrep<-25; kf<-5; bi<-br<-rep(NA_real_,nrep)
for(r in seq_len(nrep)){ set.seed(1000+r); up<-unique(pidC); fo<-sample(rep_len(seq_len(kf),length(up))); names(fo)<-up; folds<-fo[pidC]
  ti<-tp<-ri<-rp<-character(0)
  for(f in seq_len(kf)){ te<-which(folds==f); tr<-which(folds!=f); if(length(unique(Yc[tr]))<2) next
    Xtr<-lapply(X,function(m) clean_block(m[tr,,drop=FALSE],min_ev=3))
    Xte<-lapply(names(Xtr),function(b) X[[b]][te,colnames(Xtr[[b]]),drop=FALSE]); names(Xte)<-names(Xtr)
    if(any(sapply(Xtr,ncol)<2)) next
    kk<-list(Mutations=cap(kx$Mutations,Xtr$Mutations),RNA=cap(kx$RNA,Xtr$RNA),CNA=cap(kx$CNA,Xtr$CNA))
    mi<-tryCatch(block.splsda(Xtr,Yc[tr],ncomp=ncomp,keepX=kk,design=design,near.zero.var=TRUE),error=function(e)NULL)
    mr<-tryCatch(mixOmics::splsda(Xtr$RNA,Yc[tr],ncomp=ncomp,keepX=kk$RNA,near.zero.var=TRUE),error=function(e)NULL)
    if(is.null(mi)||is.null(mr)) next
    pi<-tryCatch(predict(mi,newdata=Xte)$WeightedVote[[pred_dist]][,ncomp],error=function(e)NULL)
    pq<-tryCatch(predict(mr,newdata=Xte$RNA)$class[[pred_dist]][,ncomp],error=function(e)NULL)
    if(is.null(pi)||is.null(pq)) next
    ti<-c(ti,as.character(Yc[te])); tp<-c(tp,as.character(pi)); ri<-c(ri,as.character(Yc[te])); rp<-c(rp,as.character(pq)) }
  lv<-levels(Yc); if(length(ti)>0){ bi[r]<-BERf(ti,tp,lv); br[r]<-BERf(ri,rp,lv) } }
ok<-is.finite(bi)&is.finite(br); bi<-bi[ok]; br<-br[ok]
cat(sprintf("  valid repeats %d/%d | Integrated median %.3f | RNA-only median %.3f\n",sum(ok),nrep,median(bi),median(br)))
cat(sprintf("  Paired Wilcoxon p = %.4g | median improvement = %.3f\n",wilcox.test(br,bi,paired=TRUE)$p.value,median(br-bi)))

cat("\n", strrep("=",70), "\n  ALL DONE. Figures in output/. Report B1/B3/B4/B5 above.\n", strrep("=",70), "\n")
cat("\n--- sessionInfo() ---\n"); print(sessionInfo())
################################################################################
