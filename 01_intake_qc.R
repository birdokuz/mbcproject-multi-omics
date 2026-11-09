###############################################################################
# MBCproject — VERI KABUL KAPISI (intake QC)
# -----------------------------------------------------------------------------
# AMAC: Tum pipeline'i calistirmadan ONCE, bes dosyanin AYNI export'tan gelip
#       gelmedigini ve kimliklerin tutarli oldugunu dogrulamak.
# KULLANIM: Yeni indirdiginiz tarball'i bir klasore acin, calisma dizinini
#           oraya alin ve bu script'i calistirin. Sorun varsa net mesajla durur.
###############################################################################

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

na_strings <- c("NA", "", "[Not Available]", "[Not Applicable]", "[Unknown]")

req <- c("data_clinical_patient.txt", "data_clinical_sample.txt",
         "data_mutations.txt", "data_cna.txt", "data_mrna_seq_v2_rsem.txt")

cat("== 1) Dosya varligi ==\n")
present <- file.exists(req); names(present) <- req
for (f in req) cat(sprintf("  [%s] %s\n", ifelse(present[f], "OK", "YOK"), f))
has_rna <- present["data_mrna_seq_v2_rsem.txt"]
if (!has_rna)
  cat("\n  UYARI: RNA dosyasi yok. Bu export'ta mRNA katmani bulunmuyor olabilir;\n",
      "  o zaman entegrasyonu mutation + CNA ile kurmaniz gerekir.\n")

## --- Yardimci: bir kimlik vektorunun sema dagilimini ozetle -----------------
id_scheme <- function(ids, label) {
  ids <- ids[!is.na(ids)]
  pref <- ifelse(grepl("^TCGA", ids), "TCGA",
          ifelse(grepl("^MTS",  ids), "MTS-T",
          ifelse(grepl("^MB",   ids), "MB-", "DIGER")))
  tb <- sort(table(pref), decreasing = TRUE)
  cat(sprintf("  %-10s n=%-5d | %s\n", label, length(ids),
              paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = "  ")))
  cat(sprintf("             ornek: %s\n", paste(head(ids, 3), collapse = ", ")))
  if (any(pref == "TCGA"))
    cat("   >>> DIKKAT: TCGA kimligi tespit edildi -- bu dosya MBCproject DEGIL!\n")
  invisible(ids)
}

## --- 2) Klinik ---------------------------------------------------------------
cat("\n== 2) Klinik kimlik semalari ==\n")
cs <- read.delim("data_clinical_sample.txt", skip = 4,
                 stringsAsFactors = FALSE, na.strings = na_strings)
cp <- read.delim("data_clinical_patient.txt", skip = 4,
                 stringsAsFactors = FALSE, na.strings = na_strings)
stopifnot(all(c("SAMPLE_ID","PATIENT_ID") %in% colnames(cs)))
id_scheme(cs$SAMPLE_ID,  "SAMPLE_ID")
id_scheme(cs$PATIENT_ID, "PATIENT_ID")
cat(sprintf("  Klinik: %d ornek, %d hasta\n",
            length(unique(cs$SAMPLE_ID)), length(unique(cs$PATIENT_ID))))

## --- 3) Her omik katmanin ornek kimlikleri -----------------------------------
cat("\n== 3) Omik katman kimlikleri ==\n")

mut <- read.delim("data_mutations.txt", stringsAsFactors = FALSE,
                  comment.char = "#", quote = "", fill = TRUE)
mut_ids <- unique(mut$Tumor_Sample_Barcode)
id_scheme(mut_ids, "MUT")

cna <- read.delim("data_cna.txt", stringsAsFactors = FALSE, check.names = FALSE)
# ilk 1-2 kolon gen bilgisidir (Hugo_Symbol / Entrez_Gene_Id); gerisi ornek
cna_meta <- intersect(c("Hugo_Symbol","Entrez_Gene_Id"), colnames(cna))
cna_ids  <- setdiff(colnames(cna), cna_meta)
id_scheme(cna_ids, "CNA")

if (has_rna) {
  rna <- read.delim("data_mrna_seq_v2_rsem.txt", stringsAsFactors = FALSE,
                    check.names = FALSE, nrows = 5)  # sadece basliklar yeterli
  rna_meta <- intersect(c("Hugo_Symbol","Entrez_Gene_Id"), colnames(rna))
  rna_ids  <- setdiff(colnames(rna), rna_meta)
  id_scheme(rna_ids, "RNA")
}

## --- 4) SAMPLE_ID duzeyinde ortusme -----------------------------------------
cat("\n== 4) SAMPLE_ID duzeyinde ortusme ==\n")
S <- cs$SAMPLE_ID
cat(sprintf("  MUT ∩ clinical: %d\n", length(intersect(mut_ids, S))))
cat(sprintf("  CNA ∩ clinical: %d\n", length(intersect(cna_ids, S))))
if (has_rna) cat(sprintf("  RNA ∩ clinical: %d\n", length(intersect(rna_ids, S))))
tri_sample <- if (has_rna) Reduce(intersect, list(mut_ids, cna_ids, rna_ids))
              else intersect(mut_ids, cna_ids)
cat(sprintf("  >>> UC KATMAN ortak ORNEK (SAMPLE_ID): %d\n", length(tri_sample)))

## --- 5) PATIENT_ID duzeyinde ortusme (semalar farkliysa asil koprü budur) ----
cat("\n== 5) PATIENT_ID duzeyinde ortusme ==\n")
to_patient <- function(ids) unique(na.omit(cs$PATIENT_ID[match(ids, cs$SAMPLE_ID)]))
mut_pat <- to_patient(mut_ids); cna_pat <- to_patient(cna_ids)
cat(sprintf("  MUT hasta: %d | CNA hasta: %d\n", length(mut_pat), length(cna_pat)))
if (has_rna) {
  rna_pat <- to_patient(rna_ids)
  cat(sprintf("  RNA hasta: %d\n", length(rna_pat)))
  tri_pat <- Reduce(intersect, list(mut_pat, cna_pat, rna_pat))
} else tri_pat <- intersect(mut_pat, cna_pat)
cat(sprintf("  >>> UC KATMAN ortak HASTA (PATIENT_ID): %d\n", length(tri_pat)))

## --- 6) KARAR ----------------------------------------------------------------
cat("\n== 6) Karar ==\n")
if (has_rna && any(grepl("^TCGA", rna_ids))) {
  cat("  X  RNA dosyasi TCGA -- YANLIS. Dogru MBCproject mRNA dosyasini indirin.\n")
} else if (length(tri_sample) >= 30) {
  cat(sprintf("  OK ORNEK duzeyinde entegrasyon uygun (n=%d). Pipeline'i baslatabilirsiniz.\n",
              length(tri_sample)))
} else if (length(tri_pat) >= 30) {
  cat(sprintf("  ~  ORNEK ortusmesi dusuk ama HASTA duzeyinde n=%d.\n", length(tri_pat)),
      "     Entegrasyonu PATIENT_ID uzerinden (hasta-basi tek ornek) kurun.\n")
} else {
  cat("  !  Ortusme cok dusuk. Dosyalarin ayni export'tan geldigini dogrulayin.\n")
}
cat("\n--- QC bitti ---\n")
###############################################################################
