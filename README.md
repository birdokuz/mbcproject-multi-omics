# Reproducible, genome-wide copy-number-to-expression coupling in metastatic breast cancer

Analysis code for the manuscript **"Reproducible, Genome-Wide Copy-Number-to-Expression
Coupling in Metastatic Breast Cancer: An Integrative Multi-Omics Analysis of the MBCproject."**
Secondary analysis of the Metastatic Breast Cancer Project (MBCproject), with independent
validation in TCGA-BRCA.

## What the analysis does

- Integrates somatic mutation, RNA-expression, and copy-number (CNA) data with DIABLO
  (molecular subtype as the primary axis; biopsy site as a reference).
- Quantifies copy-number-to-expression (cis) coupling (dosage sensitivity) genome-wide,
  locates driver genes within that distribution, and tests its reproducibility in TCGA-BRCA.
- Shows a gene-specific mechanism: ERBB2 is dosage-driven, whereas CDH1 is inactivated by
  (truncating) mutation (double dissociation).
- Tests pathway-level enrichment and reproducibility of dosage sensitivity (GO + Reactome).

## Data availability

Data are **not** included and must be downloaded separately, each study into its **own**
folder (cBioPortal studies share identical file names — never mix them):

- **MBCproject** — cBioPortal study `brca_mbcproject_2022`
  (https://www.cbioportal.org/study/summary?id=brca_mbcproject_2022);
  also NCI GDC (CMI-MBC) and dbGaP (phs001709).
- **TCGA-BRCA** (external validation) — cBioPortal, e.g. `brca_tcga_pan_can_atlas_2018`.

Expected files per folder: `data_clinical_patient.txt`, `data_clinical_sample.txt`,
`data_mutations.txt`, `data_cna.txt`, `data_mrna_seq_v2_rsem.txt`.

## Scripts (run in order)

| Script | Purpose |
|--------|---------|
| `01_intake_qc.R` | Data intake / identifier-consistency check (run first; catches wrong-cohort files) |
| `02_pipeline.R` | Main pipeline: preprocessing, molecular subtype, figures, DIABLO, validation (BER, confusion), cross-omics convergence, TCGA external validation, reviewer re-analysis. Produces the core objects used downstream. |
| `03_genomewide_dosage.R` | Genome-wide CNA-RNA dosage sensitivity; driver percentiles; genome-wide cross-cohort reproducibility |
| `04_cdh1_mechanism.R` | CDH1 inactivation is mutation-driven vs ERBB2 dosage-driven (double dissociation) |
| `05_pathway_gsea.R` | GO-BP + Reactome GSEA on the genome-wide dosage ranking; pathway-level cross-cohort reproducibility |
| `06_clinical_exploratory.R` | Exploratory clinical-pathological associations (reported as a negative result in the manuscript; no survival data in MBCproject) |

Set `TCGA_DIR` at the top of `02_pipeline.R` to your TCGA-BRCA folder. If TCGA is absent,
the external-validation sections self-skip, and `03`/`05` still run on MBCproject alone.

## Requirements

- R (>= 4.3)
- CRAN: `tidyverse`, `ggplot2`, `ggrepel`, `patchwork`, `caret`, `vegan`, `permute`,
  `circlize`, `RColorBrewer`
- Bioconductor: `mixOmics`, `maftools`, `ComplexHeatmap`, `clusterProfiler`,
  `org.Hs.eg.db`, `ReactomePA`

```r
install.packages(c("tidyverse","ggplot2","ggrepel","patchwork","caret","vegan",
                   "permute","circlize","RColorBrewer"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("mixOmics","maftools","ComplexHeatmap","clusterProfiler",
                       "org.Hs.eg.db","ReactomePA"))
```

## Running

```r
setwd("path/to/MBCproject_folder")   # contains the 5 data_*.txt files
source("01_intake_qc.R")             # verify data consistency
source("02_pipeline.R")              # main pipeline + validation (set TCGA_DIR inside)
source("03_genomewide_dosage.R")     # genome-wide dosage sensitivity + reproducibility
source("04_cdh1_mechanism.R")        # CDH1 mutation vs dosage (double dissociation)
source("05_pathway_gsea.R")          # pathway enrichment + reproducibility
```

All analyses use `set.seed(42)`. Figures are written to `output/`.

## Notes on reproducibility

- Sample counts differ by analysis (integrated set n=115; subtype set n=64; patient-level
  inferential analyses n=55-64; genome-wide n=103; TCGA-BRCA n=945), as described in the
  manuscript; each analysis reports its own n.
- Inferential analyses are performed at the patient level (one biopsy per patient).
- TCGA subtype is approximated from PAM50 and does not perfectly match the IHC-based
  three-class definition; the dosage-coupling reproducibility does not use subtype labels.
- Methylation is not available in MBCproject; the CDH1 analysis tests the mutation component.

## Citation

Please cite the manuscript (details on publication) and the archived release DOI (Zenodo).

## License

MIT License (see `LICENSE`).
