# Molecular signatures of placental dysfunction classify preeclampsia [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24148649.svg)](https://doi.org/10.5281/zenodo.24148649)

## Abstract
Preeclampsia is a leading cause of pregnancy-related death, accounting for over 50,000 maternal and 500,000 fetal deaths worldwide each year. Preeclampsia has been associated with confined placental mosaicism, underscoring the potential role of placental variants in driving adverse pregnancy outcomes. Here, we explored expressed rare variant burden and hypoxia response with respect to clinical maternal outcomes. We found that expressed rare variant burden was the highest in early preterm delivery preeclamptic placentas. At the single-cell level, we observed higher expressed rare variant burden in extravillous trophoblasts of preterm preeclamptic placentas. Mutational signatures associated with defective DNA mismatch repair and hypoxia response score are able to classify placental villous tissue and circulating cell-free RNA with preeclampsia. Together, these findings provide promising insights into the clinical correlates of placental expressed rare variant burden and its potential utility as a circulating biomarker for early clinical diagnosis of preeclampsia.

## Repo Contents
- rna_seq_processing: code to call and process RNA-seq variants
- wgs_validation: code to call and process WGS variants and compare to RNA-seq variants
- single_cell: code to process and visualize single-cell RNA-seq
- figure_visualization: code to create figures

## Data availability:
- Placental villous tissue RNA-seq: Raw and processed sequencing data from 76 placentas generated in this study are available from Sequence Read Archive under accession SRP466001 and BioProject under accession PRJNA1027377. Raw and processed sequencing data from 36 placentas were a reanalysis available for download from GSE186257 (Horii et al, 2023). Raw data from the validation cohort was downloaded from GSE234729 (Aisagbonhi et al, 2023).
- Circulating cell-free RNA - Raw data was downloaded from SRP352519 (Moufarrej et al, 2022).
- Deep placenta single-cell RNA-seq: Raw data was downloaded from GSE290578 (Paik et al, 2026).


## Dependencies
This code was developed using macOS Monterey/Tahoe (12.1/26.2) and Rocky Linux (8.9). Please refer to each package for its installation instructions and estimated installation time. 

- Dependencies for rna_seq_processing:
  - bcbionext-gen was installed following: https://github.com/bcbio/bcbio-nextgen

  | conda packages | Version |
  |---|---|
  | BCFtools | 1.15.1 |
  | VCFtools | 0.1.17 |
  | VEP | 107.0 |
  | vcf2maf | 1.6.21 |
  | bedtools | 2.30.0 |

  | R packages | Version |
  |---|---|
  | R-base | 4.4.3 |
  | Glimma | 2.12.0 |
  | RNASeqPower | 1.42.0 |
  | mygene | 1.38.0 |
  | GenomicFeatures | 1.54.4 |
  | AnnotationDbi | 1.64.1 |
  | Biobase | 2.62.0 |
  | GenomicRanges | 1.54.1 |
  | GenomeInfoDb | 1.38.8 |
  | IRanges | 2.36.0 |
  | S4Vectors | 0.40.2 |
  | BiocGenerics | 0.48.1 |
  | fs | 1.6.4 |
  | readr | 2.1.5 |
  | tximport | 1.30.0 |
  | edgeR | 4.0.16 |
  | limma | 3.58.1 |
  | ggplot2 | 3.5.2 |
  | ggbreak | 0.1.6 |
  | tidyr | 1.3.1 |
  | gdata | 3.0.0 |
  | data.table | 1.15.4 |
  | dplyr | 1.1.4 |
  | maftools | 2.18.0 |
  | future | 1.33.2 |
  | ADAPTS | 1.0.22 |
  | MutationalPatterns | 3.16.0 |
  | NMF | 0.28 |

- Dependencies for wgs_validation:

    | conda packages | Version |
    |---|---|
    | GATK4 | 4.2.6.1 |
    | BWA | 0.7.17 |
    | Samtools | 1.15.1 |
    | BCFtools | 1.15.1 |

- Dependencies for single_cell:

    | conda packages | Version |
    |---|---|
    | star | 2.7.11b |
    | umi_tools |1.1.2  |
    | samtools | 1.15.1 |
    | python | 3.14.3 |
    
    | R packages | Version |
    |---|---|
    | R-base | 4.4.3 |
    | R-ggplot2 | 4.0.3 |
    | R-dplyr | 1.2.1 |
    | R-data.table | 1.17.8 |
    | R-SeuratObject | 5.4.0 |
    | DoubletFinder | 2.0.6 |
    | SoupX | 1.2.1 |
    | infercnv | 1.22.0 |
    | UCell | 2.10.1 |
    | emmeans | 2.0.4 |
    | lme4 | 1.1-38 |
    | Matrix | 1.7-5 |

- Dependencies for figure_visualization:

    | python packages | Version |
    |---|---|
    | python |3.12.1|
    | pandas | 2.2.2 |
    | seaborn | 0.11.2 |
    | numpy | 2.2.6 |
    | matplotlib | 3.8.4 |
    | matplotlib-venn | 1.1.2 |
    | scipy | 1.31.1 |
    | tableone | 0.9.1 |
    | scikit-learn | 1.6.1 |
    | confidenceinterval | 0.54.1 |
    | shap | 1.0.5 |
    | statmodels | 0.14.4 |
    | yellowbrink | 1.5 |
    
    | R packages | Version |
    |---|---|
    | R-base | 4.4.3 |
    | car | 3.0-13 |
    | carData | 3.0-5 |
    | lmtest | 0.9-40 |
    | zoo | 1.8-10 |
    | emmeans | 2.0.4 |
    | tidyverse | 2.0.0 |
