# RDEB-SCC Transcriptomics Analysis

## Overview

Transcriptomic analysis of Recessive Dystrophic Epidermolysis Bullosa (RDEB), Squamous Cell Carcinoma (SCC), and healthy samples to identify differentially expressed genes, molecular pathways, immune-related changes, and potential biomarkers associated with disease progression.

## Objectives

- ## Objectives

- Identify differentially expressed genes (DEGs) and genetic variants associated with the progression of Recessive Dystrophic Epidermolysis Bullosa (RDEB) to RDEB-associated Cutaneous Squamous Cell Carcinoma (RDEB-SCC) using RNA-seq data.

- Integrate differential gene expression and variant annotation analyses to identify candidate genes involved in RDEB-SCC progression.

- Investigate the biological functions, signaling pathways, and regulatory networks of candidate genes using functional enrichment, protein–protein interaction (PPI), and transcription factor analyses.

- Identify key transcription factors regulating candidate genes and investigate their potential roles in RDEB-SCC progression.

- Characterize the immune microenvironment and evaluate gene–immune correlations to identify immune-associated molecular signatures involved in RDEB-SCC progression.

## Analysis Workflow

                         RNA-seq DATA
                             │
                             ↓
                    QUALITY CONTROL
                         FastQC
                             │
                             ↓
                     READ TRIMMING
                  fastp / Trimmomatic
                             │
                             ↓
                    GENOME ALIGNMENT
                         HISAT2
                             │
                             ↓
                   READ QUANTIFICATION
                      featureCounts
                             │
                             ↓
                DIFFERENTIAL EXPRESSION
                         DESeq2
                             │
              ┌──────────────┴──────────────┐
              ↓                             ↓
        DEG IDENTIFICATION          PCA / MA / Volcano
              │
              ↓
                                          RNA-seq DATA
                             │
                             ↓
                    QUALITY CONTROL
                         FastQC
                             │
                             ↓
                     READ TRIMMING
                  fastp / Trimmomatic
                             │
                             ↓
                    GENOME ALIGNMENT
                         HISAT2
                             │
                             ↓
                   READ QUANTIFICATION
                      featureCounts
                             │
                             ↓
                DIFFERENTIAL EXPRESSION
                         DESeq2
                             │
              ┌──────────────┴──────────────┐
              ↓                             ↓
        DEG IDENTIFICATION         VARIANT ANALYSIS
                               
              │                            │
              ↓                            ↓                      
          PCA / MA / Volcano       Variant Calling + Annotation         
              │                              │
              ↓                              ↓
      Differential
      expression gene       +       VARIANT INTEGRATION
              │
              ↓
     91 CANDIDATE GENES
              │
       ┌──────┼───────────────┐
       ↓      ↓               ↓
      GO   REACTOME      FUNCTIONAL ANALYSIS
       │      │               │
       └──────┴───────┬───────┘
                      ↓
          ┌────────────┼────────────┐
          ↓                         ↓
       WGCNA /                   IMMUNE       
    Co-expression                ANALYSIS     
          │                        │
          ↓                        ↓
    APOBEC3A                  13 Significant 
    Associated                 Immune Cell     
      Module                      Types         
          │                         │
          └────────────│────────────┘
                       ↓
              TRANSCRIPTION FACTOR
                    ANALYSIS
                       │
                       ↓
              INTEGRATED BIOMARKER
                  PRIORITIZATION
                       │
                       ↓
                 7 CANDIDATE
                 BIOMARKERS
                     
## Tools & Technologies

- FastQC
- Fastp / Trimmomatic
- HISAT2
- SAMtools
- featureCounts
- DESeq2
- BCFtools
- SnpEff / VEP
- Reactome
- xCell
- R
- WGCNA
- xcell
- gene ontology
- NetworkAnalyst

## Comparisons

- RDEB vs SCC
- SCC vs Healthy
- RDEB vs Healthy

## Project Status

🔬 Research project — analysis in progress.
