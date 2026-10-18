# RDEB–SCC Transcriptomics Analysis

## Overview

This project investigates the molecular changes associated with the progression of Recessive Dystrophic Epidermolysis Bullosa (RDEB) to RDEB-associated cutaneous Squamous Cell Carcinoma (RDEB-SCC).

Using transcriptomic data from healthy, RDEB, and SCC samples, the study aims to identify differentially expressed genes (DEGs), explore biological pathways and regulatory networks, characterize the immune microenvironment, and prioritize potential molecular biomarkers associated with disease progression.

**Project status:** Research project — analysis in progress.

## Objectives

1. Identify differentially expressed genes across healthy, RDEB, and SCC samples.
2. Integrate differential gene expression results with genetic variant annotation to identify candidate genes associated with RDEB-SCC progression.
3. Investigate the biological functions and pathways of candidate genes using Gene Ontology (GO), Reactome, and network-based analyses.
4. Explore protein–protein interaction (PPI) networks, co-expression modules, and transcription factor regulatory relationships.
5. Characterize immune cell infiltration patterns and investigate associations between gene expression and immune cell populations.
6. Prioritize candidate biomarkers using integrated gene expression, functional, regulatory, and immune-related evidence.

## Analysis Workflow

### 1. RNA-seq Preprocessing and Quality Control
- **FastQC:** Assess raw sequencing read quality.
- **fastp / Trimmomatic:** Remove adapters and low-quality reads.

### 2. Genome Alignment and Read Quantification
- **HISAT2:** Align sequencing reads to the reference genome.
- **SAMtools:** Process and inspect alignment files.
- **featureCounts:** Generate gene-level read counts.

### 3. Differential Gene Expression Analysis
- **DESeq2:** Identify differentially expressed genes.
- **Visualization:** Generate PCA plots, MA plots, and volcano plots.
- **Comparisons:**
  - RDEB vs SCC
  - SCC vs Healthy
  - RDEB vs Healthy

### 4. Variant Analysis and Integration
- **BCFtools:** Process and analyze variant calls.
- **SnpEff / VEP:** Annotate genetic variants and their predicted effects.
- Integrate variant annotation with differential expression results to prioritize candidate genes.

*Note: RNA-seq variant analysis identifies variants detectable in the RNA-seq data; it does not replace comprehensive DNA-based variant analysis.*

### 5. Candidate Gene and Functional Analysis
- Investigate the 91 candidate genes identified during the study.
- Perform Gene Ontology and Reactome pathway enrichment.
- Explore relevant signaling pathways and biological processes.
- Construct protein–protein interaction networks using appropriate network analysis tools.

### 6. Co-expression and Immune Microenvironment Analysis
- **WGCNA:** Identify gene co-expression modules and investigate their relationships with disease groups.
- **xCell:** Estimate relative immune and stromal cell enrichment from transcriptomic data.
- Investigate immune-associated genes, immune cell differences, and gene–immune correlations.
- Examine the reported APOBEC3A-associated module and 13 significant immune cell types.

### 7. Transcription Factor and Biomarker Prioritization
- Investigate transcription factors and their potential regulatory relationships with candidate genes.
- Integrate differential expression, variant annotation, pathway, network, and immune analysis results.
- Prioritize seven candidate biomarkers for further investigation.

## Tools and Technologies

| Analysis | Tools |
|---|---|
| Quality control | FastQC |
| Read trimming | fastp / Trimmomatic |
| Genome alignment | HISAT2 |
| Alignment processing | SAMtools |
| Gene quantification | featureCounts |
| Differential expression | DESeq2 |
| Variant processing | BCFtools |
| Variant annotation | SnpEff / VEP |
| Functional enrichment | Gene Ontology, Reactome |
| Co-expression analysis | WGCNA |
| Immune microenvironment | xCell |
| Network and regulatory analysis | NetworkAnalyst and appropriate PPI/TF resources |
| Programming and analysis | R, Linux |

## Key Study Outputs

- Differentially expressed genes across three sample-group comparisons.
- A set of 91 candidate genes for integrated analysis.
- Functional enrichment and pathway analysis results.
- Gene co-expression modules, including an APOBEC3A-associated module under investigation.
- Immune microenvironment analysis, including 13 reported significant immune cell types.
- A prioritized set of seven candidate biomarkers.

*These outputs are subject to validation and may be updated as the analysis progresses.*

## Repository Structure

```text
RDEB-SCC-Transcriptomics-Analysis/
├── README.md
├── data/
│   └── README.md
├── scripts/
│   └── README.md
├── results/
│   └── README.md
├── figures/
│   └── README.md
├── docs/
│   └── workflow.md
└── .gitignore
```

The `data/` directory should contain documentation and, where permitted, data links rather than private or restricted datasets. Analysis scripts, result tables, and figures should be added as they become available.

## Reproducibility and Limitations

- Document the reference genome, annotation release, software versions, and analysis parameters.
- Record sample metadata and group assignments accurately.
- Apply appropriate multiple-testing correction in differential expression and enrichment analyses.
- Validate candidate biomarkers using independent datasets or experimental evidence where possible.
- Interpret immune cell enrichment and variant effects according to the limitations of the underlying methods.

## Research Status

This project is ongoing. The final conclusions will depend on completion of the analyses, statistical validation, and biological interpretation of the results.

