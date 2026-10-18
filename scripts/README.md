# Analysis Scripts

This directory contains the scripts used in the RDEB–SCC transcriptomics analysis workflow.

## Planned workflow

1. RNA-seq quality control using FastQC.
2. Adapter and quality trimming using fastp or Trimmomatic.
3. Genome alignment using HISAT2.
4. Alignment processing using SAMtools.
5. Gene-level read quantification using featureCounts.
6. Differential gene expression analysis using DESeq2.
7. Variant processing and annotation using BCFtools and SnpEff or VEP.
8. Functional enrichment, co-expression, immune microenvironment, and regulatory network analyses.

Scripts will be added and documented as each analysis step is implemented and tested.
