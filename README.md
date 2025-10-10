# ATAC-seq Data Processing, Peak Calling, Motif, and Differential Accessibility Pipeline

This repository provides end-to-end R scripts for bulk ATAC-seq data analysis, including: quality control, alignment, peak calling, peak annotation, motif discovery, and differential accessibility analysis. The provided scripts are modular and can be adapted to both human and mouse datasets.

***

## Repository Structure

| File Name                                                 | Description                                                    |
|-----------------------------------------------------------|----------------------------------------------------------------|
| `Script-1-file-processing-and-peak-files-from-BAM-files.R`| Raw data QC, alignment, BAM processing, peak calling, and annotation for a single ATAC-seq dataset.|
| `Script-2-differential-ATAC-Motifs.R`                     | Consensus peakset generation, differential analysis across samples, and motif discovery (with chromVAR). |

***

## Workflow Overview

1. **Script-1-file-processing-and-peak-files-from-BAM-files.R**
   - Performs quality control on raw FASTQ files.
   - Aligns to the reference genome (using Rsubread or Rbowtie2).
   - Sorts and indexes BAM files.
   - Calls peaks (MACS2) and removes ENCODE blacklisted regions.
   - Annotates peaks to genes using ChIPseeker.
   - Generates QC reports (ChIPQC, ATACseqQC).
   - Prepares nucleosome-free/mono-/di-nucleosome region BAMs.
   - Outputs: sorted BAM, peaks (narrowPeak), annotated peaks (Excel), QC reports.

2. **Script-2-differential-ATAC-Motifs.R**
   - Suitable for datasets with multiple conditions (e.g., tissue comparisons).
   - Defines consensus and non-redundant peaks from multiple samples.
   - Generates a count matrix for DESeq2 or other count-based methods.
   - Differential analysis of chromatin accessibility.
   - Annotates differentially accessible regions and performs GO enrichment.
   - Motif discovery using MotifDb and JASPAR2020.
   - Integrates chromVAR for motif deviation analysis.
   - Outputs: normalized counts, annotation tables, motif analysis results, chromVAR variability statistics.

***

## Requirements

- **R version 4.0+**
- Key Bioconductor packages:
  - ShortRead, GenomicAlignments, Rsubread, Rsamtools, rtracklayer, BSgenome, ChIPQC, ChIPseeker, ATACseqQC, soGGi, rGREAT, clusterProfiler, MotifDb, JASPAR2020, TFBSTools, Biostrings, chromVAR, openxlsx, dplyr, ggplot2, pheatmap
- **MACS2** and **python** in your system path (for peak calling).
- Human (`BSgenome.Hsapiens.UCSC.hg38`) or mouse (`BSgenome.Mmusculus.UCSC.mm10`) reference genome.

Install missing Bioconductor packages using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
    "ShortRead", "GenomicAlignments", "Rsubread", "Rsamtools", "rtracklayer",
    "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10", "ChIPQC", "ChIPseeker", "ATACseqQC", "soGGi",
    "rGREAT", "clusterProfiler", "MotifDb", "JASPAR2020", "TFBSTools", "Biostrings", "chromVAR"
))
```

***

## Usage

1. Adjust file paths, sample names, and reference genome as needed at the start of each script.
2. Run `Script-1-file-processing-and-peak-files-from-BAM-files.R` from top to bottom for initial QC, alignment, peak calling, and annotation of a single dataset.
3. To perform comparative/differential accessibility and motif analysis on multiple datasets, run `Script-2-differential-ATAC-Motifs.R`.
4. Outputs generated at each step (e.g., BAMs, peak files, annotations) will be saved in user-specified directories.

***

## Notes

- Results (e.g., peaks, annotations, plots) depend on the sample type, sequencing depth, and reference genome version selected.
- The scripts provide points where user inspection is recommended (e.g., fragment size distribution, QC metrics).
- Downstream visualization and interaction with results (e.g., in IGV or R) are possible using the output bigWig and Excel files.
- For non-model organisms, reference genomes and annotation steps will need to be adjusted accordingly.
