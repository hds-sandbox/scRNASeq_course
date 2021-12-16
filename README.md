# single-cell RNAseq workshop
Repository for the scRNA-Seq course of the Danish Health Data Science Sandbox project.

Created: November 2021

NOTE: workshop still under construction

This workshop contains a basic tutorial on how to approach scRNAseq experiments, starting from fastq files out of the sequencer. Thus, the workshop does not include any information of laboratory protocols, library preparation or any wet-lab related procedures. This workshop is based on a collection of modified tutorials from the [nf-core](https://nf-co.re/scrnaseq) pipeline, the [Seurat](https://satijalab.org/seurat/) and [gProfiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) vignettes.

## Syllabus:
1. Introduction to single cell RNA-Seq
2. Preprocessing of RNAseq reads (fastq). 
	- Trimming and filtering (TrimGalore)
	- Alignment (STARsolo)
	- Demultiplexing (STARsolo)
	- Feature count (STARsolo)
	- QC (FastQC and MultiQC)

3. Read normalization and QC (Seurat)
4. Exploratory analysis and clustering (Seurat)
5. Differential Expression Analysis (Seurat)
6. Pseudotime and trajectories (Seurat)
7. Functional Analysis (gprofiler2)

## Workshop requirements:
- Knowledge of R, Rstudio and Rmarkdown
- Basic knowledge of scRNAseq techonology
- Basic knowledge of data science and statistics such as PCA, clustering and statistical testing

## Intended use
The aim of this repository is to run a comprehensive but introductory workshop on sc-RNAseq bioinformatic analyses. Each of the modules of this workshop is accompanied by a powerpoint slideshow explaining the steps and the theory behind a typical bioinformatics analysis (ideally with a teacher). Many of the slides are annotated with extra information and/or point to original sources for extra reading material. The [example Rmarkdown](./Notebooks/slides/RNAseq_analysis_basics.Rmd) compiles modules 3-7 and can be used as a stand-alone template for a standard scRNA-Seq analysis.

## Acknowledgements:
- [Center for Health Data Science](https://heads.ku.dk/), University of Copenhagen.
