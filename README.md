# single-cell RNAseq workshop
Repository for the scRNA-Seq course of the Danish Health Data Science Sandbox project.

Created: November 2021

NOTE: workshop still under construction

This workshop contains a basic tutorial on how to approach scRNAseq experiments, starting from fastq files out of the sequencer. Thus, the workshop does not include any information of laboratory protocols, library preparation or any wet-lab related procedures. 

## Syllabus:
1. Data explanation
2. Read normalization and QC (python: scanpy | R: Seurat)
3. Exploratory analysis and clustering (python: scanpy | R: Seurat)
4. Differential Expression Analysis (python: scanpy | R: Seurat)
5. Pseudotime and trajectories (python: scanpy/Palantir | R: Seurat/Slingshot)
6. Functional Analysis (python: scanpy/gseapy) 

## Workshop requirements:
- Knowledge of R, Rstudio and Rmarkdown or Python and Jupyter Notebooks.
- Basic knowledge of scRNAseq technology
- Basic knowledge of data science and statistics such as PCA, clustering and statistical testing

## Intended use
The aim of this repository is to run a comprehensive but introductory workshop on sc-RNAseq bioinformatic analyses. Each of the modules of this workshop is accompanied by a powerpoint slideshow explaining the steps and the theory behind a typical bioinformatics analysis (ideally with a teacher). Many of the slides are annotated with extra information and/or point to original sources for extra reading material. 

The [example Rmarkdown](./Notebooks/R/scRNAseq_Seurat.html) compiles modules 3-7 and can be used as a stand-alone template for a standard scRNA-Seq analysis. The analysis is in R is based on a collection of modified tutorials from the [nf-core](https://nf-co.re/scrnaseq) pipeline, the [Seurat](https://satijalab.org/seurat/), [gProfiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) and [slingshot](https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html#constructing-smooth-curves-and-ordering-cells) vignettes.

## Collaborators
- Samuele Soraggi. Python Data Scientist. University of Aarhus
- Jose Alejandro Romero Herrera. R Data Scientist. University of Copenhagen

## Acknowledgements:
- [Center for Health Data Science](https://heads.ku.dk/), University of Copenhagen.
