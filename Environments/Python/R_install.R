#remotes::install_github("theislab/kBET", lib='sandbox_scRNA_testAndFeedback/scrna-environment/lib/R/library/', force=TRUE)

#r = getOption("repos")
#r["CRAN"] = "http://cran.us.r-project.org "
#options(repos = r)
#install.packages("rlang", lib='sandbox_scRNA_testAndFeedback/scrna-environment/lib/R/library/')
remotes::install_github("r-lib/rlang", lib='sandbox_scRNA_testAndFeedback/scrna-environment/lib/R/library/', force=TRUE)
remotes::install_github("satijalab/sctransform", ref="develop", lib='sandbox_scRNA_testAndFeedback/scrna-environment/lib/R/library/', force=TRUE)
#needs installation of r-dtw and r-pheatmap through conda
#remotes::install_github("shenorrLab/cellAlign", lib='sandbox_scRNA_testAndFeedback/scrna-environment/lib/R/library/')
