# install required packages

# It is important that you have R version 4.0.2
if (paste0(R.version$major, ".", R.version$minor) != "4.0.4") {
  stop("Please use R version 4.0.3")
}


# Install CRAN packages
cran_pkgs <- c("Rcpp", "openssl", "CpGassoc", "rmarkdown", "knitr",
               "matrixStats", "reshape", "glmnet", "statmod", "XML",
               "BiocManager", "pryr", "data.table", "qqman", "RPMM",
               "MASS", "sandwich", "lmtest", "foreach", "stringi",
               "doParallel", "magrittr", "purrr")
install.packages(methpackagesCRAN)

# Install Bioc pacages
bioc_pkgs <- c("minfi", "FlowSorted.Blood.450k", "missMethyl", "ENmix",
               "IlluminaHumanMethylation450kanno.ilmn12.hg19",
               "IlluminaHumanMethylation450kmanifest",
               "IlluminaHumanMethylationEPICmanifest", "sva",
               "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
               "DMRcate", "shinyMethyl", "bumphunter", "wateRmelon",
               "FDb.InfiniumMethylation.hg19")
BiocManager::install(methpackagesBioC)

# update all pacages
update.packages(ask = FALSE)
