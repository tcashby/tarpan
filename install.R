#install dependiences
list.of.packages <- c("shiny","shinyBS","shinythemes","DT","DBI","RSQLite","RColorBrewer",
                      "ggpubr","RCircos","config")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(list.of.packages)
rm(new.packages)
#install VariantAnnotation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

#older versions of R
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("VariantAnnotation")
