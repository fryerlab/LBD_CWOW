#----------------------
# Load libraries
#     Libraries may by running `install.packages("packageName")` or  BiocManager::install("packageName")
#     Some libraries were installed from git repos:
#        install.packages("remotes")
#        remotes::install_github("GabrielHoffman/mvIC")
#-----------------------
library(ComplexUpset)
library("llapack")
library(gprofiler2)
library(airway)
library(enrichplot)
library(DOSE) 
library(plyr)
library(scales)
library(forcats)
library(rmarkdown)
library(BiocParallel)
library(dplyr)
library(edgeR)
library(limma)
library(ggrepel)
library(ggplot2)
library(gplots)
library(grDevices)
require(philentropy) 
library(rtracklayer)
library(stringr)
require(variancePartition) 
library(reshape)
library(Glimma)
library(plyr)
library(corrplot)
library(ggpubr)
library(tidyverse)
library(caret)
library(glmnet)
library(vroom)
library(matrixStats)
library("data.table")
library(DESeq2)
library(dittoSeq) 
library(Hmisc)
library(tidyr)
library(gridExtra)
library(grid)
require(openxlsx)
library(UpSetR)
library(mvIC) 
library(RColorBrewer)
library(vctrs)

#-----------------------
# paths, colors, shapes
#-----------------------
LBD <- "LBD"
AD <- "AD"
PA <- "PA"
CONTROL <- "CONTROL"
control_color <- "#4682B4" # gray
AD_color <- "#B4464B" # yellow gold
PA_color <- "#B4AF46" # brown gold
LBD_color <- "gray35" # green
control_shape <- c(15) # square
AD_shape <- c(16) # circle
PA_shape <- c(17) # triangle
LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B4", "#B4AF46","#B4464B", "gray35") # disease type colors 
SexColors <- c("#490092", "#D55E00") # orange and purple
colorbindColors <- dittoColors()
correlationColors <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

pathToRef = c("/LBD_CWOW/references/")
pathToRawData = c("/LBD_CWOW/raw_fastq/")

# Save files as PDFs
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}

# read in metadata
metadata <- read.delim("scripts/RNA_metadata.tsv")

# set factor levels
metadata$TYPE <- factor(metadata$TYPE, levels = c("CONTROL", "PA", "AD", "LBD"))
