pathToRef = c("/LBD_CWOW/references/human/")

# read in annotation file
gtf.file <- paste0(pathToRef, "gencode.v38.annotation.gtf")
genes.gtf <- rtracklayer::import(gtf.file)
genes.gtf <- as.data.frame(genes.gtf)
genes.gtf <- genes.gtf[genes.gtf$type == "gene",]
