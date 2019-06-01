# Author: Benjamin Frot
rm(list=ls())

transcripts <- read.csv("./gencode.v19.genes.v7.patched_contigs.gtf", header=F, sep="\t", 
                        nrows=-1, skip=6)
transcripts <- transcripts[transcripts$V3 == 'gene',]
transcripts$V9 <- as.character(transcripts$V9)
# Keep only protein coding genes
idx <- grep("protein_coding", transcripts$V9)
transcripts <- transcripts[idx, ]
# Get their names
valid_genes <- sapply(strsplit(sapply(strsplit(transcripts$V9,";"), function(x) x[1]), " "), function(x) x[2])
gene_names <- sapply(strsplit(sapply(strsplit(transcripts$V9,";"), function(x) x[5]), " "), function(x) x[3])

library(hash)
ens_to_symbol <- hash()
for(i in 1:length(valid_genes)) {
  ens_to_symbol[[valid_genes[i]]] <- gene_names[i]
}

# Now go through all tissues and see whether they have enough samples
fns <- list.files("./GTEx_Analysis_v7_eQTL_covariates/", full.names = T)
datasets <- list()
for(fn in fns) {
  C <- read.csv(fn, header=T, sep="\t")
  if(dim(C)[2] < 300) {
    next()
  }
  tissue_name <- strsplit(strsplit(fn, ".", fixed = T)[[1]][2], "/")[[1]][4]
  d_fn <- paste("./GTEx_Analysis_v7_eQTL_expression_matrices/", tissue_name, ".v7.normalized_expression.bed.gz", sep="")
  all_data <- read.csv(d_fn, header = T, sep="\t")
  
  to_keep <- which(all_data$gene_id %in% valid_genes)
  all_data <- all_data[to_keep,]
  
  # Get rid of X
  all_data <- all_data[all_data$X.chr != 'X', ]
  
  expr <- (all_data[,5:ncol(all_data)])
  rownames(expr) <- all_data$gene_id
  chrs <- all_data$X.chr
  expr <- t(expr)
  
  covariates <- C[,2:ncol(C)]
  l <- lm(expr ~ t(covariates))
  expr2 <- l$residuals
  
  L <- list()
  L$covariates <- C
  L$expr.normal <- expr
  L$expr.unconfounded <- expr2
  L$chromosomes <- chrs
  datasets[[tissue_name]] <- L
  print(tissue_name)
  print(dim(expr))
  print(dim(expr2))
  print(dim(C))
  print(dim(covariates))
}

save(file = "all_tissues.preprocessed.data.RData.gz", compress = T, ens_to_symbol, datasets)

# Save the list of ensembl and gene symbols
e_to_g <- matrix(NA, nrow=0, ncol=2)
for (k in keys(ens_to_symbol)) {
  e_to_g <- rbind(e_to_g, c(k , ens_to_symbol[[k]]))
}
write.csv(file="gene_name_mapping.csv", e_to_g)
