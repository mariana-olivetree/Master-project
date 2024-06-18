library(gamlss)

#count data
counts = seqdata

#creating gene indexes
genes.index = screen.gene(seqdata, gene.col = NULL)$feature$index

#Summing the samples
sampleSum = apply(seqdata, 2, sum)

#removal of NAs
genes.index = genes.index[!is.na(genes.index)]

j <<- 1; cat("j=")

gene.BEZI <- sapply(genes.index, function(i) {
  cat(j,".." ); j <<- j+1
  tmp <- data.frame(y = seqdata[i,] %>% as.numeric,
                    sampleSum = sampleSum)
  print(dim(tmp))
  LB_nobatch(tmp) 
})

#Transpose matrix
genes = t(gene.BEZI)

#Turn into dataframe
genes = as.data.frame(genes)

#Get significant genes
sig_genes = genes[genes$pval < 0.05, ]
