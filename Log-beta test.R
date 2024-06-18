library(gamlss)

#count data
counts = seqdata

#metadata
metadata = metadata

#Turning the variable of interest into a numerical variable
metadata$gender_num = factor(ifelse(metadata$gender== "female",0,1))

#creating gene indexes
genes.index = screen.gene(seqdata, gene.col = NULL)$feature$index

#creating gene list
#gene.list = counts[genes.index,1]

#Summing the samples
sampleSum = apply(seqdata, 2, sum)

#removal of NAs
genes.index = genes.index[!is.na(genes.index)]

sample_id = metadata$barcode

variable =  metadata$gender_num

j <<- 1; cat("j=")

gene.BEZI <- sapply(genes.index, function(i) {
  cat(j,".." ); j <<- j+1
  tmp <- data.frame(y = seqdata[i,] %>% as.numeric, id = sample_id, 
                    phenotype = variable,
                    sampleSum = sampleSum)
  print(dim(tmp))
  LB(tmp) 
})

#Transpose matrix
genes = t(gene.BEZI)

#Turn into dataframe
genes = as.data.frame(genes)

#Get significant genes
sig_genes = genes[genes$pval < 0.05, ]



