source("C:\\Users\\olive\\Desktop\\Projeto\\LB functions.R")
library(tidyverse)
library(readr)

metadata<- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/metadata.rds")

data <-readr::read_tsv("C:\\Users\\olive\\Downloads\\ecs_3 (1).tsv\\ecs_3 (1).tsv")

counts = data

counts = as.data.frame(counts)

############### Fixing dataframe

rownames(counts) = counts[,1]

counts = counts[,-1]

column_names <- colnames(counts)

col_names <- gsub("_level4ec", "", column_names)

# Atualizar os nomes das colunas no dataframe
colnames(counts) <- col_names


############### Filtering metadata

column_check <- col_names %in% metadata$External.ID

print(column_check)

filtered_metadata <- metadata[metadata$External.ID %in% col_names, ]

metadata_filtrado_2 <- filtered_metadata[filtered_metadata$data_type == "metatranscriptomics", ]

# Obter os IDs externos filtrados
filtered_external_ids <- metadata_filtrado_2$External.ID

# Filtrar o dataframe counts para manter apenas as colunas correspondentes aos IDs externos filtrados
counts_filtrado <- counts[, colnames(counts) %in% filtered_external_ids]

#Saving file for posterior analyses
#saveRDS(counts_filtrado, "C:\\Users\\olive\\Desktop\\Projeto\\Data ecs\\counts_filtrado.rds")


###############


#Turning the variable of interest into a numerical variable
metadata_filtrado_2$health = factor(ifelse(metadata_filtrado_2$diagnosis== "nonIBD",0,1))

#Saving file for posterior analyses
#saveRDS(metadata_filtrado_2, "C:\\Users\\olive\\Desktop\\Projeto\\Data ecs\\metadata_numvariable.rds")

#creating gene indexes
genes.index = screen.gene(counts_filtrado, gene.col = NULL)$feature$index

#Summing the samples
#sampleSum <- apply(counts_filtrado, 2, sum)

sampleSum_df <- counts_filtrado %>% 
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))

sampleSum <- unlist(sampleSum_df)


genes.index = genes.index[!is.na(genes.index)]

sample_id = metadata_filtrado_2$External.ID

variable =  metadata_filtrado_2$health

j <<- 1; cat("j=")

gene.BEZI <- sapply(genes.index, function(i) {
  cat(j,".." ); j <<- j+1
  tmp <- data.frame(y = counts_filtrado[i,] %>% as.numeric, id = sample_id, 
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

#saving files
saveRDS(genes, "C:\\Users\\olive\\Desktop\\Projeto\\Data ecs\\dge_genes_ecs.rds")
saveRDS(sig_genes, "C:\\Users\\olive\\Desktop\\Projeto\\Data ecs\\sig_genes_ecs.rds")




