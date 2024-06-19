# Differential expression analysis
# By João Sequeira
# Sep 2017

library(readxl)

setwd("C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs")

# Upload count data and metadata
counts_filtrado = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\Data ecs\\counts_filtrado.rds")

metadata_filtrado_2 = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\Data ecs\\metadata_numvariable.rds")

# Renaming conditions
conditions <- metadata_filtrado_2$health

# Rename counts
total <- counts_filtrado

#total <- read.table("simulated_de_input.txt", h=T, row.names=1, sep = '\t')
foldchange <- 2
output <- "."

# RNA-Seq differential analysis
for (package in c("DESeq2", "pheatmap", "RColorBrewer")){
  eval(bquote(library(.(package))))
}
conditions <- factor(conditions)
total[is.na(total)] <- 0
total <- total[ rowSums(total) > 1, ]
cd <- data.frame(conditions)
colnames(cd)[1] <- "condition"
rownames(cd) <- colnames(total)

total <- round(total)

dds <- DESeqDataSetFromMatrix(countData = total, colData = cd, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, lfcThreshold = log2(foldchange))
data <- counts(estimateSizeFactors(dds), normalized=TRUE)
write.table(
  data, file=paste0(file=output, "/normalized_counts_fold2.tsv"), sep='\t', col.names = NA, quote=FALSE)

# Saving files
saveRDS(dds, file = "C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\dds_fold2.rds")
saveRDS(res, file = "C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\res_fold2.rds")


# Bland–Altman plot
jpeg(paste0(output, "/ma.jpeg"))
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

# Normalized counts
jpeg(paste0(output, "/counts.jpeg"))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

# Protein expressions differential analysis
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), file=paste0(
  output, "/condition_treated_results.tsv"), sep='\t', col.names = NA, quote=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
select <- rownames(head(resOrdered,30))
vsd.counts <- assay(vsd)[select,]
jpeg(paste0(output, "/gene_expression.jpeg"))
pheatmap(vsd.counts)
dev.off()

#tsv = read_tsv("C:\\Users\\olive\\Downloads\\condition_treated_results.tsv")

#install.packages("br")
#library(br)
# Sample expressions differential analysis
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(total)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(paste0(output, "/sample_distances.jpeg"))
pheatmap(
  sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
dev.off()

# Principal components analysis
jpeg(paste0(output, "/pca.jpeg"))
plotPCA(vsd, intgroup = c("condition"))
dev.off()

print("DE analysis finished.")