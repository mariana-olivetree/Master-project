

#LB ecs


siggenes_ecs_lb = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\Results\\Data ecs\\sig_genes_ecs.rds")


library(gridExtra)
library(grid)
d <- head(siggenes_ecs_lb[c(3:5,7),])
grid.table(d)


#DESeq2 ecs

dds_ecs = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\dds.rds")
res_ecs = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\res.rds")

p_values_ecs_deseq = res_ecs$pvalue

genes_ecs <- rownames(res_ecs)

results_pval_ecs <- data.frame(Gene = genes_ecs, P_Value = p_values_ecs_deseq)

sig_ecs_deseq <- results_pval_ecs[results_pval_ecs$P_Value < 0.05, ]

#Saving files
saveRDS(sig_ecs_deseq, "C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\sig_genes.rds" )

siggenes_ecs_deseq = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\sig_genes.rds" )



# Comparison of tools


siggenes_ecs_deseq

no_na_siggenes_ecs_deseq = na.omit(siggenes_ecs_deseq)


matches_ecs = intersect(rownames(siggenes_ecs_lb), no_na_siggenes_ecs_deseq$Gene)


# Instalar o pacote VennDiagram se ainda não estiver instalado
install.packages("VennDiagram")

# Carregar o pacote VennDiagram
library(VennDiagram)

# Conjuntos de genes
genes_lb <- rownames(siggenes_ecs_lb)
genes_deseq <- no_na_siggenes_ecs_deseq$Gene

# Encontrar genes comuns
matches_ecs <- intersect(genes_lb, genes_deseq)

library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")

# Criar o diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    "LB test" = genes_lb,
    "DESeq2" = genes_deseq
  ),
  filename = NULL,
  lwd = 2,
  lty = 'blank',
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = 0,
  cat.dist = 0.02,
  # Numbers
  fontface = "bold",
  fontfamily = "serif",
)

# Mostrar o diagrama de Venn
grid.draw(venn.plot)



#fold change 2

genes_ecs_2_deseq = readRDS("C:\\Users\\olive\\Desktop\\Projeto\\DESeq MOSCA\\ecs\\res_fold2.rds")


p_values_ecs_2_deseq = genes_ecs_2_deseq$pvalue


genes_ecs_2 <- rownames(genes_ecs_2_deseq)

results_pval_ecs_2 <- data.frame(Gene = genes_ecs_2, P_Value = p_values_ecs_2_deseq)

sig_ecs_2_deseq <- results_pval_ecs_2[results_pval_ecs_2$P_Value < 0.05, ]

no_na_siggenes_2_ecs_deseq = na.omit(sig_ecs_2_deseq)

# Conjuntos de genes
genes_2_deseq <- no_na_siggenes_2_ecs_deseq$Gene

# Encontrar genes comuns
matches_ecs_fold_2 <- intersect(genes_lb, genes_2_deseq)


# Conjuntos de genes
genes_lb <- rownames(siggenes_ecs_lb)
genes_deseq <- no_na_siggenes_ecs_deseq$Gene


# Criar o diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    "LB test" = genes_lb,
    "DESeq2" = genes_2_deseq
  ),
  filename = NULL,
  lwd = 2,
  lty = 'blank',
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = 0,
  cat.dist = 0.04,
  # Numbers
  fontface = "bold",
  fontfamily = "serif",
)

# Mostrar o diagrama de Venn
grid.draw(venn.plot)





