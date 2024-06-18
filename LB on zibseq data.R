gene.marginal.RPK.RNA = testdata[, c(1, 9:248)]
outcome.118 = as.numeric(as.data.frame(testdata[1:240,c(2,3)]))
#outcome.118 = as.numeric(testdata[1:240, c(2)]) #coloquei só as 240 linhas para dar dps match



library(tidyverse)

(gene.marginal.RPK.RNA[,-1] %>% as.matrix %>% as.vector) -> a
a[a>0] %>% sort %>% head

sampleSum = apply(gene.marginal.RPK.RNA[,-1], 2, sum)

##########LB

screen.gene <- function(data, avg.detect = 2, gene.col = NULL) {
  # data: each column = sample, each row = genes
  # gene.col: the column where gene names are provided. If null, ignored.
  gene.name = data[,gene.col]
  if (!is.null(gene.col)) data = data[,-gene.col]
  n.genes = dim(data)[1]
  n.sample = dim(data)[2]
  
  useF = which(rowSums(data) > avg.detect * n.sample)
  print(paste(length(useF), "out of ", n.genes, " will be used."))
  list(sample.size = n.sample, 
       stat = c(use = length(useF), out.of = n.genes, percentage = round(length(useF)/n.genes,2)*100),
       feature = list(index = useF, names = if (is.null(gene.col)) NA else gene.name[useF]))
}

gene.index = screen.gene(gene.marginal.RPK.RNA)

args = commandArgs(trailingOnly=TRUE)  # passed from Script_Boyang
method = args[1] %>% as.numeric  # 1..5
k = args[2] %>% as.numeric

gene.index.k = gene.index[(k-1)*5000 + (1:5000)]

gene.index.k = gene.index.k[!is.na(gene.index.k)]

j <<- 1
gene.BEZI <- sapply(gene.index.k, function(i) {
  cat(j, ".."); j <<- j + 1
  tmp <- data.frame(y = as.numeric(gene.marginal.RPK.RNA[i, -1]), 
                    outcome = outcome.118,
                    sampleSum = sampleSum)
  # Chamando LB() com tmp estruturado corretamente
  LB(tmp$outcome, tmp$y, tmp$sampleSum)
})

tmp <- data.frame(y = as.numeric(gene.marginal.RPK.RNA[, -1]), 
                  outcome = outcome.118,
                  sampleSum = sampleSum)
