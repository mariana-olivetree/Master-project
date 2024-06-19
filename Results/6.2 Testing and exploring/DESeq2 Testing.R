#Load the package
library("ZIBseq")
library(DESeq2)

#Using ZIBseq testdata
data(testdata)

#Metatranscriptomic data
data_mt = testdata[,9:248]

#Number of columns of the metatranscriptomic data
columns = dim(data_mt)[2]

#Converts the data to numeric format
for (col in 1:columns){data_mt[,col]=as.numeric(as.character(data_mt[,col]))}

# Calculate the average reads for each column
average_reads = colMeans(data_mt)

# Find the columns with average reads below 2
columns_to_keep = average_reads >= 2

# Subset the dataframe to keep only the columns with average reads >= 2
data_mt <- data_mt[, columns_to_keep]

#Deseq2 must use this form of countData
data_mt = t(data_mt)

#Meta data
meta_data = testdata[,1:8]

#Preparing dataframe to compare (design)
group = testdata[,2]

#Converts the group vector to numeric format
group = as.numeric(group)

group

#Sets all values in group that are less than 4 (N, O and OB) to 0
group[which(group<4)]="other"

#Sets all values in gr that are equal to 4 (OM) to 1
group[which(group==4)]="OM"

#Vector becomes binary
factor(group)

#Making metadata correspond to the design
meta_data$group = ifelse(meta_data$group == "OM", "OM", "other")

## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data_mt, colData = meta_data , design = ~ group)

## Run analysis
dds <- DESeq(dds)

#Obtain results
res = results(dds, contrast = c("group", "OM", "other"),alpha = 0.05)

#Quantity of differentially expressed taxa according to pvalues
sum(res$padj < 0.05, na.rm=TRUE)

#calculate q value
library(qvalue)

q_value = qvalue(res$pvalue)

q_value$qvalues

#Quantity of significantly expressed taxa according to qvalues
sum(q_value$qvalues < 0.05, na.rm=TRUE)

row_names <- rownames(data_mt)

row_names

named_q_values <- setNames(q_value$qvalues, row_names)

named_q_values

# Assuming named_q_values contains the named vector of q-values

# Select names where q-values are less than 0.05
significant_names <- names(named_q_values[named_q_values < 0.05])

# Now, significant_names contains the names for which the q-values are less than 0.05

significant_names
