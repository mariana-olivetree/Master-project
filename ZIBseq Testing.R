
#Instalation of package
install.packages("C:/Users/olive/Desktop/ZIBseq_1.2/ZIBseq", repos = NULL, type = "source")

#Load the package
library("ZIBseq")

#Testing ZIBseq
data(testdata)

#Metatranscriptomic data
data_mt = testdata[,9:248]

#Number of columns of the metatranscriptomic data
columns = dim(data_mt)[2]

#Converts the data to numeric format
for (col in 1:columns){data_mt[,col]=as.numeric(as.character(data_mt[,col]))}

#Vector with class labels for the samples
group = testdata[,2]

#Converts the group vector to numeric format
group = as.numeric(group)

#Create a mapping dictionary
group_mapping = c("N" = 1, "O" = 2, "OM" = 3, "OB" = 4)

#Get original groups corresponding to value 1
original_groups_1 <- names(group_mapping[group_mapping == 1])
original_groups_2 <- names(group_mapping[group_mapping == 2])
original_groups_3 <- names(group_mapping[group_mapping == 3])
original_groups_4 <- names(group_mapping[group_mapping == 4])

#Display the original groups for value 1
cat("Mapped value of 1 corresponds to original groups:", original_groups_1, "\n")
cat("Mapped value of 2 corresponds to original groups:", original_groups_2, "\n")
cat("Mapped value of 3 corresponds to original groups:", original_groups_3, "\n")
cat("Mapped value of 4 corresponds to original groups:", original_groups_4, "\n")

#Sets all values in group that are less than 4 (N, O and OM) to 0
group[which(group<4)]=0

#Sets all values in gr that are equal to 4 (OB) to 1
group[which(group==4)]=1

#Vector becomes binary
group

#Perfoms DGE
result = ZIBseq(data = data_mt, outcome = group)

result$qvalues

result$pvalues

result$sigFeature

result$useFeature

#Performs DGE with the square-root transformation of the data
result_transformed = ZIBseq(data = data_mt, outcome = group, transform = TRUE)

result_transformed$sigFeature

result_transformed$useFeature

result_transformed$qvalues

result_transformed$pvalues

