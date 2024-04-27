
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

#group -> N, O, OM, OB, N

#Converts the group vector to numeric format
group = as.numeric(group)

#Sets all values in group that are less than 4 (N, O and OB) to 0
group[which(group<4)]=0

#Sets all values in gr that are equal to 4 (OM) to 1
group[which(group==4)]=1

#Vector becomes binary
class(group)

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

##########Perfoming DGE with another group#################

#Vector with class labels for the samples
group_2 = testdata[,2]

#Converts the group vector to numeric format
group_2 = as.numeric(group_2)

#Selecting OB (Obese) -> number 3

group_2[which(group_2<3)]=0

group_2[which(group_2==4)]=0

group_2[which(group_2==3)] = 1

class(group_2)

#Perfoms DGE
result_OB = ZIBseq(data = data_mt, outcome = group_2, transform = TRUE)

result_OB$qvalues

result_OB$pvalues

result_OB$sigFeature

result_OB$useFeature

##########Perfoming DGE with 2 groups#################

#Vector with class labels for the samples
group_3 = testdata[,2]

#Converts the group vector to numeric format
group_3 = as.numeric(group_3)

group_3

#Selecting OB (Obese) -> number 3

group_3[which(group_3<3)] = 0

group_3[which(group_3==4)] = 1

group_3[which(group_3==3)] = 1

group_3

class(group_3)

#Perfoms DGE
result_3 = ZIBseq(data = data_mt, outcome = group_3)

result_3$qvalues

result_3$pvalues

result_3$sigFeature

result_3$useFeature

##########Perfoming with categorical groups#################

#Vector with class labels for the samples
group_4 = testdata[,2]

group_4

#Perfoms DGE
result_4 = ZIBseq(data = data_mt, outcome = group_4)

result_4$qvalues

result_4$pvalues

result_4$sigFeature

result_4$useFeature

############Using another dataset###########################



