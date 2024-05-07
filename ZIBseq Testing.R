
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

group = testdata[,2]

group
##########Using categorical values#########################


# Recode the "group" variable
group_1 = ifelse(group == "OM", "OM", "N_OW_OB")

# Run the ZIBseq analysis
result = ZIBseq(data = data_mt, outcome = group_1)

result$qvalues

result$pvalues

result$sigFeature

result$useFeature


##########Using categorical + sqrt#######################


# Recode the "group" variable
group_sqrt = ifelse(group == "OM", "OM", "N_OW_OB")

# Run the ZIBseq analysis
result_sqrt = ZIBseq(data = data_mt, outcome = group_sqrt, transform = TRUE)

result_sqrt$qvalues

result_sqrt$pvalues

result_sqrt$sigFeature

result_sqrt$useFeature


##########Using ordinal#################################


#Converts the group vector to numeric format
group_num = as.numeric(group)

group_num

#Sets all values in group that are less than 4 (N, O and OB) to 0
group_num[which(group_num<4)]=0

#Sets all values in gr that are equal to 4 (OM) to 1
group_num[which(group_num==4)]=1

#Vector becomes binary
group_num

#Perfoms DGE
result_num = ZIBseq(data = data_mt, outcome = group_num)

result_num$qvalues

result_num$pvalues

result_num$sigFeature

result_num$useFeature


##########Using ordinal + sqrt##########################

#Perfoms DGE
result_num_sqrt_2 = ZIBseq(data = data_mt, outcome = group_num, transform = TRUE)

result_num_sqrt_2$qvalues

result_num_sqrt_2$pvalues

sum(resut_num_sqrt$padj < 0.05)

result_num_sqrt_2$sigFeature

result_num_sqrt$useFeature

result_num_sqrt$qvalues

result_num_sqrt$pvalues

####################################

#Converts the group vector to numeric format
# Define as categorias do grupo BMI na ordem desejada (N, O, OM, OB)
categories <- c("N", "O", "OB", "OM")

# Converte o grupo BMI para fatores ordenados
group_ord <- factor(group, levels = categories, ordered = TRUE)

# Executa o DGE considerando o grupo BMI como variável ordinal
result_ord <- ZIBseq(data = data_mt, outcome = group_ord)

result_ord

##################################


#Other TESTING

test_data2 = testdata

test_data2$Acidovorax[test_data2$group == "N"] <- 400

#Number of columns of the metatranscriptomic data
columns = dim(test_data2)[2]

#Converts the data to numeric format
for (col in 1:columns){test_data2[,col]=as.numeric(as.character(test_data2[,col]))}

group_acido = testdata[,2]

group_acido = as.numeric(group_acido)

group_acido

results_acido = ZIBseq(data = test_data2, outcome = group_acido )

results_acido$sigFeature

results_acido = ZIBseq(data = test_data2, outcome = group_acido, transform = TRUE )

results_acido

####################3

test_data3 = testdata

test_data3$Ahrensia[test_data3$group == "OM"] <- 10

#Number of columns of the metatranscriptomic data
columns = dim(test_data3)[2]

#Converts the data to numeric format
for (col in 1:columns){test_data3[,col]=as.numeric(as.character(test_data3[,col]))}

group_acido = testdata[,2]

group_acido = as.numeric(group_acido)

group_acido

results_acido = ZIBseq(data = test_data3, outcome = group_acido )

results_acido$sigFeature

results_acido = ZIBseq(data = test_data2, outcome = group_acido, transform = TRUE )

results_acido


bmi = testdata[,3]

bmi = as.numeric(bmi)




results_ordinal = ZIBseq(data = data_mt, outcome = bmi)

results_ordinal
