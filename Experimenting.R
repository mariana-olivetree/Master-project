library(magrittr); 

############### File Upload ####################################################

# Generated taxonomic and functional profiles
#metadata_generated = read.delim("C:/Users/olive/Desktop/Projeto PC/Fase 2/Ficheiros IBD/HMP2_metadata_subset_multivariable.tsv") 
#taxa_profiles_generated = read.delim("C:/Users/olive/Desktop/Projeto PC/Fase 2/Ficheiros IBD/metaphlan_taxonomic_profiles.tsv")
#path_abundance_generated = read.delim("C:/Users/olive/Desktop/Projeto PC/Fase 2/Ficheiros IBD/pathabundance_relab.tsv")

#saveRDS(metadata_generated, file = "C:/Users/olive/Desktop/Projeto PC/ficheiros/metadata_generated.rds")
#saveRDS(taxa_profiles_generated, file = "C:/Users/olive/Desktop/Projeto PC/ficheiros/taxa_profiles_generated.rds")
#saveRDS(path_abundance_generated, file = "C:/Users/olive/Desktop/Projeto PC/ficheiros/path_abundance_generated.rds")

metadata_generated <- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/metadata_generated.rds")
taxa_profiles_generated <- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/taxa_profiles_generated.rds")
path_abundance_generated <- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/path_abundance_generated.rds")

# IBD Metadata
#metadata = read.csv("C:/Users/olive/Desktop/Projeto PC/Fase 2/Ficheiros IBD/metadata.csv")

saveRDS(metadata, file = "C:/Users/olive/Desktop/Projeto PC/ficheiros/metadata.rds")
#metadata<- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/metadata.rds")


# MTX Dataset
#data_mtx = read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtx/ecs_3.tsv")
#818 samples, same as the authors (including the column "Feature.Sample")

#saveRDS(data_mtx, file = "C:/Users/olive/Desktop/Projeto PC/ficheiros/data_mtx.rds")
data_mtx<- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/data_mtx.rds")


# MTG Dataset
#hmp_mtg_ecs3 =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtg/ecs_3.tsv") 

#saveRDS(data_mtx, file = "C:/Users/olive/Desktop/Projeto PC/ficheiros/hmp_mtg_ecs3.rds")
hmp_mtg_ecs3<- readRDS("C:/Users/olive/Desktop/Projeto PC/ficheiros/hmp_mtg_ecs3.rds")

############### DataFrame Fixing ###############################################


                               # GENERATED METADATA

colnames(metadata_generated) = gsub("_.*", "", colnames(metadata_generated))

                               # MTX COLUMNS

# Extract column names without the "_level..." part
new_colnames <- gsub("_.*", "", colnames(data_mtx))

# Assign new column names to the dataset
colnames(data_mtx) <- new_colnames

                               # MTG COLUMNS

# Extract column names without the "_level..." part
new_colnames_1 <- gsub("_.*", "", colnames(hmp_mtg_ecs3))

# Assign new column names to the dataset
colnames(hmp_mtg_ecs3) <- new_colnames_1


############### Analysis of Generated taxonomic and functional profiles ########

# Transpose of dataframe
meta_transposta <- t(metadata_generated)

#Fixing column names
colnames(meta_transposta) = c("site", "age", "diagnosis", "antibiotics", 
                              "dysbiosis")

meta_transposta = meta_transposta[2:130,]

meta_transposta = as.data.frame(meta_transposta)

table(meta_transposta$diagnosis)

#Quantity of CD and UC not as described. 
# But if 13 of CD and 12 of UC are removed --> 25 samples removed
# We get 129-25 = 104, number of patients Cho et al. mention

#Let's check if we have, at least, 25 duplicated samples
metadata_dupl = meta_transposta[2:4]
rownames(metadata_dupl) = NULL
sum(duplicated(metadata_dupl))
#Yes! This may, then, correspond to our dataset


############### Analysis of IBD Metadada #######################################

colnames(metadata)

#mtx samples
table(metadata$data_type)["metatranscriptomics"]
#835 samples

#mtg samples
table(metadata$data_type)["metagenomics"]

#Number of participants
length(unique(metadata$Participant.ID))

anyNA(metadata$Participant.ID)
#There is no Missing value

#Now, we want to know which 104 Partipants IDs were selected


                            # MTX DATATYPE

#Subdataframe of mtx datatype
metadata_mtx = subset(metadata, data_type == "metatranscriptomics")

length(unique(metadata_mtx$Participant.ID))
#109, not 132
 

                            #MTG DATATYPE

#Subdataframe of mtg datatype
metadata_mtg = subset(metadata, data_type == "metagenomics")

length(unique(metadata_mtg$Participant.ID))
#130, not 132


############### Is There a Match between MTX and IBD METADADA? #################

#Removing "Feature.Sample" from vector - MTX 
col_data_mtx = colnames(data_mtx)[-1]

# METADADA
external_id_metadata_tx = metadata_mtx$External.ID

# Match?
match_samples = intersect(col_data_mtx, external_id_metadata_tx)
#716, not 818

# Metadata with MATCHES
metadata_mtx_fil<- metadata_mtx[metadata_mtx$External.ID %in% match_samples, ]
#716, not 818

#Let's obtain the MTX Filtered according to MATCHES
external_id_fil = metadata_mtx_fil$External.ID

data_mtx_zero = data_mtx[-1]

#FINAL DATASET OF MTX -> MTX with MATCHES
data_mtx_fil = data_mtx_zero[names(data_mtx_zero) %in% external_id_fil ]

################### Checking zero proportion

#without filtering
zero.proportion.mtx =
  data_mtx_fil %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx %>% mean #94.2% 

################### Female vs Male

participant_id_meta = metadata_mtx_fil$Participant.ID

reads_raw = metadata_mtx_fil$reads_raw

filtering = data.frame(Participant.ID = participant_id_meta, reads_raw = reads_raw)

filtering

filtering_unique <- filtering[!duplicated(filtering$Participant.ID), ]

reads_raw_fil = filtering_unique$reads_raw

metadata_fem_mas = metadata_mtx_fil[metadata_mtx_fil$reads_raw %in% reads_raw_fil, ]

table(metadata_fem_mas$sex)
#52 female + 52 male

############### PARTICIPANT ID FROM MTX ########################################

p.id = metadata_fem_mas$Participant.ID

metadata_pid = metadata[metadata$Participant.ID %in% p.id,]

metadata_pid_mtx = subset(metadata_pid, data_type == "metatranscriptomics")

metadata_pid_mtx$External.ID = gsub("_.*", "", metadata_pid_mtx$External.ID)

data_mtx_fil_2 = data_mtx[names(data_mtx_zero) %in% metadata_pid_mtx$External.ID ]

zero.proportion.mtx_fil_2 = 
  data_mtx_fil_2 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx_fil_2 %>% mean #94.2%


############### Is There a Match between MTG and IBD METADADA? #################

#Subdataframe of mtx datatype
metadata_mtg = subset(metadata, data_type == "metagenomics")

length(unique(metadata_mtg$Participant.ID))
#130, not 132

#Cleaning of external.id column
metadata_mtg$External.ID = gsub("_.*", "", metadata_mtg$External.ID)

#Column names from mtg dataset
colnames_mtg = colnames(hmp_mtg_ecs3)[-1]

external_id_metadata_mtg = metadata_mtg$External.ID

match_samples_mtg = intersect(colnames_mtg, external_id_metadata_mtg)
#1593, close to 1595 like the authors stated

metadata_mtg_fil = metadata_mtg[metadata_mtg$External.ID %in% match_samples_mtg, ]

length(unique(metadata_mtg_fil$Participant.ID))

length(unique(metadata_mtg_fil$External.ID)) #there are 1593 samples still

mtg_dataset_fil = hmp_mtg_ecs3[names(hmp_mtg_ecs3) %in% match_samples_mtg]


zero.proportion.mtg_fil_ =
  mtg_dataset_fil %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_fil_ %>% mean #91.7%

############### PARTICIPANT ID FROM MTX ########################################

p.id = metadata_fem_mas$Participant.ID

metadata_pid = metadata[metadata$Participant.ID %in% p.id,]

metadata_pid_mtg = subset(metadata_pid, data_type == "metagenomics")

metadata_pid_mtg$External.ID = gsub("_.*", "", metadata_pid_mtg$External.ID)

# THIS IS looks like the FINAL MTG DATASET !!!!!!!!!!!!!!!!! 1585 metagenomics
data_mtg_fil_2 = hmp_mtg_ecs3[names(hmp_mtg_ecs3) %in% metadata_pid_mtg$External.ID ]


zero.proportion.mtg_fil_2 = 
  data_mtg_fil_2 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_fil_2 %>% mean #91.7%

############## Is There a Match between MTX/MTG and GENERATED METADADA? ########


colnames_metadata_gen = colnames(metadata_generated)[-1]

                              # MTX Dataset

intersect(names(data_mtx), colnames_metadata_gen)
#Very few matches, something is wrong. Taxa profiles may come from mtg samples.

                              # MTG Dataset

colnames_mtg = colnames(hmp_mtg_ecs3)[-1]

length(intersect(colnames_metadata_gen, colnames_mtg))
#124 matches, should be 130. 



################# Filtering IBD METADATA visit number = 1 ######################

metadata_filtrado <- subset(metadata, visit_num == 1)

length(unique(metadata_filtrado$Participant.ID))

sample_names = metadata_filtrado$External.ID

as.vector(sample_names)

# Analyzing correspondence with data_mtx
colunas_corr = intersect(names(data_mtx), sample_names)
#No matching








