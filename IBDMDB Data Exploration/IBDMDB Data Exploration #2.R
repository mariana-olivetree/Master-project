

################Analysis of dataframes from IBD Database########################

###############  HMP2 PILOT

#MTX

#genefamilies
pilot_mtx_genefamilies =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtx/genefamilies.tsv") 

zero.proportion.pilot_mtx_genefamilies =
  pilot_mtx_genefamilies %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtx_genefamilies %>% mean #0



#taxonomic profiles

pilot_mtx_taxa =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtx/taxonomic_profiles.tsv") 

zero.proportion.pilot_mtx_taxa =
  pilot_mtx_taxa %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtx_taxa %>% mean #0

#pathabundance

pilot_mtx_pa =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtx/pathabundance.tsv") 

zero.proportion.pilot_mtx_pa =
  pilot_mtx_pa %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtx_pa %>% mean #0



#pathcoverage

pilot_mtx_pc =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtx/pathcoverage.tsv") 

zero.proportion.pilot_mtx_pc =
  pilot_mtx_pc %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtx_pc %>% mean #0

#Não parece que é o HMP2 Pilot pq nenhum dos mtx dá bem



#MTG

#genefamilies

pilot_mtg_genefamilies =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtg/genefamilies.tsv") 

saveRDS(pilot_mtg_genefamilies, file = "C:/Users/olive/Desktop/Projeto PC/data_mtg_big.rds")

zero.proportion.pilot_mtg_genefamilies =
  pilot_mtg_genefamilies %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtg_genefamilies %>% mean #este é aquele de 1 milhão de samples


#taxonomic profiles

pilot_mtg_taxa =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtg/taxonomic_profiles.tsv") 

zero.proportion.pilot_mtg_taxa =
  pilot_mtg_taxa %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtg_taxa %>% mean #0



#pathabundance

pilot_mtg_pa =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtg/pathabundance.tsv") 

zero.proportion.pilot_mtg_pa =
  pilot_mtg_pa %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtg_pa %>% mean #0.003154574



#pathcoverage

pilot_mtg_pc =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2 Pilot/mtg/pathcoverage.tsv") 

zero.proportion.pilot_mtg_pc =
  pilot_mtg_pc %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.pilot_mtg_pc %>% mean #0.003153705


################## HMP2

#MTX

#ecs_3

hmp_mtx_ecs_3 =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtx/ecs_3.tsv") 

zero.proportion.mtx_ecs3 =
  hmp_mtx_ecs_3 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx_ecs3 %>% mean #0.9423719
#este tem 818 amostras como o artigo descreve

#ecs_relab

hmp_mtx_ecs_relab =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtx/ecs_relab.tsv") 

zero.proportion.mtx_ecs_relab =
  hmp_mtx_ecs_relab %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx_ecs_relab %>% mean #0

#ecs_relab_1

hmp_mtx_ecs_relab_1 =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtx/ecs_relab_1.tsv") 

zero.proportion.mtx_ecs_relab_1 =
  hmp_mtx_ecs_relab_1 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx_ecs_relab_1 %>% mean #0

#ecs_relab_2

hmp_mtx_ecs_relab_2 =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtx/ecs_relab_2.tsv") 

zero.proportion.mtx_ecs_relab_2 =
  hmp_mtx_ecs_relab_2 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx_ecs_relab_2 %>% mean #0

#ecs_relab_3

hmp_mtx_ecs_relab_3 =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtx/ecs_relab_3.tsv") 

zero.proportion.mtx_ecs_relab_3 =
  hmp_mtx_ecs_relab_3 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtx_ecs_relab_3 %>% mean #0


#MTG

#ecs_3

hmp_mtg_ecs3 =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtg/ecs_3.tsv") 

zero.proportion.mtg_ecs3 =
  hmp_mtg_ecs3 %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_ecs3 %>% mean #0.9165769

#Filtering data_mtx with External.ID of metadata_mtx_fil
intersect(external_id_fil, gsub("_.*", "", colnames(hmp_mtg_ecs3)))

mtg_fil = hmp_mtg_ecs3[, external_id_fil]

zero.proportion.mtg =
  mtg_fil %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg %>% mean #91.2%

#ecs_relab

hmp_mtg_ecs_relab =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtg/ecs_relab.tsv") 

zero.proportion.mtg_ecs_relab =
  hmp_mtg_ecs_relab %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_ecs_relab %>% mean #0.003660724

#taxonomic profile

hmp_mtg_taxa =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtg/taxa.tsv") 

zero.proportion.mtg_taxa =
  hmp_mtg_taxa %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_taxa %>% mean #0

#path abundance relab

hmp_mtg_pa_relab =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtg/pa_relab.tsv") 

zero.proportion.mtg_pa_relab =
  hmp_mtg_pa_relab %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_pa_relab %>% mean #0.004880969

#path abundance

hmp_mtg_pa =  read.delim("C:/Users/olive/Desktop/Projeto PC/HMP2/mtg/pa.tsv") 

zero.proportion.mtg_pa =
  hmp_mtg_pa %>% 
  apply(1, function(x) mean(x == 0)) 

zero.proportion.mtg_pa %>% mean #0.008541794








