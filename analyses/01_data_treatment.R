# Script parameters ----

## Library ----

library("data.table")
library("dplyr")
library("ggplot2")

## List of rodents identified putative pathogens ----
#Beware, using this list should replace prior exploration of the cluster present in the metabarcoding dataset
putative_pathos_genra <- c("Bartonella", "Neoehrlichia", "Mycoplasma", "Borrelia", "Borreliella", "Streptobacillus", "Ehrlichia", "Spiroplasma","Brevinema",
                           "Orientia", "Rickettsia", "Leptospira", "Yersinia", "Actinobacillus", "Treponema", "Chlamydia", "Neisseria", "Pasteurella",
                           "Francisella", "Brucella", "Coxiella")
putative_pathos_family <- c("Sarcocystidae")


# Import data ----

#Import hosts and line modalities file
d_host <- readr::read_csv(here::here("data/", "raw-data", "host_data", "20250910_SM_beprep.csv") )

# Import 16S filtered file
file16s <- data.table::fread(file = here::here( "data", "raw-data","16s","2025.09.08-16S.mergedfinal.csv"))

# Import rpoB filtered file
filerpoB_run01_05 <- data.table::fread(file = here::here( "data","raw-data","rpoB_run01-05","Run01-05-rpoB_filtered-merged.txt"))

# Import lipl32 file
filelipl32 <- readxl::read_excel(path = here::here( "data","raw-data", "lepto","20250116_data_qPCR_lepto_BePrep.xlsx"))

# Import helminths file
file_helm <- readxl::read_excel(here::here("data", "raw-data", "helminths", "20250317_fiche_dissection.xlsx"),  col_types = c("text", "numeric", "numeric", 
              "numeric", "numeric", "numeric", 
              "numeric", "numeric", "skip", "skip"))

# Import rodent macroparasite
d_macroparasite <- data.table::fread(file = here::here("data", "derived-data", "ticks", "rodents_tick", "20240731_macroparasite.csv") )


# Metabarcoding data treatment ----

## Generate files containing samples of interest ----

#Extract id of small mammals caught and dissected in Beprep
d_host <- d_host %>% 
  filter(!is.na(numero_centre))%>%
  filter(stringr::str_detect(numero_centre, pattern = "NCHA"))

sm_id <- unique(d_host %>% 
                  pull(numero_centre))

#Identify taxonomy columns of 16S and rpoB files 
taxo_name_16s <- colnames(file16s [, 1:which(colnames(file16s) == "observation_sum") ])

taxo_name_rpoB <- colnames(filerpoB_run01_05 [, 1:which(colnames(filerpoB_run01_05) == "observation_sum") ])

#Calculate total run read number per cluster 
#(considering we added filtering steps prior import of data)
file16s <- file16s %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name_16s))) %>%
  relocate(clean_totalrunreads, .after = names(file16s)[which(colnames(file16s) == "observation_sum")])

filerpoB_run01_05 <- filerpoB_run01_05 %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name_rpoB))) %>%
  relocate(clean_totalrunreads, .after = names(filerpoB_run01_05)[which(colnames(filerpoB_run01_05) == "observation_sum")])

#Add new column name to names of taxonomy columns
taxo_name_16s <- c(taxo_name_16s, "clean_totalrunreads")

taxo_name_rpoB <- c(taxo_name_rpoB, "clean_totalrunreads")

# #Take away problematic samples - if needed
# file16s <- file16s |>

#filerpoB_run01_05 <- filerpoB_run01_05 |>


#Generate new file containing only Spleen BePrep's samples

#by selecting those samples (for 16S)
beprep_sample_16s_sp <- file16s |>
  select( (contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
beprep_16s_sp <- file16s |>
  select(taxo_name_16s) |>
  cbind(beprep_sample_16s_sp)

#by selecting those samples (for rpoB)
beprep_sample_rpoB_sp <- filerpoB_run01_05 |>
  select( (contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
beprep_rpoB_sp <- filerpoB_run01_05 |>
  select(taxo_name_rpoB) |>
  cbind(beprep_sample_rpoB_sp)

#Delete the PCZymo positive in 16S after visual control
beprep_16s_sp <- beprep_16s_sp |>
  select(!contains("PCzymo"))

#Verify absence of PCZymo positive in rpoB data
beprep_rpoB_sp |>
  select(contains("PCzymo"))

#Calculate total number of reads for each cluster considering chosen samples
beprep_16s_sp <- beprep_16s_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name_16s)))%>%
  relocate(totalreads, .after = names(beprep_16s_sp)[which(colnames(file16s) == "clean_totalrunreads")])

beprep_rpoB_sp <- beprep_rpoB_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name_rpoB)))%>%
  relocate(totalreads, .after = names(beprep_rpoB_sp)[which(colnames(filerpoB_run01_05) == "clean_totalrunreads")])

#Order rows by new total read number
beprep_16s_sp <- data.table::setorder(beprep_16s_sp, -totalreads)

beprep_rpoB_sp <- data.table::setorder(beprep_rpoB_sp, -totalreads)


#Delete cluster with no reads for the selected individuals
beprep_16s_sp <- beprep_16s_sp |> 
  filter(totalreads >0)

beprep_rpoB_sp <- beprep_rpoB_sp %>% 
  filter(totalreads >0)

#Write a file containing Spleen OTUs containing our samples only :
data.table::fwrite(beprep_16s_sp, here::here("data", "derived-data","16S", "20250910_beprep_16s_sp.txt") )

data.table::fwrite(beprep_rpoB_sp, here::here("data", "derived-data","rpoB", "20250123_beprep_rpoB_sp.txt") )


## Split taxonomy and filter affiliation failures ----

#Split the taxonomy to access taxon
pathos_16s <- beprep_16s_sp |>
  tidyr::separate_wider_delim(cols = blast_taxonomy,
                              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                              delim = ";")

pathos_rpoB <- beprep_rpoB_sp |>
  tidyr::separate_wider_delim(cols = blast_taxonomy,
                              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                              delim = ";")

#Delete cluster without real affiliation
pathos_16s <- pathos_16s %>%
  filter(!(domain == "no" | domain == "Plantae"))

pathos_rpoB <- pathos_rpoB %>%
  filter(!domain == "no")


## Generate putative pathogens data ----

#For 16S Keep only OTUs with known putative pathogens 
#(you should look at raw data before this step, ensure it is ok and you're not missing pathogens taxa as the list may not be complete)
pathos_16s <- pathos_16s |>
  filter(genus %in% putative_pathos_genra 
         | family %in% putative_pathos_family)

#No need to apply such filter for rpoB, marker is specific enough to only amplify pathogenic species


#Write files containing Spleen OTUs with BePrep samples only AND only putative pathogens :
data.table::fwrite(pathos_16s, here::here("data/", "derived-data/","16S", "20250123_beprep_16s_sp_pathos.txt") )

data.table::fwrite(pathos_rpoB, here::here("data", "derived-data","rpoB", "20250123_beprep_rpoB_sp_pathos.txt") )

## External manipulation of pathogen data ----

#Seek for species level identification for specific taxa (e.g. Mycoplasma) throught external manipulation (e.g phylogenic tree and Genbank)
#A new datafile is created BY HAND to add the checked species information (in new column "species_adjusted") and then imported
#If probable pseudo-gene affiliation is observed in data file they may also get deleted ny operator at this step
pathos_16s_checked <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01-04-05","20241204_beprep_16s_sp_pathos-adjusted.txt"))
pathos_16s_checked$species_adjusted[pathos_16s_checked$species_adjusted == ""] <- NA

pathos_rpoB_checked <- data.table::fread(file = here::here( "data","raw-data/","rpoB_run01-05","20241204_beprep_rpoB_sp_pathos-adjusted.txt"))
pathos_rpoB_checked$species_adjusted[pathos_rpoB_checked$species_adjusted == ""] <- NA


# Change species value when there is a species adjusted value (generated by operator manually) and drop the species adjusted column
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
  select(!species_adjusted)

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
  select(!species_adjusted)


## Additionnal filtering and correction ----
#Rename inappropriate affiliation 
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(species = if_else(species == "Haemobartonella muris", "Mycoplasma haemomuris", species))

#FOR rpoB : reads based filter
pathos_rpoB_checked <- pathos_rpoB_checked %>%
  filter(totalreads >= 100)


## Best identified taxa column ----

#If there are unknown or Multi-affiliation in species column, change it to NA
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species ))

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species )) 

#Creation of a new column containing best available identity (between species or genus here)
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(
    identification16s = if_else(!is.na(species), species, genus),
    best_identification = if_else(!is.na(species), "species", "genus")
  ) %>%
  relocate(identification16s, .after = species) %>%
  relocate(best_identification, .after = identification16s)

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(
    identificationrpoB = if_else(!is.na(species), species, genus),
    best_identificationrpoB = if_else(!is.na(species), "species", "genus")
  ) %>%
  relocate(identificationrpoB, .after = species) %>%
  relocate(best_identificationrpoB, .after = identificationrpoB)

#Replace space in taxa name 
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(identification16s = gsub(" ", "_", identification16s ))

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(identificationrpoB = gsub(" ", "_", identificationrpoB ))


## Compress metabarcoding affiliation ----
#by best identified taxa

#Keep only best taxa and samples
taxa_pathos_16s <- pathos_16s_checked %>%
  select(c("identification16s", names(pathos_16s_checked)[(which(names(pathos_16s_checked) == "totalreads") + 1):ncol(pathos_16s_checked)]))

taxa_pathos_rpoB <- pathos_rpoB_checked %>%
  select(c("identificationrpoB", names(beprep_rpoB_sp)[(which(names(beprep_rpoB_sp) == "totalreads") + 1):ncol(beprep_rpoB_sp)]))


#Sum rows by best taxa
taxa_pathos_16s <- taxa_pathos_16s %>%
  group_by(identification16s) %>%
  summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop")

taxa_pathos_rpoB <- taxa_pathos_rpoB %>%
  group_by(identificationrpoB) %>%
  summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop")

#Transpose the table
taxa_transposed_16S <- taxa_pathos_16s |>
  data.table::transpose(
    keep.names = "numero_centre_16s",
    make.names = 1
  )

taxa_transposed_rpoB <- taxa_pathos_rpoB |>
  data.table::transpose(
    keep.names = "numero_centre_rpoB",
    make.names = 1
  )

#Finish sample names transformation to join to host data
taxa_transposed_16S <- taxa_transposed_16S |>
  mutate(numero_centre_16s = gsub(".SP", "", numero_centre_16s ))

taxa_transposed_rpoB <- taxa_transposed_rpoB |>
  mutate(numero_centre_rpoB = gsub(".SP", "", numero_centre_rpoB ))


## Join 16S and rpoB  data ----
taxa_transposed_16S_rpoB <- left_join(taxa_transposed_16S, 
                                      taxa_transposed_rpoB,
                                      by = c("numero_centre_16s" = "numero_centre_rpoB")) %>%
  rename(numero_centre_combined = numero_centre_16s)

#Replace NA by 0 in new data
taxa_transposed_16S_rpoB <- taxa_transposed_16S_rpoB %>%
  mutate(across(-numero_centre_combined, ~ if_else(is.na(.), 0, .)))

#See benefit of Neoehrlichia_mikurensis columns combining
taxa_transposed_16S_rpoB %>%
  filter(Neoehrlichia_mikurensis.y > 0) %>%
  filter(Neoehrlichia_mikurensis.x == 0) %>%
  select(numero_centre_combined) %>%
  pull()

#Combine Neoehrlichia_mikurensis columns 
taxa_transposed_16S_rpoB <- taxa_transposed_16S_rpoB %>%
  mutate(Neoehrlichia_mikurensis = Neoehrlichia_mikurensis.x + Neoehrlichia_mikurensis.y) %>%
  select(- c("Neoehrlichia_mikurensis.x", "Neoehrlichia_mikurensis.y") ) 

#Control rpoB metabarcoding - identify failed sample
#To identify samples which were positives for Bartonella in 16S but are not with rpoB
taxa_transposed_16S_rpoB %>%
  filter(rowSums(across( names(taxa_transposed_rpoB[,-1]))) == 0) %>%
  filter(Bartonella > 0 ) %>%
  select(numero_centre_combined) %>%
  pull()


# Lipl32 data treatment ----

#Transform lipl32 qualitative column
filelipl32 <- filelipl32 %>%
  mutate(lipL32_result = if_else(lipL32_result == "-", "0", lipL32_result) |> as.numeric()) %>%
  select(ID_rodent, lipL32_result)


## Join lipl32 to metabarcoding data ----

#Join data
taxa_transposed_16S_rpoB_lipl32 <- left_join(taxa_transposed_16S_rpoB, 
                                             filelipl32,
                                      by = c("numero_centre_combined" = "ID_rodent"))

#Identify sample without lipl32 results (for them, metabarcoding results will be used for Leptospira)
taxa_transposed_16S_rpoB_lipl32 %>%
  filter(is.na(lipL32_result)) %>%
  select(numero_centre_combined)

#See benefit of Lepto and lipl32 columns combining
taxa_transposed_16S_rpoB_lipl32 %>%
  filter(lipL32_result > 0) %>%
  filter(Leptospira == 0) %>%
  select(numero_centre_combined) %>%
  pull()

taxa_transposed_16S_rpoB_lipl32 %>%
  filter(Leptospira > 0) %>%
  filter(lipL32_result == 0) %>%
  select(numero_centre_combined) %>%
  pull()

#Combine Leptospira columns 
taxa_transposed_16S_rpoB_lipl32 <- taxa_transposed_16S_rpoB_lipl32 %>%
  mutate(Leptospira = case_when(
    (Leptospira == 0) & !is.na(lipL32_result) & lipL32_result > 0 ~ lipL32_result,
    Leptospira > 0 & !is.na(lipL32_result) & lipL32_result > 0 ~ lipL32_result,
    TRUE ~ Leptospira  
    )) %>%
  select(-lipL32_result)



# jsp ce que je voulais faire en bas, mais renommer numero_centre peut-être pas mal (genre rpoB-16s?)
# genre controle quelles echantillons de l'un passé dans l'autre (tous les barto + 16s passé rpoB?)
# peut-être faisable avant ça deja ?
# + ça a changé car j'ai modif pour ajout de lipl32



# ## TO REWORK Join to macroparasite for complete pathogen data frame (hasbeen, and was only related to ticks) ----
# 
# # Select columns for macroparasite of interest (data from temp_tique script) and transform NA to 0 i effectif_tick
# d_macroparasite <- d_macroparasite %>%
#   select( c("effectif_tick", "numero_centre", "Ixodida", "Siphonaptera"))
# 
# # Join macroparasite on host trapping data, keep every lines from trapping data
# d_macroparasite <- left_join(d_host, d_macroparasite, by = "numero_centre" )
# 
# # Transform NA to 0
# d_macroparasite <- d_macroparasite %>%
#   mutate( across(c("Ixodida", "Siphonaptera"), ~if_else(. == 0, NA, .) ) ) %>%
#   mutate( across(c("Ixodida", "Siphonaptera"), ~if_else( is.na(.), 0, 1) ) ) %>%
#   mutate( effectif_tick = if_else( is.na(effectif_tick), 0, effectif_tick) )
# 
# # Join 16s pathogens on host trapping/macroparasite data - complete join, keep only 16s screened individuals /!\
# rodent_pathos <- merge(x = d_macroparasite, y = taxa_transposed_16S_rpoB, by.x = "numero_centre",  by.y = "numero_centre_combined")


# # Extract Pathogens name
# pathos_name <- rodent_pathos |>
#   select(names(rodent_pathos)[(which(names(rodent_pathos) == "effectif_tick") + 1):ncol(rodent_pathos)]) %>%
#   colnames()
# pathos_name



# Join rodent and line data ----
rodent_pathos <- left_join(d_host,
                           taxa_transposed_16S_rpoB_lipl32,
                           by = c("numero_centre" = "numero_centre_combined"))

# Vector containing pathogen names
pathos_name <- rodent_pathos %>%
  select(names(rodent_pathos)[(which(names(rodent_pathos) == tail(colnames(d_host), 1) ) + 1):ncol(rodent_pathos)]) %>%
  colnames()

#Transform reads into presence/absence data
rodent_pathos <- rodent_pathos %>%
  mutate(across(all_of(pathos_name), ~ replace(., . > 0, 1)))

# Beware if some pathogen for which abundance is important are to be added before this step (ticks for example)
#please apply function not on all pathos name but pathos name minus those specific names


# Generate sub-datasets ----

## Apply A. sylvaticus only ----

#Keep only Apodemus sylvaticus for analysis
d_apo_pathos <- rodent_pathos %>%
  filter(taxon_mamm == "Apodemus sylvaticus") %>%
  filter(!category == "broadleaved_forest")                                          #ATTENTION temp pour evmc, retirer ensuite----

d_apo_pathos %>%
  select(all_of(c(
    "Borreliella",
    "Neoehrlichia_mikurensis",
    "Streptobacillus",
    "Mycoplasma_haemomuris",
    "Mycoplasma_coccoides",
    "Mycoplasma_penetrans",
    "Ehrlichia",
    "Neisseria",
    "Spiroplasma",
    "Treponema",
    "Bartonella_taylorii",
    "Bartonella_grahamii",
    "Bartonella_doshiae",
    "Bartonella_birtlesii",
    "Bartonella_elizabethae",
    "Leptospira"
  ))) %>%
  summarise(across(everything(), ~ sum(. > 0))) %>%
  tidyr::pivot_longer(cols = everything(),
               names_to = "pathogen",
               values_to = "n_positive") %>%
  fwrite(here::here("nb_indiv_evmc.csv") )




# Filter empty pathogen for Apodemus sylvaticus
d_apo_pathos <- d_apo_pathos %>%
  select(-names(which(colSums(d_apo_pathos[, pathos_name]) == 0)))

# # Regroup some specific taxa for analysis purposes 
# d_apo_pathos <- d_apo_pathos %>%
#   mutate(Borrelia = Borrelia + Borreliella_afzelii) %>%
#   select(- Borreliella_afzelii)

#Generate apo pathos name vector
pathos_name_apo <- pathos_name[pathos_name %in% names(d_apo_pathos)]
setdiff(pathos_name, pathos_name_apo)

# Identify pathos with GLOBAL prevalence >=10 for Apodemus data
patho10_apo <- d_apo_pathos %>%
  summarise(across(all_of(pathos_name_apo), ~ sum(. > 0) / n() ))  %>%
  unlist()
patho10_apo <- names(patho10_apo[patho10_apo >= 0.10])
patho10_apo


## More sub datasets based on filters ----

#Add new variable : pathogen richness
#But before choose if certain taxa should be excluded from pathogen richness calculation
candidate_richness_pathos <- setdiff(pathos_name_apo, "Bartonella")

candidate_richness_pathos <- c(                                     # ATTENTION TEMP EVMC ---- 
  "Borreliella",
  "Neoehrlichia_mikurensis",
  "Streptobacillus",
  "Mycoplasma_haemomuris",
  "Mycoplasma_coccoides",
  "Mycoplasma_penetrans",
  "Ehrlichia",
  "Neisseria",
  "Spiroplasma",
  "Treponema",
  "Bartonella_taylorii",
  "Bartonella_grahamii",
  "Bartonella_doshiae",
  "Bartonella_birtlesii",
  "Bartonella_elizabethae",
  "Leptospira"
)

data_for_m <- d_apo_pathos %>%
  mutate(number_pathos = rowSums(across(all_of(candidate_richness_pathos), ~ (. > 0))))  

#Take away rows with NA for some of our factor
data_for_m %>%
  filter(is.na(sexe))

data_for_m <- data_for_m %>%
  filter(!is.na(sexe))


#Generate dataset without broadleaved forest analysis
data_for_m_noforests <- data_for_m %>%
  filter(category != "broadleaved_forest")


#Generate dataset without broadleaved forest with lines containing PCA values
data_for_m_noforests_pca <- data_for_m_noforests %>%
  filter(!is.na(PCA_lines_local_axis1))

  
#Generate other dataset without broadleaved forest, without year that do not contains them
data_for_m_forestsyear <- data_for_m %>%
  filter(code_mission %in% (data_for_m %>% 
           filter(category == "broadleaved_forest") %>%
           select(code_mission) %>% 
           pull() %>%
           unique()))





## Parallel dataset containing helminths abundance data ----

#Extract helminths names
pathos_name_apo_helm <- file_helm %>%
  select(where(is.numeric)) %>%
  names()

#Join with data_for_m (apodemus data with some conversion)
data_for_m_helm <- data_for_m %>%
  left_join(file_helm,
            by = c("numero_centre" = "code_rongeur"))

#Only keep individuals dissected for helminths
data_for_m_helm <- data_for_m_helm %>%
  filter(!if_any(all_of(pathos_name_apo_helm), is.na))

#Check that every individual in helm file is retrieved in new file
setdiff(data_for_m_helm$numero_centre, file_helm$code_rongeur)
setdiff(file_helm$code_rongeur, data_for_m_helm$numero_centre)

#Regenerate and relocate pathogen richness
data_for_m_helm <- data_for_m_helm %>%
  mutate(number_helm_sp = rowSums(across(all_of( c(pathos_name_apo_helm)), ~ . > 0))) %>%
  mutate(number_pathos = rowSums(across(all_of( c(candidate_richness_pathos, pathos_name_apo_helm)), ~ . > 0))) %>%
  relocate(number_helm_sp, .after = last_col()) %>%
  relocate(number_pathos, .after = last_col())

hist(data_for_m_helm$number_pathos)
hist(data_for_m_helm$number_helm_sp)


# Identify pathos with GLOBAL prevalence >=10 for this subset of apodemus data
patho10_apo_helm <- data_for_m_helm %>%
  summarise( across(all_of(c(pathos_name_apo, pathos_name_apo_helm)), ~ sum(. > 0) / n() ))  %>%
  unlist()
patho10_apo_helm <- names(patho10_apo_helm[patho10_apo_helm >= 0.10])
patho10_apo_helm


#Generate helm data from lines containing PCA values
data_for_m_helm_pca <- data_for_m_helm %>%
  filter(!is.na(PCA_lines_local_axis1))
