
# /!\ WARNING ----
#THIS SCRIPT IS TEMPORARY AND SHOULD NOT BE CONSERVED IN THIS R PROJECT LATER ----

#the script extracting the host data from bpm has to be run before this script
#this script should only be used to analyse 16S or rpob/glta data from spleen, in order to identify pathogenic bacteria
#this script is dependent on the generation of d_host dataframe from temp_import_survey and d_macroparasite from rodent_macroparasite to join dataframes

# /!\ WARNING ----



# Some parameters ----

library("data.table")
library("dplyr")
library("ggplot2")

# Graphical parameters 

category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
category_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
type_order <- c("CT_LB", "CT_HB", "NC_LB", "NC_HB", "C_LB", "C_HB", "B" )
brd_palette <- c("LB" = "#FFB3BA", "HB" = "#B3E2B3")
mission_order <- c("Juin 2023", "Octobre 2023", "Juin 2024", "Septembre 2024")
mission_color <- c("Juin 2023" = "#66c2a5", "Octobre 2023" = "#fc8d62", "Juin 2024" = "#8da0cb", "Septembre 2024" = "#ab7a82")


# Function to count individuals number in graphs
count_summary <- function(x, y_position = max(x) + 0.5) {
  return(data.frame(y = y_position, label = paste0("n = ", length(x))))
}

# Function to darken color of a palette
darken_color <- function(color, amount = 0.2) {
  color_rgb <- col2rgb(color) / 255
  darkened_rgb <- color_rgb * (1 - amount)
  rgb(darkened_rgb[1], darkened_rgb[2], darkened_rgb[3])
}


# GENERATE DATA ----

## Import data ----
# Import hosts and line modalities file
d_host <- readr::read_csv( here::here("data/", "raw-data/", "host_data", "20241203_bpm_modalities.csv") )

# Import 16S filtered file
file16s_run00_01_04_05 <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01-04-05","Run00-01-04-05_16S_filtered-merged_postfrogs.txt"))

# Import 16S filtered file
filerpoB_run01_05 <- data.table::fread(file = here::here( "data","raw-data/","rpoB_run01-05","Run01-05-rpoB_filtered-merged.txt"))

# Import rodent macroparasite
d_macroparasite <- data.table::fread(file = here::here("data", "derived-data", "ticks", "rodents_tick", "20240731_macroparasite.csv") )


## Generate a file containing samples of interest ----

# Extract id of small mammals caught and dissected in Beprep
sm_id <- unique(d_host %>% 
                  filter(!is.na(numero_centre))%>%
                  filter(stringr::str_detect(numero_centre, pattern = "NCHA")) %>%
                  pull(numero_centre)
)

# Identifiy taxonomy columns
taxo_name <- colnames( file16s_run00_01_04_05 [, 1:which(colnames(file16s_run00_01_04_05) == "observation_sum") ])

# Calculate total run read number per cluster (after our filtering)
file16s_run00_01_04_05 <- file16s_run00_01_04_05 %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name))) %>%
  relocate(clean_totalrunreads, .after = names(file16s_run00_01_04_05)[which(colnames(file16s_run00_01_04_05) == "observation_sum")])

# Reattribute taxonomy columns
taxo_name <- colnames( file16s_run00_01_04_05 [, 1:which(colnames(file16s_run00_01_04_05) == "clean_totalrunreads") ])

# Take away problematic samples - if needed
#file16s_run00_01_04_05 <- file16s_run00_01_04_05 |>

# Generate new file containing only Spleen BePrep's samples
#by selecting those samples
beprep_sample_16s_sp <- file16s_run00_01_04_05 |>
  select( (contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
beprep16s_sp <- file16s_run00_01_04_05 |>
  select( taxo_name ) |>
  cbind(beprep_sample_16s_sp)

# Delete the PCZymo positive after visual control
beprep16s_sp <- beprep16s_sp |>
  select(!contains("PCzymo"))

# Calculate total number of reads for each cluster considering chosen samples
beprep16s_sp <- beprep16s_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name)))%>%
  relocate(totalreads, .after = names(beprep16s_sp)[which(colnames(file16s_run00_01_04_05) == "clean_totalrunreads")])

# Order rows by new total read number
beprep16s_sp <- data.table::setorder(beprep16s_sp, -totalreads)

# Delete cluster with no reads for the selected individuals
beprep16s_sp <- beprep16s_sp |> 
  filter(totalreads >0)

# Write a file containing 16S Spleen OTUs containing our samples only :
data.table::fwrite(beprep16s_sp, here::here("data", "derived-data","16S", "20241204_beprep16s_sp.txt") )



## Generate 16S spleen putative pathogen data ----

# Split the taxonomy to access taxon
pathos16s <- beprep16s_sp |>
  tidyr::separate_wider_delim(cols = blast_taxonomy,
                              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                              delim = ";")

# List of rodents identified putative pathogens
putative_pathos_genra <- c("Bartonella", "Neoehrlichia", "Mycoplasma", "Borrelia", "Borreliella", "Streptobacillus", "Ehrlichia", "Spiroplasma","Brevinema",
                           "Orientia", "Rickettsia", "Leptospira", "Yersinia", "Actinobacillus", "Treponema", "Chlamydia", "Neisseria", "Pasteurella",
                           "Francisella", "Brucella", "Coxiella")
putative_pathos_family <- c("Sarcocystidae")

# Keep only OTUs with known putative pathogens (you should look at raw data before this step, ensure it is ok and you're not missing pathogens taxa as the list is not complete)
pathos16s <- pathos16s |>
  filter(genus %in% putative_pathos_genra 
         | family %in% putative_pathos_family)

# Write a file containing 16S Spleen OTUs containing our samples only AND only putative pathogens :
data.table::fwrite(pathos16s, here::here("data/", "derived-data/","16S", "20241204_beprep16s_sp_pathos.txt") )


## Seek for species level identification for specific taxa (e.g. Mycoplasma) throught external manipulation (phylogenic tree) ----
# A new datafile is created BY HAND to add the checked species information (in new column "species_adjusted") and then imported
# If probable pseudo-gene affiliation is observed in data file they may also get deleted ny operator at this step
pathos16s_checked <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01-04-05","20241204_beprep16s_sp_pathos-adjusted.txt"))
pathos16s_checked$species_adjusted[pathos16s_checked$species_adjusted == ""] <- NA

# Change species value when there is a species adjusted value (generated by operator manually ) and drop the species adjusted column
pathos16s_checked <- pathos16s_checked %>%
  mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
  select(!species_adjusted)

# Rename inappropriate affiliation 
pathos16s_checked <- pathos16s_checked %>%
  mutate(species = if_else(species == "Haemobartonella muris", "Mycoplasma haemomuris", species))


## Create new column for best identified taxa ----

# If there are unknown or Multi-affiliation in species column, change it to NA
pathos16s_checked <- pathos16s_checked %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species )) 

# Creation of a new column containing best available identity (between species or genus here)
pathos16s_checked <- pathos16s_checked %>%
  mutate(
    identification16s = if_else(!is.na(species), species, genus),
    best_identification = if_else(!is.na(species), "species", "genus")
  ) %>%
  relocate(identification16s, .after = species) %>%
  relocate(best_identification, .after = identification16s)

# Replace space in taxa name 
pathos16s_checked <- pathos16s_checked %>%
  mutate(identification16s = gsub(" ", "_", identification16s ))


## Simplification, compress 16s information by best identified taxa ----

# Keep only best taxa and samples
taxa_pathos16s <- pathos16s_checked %>%
  select(c("identification16s", names(pathos16s_checked)[(which(names(pathos16s_checked) == "totalreads") + 1):ncol(pathos16s_checked)]))

# Sum rows by best taxa
taxa_pathos16s <- taxa_pathos16s %>%
  group_by(identification16s) %>%
  summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop")

# Transpose the table
taxa_transposed_16S <- taxa_pathos16s |>
  data.table::transpose(
    keep.names = "numero_centre_16s",
    make.names = 1
  )

# Finish sample name transformation to join to host data
taxa_transposed_16S <- taxa_transposed_16S |>
  mutate(numero_centre_16s = gsub(".SP", "", numero_centre_16s ))







# Repeat process for rpoB data ----BEWARRRREEE ----


## Generate a file containing samples of interest ----

# Identifiy taxonomy columns
taxo_name <- colnames( filerpoB_run01_05 [, 1:which(colnames(filerpoB_run01_05) == "observation_sum") ])

# Calculate total run read number per cluster (after our filtering)
filerpoB_run01_05 <- filerpoB_run01_05 %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name))) %>%
  relocate(clean_totalrunreads, .after = names(filerpoB_run01_05)[which(colnames(filerpoB_run01_05) == "observation_sum")])

# Reattribute taxonomy columns
taxo_name <- colnames( filerpoB_run01_05 [, 1:which(colnames(filerpoB_run01_05) == "clean_totalrunreads") ])

# Take away problematic samples - if needed
#filerpoB_run01_05 <- filerpoB_run01_05 |>

# Generate new file containing only Spleen BePrep's samples
#by selecting those samples
beprep_sample_rpoB_sp <- filerpoB_run01_05 |>
  select( (contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
bepreprpoB_sp <- filerpoB_run01_05 |>
  select( taxo_name ) |>
  cbind(beprep_sample_rpoB_sp)

# Verify absence of PCZymo positive in rpoB data
bepreprpoB_sp |>
  select(contains("PCzymo"))

# Calculate total number of reads for each cluster considering chosen samples
bepreprpoB_sp <- bepreprpoB_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name)))%>%
  relocate(totalreads, .after = names(bepreprpoB_sp)[which(colnames(filerpoB_run01_05) == "clean_totalrunreads")])

# Order rows by new total read number
bepreprpoB_sp <- data.table::setorder(bepreprpoB_sp, -totalreads)

# Delete cluster with no reads for the selected individuals
bepreprpoB_sp <- bepreprpoB_sp %>% 
  filter(totalreads >0)




##Generate rpoB spleen  data ----
  
# Split the taxonomy to access taxon
bepreprpoB_sp <- bepreprpoB_sp |>
  tidyr::separate_wider_delim(cols = blast_taxonomy,
                              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                              delim = ";")

# Delete cluster without real affiliation
bepreprpoB_sp <- bepreprpoB_sp %>%
  filter(!domain == "no")

# Write a file containing rpoB Spleen OTUs containing our samples only :
data.table::fwrite(bepreprpoB_sp, here::here("data", "derived-data","rpoB", "20241204_bepreprpoB_sp.txt") )


# Temporary filter for rpoB #################
bepreprpoB_sp <- bepreprpoB_sp %>%
  filter(totalreads > 100 )


# ## Seek for species level identification for specific taxa (e.g. Mycoplasma) throught external manipulation (phylogenic tree) ----
# # A new datafile is created BY HAND to add the checked species information (in new column "species_adjusted") and then imported
# # If probable pseudo-gene affiliation is observed in data file they may also get deleted ny operator at this step
# pathos16s_checked <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01-04-05","20241204_beprep16s_sp_pathos-adjusted.txt"))
# pathos16s_checked$species_adjusted[pathos16s_checked$species_adjusted == ""] <- NA
# 
# # Change species value when there is a species adjusted value (generated by operator manually ) and drop the species adjusted column
# pathos16s_checked <- pathos16s_checked %>%
#   mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
#   select(!species_adjusted)
# 
# # Rename inappropriate affiliation 
# pathos16s_checked <- pathos16s_checked %>%
#   mutate(species = if_else(species == "Haemobartonella muris", "Mycoplasma haemomuris", species))


## Create new column for best identified taxa ----

# If there are unknown or Multi-affiliation in species column, change it to NA
bepreprpoB_sp <- bepreprpoB_sp %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species )) 

# Creation of a new column containing best available identity (between species or genus here)
bepreprpoB_sp <- bepreprpoB_sp %>%
  mutate(
    identificationrpoB = if_else(!is.na(species), species, genus),
    best_identificationrpoB = if_else(!is.na(species), "species", "genus")
  ) %>%
  relocate(identificationrpoB, .after = species) %>%
  relocate(best_identificationrpoB, .after = identificationrpoB)

# Replace space in taxa name 
bepreprpoB_sp <- bepreprpoB_sp %>%
  mutate(identificationrpoB = gsub(" ", "_", identificationrpoB ))


## Simplification, compress 16s information by best identified taxa ----

# Keep only best taxa and samples
taxa_pathosrpoB <- bepreprpoB_sp %>%
  select(c("identificationrpoB", names(bepreprpoB_sp)[(which(names(bepreprpoB_sp) == "totalreads") + 1):ncol(bepreprpoB_sp)]))

# Sum rows by best taxa
taxa_pathosrpoB <- taxa_pathosrpoB %>%
  group_by(identificationrpoB) %>%
  summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop")

# Transpose the table
taxa_transposed_rpoB <- taxa_pathosrpoB |>
  data.table::transpose(
    keep.names = "numero_centre_rpoB",
    make.names = 1
  )

# Finish sample name transformation to join to host data
taxa_transposed_rpoB <- taxa_transposed_rpoB |>
  mutate(numero_centre_rpoB = gsub(".SP", "", numero_centre_rpoB ))





# Join rpoB and 16S data ----
taxa_transposed_16S_rpoB <- left_join(taxa_transposed_16S, 
                                      taxa_transposed_rpoB,
                                      by = c("numero_centre_16s" = "numero_centre_rpoB"))

# Replace NA by 0 in new data
taxa_transposed_16S_rpoB <- taxa_transposed_16S_rpoB %>%
  mutate(across(-numero_centre_16s, ~ if_else(is.na(.), 0, .)))

# Combine relevant taxa
taxa_transposed_16S_rpoB <- taxa_transposed_16S_rpoB %>%
  mutate(Neoehrlichia_mikurensis = Neoehrlichia_mikurensis + Neoehrlichia_lotoris) %>%
  select(-Neoehrlichia_lotoris)



# jsp ce que je voulais faire en bas, mais renommer numero_centre peut-être pas mal (genre rpoB-16s?)
# genre controle quelles echantillons de l'un passé dans l'autre (tous les barto + 16s passé rpoB?)
# peut-être faisable avant ça deja ?

# taxa_transposed_16S_rpoB %>%
#   mutate(numero_centre_sp = numero_centre_16s) %>%
#   filter(-) 










# ## Join to macroparasite for complete pathogen data frame ----
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
# rodent_pathos <- merge(x = d_macroparasite, y = taxa_transposed_16S_rpoB, by.x = "numero_centre",  by.y = "numero_centre_16s")


# # Extract Pathogens name
# pathos_name <- rodent_pathos |>
#   select(names(rodent_pathos)[(which(names(rodent_pathos) == "effectif_tick") + 1):ncol(rodent_pathos)]) %>%
#   colnames()
# pathos_name



# Temporary solution : -attention-
rodent_pathos <- left_join(d_host, taxa_transposed_16S_rpoB,
                           by = c("numero_centre" = "numero_centre_16s"))

# devait etre fait avant surement mais bon :
rodent_pathos <- rodent_pathos %>%
  filter(stringr::str_detect(numero_centre, pattern = "NCHA") ) 


pathos_name <- rodent_pathos %>%
  select(names(rodent_pathos)[(which(names(rodent_pathos) == "broadleaved_class") + 1):ncol(rodent_pathos)]) %>%
  colnames()



# Global data exploration ----

# Number infected per taxa
rodent_pathos %>%
  mutate(across(all_of(pathos_name), ~ replace(., . > 0, 1))) %>%
  group_by(taxon_mamm) %>%
  summarise(across(all_of(pathos_name), sum)) 



# List of pathogen per species
list_pathos_per_species <- list()

for (i in unique(rodent_pathos$taxon_mamm)) {
  list_pathos_per_species[[i]] <- rodent_pathos %>%
    filter(taxon_mamm == i) %>%
    select(names(which(colSums(.[, pathos_name]) > 0))) %>%
    colnames()
}
list_pathos_per_species




# ANALYSIS (ONLY FOR APODEMUS FOR NOW) ----

## Generate data ----

# Keep only Apodemus for analysis
d_apo_pathos_glm <- rodent_pathos %>%
  filter(taxon_mamm == "Apodemus sylvaticus")

# Filter empty pathogen for Apodemus 
d_apo_pathos_glm <- d_apo_pathos_glm %>%
  select(-names(which(colSums(d_apo_pathos_glm[, pathos_name]) == 0)))

# # Regroup some specific taxa for analysis purposes 
# d_apo_pathos_glm <- d_apo_pathos_glm %>%
#   mutate(Borrelia = Borrelia + Borreliella_afzelii) %>%
#   select(- Borreliella_afzelii)

# # Generate apo pathos name vector
# pathos_name_apo <- d_apo_pathos_glm |>
#   select(names(d_apo_pathos_glm)[(which(names(d_apo_pathos_glm) == "effectif_tick") + 1):ncol(d_apo_pathos_glm)]) %>%
#   colnames()
# setdiff(pathos_name, pathos_name_apo)

# Generate apo pathos name vector
pathos_name_apo <- d_apo_pathos_glm |>
  select(names(d_apo_pathos_glm)[(which(names(d_apo_pathos_glm) == "broadleaved_class") + 1):ncol(d_apo_pathos_glm)]) %>%
  colnames()
setdiff(pathos_name, pathos_name_apo)


# Generate the dataframe for logistic analysis (0/1)
data_for_m <- d_apo_pathos_glm %>%
  mutate(across(all_of(pathos_name_apo), ~ replace(., . > 0, 1)))

# Add new variable : pathogen richness
data_for_m <- data_for_m %>%
  mutate(number_pathos = rowSums(across(all_of(pathos_name_apo))) ) 

# Identify pathos with GLOBAL prevalence >=10 for whole data
patho10_apo <- d_apo_pathos_glm %>%
  summarise( across(all_of(pathos_name_apo), ~ sum(. > 0) / n() ))  %>%
  unlist()
patho10_apo <- names(patho10_apo[patho10_apo >= 0.10])
patho10_apo




# Take away rows with NA for some of our factor
data_for_m <- data_for_m %>%
  filter(!is.na(sexe))

# Generate dataset without broadleaved forest for different analysis
data_for_m_noforests <- data_for_m %>%
  filter(category != "broadleaved_forest")
  





## Exploration (for apodemus only) ----

# PREVALENCE ALL MISSION
data_for_m %>%
  summarise(across(all_of(pathos_name_apo), ~ sum(.) * 100 / n()), .groups = "drop") %>%
  tidyr::pivot_longer(
    cols = all_of(pathos_name_apo),
    names_to = "pathogens",
    values_to = "prevalence"
  )


# PREVALENCE PER MISSION
data_for_m %>%
  group_by(code_mission) %>%
  summarise(across(all_of(pathos_name_apo), ~ sum(.) * 100 / n()), .groups = "drop") 


# PREVALENCES PER SITE (beware, sometiques really few individuals per site rise prevalence + line with no individuals are not integrated)
# Init pathogen prevalence per grouping factor 
grouping_prevalence_factor <- c("numero_ligne")


# Number of positive individuals per grouping factor
pp <- data_for_m %>%
  group_by(across(all_of(grouping_prevalence_factor))) %>%
  summarise(
    effectif = n(),
    type = unique(type),
    broadleaved_class = unique(broadleaved_class),
    treatment = unique(treatment),
    category = unique(category),
    across(all_of(pathos_name_apo), ~ sum(. > 0)),  
    .groups = "drop"
  )

# Calculate the total number of individuals and prevalence for each pathogen
pp_prevalence <- pp %>%
  group_by(across(all_of(grouping_prevalence_factor))) %>%
  summarise(
    type = unique(type),
    broadleaved_class = unique(broadleaved_class),
    treatment = unique(treatment),
    category = unique(category),
    effectif = sum(effectif),  # Total number of individuals
    across(all_of(pathos_name_apo), sum),  # Sum of positive cases for each pathogen
    .groups = "drop"
  ) %>%
  mutate(across(all_of(pathos_name_apo), ~ . / effectif))

# Pivot prevalence table
pp_prevalence <- pp_prevalence %>%
  tidyr::pivot_longer ( cols = all_of(pathos_name_apo),
                        names_to = "pathos",
                        values_to = "prevalence")

# Prevalence watching only for the one above 10 % globally
pp_prevalence <- pp_prevalence %>%
  filter(pathos %in% patho10_apo)

# Watch prevalences with plot for all pathos
pp_prevalence %>%  
  filter(!pathos %in% c("Siphonaptera", "Ixodida")) %>%
  ggplot(aes(x = prevalence, y = pathos, fill = category)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.6) +
  geom_jitter(aes(color = category), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = category_palette) +
  scale_color_manual(values = category_palette) + 
  labs(x = "Prevalence par ligne de piegeage", y = "Pathogènes") +
  guides(fill = guide_legend(title = "Type de ligne"), 
         color = guide_legend(title = "Type de ligne")) +
  theme_minimal()

pp_prevalence %>%  
  filter(!pathos %in% c("Siphonaptera", "Ixodida")) %>%
  ggplot(aes(x = prevalence, y = pathos, fill = category)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.6) +
  geom_jitter(aes(color = category), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = category_palette) +
  scale_color_manual(values = category_palette) + 
  labs(x = "Prevalence par ligne de piegeage", y = "Pathogènes") +
  guides(fill = guide_legend(title = "Type de ligne"), 
         color = guide_legend(title = "Type de ligne")) +
  theme_minimal() +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = category),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)


# test neoehrlichia
pp_prevalence %>%  
  filter(pathos == "Neoehrlichia_mikurensis") %>%
  ggplot(aes(x = prevalence, y = broadleaved_class, fill = broadleaved_class)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = brd_palette) +
  labs(x = "Prevalence par ligne de piegeage", y = "Broadleaved Status") +
  guides(fill = guide_legend(title = "Broadleaved Status")) +
  theme_minimal() +
  stat_summary(fun.data = count_summary, 
               aes(group = category),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)


data_for_m %>% count(Neoehrlichia_mikurensis, type)

ggplot(data_for_m, aes(x = as.factor(Mycoplasma_haemomuris), y = poids)) +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = sexe),
               dotsize = 0.7 ,
               position=position_dodge(0.4)) + 
  coord_flip() +
  stat_summary(fun = mean, geom="point", shape=18,
               size=5, color="purple") +
  theme_minimal()

hist(data_for_m$number_pathos)


## GLMMs ----

### Model : Neoehrlichia_mikurensis  ----
rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ treatment * broadleaved_class + scale(poids) + code_mission  + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ treatment + broadleaved_class  + code_mission + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

DHARMa::simulateResiduals(m_neoeh_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_neoeh_r,.~.,test="Chisq")
summary(m_neoeh_r)

em <- emmeans::emmeans(m_neoeh_r, specs = pairwise ~ code_mission, adjust = "Tukey", type = "response" )
em$contrasts

plot(em, comparisons = TRUE)


ggstats::ggcoef_model(m_neoeh_r)

gtsummary::tbl_regression(m_neoeh_r)



# Alternative - just test category

rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_neoeh_r)

gtsummary::tbl_regression(m_neoeh_r)





# 
 #### Test graph ----
# 
# # Extract fixed effects coefficients
# fixed_effects <- coef(summary(m_neoeh_r))[, "Estimate"] # ou : fixed_effects <- lme4::fixef(m_neoeh_r)
# 
# # Create data for plotting
# plot_data <- expand.grid(
#   broadleaved_class = unique(data_for_m$broadleaved_class),
#   code_mission = unique(data_for_m$code_mission),
#   poids = seq(min(data_for_m$poids), max(data_for_m$poids), length.out = 100)  
# )
# 
# # Predict response variable
# plot_data$predicted_prob <- predict(m_neoeh_r, newdata = plot_data, type = "response", re.form = NA)
# 
# # Plot
# library(ggplot2)
# ggplot(plot_data, aes(x = poids, y = predicted_prob, color = broadleaved_class)) +
#   geom_line() +
#   facet_wrap(~ code_mission) +
#   labs(x = "Poids", y = "Predicted Probability", color = "Broadleaved Status")
# 



### Model : Mycoplasma haemomuris ----
rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_mycoplasma_r_best <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ scale(poids) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

DHARMa::simulateResiduals(m_mycoplasma_r_best, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_mycoplasma_r_best,.~.,test="Chisq")
summary(m_mycoplasma_r_best)

ggstats::ggcoef_model(m_mycoplasma_r_best)

gtsummary::tbl_regression(m_mycoplasma_r_best)



#Alternative

rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


ggstats::ggcoef_model(m_mycoplasma_r_best)

gtsummary::tbl_regression(m_mycoplasma_r_best)








### Model : Mycoplasma coccoides ----

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ poids +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(m_mycoplasmacoco_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_mycoplasmacoco_r,.~.,test="Chisq")
summary(m_mycoplasmacoco_r)


ggstats::ggcoef_model(m_mycoplasmacoco_r)

gtsummary::tbl_regression(m_mycoplasmacoco_r)




# Alternative
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels









### Model : Bartonella ----
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_r_best <- lme4::glmer(
  formula = Bartonella ~ scale(poids) + code_mission + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

DHARMa::simulateResiduals(m_barto_r_best, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_barto_r_best,.~.,test="Chisq")
summary(m_barto_r_best)


ggstats::ggcoef_model(m_barto_r_best)

gtsummary::tbl_regression(m_barto_r_best)



#Alternative
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels






### Model : Bartonella_taylorii ----
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_t_r)

gtsummary::tbl_regression(m_barto_t_r)


#Alternative
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels




### Model : Bartonella_grahamii ----
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_g_r)

gtsummary::tbl_regression(m_barto_g_r)


#Alternative
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels



### Model : Bartonella_birtlesii ----

m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_b_r)

gtsummary::tbl_regression(m_barto_b_r)



#Alternative
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels





















### Model : Ixodida  ----
m_ixod_r <- lme4::glmer(
  formula = Ixodida ~ broadleaved_class * treatment + code_mission + poids + sexe +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

SelectionModels<- MuMIn::dredge(m_ixod_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


### Model : Siphonaptera  ----
m_siphon_r <- lme4::glmer(
  formula = Siphonaptera ~ broadleaved_class * treatment + code_mission + poids + sexe +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

SelectionModels<- MuMIn::dredge(m_siphon_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_siphon_r1 <- lme4::glmer(
  formula = Siphonaptera ~ code_mission  +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
DHARMa::simulateResiduals(m_siphon_r1, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_siphon_r1,.~.,test="Chisq")
summary(m_siphon_r1)


m_siphon_r2 <- lme4::glmer(
  formula = Siphonaptera ~ sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(m_siphon_r2, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_siphon_r2,.~.,test="Chisq")
summary(m_siphon_r2)

em <- emmeans::emmeans(m_siphon_r1, specs = pairwise ~ code_mission, adjust = "Tukey", type = "response" )
em$contrasts

plot(em, comparisons = TRUE)

###Model : Pathogen Richnesse  ----
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ treatment * broadleaved_class + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

DHARMa::simulateResiduals(m_pathnumber_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~  code_mission + poids  +(1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

m_pathnumber_r %>%
  RVAideMemoire::overdisp.glmer()

DHARMa::simulateResiduals(m_pathnumber_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_pathnumber_r,.~.,test="Chisq")
summary(m_pathnumber_r)





# Alternative
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


# ANALYSIS (NOT ONLY ON APODEMUS) ----

### Model : Parasite Diversity for all rodents  ----

#### Generate alpha diversity index ----

alpha_div_pathos <- rodent_pathos %>%
  mutate(across(all_of(pathos_name), ~ replace(., . > 0, 1))) %>%
  group_by(numero_ligne, code_mission) %>%
  summarise(treatment = unique(treatment),
            broadleaved_class = unique(broadleaved_class),
            across(all_of(pathos_name), ~ sum(. > 0))
  ) %>%
  mutate( richness = rowSums(across(all_of(pathos_name)) > 0)) %>%
  mutate (shannon = vegan::diversity(across(all_of(pathos_name)) ))

hist(alpha_div_pathos$richness)
hist(alpha_div_pathos$shannon)



#### Poisson regression ----

m_rich_pathos <- stats::glm(
  formula = richness  ~ broadleaved_class * treatment ,
  family = poisson(link="log"),
  data = alpha_div_pathos,
  na.action="na.fail"
)

DHARMa::simulateResiduals(m_rich_pathos, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

summary(m_rich_pathos)

performance::check_zeroinflation(m_rich_pathos) 


#### zero_inflated poisson regression ----
m_richness_zip <-pscl::zeroinfl(richness ~ treatment + broadleaved_class + code_mission, 
                                data = alpha_div_pathos, dist = "poisson", na.action = "na.fail")

summary(m_richness_zip)
drop1(m_richness_zip,.~.,test="Chisq")


SelectionModels<- MuMIn::dredge(m_richness_zip, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_richness_zip <-pscl::zeroinfl(richness ~ treatment, 
                                data = alpha_div_pathos, dist = "poisson", na.action = "na.fail")

summary(m_richness_zip)
drop1(m_richness_zip,.~.,test="Chisq")

## Heatmap matrix ----

# Prevalence matrix calculation (all rodents)
#for type
graph_sequence <- c("CT_LB", "CT_HB", "NC_LB", "NC_HB", "C_LB", "C_HB" )

matrix_pathoss <- d_apo_pathos_glm %>%
  group_by(treatment, broadleaved_class) %>%
  summarise(
    across(all_of(pathos_name_apo), 
           ~ round(sum(. > 0) / n(), digits = 2)) ) |>
  mutate(combined_name = paste( treatment, broadleaved_class, sep = "_")) |>
  arrange(factor(combined_name, levels = graph_sequence)) |>
  tibble::column_to_rownames("combined_name") |>
  select(all_of(pathos_name_apo)) |>
  as.matrix() 

effectif_df <- d_apo_pathos_glm %>%
  mutate(combined_name = paste(treatment, broadleaved_class, sep = "_")) %>%
  group_by(treatment, broadleaved_class) %>%
  summarise(
    effectif = n(),
    combined_name = unique(combined_name)
  ) %>%
  arrange(factor(combined_name, levels = graph_sequence))

# Check the alignment by comparing combined names
if(!all(rownames(matrix_pathoss) == effectif_df$combined_name)) {
  stop("Mismatch between matrix_pathoss row names and effectif combined names")
}

effectif <- effectif_df %>%
  pull(effectif) 

cividis_palette <- cividis::cividis(2)

superheat::superheat(X = matrix_pathoss,
                     X.text.size = 4,
                     left.label.col = "white",
                     left.label.text.size = 3,
                     bottom.label.text.size = 2,
                     bottom.label.col = "white",
                     bottom.label.text.angle = 45,
                     heat.lim = c(0,1),
                     grid.hline = FALSE,
                     grid.vline = FALSE,
                     heat.na.col = "white",
                     yr = effectif,
                     yr.axis.name = "Headcount",
                     yr.plot.type = "bar",
                     heat.pal = cividis_palette,
                     pretty.order.cols = TRUE)




library(ComplexHeatmap)
row_ha = ComplexHeatmap::rowAnnotation(effectif = anno_barplot(
  effectif, 
  bar_width = 1, 
  baseline = 0,
  add_numbers = TRUE,
  gp = gpar(fill = cividis::cividis(2, begin = 0.2, en = 0.9)),
  width = unit(3, "cm"),
  border = FALSE
))


# Créer l'objet Heatmap
ComplexHeatmap::Heatmap(matrix_pathoss,
                        name = "Prévalence",
                        row_names_side = "left",
                        col = cividis_palette, 
                        show_row_names = TRUE, 
                        show_column_names = TRUE, 
                        cluster_columns = TRUE, 
                        show_column_dend = FALSE,
                        cluster_rows = FALSE,
                        row_order = graph_sequence,
                        right_annotation = row_ha,  
) 




# -------------BELOW is TEST GRAPH FOR EWDA POSTER -------

category_palette <- c("pine_edge" = "#8FC08C", "hedgerows" = "#A95738", "broadleaved_forest" = "#B3E2B3")


# Generate a darker version of the fill colors
darker_palette <- sapply(category_palette, darken_color, amount = 0.25)

## Neoehrlichia solo ----

# Neoehrlichia line_ype x brdld detail
pp_prevalence %>%  
  filter(pathos == "Neoehrlichia_mikurensis") %>%
  ggplot(aes(x = prevalence, y = broadleaved_class, fill = broadleaved_class)) +
  geom_boxplot(aes(fill = category,
                   color = category), 
               position = position_dodge(0.75), 
               size = 1.5, 
               width = 0.5,
               alpha = 0.7,
               outlier.shape = NA) +  
  geom_dotplot(aes(fill = category), 
               binaxis = 'x', 
               stackdir = 'center', 
               position = position_dodge(0.75), 
               dotsize = 0.5, 
               color = "black",  
               alpha = 0.8) +  
  scale_fill_manual(values = c(category_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall prevalence per line", y = "Broadeleaved_status") +
  theme_classic()  +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = interaction(broadleaved_class, category)),  
               geom = "text", 
               position = position_dodge(0.75),
               color = "black", 
               size = 3)

# Neoehrlichia line_ype x brdld sans details
pp_prevalence %>%  
  filter(pathos == "Neoehrlichia_mikurensis") %>%
  ggplot(aes(x = prevalence, y = broadleaved_class)) +
  geom_boxplot(aes(fill = category,
                   color = category), 
               position = position_dodge(0.75), 
               size = 1.5, 
               width = 0.5,
               alpha = 0.7) +  
  scale_fill_manual(values = c(category_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall prevalence per line", y = "broadleaved_class") +
  theme_classic()  +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = broadleaved_class),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)

# Neoehrlichia brdld sans details

brd_palette <- c("LB" = "orange", "HB" = "blue")
darker_palette2 <- sapply(brd_palette, darken_color, amount = 0.25)

pp_prevalence %>%  
  filter(pathos == "Neoehrlichia_mikurensis") %>%
  ggplot(aes(x = prevalence, y = broadleaved_class, fill = broadleaved_class, color = broadleaved_class)) +
  geom_boxplot(position = position_dodge(0.75), 
               size = 1.5, 
               width = 0.5,
               alpha = 0.7) +  
  scale_fill_manual(values = brd_palette) + 
  scale_color_manual(values = darker_palette2) + 
  labs(x = "Overall prevalence per line", y = "Broadeleaved_status") +
  theme_classic()  +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = broadleaved_class),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)



## Mycoplasma haemomuris solo ----
pp_prevalence %>%  
  filter(pathos == "Haemobartonella_muris") %>%
  ggplot(aes(x = prevalence, y = category, fill = category)) +
  geom_boxplot(aes(fill = category,
                   color = category), 
               position = position_dodge(0.75), 
               size = 1.2, 
               width = 0.5,
               alpha = 0.7) +  
  scale_fill_manual(values = c(category_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall prevalence per line", y = "Line type") +
  theme_classic()  +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = category),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)


## combined plot :  ----

plot1 <-pp_prevalence %>%  
  filter(pathos == "Neoehrlichia_mikurensis") %>%
  ggplot(aes(x = prevalence, y = broadleaved_class, fill = broadleaved_class, color = broadleaved_class)) +
  geom_boxplot(position = position_dodge(0.75), 
               size = 1.5, 
               width = 0.5,
               alpha = 0.7) +  
  scale_fill_manual(values = brd_palette) + 
  scale_color_manual(values = darker_palette2) + 
  labs(x = "Overall prevalence per line", y = "Broadeleaved_status") +
  theme_classic()  +
  theme(axis.title.x = element_blank(),   
        axis.text.x = element_blank(),    
        axis.ticks.x = element_blank(),
        axis.line.x.bottom = element_blank()) +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = broadleaved_class),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)

plot2 <- pp_prevalence %>%  
  filter(pathos == "Haemobartonella_muris") %>%
  ggplot(aes(x = prevalence, y = category, fill = category)) +
  geom_boxplot(aes(fill = category,
                   color = category), 
               position = position_dodge(0.75), 
               size = 1.2, 
               width = 0.5,
               alpha = 0.7) +  
  scale_fill_manual(values = c(category_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall prevalence per line", y = "Line type") +
  theme_classic()  +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = category),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)

library(patchwork)

combined_plot <- plot1 / plot2 + plot_layout(heights = c(1, 1))
combined_plot

## IN POSTER ----

plot <- pp_prevalence %>%  
  filter(pathos == "Neoehrlichia_mikurensis") %>%
  ggplot(aes(x = prevalence, y = broadleaved_class)) +
  geom_boxplot(aes(fill = category, color = category), 
               position = position_dodge(0.75), 
               size = 2.5,       
               width = 0.6,      
               alpha = 0.8,
               outlier.size = 5) +    
  scale_fill_manual(values = c(category_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall Prevalence per Line", y = "Broadleaved Status") +
  theme_classic(base_size = 20) + 
  theme(
    axis.title = element_text(size = 24, face = "bold"),  
    axis.text = element_text(size = 20),                  
    legend.title = element_text(size = 22),               
    legend.text = element_text(size = 20),                
    plot.title = element_text(size = 28, face = "bold")   
  ) +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = broadleaved_class),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 6) 
plot

ggsave(filename = here::here("figures","neoehrlichia_plot_ewda.pdf"), plot = plot, 
       width = 10.5, height = 7.5, units = "in", device = "pdf")
=