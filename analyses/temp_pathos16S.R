
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


# List of rodents identified putative pathogens
putative_pathos_genra <- c("Bartonella", "Neoehrlichia", "Mycoplasma", "Borrelia", "Borreliella", "Streptobacillus", "Ehrlichia", "Spiroplasma","Brevinema",
                           "Orientia", "Rickettsia", "Leptospira", "Yersinia", "Actinobacillus", "Treponema", "Chlamydia", "Neisseria", "Pasteurella",
                           "Francisella", "Brucella", "Coxiella")
putative_pathos_family <- c("Sarcocystidae")


# GENERATE DATA ----

## Import data ----
# Import hosts and line modalities file
d_host <- readr::read_csv( here::here("data/", "raw-data/", "host_data", "20250123_bpm_modalities.csv") )

# Import 16S filtered file
file16s_run00_01_04_05 <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01-04-05","Run00-01-04-05_16S_filtered-merged_postfrogs.txt"))

# Import rpoB filtered file
filerpoB_run01_05 <- data.table::fread(file = here::here( "data","raw-data/","rpoB_run01-05","Run01-05-rpoB_filtered-merged.txt"))

# Import rodent macroparasite
d_macroparasite <- data.table::fread(file = here::here("data", "derived-data", "ticks", "rodents_tick", "20240731_macroparasite.csv") )

# Import lipl32 file
filelipl32 <- readxl::read_excel(path = here::here( "data","raw-data/","lepto","20250116_data_qPCR_lepto_BePrep.xlsx"))


## Generate a file containing samples of interest ----

# Extract id of small mammals caught and dissected in Beprep
sm_id <- unique(d_host %>% 
                  filter(!is.na(numero_centre))%>%
                  filter(stringr::str_detect(numero_centre, pattern = "NCHA")) %>%
                  pull(numero_centre)
)

# Identify taxonomy 16S and rpoB columns
taxo_name_16s <- colnames(file16s_run00_01_04_05 [, 1:which(colnames(file16s_run00_01_04_05) == "observation_sum") ])

taxo_name_rpoB <- colnames(filerpoB_run01_05 [, 1:which(colnames(filerpoB_run01_05) == "observation_sum") ])

# Calculate total run read number per cluster (after our filtering)
file16s_run00_01_04_05 <- file16s_run00_01_04_05 %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name_16s))) %>%
  relocate(clean_totalrunreads, .after = names(file16s_run00_01_04_05)[which(colnames(file16s_run00_01_04_05) == "observation_sum")])

filerpoB_run01_05 <- filerpoB_run01_05 %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name_rpoB))) %>%
  relocate(clean_totalrunreads, .after = names(filerpoB_run01_05)[which(colnames(filerpoB_run01_05) == "observation_sum")])

# Reattribute taxonomy columns
taxo_name_16s <- colnames( file16s_run00_01_04_05 [, 1:which(colnames(file16s_run00_01_04_05) == "clean_totalrunreads") ])

taxo_name_rpoB <- colnames(filerpoB_run01_05 [, 1:which(colnames(filerpoB_run01_05) == "clean_totalrunreads") ])

# Take away problematic samples - if needed
#file16s_run00_01_04_05 <- file16s_run00_01_04_05 |>

#filerpoB_run01_05 <- filerpoB_run01_05 |>

# Generate new file containing only Spleen BePrep's samples
#by selecting those samples
beprep_sample_16s_sp <- file16s_run00_01_04_05 |>
  select( (contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
beprep_16s_sp <- file16s_run00_01_04_05 |>
  select(taxo_name_16s) |>
  cbind(beprep_sample_16s_sp)

#by selecting those samples
beprep_sample_rpoB_sp <- filerpoB_run01_05 |>
  select( (contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
beprep_rpoB_sp <- filerpoB_run01_05 |>
  select(taxo_name_rpoB) |>
  cbind(beprep_sample_rpoB_sp)

# Delete the PCZymo positive in 16S after visual control
beprep_16s_sp <- beprep_16s_sp |>
  select(!contains("PCzymo"))

# Verify absence of PCZymo positive in rpoB data
beprep_rpoB_sp |>
  select(contains("PCzymo"))

# Calculate total number of reads for each cluster considering chosen samples
beprep_16s_sp <- beprep_16s_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name_16s)))%>%
  relocate(totalreads, .after = names(beprep_16s_sp)[which(colnames(file16s_run00_01_04_05) == "clean_totalrunreads")])

beprep_rpoB_sp <- beprep_rpoB_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name_rpoB)))%>%
  relocate(totalreads, .after = names(beprep_rpoB_sp)[which(colnames(filerpoB_run01_05) == "clean_totalrunreads")])

# Order rows by new total read number
beprep_16s_sp <- data.table::setorder(beprep_16s_sp, -totalreads)

beprep_rpoB_sp <- data.table::setorder(beprep_rpoB_sp, -totalreads)

# Delete cluster with no reads for the selected individuals
beprep_16s_sp <- beprep_16s_sp |> 
  filter(totalreads >0)

beprep_rpoB_sp <- beprep_rpoB_sp %>% 
  filter(totalreads >0)

# Write a file containing Spleen OTUs containing our samples only :
data.table::fwrite(beprep_16s_sp, here::here("data", "derived-data","16S", "20250123_beprep_16s_sp.txt") )

data.table::fwrite(beprep_rpoB_sp, here::here("data", "derived-data","rpoB", "20250123_beprep_rpoB_sp.txt") )


## Generate spleen putative pathogen data ----

# Split the taxonomy to access taxon
pathos_16s <- beprep_16s_sp |>
  tidyr::separate_wider_delim(cols = blast_taxonomy,
                              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                              delim = ";")

pathos_rpoB <- beprep_rpoB_sp |>
  tidyr::separate_wider_delim(cols = blast_taxonomy,
                              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                              delim = ";")

# Delete cluster without real affiliation
pathos_16s <- pathos_16s %>%
  filter(!(domain == "no" | domain == "Plantae"))

pathos_rpoB <- pathos_rpoB %>%
  filter(!domain == "no")

# FOR 16S Keep only OTUs with known putative pathogens (you should look at raw data before this step, ensure it is ok and you're not missing pathogens taxa as the list may not be complete)
pathos_16s <- pathos_16s |>
  filter(genus %in% putative_pathos_genra 
         | family %in% putative_pathos_family)

# Write a file containing 16S Spleen OTUs containing our samples only AND only putative pathogens :
data.table::fwrite(pathos_16s, here::here("data/", "derived-data/","16S", "20250123_beprep_16s_sp_pathos.txt") )

data.table::fwrite(pathos_rpoB, here::here("data", "derived-data","rpoB", "20250123_beprep_rpoB_sp_pathos.txt") )


## Seek for species level identification for specific taxa (e.g. Mycoplasma) throught external manipulation (e.g phylogenic tree and Genbank) ----
# A new datafile is created BY HAND to add the checked species information (in new column "species_adjusted") and then imported
# If probable pseudo-gene affiliation is observed in data file they may also get deleted ny operator at this step
pathos_16s_checked <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01-04-05","20241204_beprep_16s_sp_pathos-adjusted.txt"))
pathos_16s_checked$species_adjusted[pathos_16s_checked$species_adjusted == ""] <- NA

pathos_rpoB_checked <- data.table::fread(file = here::here( "data","raw-data/","rpoB_run01-05","20241204_beprep_rpoB_sp_pathos-adjusted.txt"))
pathos_rpoB_checked$species_adjusted[pathos_rpoB_checked$species_adjusted == ""] <- NA


# Change species value when there is a species adjusted value (generated by operator manually ) and drop the species adjusted column
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
  select(!species_adjusted)

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
  select(!species_adjusted)

# Rename inappropriate affiliation 
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(species = if_else(species == "Haemobartonella muris", "Mycoplasma haemomuris", species))

# FOR RPOB : reads based filter
pathos_rpoB_checked <- pathos_rpoB_checked %>%
  filter(totalreads >= 100)


## Create new column for best identified taxa ----

# If there are unknown or Multi-affiliation in species column, change it to NA
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species ))

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species )) 

# Creation of a new column containing best available identity (between species or genus here)
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

# Replace space in taxa name 
pathos_16s_checked <- pathos_16s_checked %>%
  mutate(identification16s = gsub(" ", "_", identification16s ))

pathos_rpoB_checked <- pathos_rpoB_checked %>%
  mutate(identificationrpoB = gsub(" ", "_", identificationrpoB ))


## Simplification, compress spleen pathogens information by best identified taxa ----

# Keep only best taxa and samples
taxa_pathos_16s <- pathos_16s_checked %>%
  select(c("identification16s", names(pathos_16s_checked)[(which(names(pathos_16s_checked) == "totalreads") + 1):ncol(pathos_16s_checked)]))

taxa_pathos_rpoB <- pathos_rpoB_checked %>%
  select(c("identificationrpoB", names(beprep_rpoB_sp)[(which(names(beprep_rpoB_sp) == "totalreads") + 1):ncol(beprep_rpoB_sp)]))


# Sum rows by best taxa
taxa_pathos_16s <- taxa_pathos_16s %>%
  group_by(identification16s) %>%
  summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop")

taxa_pathos_rpoB <- taxa_pathos_rpoB %>%
  group_by(identificationrpoB) %>%
  summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop")

# Transpose the table
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

# Finish sample name transformation to join to host data
taxa_transposed_16S <- taxa_transposed_16S |>
  mutate(numero_centre_16s = gsub(".SP", "", numero_centre_16s ))

taxa_transposed_rpoB <- taxa_transposed_rpoB |>
  mutate(numero_centre_rpoB = gsub(".SP", "", numero_centre_rpoB ))


## Join 16S and rpoB  data ----
taxa_transposed_16S_rpoB <- left_join(taxa_transposed_16S, 
                                      taxa_transposed_rpoB,
                                      by = c("numero_centre_16s" = "numero_centre_rpoB")) %>%
  rename(numero_centre_combined = numero_centre_16s)

# Replace NA by 0 in new data
taxa_transposed_16S_rpoB <- taxa_transposed_16S_rpoB %>%
  mutate(across(-numero_centre_combined, ~ if_else(is.na(.), 0, .)))

# See benefit of Neoehrlichia_mikurensis columns combining
taxa_transposed_16S_rpoB %>%
  filter(Neoehrlichia_mikurensis.y > 0) %>%
  filter(Neoehrlichia_mikurensis.x == 0) %>%
  select(numero_centre_combined) %>%
  pull()

# Combine Neoehrlichia_mikurensis columns 
taxa_transposed_16S_rpoB <- taxa_transposed_16S_rpoB %>%
  mutate(Neoehrlichia_mikurensis = Neoehrlichia_mikurensis.x + Neoehrlichia_mikurensis.y) %>%
  select(- c("Neoehrlichia_mikurensis.x", "Neoehrlichia_mikurensis.y") ) 

# Control rpoB metabarcoding - identify failed sample
taxa_transposed_16S_rpoB %>%
  filter(rowSums(across( names(taxa_transposed_rpoB[,-1]))) == 0) %>%
  filter(Bartonella > 0 ) %>%
  select(numero_centre_combined) %>%
  pull()
  


# jsp ce que je voulais faire en bas, mais renommer numero_centre peut-être pas mal (genre rpoB-16s?)
# genre controle quelles echantillons de l'un passé dans l'autre (tous les barto + 16s passé rpoB?)
# peut-être faisable avant ça deja ?



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
# rodent_pathos <- merge(x = d_macroparasite, y = taxa_transposed_16S_rpoB, by.x = "numero_centre",  by.y = "numero_centre_combined")


# # Extract Pathogens name
# pathos_name <- rodent_pathos |>
#   select(names(rodent_pathos)[(which(names(rodent_pathos) == "effectif_tick") + 1):ncol(rodent_pathos)]) %>%
#   colnames()
# pathos_name



# Temporary solution : -attention-------------
rodent_pathos <- left_join(d_host %>%
                             filter(stringr::str_detect(numero_centre, pattern = "NCHA")),
                           taxa_transposed_16S_rpoB,
                           by = c("numero_centre" = "numero_centre_combined"))

rodent_pathos <- rodent_pathos %>%
  filter(stringr::str_detect(numero_centre, pattern = "NCHA"))

# Vector containing pathogen names
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
candidate_richness_pathos <- setdiff(pathos_name_apo, "Bartonella")

data_for_m <- data_for_m %>%
  mutate(number_pathos = rowSums(across(all_of(candidate_richness_pathos))) ) 

# Identify pathos with GLOBAL prevalence >=10 for whole data
patho10_apo <- d_apo_pathos_glm %>%
  summarise( across(all_of(pathos_name_apo), ~ sum(. > 0) / n() ))  %>%
  unlist()
patho10_apo <- names(patho10_apo[patho10_apo >= 0.10])
patho10_apo


# Take away rows with NA for some of our factor
data_for_m <- data_for_m %>%
  filter(!is.na(sexe))

# Generate dataset without broadleaved forest analysis
data_for_m_noforests <- data_for_m %>%
  filter(category != "broadleaved_forest")

# Generate dataset without broadleaved forest with lines containing PCA values
data_for_m_noforests_pca <- data_for_m_noforests %>%
  filter(!is.na(PCA_axis1))

  
# Generate other dataset without broadleaved forest, without year that do not contains them
data_for_m_forestsyear <- data_for_m %>%
  filter(code_mission %in% (data_for_m %>% 
           filter(category == "broadleaved_forest") %>%
           select(code_mission) %>% 
           pull() %>%
           unique()))


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

table(data_for_m_noforests$Neoehrlichia_mikurensis, data_for_m_noforests$season)

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
  filter(pathos == "Mycoplasma_haemomuris") %>%
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

