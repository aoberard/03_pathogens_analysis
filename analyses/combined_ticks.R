# Script parameters ----

set.seed(123)

## Library ----
library("data.table")
library("dplyr")
library("ggplot2")


## Filtering parameters ----
# mission for which tick were passed individually in microfluidic test
microfluidic_mission <- c("Juin 2023")

# mission for which macroparasite data were written in BPM
macroparasite_mission <- c("Juin 2023", "Octobre 2023")

# Species of host which interest us
microfluidic_hostspecies <- c("Apodemus sylvaticus")


## Graphical parameters ----
type_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
type_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
treatment_order <- c("CT-LB", "CT-HB", "NC-LB", "NC-HB", "C-LB", "C-HB", "B" )
mission_order <- c("Juin 2023", "Octobre 2023", "Juin 2024")
mission_color <- c("Juin 2023" = "#66c2a5", "Octobre 2023" = "#fc8d62", "Juin 2024" = "#8da0cb")

#new g parameters
category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
category_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
type_order <- c("CT_LB", "CT_HB", "NC_LB", "NC_HB", "C_LB", "C_HB", "B" )
type_palette <- c("CT_LB" = "#FFB3BA", "CT_HB" = "#FFB3BA" , "NC_LB" = "#FFDFBA", "NC_HB" = "#FFDFBA" , "C_LB" = "#FFDFBA" , "C_HB"= "#FFDFBA" , "B" = "#B3E2B3")
brd_palette <- c("LB" = "#FFB3BA", "HB" = "#B3E2B3")
mission_order <- c("Juin_2023", "Octobre_2023", "Juin_2024", "Septembre_2024")
mission_color <- c("Juin_2023" = "#66c2a5", "Octobre_2023" = "#fc8d62", "Juin_2024" = "#8da0cb", "Septembre_2024" = "#ab7a82")


## Functions ----
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


# Import data ----

# Import hosts and line modalities file
hosts <- readr::read_csv( here::here("data", "raw-data", "host_data", "20250123_bpm_modalities.csv") )

# Import macroparasites from hosts file
macroparasite <- fread(file = here::here("data/", "derived-data/", "ticks", "rodents_tick", "20240731_macroparasite.csv") )

# Import ticks collected from the environment data
envticks <- readxl::read_excel(here::here("data", "raw-data", "raw-ticks", "collect", "20240730_collect_tick.xlsx"))

# Import microfluidic information
microfluidic <- data.table::fread(file = here::here("data", "raw-data", "raw-ticks", "microfluidic", "20240411_microfluidic_tick.txt") )



# Data management  ----

## Host data treatment ----

# Take away not dissected individuals from hosts data
hosts <- hosts %>%
  filter(stringr::str_detect(numero_centre, pattern = "NCHA"))

# Convert columns to character type (if necessary) and filter hosts data
hosts_processed <- hosts %>%
  filter(code_resultat == 1) %>%
  select(numero_centre, taxon_mamm, numero_ligne, code_mission, year, season, connectivity, broadleaved_status, line_treatment, line_type) %>%
  mutate(across(everything(), as.character))


## Macroparasite - host data ----

# Convert columns to character type (if necessary) for macroparasite data
macroparasite_processed <- macroparasite %>%
  mutate(across(everything(), as.character))

# Perform the full join
host_macro <- full_join(macroparasite_processed, hosts_processed)

# Transform NA effectif_tick into 0, and make it numeric again
host_macro$effectif_tick[ is.na(host_macro$effectif_tick)] <- 0
host_macro$effectif_tick <- as.numeric(host_macro$effectif_tick)

# Generate specific data for macroparasite analysis by taking away rows for each data were not written in BPM
host_macro_uptodate <- host_macro %>%
  filter(code_mission %in% macroparasite_mission)


## Environment ticks ----

# Join environment ticks to line modalities from hosts data
envticks_processed <- envticks %>%
  left_join(hosts %>%
              select(numero_ligne, connectivity, broadleaved_status, line_treatment, line_type) %>%
              mutate(numero_ligne = as.character(numero_ligne)) %>%
              distinct() )


# Generate year and month variable
envticks_processed <- envticks_processed %>%
  mutate(
    year = stringr::str_extract(string = code_mission, pattern = "\\d*$"),
    month = stringr::str_extract(string = code_mission, pattern = "(?<=\\s-\\sBePrep\\s-\\s)\\w+")
  )

#choose which month correspond to each season
Spring <- c("Juin")
Autumn <- c("Octobre")

envticks_processed <- envticks_processed %>%
  mutate(
    season = case_when(
      month %in% Spring ~ "Spring",
      month %in% Autumn ~ "Autumn",
      TRUE ~ "unknown"
    )
  ) %>%
  relocate(
    year, .after = which(names(envticks_processed) == "code_mission") 
  ) %>%
  relocate(
    season, .after = which(names(envticks_processed) == "code_mission") +1
  ) %>%
  mutate(
    code_mission = paste(month, year)
  ) %>%
  select(!month)

# Error check for modality attribution
unknown_rows <- envticks_processed %>% filter(is.na(line_treatment) | line_type == "unknown" | season == "unknown" )
if (nrow(unknown_rows) > 0) {
  warning("There are rows with unknown values. Check 'unknown_rows' for details.")
}


## Environment and hosts tick combination ----

# Add a new column to each table to indicate the source of ticks
host_macro_combined <- host_macro %>%
  mutate(source = "rodent")

envticks_processed_combined <- envticks_processed %>%
  mutate(source = "environment")

# Combine the tables into one data frame
combined_table <- dplyr::bind_rows(host_macro_combined, envticks_processed_combined)


## Microfluidic data ----

# Keep only data for ticks found on rodents (no human)
microfluidic_processed <- microfluidic %>%
  dplyr::filter(sample_nature == "tick") %>%
  dplyr::filter(host_nature =="small_mammal")
  
# Control if the positive control is positive
failed_samples <- microfluidic_processed %>%
  filter(`E..coli.(temoin.PCR)` == 0 | is.na(`E..coli.(temoin.PCR)`))

if (nrow(failed_samples) > 0) {
  message(paste0("There are ", nrow(failed_samples), 
                 " samples with apparent failed PCR. They are taken away from data"))
  
  microfluidic_processed <- microfluidic_processed %>%
    filter(`E..coli.(temoin.PCR)` != 0 | !is.na(`E..coli.(temoin.PCR)`) ) %>%
    select(-`E..coli.(temoin.PCR)`)
} else {
  microfluidic_processed <- microfluidic_processed %>%
    select(-`E..coli.(temoin.PCR)`)
  message("Every sample has its PCR control positive")
}

# Generate a tick species column
microfluidic_processed <- microfluidic_processed %>%
  dplyr::mutate(tick_species = dplyr::case_when(
    sample_nature == "tick" & Dermacentor.reticulatus_ITS2 == 1 & Ixodes.ricinus_ITS2 == 0 & Dermacentor.marginatus_ITS2 == 0 ~ "Dermacentor_reticulatus",
    sample_nature == "tick" & Dermacentor.reticulatus_ITS2 == 0 & Ixodes.ricinus_ITS2 == 1 & Dermacentor.marginatus_ITS2 == 0 ~ "Ixodes_ricinus",
    sample_nature == "tick" & Dermacentor.reticulatus_ITS2 == 0 & Ixodes.ricinus_ITS2 == 0 & Dermacentor.marginatus_ITS2 == 1 ~ "Dermacentor_marginatus",
    sample_nature == "tick" & Dermacentor.reticulatus_ITS2 == 0 & Ixodes.ricinus_ITS2 == 0 & Dermacentor.marginatus_ITS2 == 0 & Tique.spp_16S == 1 ~ "unknown_tick",
    sample_nature == "tick" & rowSums(select(., c("Dermacentor.reticulatus_ITS2", "Ixodes.ricinus_ITS2", "Dermacentor.marginatus_ITS2"))) >1 ~ "multiple affiliation",
    TRUE ~ "problem"
  )) 

# Check presence of multiple affiliations
multiple_affiliations <- microfluidic_processed %>%
  filter(tick_species == "multiple affiliation")

if (nrow(multiple_affiliations) > 0) {
  message(paste0("There are ", nrow(multiple_affiliations), " samples with multiple affiliations. They will be converted to 'unknown_tick'"))
  
    microfluidic_processed <- microfluidic_processed %>%
    mutate(tick_species = if_else(tick_species == "multiple affiliation", "unknown_tick", tick_species))
} else {
  message("No samples with multiple affiliations were found")
}

# Check presence of problem tick species
problem_entries <- microfluidic_processed %>%
  filter(tick_species == "problem")

if (nrow(problem_entries) > 0) {
  message(paste0("There are ", nrow(problem_entries), " samples with 'problem' species identification"))
} else {
  message("No samples with 'problem' species identification were found.")
}

# Check presence of unknown tick species
unknown_tick_entries <- microfluidic_processed %>%
  filter(tick_species == "unknown_tick")

if (nrow(unknown_tick_entries) > 0) {
  message(paste0("There are ", nrow(unknown_tick_entries), " samples with 'unknown_tick' species identification."))
} else {
  message("No samples with 'unknown_tick' species identification were found.")
}

# Drop old tick id columns
microfluidic_processed <- microfluidic_processed %>%
  select(-c("Ixodes.ricinus_ITS2","Tique.spp_16S","Dermacentor.reticulatus_ITS2","Dermacentor.marginatus_ITS2"))

# Delete column without any detection
microfluidic_processed <- microfluidic_processed %>%
  select(where(~ any(. != 0)))

# Due to potential symbiont amplification, delete some pathogen column
microfluidic_processed <- microfluidic_processed %>%
  select(- c("Rickettsia.spp._gltA", "Francisella.tularensis_fopA." ) ) 


### Microfluidic prevalence data ----

# Identify column of pathogen (that were not purposely taken away and that are not empty)
microfluidic_pathos_name <- microfluidic_processed %>%
  select(where(is.integer)) %>%
  colnames()

# Join microfluidic to host data through host_macro
microfluidic_prevalence <- dplyr::left_join(microfluidic_processed, host_macro, by = c("host_name" = "numero_centre"))

# Here the unknown ticks and problem ticks are kept for now !!!



### Microfluidic community data ----

# Delete rows with "problem" or "unknown" in tick_species
microfluidic_community <- microfluidic_processed %>%
  filter(tick_species != "problem") %>%
  filter(tick_species != "unknown_tick") 

# Generate community data for hosts of tick in microfluidic data
microfluidic_community <- microfluidic_community %>%
  group_by(host_name, tick_species) %>%
  summarise(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = tick_species, values_from = count, values_fill = list(count = 0))

# Join to host macro_file for the pertinent subset of individuals
microfluidic_community <- left_join(x = host_macro %>%
                                      filter(code_mission %in% microfluidic_mission) %>%
                                      filter(taxon_mamm %in% microfluidic_hostspecies),
                                    y = microfluidic_community, by = c("numero_centre"= "host_name" ))
                             
# Remplace NA by 0 for individuals without ticks
microfluidic_community <- microfluidic_community %>%
  mutate(across(c(Dermacentor_reticulatus, Dermacentor_marginatus, Ixodes_ricinus),
               ~ ifelse(is.na(.), 0, .)))
  
# Calculate tick species richness
microfluidic_community <- microfluidic_community %>%
  mutate(tick_richness = rowSums(across(c(Dermacentor_reticulatus, Dermacentor_marginatus, Ixodes_ricinus), ~ . > 0))) %>%
  relocate(tick_richness, .after = Dermacentor_reticulatus)

# Take off not parasited individuals for community matrix
microfluidic_matrix_commu <- microfluidic_community %>%
  rowwise() %>%
  filter(sum(c_across(c(Dermacentor_reticulatus, Dermacentor_marginatus, Ixodes_ricinus))) > 0) %>%
  ungroup()

# Extract community matrix for parasited samples
microfluidic_matrix_commu <- microfluidic_matrix_commu %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(c(Dermacentor_reticulatus, Dermacentor_marginatus, Ixodes_ricinus) ) %>%
  as.matrix()
  
# Calculate beta diversity matrix
#Bray-Curtis : 
microfluidic_matrix_bray <- microfluidic_matrix_commu %>%
  vegan::vegdist(, method="bray", binary = FALSE)





# Analysis ----

## Tick number on host ----
# Beware that the data used take only the samples for wich data were put in BPM (go see parameters)
# Beware of different host species
# This places is to study impact of things on tick number on host (not occurence)

### Exploration
host_macro_uptodate %>%
  group_by(taxon_mamm) %>%
  summarise(
    effectif = n(),
    mean_tick = mean(effectif_tick))

hist(host_macro_uptodate$effectif_tick)


### GLMMs

# Making data frame for glm model /!\ spring 2023 only for poster EWDA
dglm_host_effectif_tick <- host_macro_uptodate %>%
  filter(taxon_mamm %in% c("Apodemus sylvaticus")) 

# GLM making (poisson) 
glm_host_effectif_tick <- lme4::glmer.nb(
  formula = effectif_tick ~ broadleaved_status * connectivity + (1|numero_ligne),
  data = dglm_host_effectif_tick,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(glm_host_effectif_tick, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

glm_host_effectif_tick |>
  RVAideMemoire::overdisp.glmer()

SelectionModels<- MuMIn::dredge(glm_host_effectif_tick, rank = "AICc")              
TopModels<-subset(SelectionModels, delta < 2)
TopModels

glm_host_effectif_tick <- lme4::glmer.nb(
  formula = effectif_tick ~ broadleaved_status + (1|numero_ligne),
  data = dglm_host_effectif_tick,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(glm_host_effectif_tick, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(glm_host_effectif_tick,.~.,test="Chisq")
summary(glm_host_effectif_tick)

em <- emmeans::emmeans(glm_host_effectif_tick, "broadleaved_status", type = "response")
em
plot(em, comparisons = TRUE)

emmeans::contrast(em, "pairwise", adjust = "Tukey")


## Environment ticks ----

### Exploration 

envticks_processed %>%
  group_by(code_mission) %>%
  summarise(
    effectif_tick = sum(effectif_tick)
  )

envticks_processed %>%
  group_by(code_mission, line_type) %>%
  summarise(
    effectif_tick = sum(effectif_tick)
  )

hist(envticks_processed$effectif_tick)

envticks_processed %>%
  ggplot(aes (x = factor(line_type, levels = type_order), y = effectif_tick, , fill = line_type)) + 
  facet_grid(~ factor(code_mission, levels = mission_order) ) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_jitter(width = 0.05), dotsize = 0.5) +
  scale_fill_manual(values = type_palette) +
  labs(x = "Type de ligne", y = "Nbr de tiques collectées au drap par ligne") +
  guides(fill = guide_legend(title = "Type de ligne")) +
  theme_minimal()


### GLMMs 

# Making data frame for glm model /!\ spring 2023 only for poster EWDA -normally should put interaction between type and code_mission
dglm_environment_tick <- envticks_processed %>%
  filter(code_mission == "Juin 2023")

# GLM making (poisson) 
rm(glm_environment_tick)
glm_environment_tick <- lme4::glmer(
  formula = effectif_tick ~ line_type + (1|numero_ligne),
  data = dglm_environment_tick,
  family = poisson(link = "log"),  
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)

DHARMa::simulateResiduals(glm_environment_tick, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

glm_environment_tick |>
  RVAideMemoire::overdisp.glmer()

SelectionModels<- MuMIn::dredge(glm_environment_tick, rank = "AICc")              
TopModels<-subset(SelectionModels, delta < 2)
TopModels

rm(glm_environment_tick)
glm_environment_tick <- lme4::glmer(
  formula = effectif_tick ~ line_type + (1|numero_ligne),
  data = dglm_environment_tick,
  family = poisson(link = "log"),  
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)

DHARMa::simulateResiduals(glm_environment_tick, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(glm_environment_tick,.~.,test="Chisq")
summary(glm_environment_tick)

em <- emmeans::emmeans(glm_environment_tick, "line_type", type = "response")
em
plot(em, comparisons = TRUE)

emmeans::contrast(em, "pairwise", adjust = "Tukey")


## Environment and hosts tick combination ----

# Graphical view next to each other of tick on hosts and caught in environment from combined table
ggplot(combined_table, aes(x = factor(line_type, levels = type_order), y = effectif_tick, fill = line_type) ) +
  facet_grid(source ~ factor(code_mission, levels = mission_order)) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  stat_summary(fun.data = count_summary, geom = "text", vjust = 0) + 
  scale_fill_manual(values = type_palette) +
  labs(title = "Boxplot of Tick Numbers rodents and the environment",
       x = "Linetype",
       y = "Tick Number") +
  theme_minimal()

## Microfluidic pathogens prevalence ----

### Exploration 

# GLOBAL PREVALENCE
# Calculate GLOBAL prevalence of left pathogens
prevalence_global_microfluidic <- microfluidic_prevalence %>%
  summarise( across(all_of(microfluidic_pathos_name), ~ sum(. > 0) / n() ))  %>%
  unlist() 

# Pathogens with more than 10% prevalence GLOBALLY
names(prevalence_global_microfluidic[prevalence_global_microfluidic >= 0.10])


# PREVALENCE PER TICK SPECIES
# Calculate prevalence of left pathogens per TICK SPECIES
prevalence_per_species_microfluidic <- microfluidic_prevalence %>%
  group_by(tick_species) %>%
  summarise(across(all_of(microfluidic_pathos_name), ~ sum(. > 0) / n()), .groups = "drop")

# Pathogens with more than 10% prevalence per tick species
prevalence_per_species_microfluidic %>%
  tidyr::pivot_longer(cols = -tick_species, names_to = "pathogen", values_to = "prevalence") %>%
  filter(prevalence >= 0.10)

# Ticks per mission in data
microfluidic_prevalence %>%
  count(mission)

# Ticks stages and species in data
microfluidic_prevalence %>%
  count(tick_species, stage)


### GLMMs 

dglm_microfluidic_pathos <- microfluidic_prevalence 

#### Model : Borrelia.spp._23S
glm_microfluidic_borrelia <- lme4::glmer(
  formula = Borrelia.spp._23S ~ broadleaved_status * connectivity + tick_species + (1|numero_ligne) + (1|host_name),
  family = binomial(link = "logit"),
  data = dglm_microfluidic_pathos,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

SelectionModels<- MuMIn::dredge(glm_microfluidic_borrelia, rank = "AICc")              
TopModels <- subset(SelectionModels, delta < 2)
TopModels


#### Model : Bartonella.spp._ssrA
glm_microfluidic_bartonella <- lme4::glmer(
  formula = Bartonella.spp._ssrA ~ broadleaved_status * connectivity  +  (1|numero_ligne)+ (1|host_name),
  family = binomial(link = "logit"),
  data = dglm_microfluidic_pathos,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

SelectionModels<- MuMIn::dredge(glm_microfluidic_bartonella, rank = "AICc")              
TopModels <- subset(SelectionModels, delta < 2)
TopModels



## Microfluidic ticks community ----

### Alpha diversity ----

#### Richness
##### Exploration
hist(microfluidic_community$tick_richness)

microfluidic_community %>%
  ggplot(aes(x = factor(line_type, levels = type_order), y = tick_richness, fill = line_type)) +
  facet_grid(~ factor(code_mission, levels = mission_order)) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  stat_summary(fun.data = count_summary, geom = "text", vjust = 0) + 
  scale_fill_manual(values = type_palette) +
  labs(x = "Type de ligne", y = "Tick richness") +
  guides(fill = guide_legend(title = "Type de ligne")) +
  theme_minimal()

##### GLMMs
glm_microfluidic_tick_rich <- lme4::glmer(
  formula = tick_richness ~ connectivity * broadleaved_status + (1|numero_ligne),
  data = microfluidic_community,
  family = poisson(link = "log"),  
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)

DHARMa::simulateResiduals(glm_microfluidic_tick_rich, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

SelectionModels<- MuMIn::dredge(glm_microfluidic_tick_rich, rank = "AICc")              
TopModels<-subset(SelectionModels, delta < 2)
TopModels


### Beta diversity ----

#### Bray-Curtis ----

#### NMDS BRAY CURTIS
# Calculate NMDS
nmds_bray <- vegan::metaMDS(
  comm = microfluidic_matrix_commu,
  distance = "bray",
  k = 2,
  autotransform = FALSE
)
cat(paste0("Final NMDS has a stress of ", round(nmds_bray$stress, 3), "\n"))
vegan::stressplot(nmds_bray)
# Extract scores
nmdspoint <- vegan::scores(nmds_bray)$sites %>%
  as_tibble(rownames = "numero_centre")
nmdsvariable <- vegan::scores(nmds_bray)$species %>%
  as_tibble(rownames = "Species")

nmdspoint %>% 
  left_join(microfluidic_community) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = as.factor(line_type))) +
  geom_point() +
  stat_ellipse(show.legend = FALSE) +
  geom_text(data = nmdsvariable, aes(x = NMDS1, y = NMDS2, label = Species), colour = "grey20") +
  theme_minimal()

#### PERMANOVA Bray-Curtis
# Check and align data without modifying the original data frame
if (class(microfluidic_matrix_bray) == "dist") {
  
  # Extract row names of the distance matrix
  matrix_row_names <- rownames(as.matrix(microfluidic_matrix_bray))
  
  # Ensure the 'numero_centre' column in microfluidic_community matches the row names
  if (all(matrix_row_names %in% microfluidic_community$numero_centre)) {
    
    # Filter and arrange data frame to ensure it matches the order of distance matrix row names
    filtered_community <- microfluidic_community %>%
      filter(numero_centre %in% matrix_row_names) %>%
      arrange(match(numero_centre, matrix_row_names))
    
    # Check if the row order in the filtered data frame matches the distance matrix
    if (all(matrix_row_names == filtered_community$numero_centre)) {
      set.seed(123)
      
      # PERMANOVA
      permanova_bray <- vegan::adonis2(microfluidic_matrix_bray ~ line_type,
              data = filtered_community,
              permutations = 9999,
              na.action = "na.fail")
      
      # PERMDISP 
      permdisp_bray <- vegan::betadisper(microfluidic_matrix_bray, 
                        filtered_community$line_type,
                        add = TRUE) %>% # To avoid negative eigenvalues
        vegan::permutest( permutations = 9999)
      
    } else {
      stop("The rows in the distance matrix and the filtered data frame are not in the same order.")
    }
    
  } else {
    stop("Some row names in the distance matrix are not present in the 'numero_centre' column of the data frame.")
  }
  
} else {
  stop("The provided matrix is not a distance matrix.")
}

permanova_bray
permdisp_bray































# -------------BELOW is TEST GRAPH FOR EWDA POSTER -------
type_palette <- c("pine_edge" = "#8FC08C", "hedgerows" = "#A95738", "broadleaved_forest" = "#B3E2B3")


# Generate a darker version of the fill colors
darker_palette <- sapply(type_palette, darken_color, amount = 0.25)


host_macro %>%
  filter(code_mission %in% c("Juin 2023", "Octobre 2023") )  %>%
  ggplot( aes (x = factor(line_type, levels = type_order), y = effectif_tick, fill = line_type)) +
  facet_grid(~ factor(code_mission, levels = mission_order) ) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_jitter(width = 0.0001), dotsize = 0.5) +
  scale_fill_manual(values = type_palette) +
  labs(x = "Type de ligne", y = "Nbr de tiques collectées sur les hôtes") +
  guides(fill = guide_legend(title = "Type de ligne")) +
  theme_minimal()

combined_table %>%
  filter(code_mission == "Juin 2023") %>%
  ggplot(aes(x = effectif_tick, y = factor(line_type, levels = type_order), fill = line_type) ) +
  facet_grid(source ~ factor(code_mission, levels = mission_order)) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'x', stackdir = 'center', dotsize = 0.5) +
  scale_fill_manual(values = type_palette) +
  labs(title = "Boxplot of tick numbers from rodents and the environment",
       x = "Tick Number",
       y = "Linetype") +
  theme_minimal()+
  stat_summary(fun.data = function(x) count_summary(x, y_position = 25), 
               aes(group = line_type),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)


combined_table %>%
  filter(code_mission %in% c("Juin 2023") )%>%
  ggplot(aes(y = effectif_tick, x = factor(line_type, levels = type_order), fill = line_type)) +
  facet_grid(~ source ) +
  geom_boxplot(aes(fill = line_type,
                   color = line_type), 
               position = position_dodge(0.75), 
               size = 1.2, 
               width = 0.5,
               alpha = 0.7) +  
  scale_fill_manual(values = c(type_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Type de ligne", y = "Nbr de tiques collectées au drap par ligne") +
  guides(fill = guide_legend(title = "Type de ligne")) +
  theme_minimal()  +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 25), 
               aes(group = line_type),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 3)


(code_mission %in% c("Juin 2023")) %>%
  ggplot(aes(x = effectif_tick, y = factor(line_type, levels = type_order))) + 
  geom_boxplot(aes(fill = line_type, color = line_type), 
               position = position_dodge(0.75), 
               size = 2.5,       
               width = 0.6,      
               alpha = 0.8,
               outlier.size = 5) + 
  scale_fill_manual(values = c(type_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) +  
  labs(x = "Type de ligne", y = "Nbr de tiques collectées au drap par ligne") +
  theme_classic(base_size = 20) + 
  theme(
    axis.title = element_text(size = 24, face = "bold"),  
    axis.text = element_text(size = 20),                  
    legend.title = element_text(size = 22),               
    legend.text = element_text(size = 20),                
    plot.title = element_text(size = 28, face = "bold")   
  ) +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 26), 
               aes(group = broadleaved_status),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 6) 




## IN POSTER ----
plot <- envticks_processed %>%
  filter(code_mission %in% c("Juin 2023")) %>%
  ggplot(aes(x = effectif_tick, y = factor(line_type))) + 
  geom_boxplot(aes(fill = line_type, color = line_type), 
               position = position_dodge(0.75), 
               size = 2.5,       
               width = 0.6,      
               alpha = 0.8,
               outlier.size = 5) + 
  scale_fill_manual(values = c(type_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) +  
  labs(x = "Number of collected ticks", y = "Line type") +
  theme_classic(base_size = 20) + 
  theme(
    axis.title = element_text(size = 24, face = "bold"),  
    axis.text = element_text(size = 20),                  
    legend.title = element_text(size = 22),               
    legend.text = element_text(size = 20),                
    plot.title = element_text(size = 28, face = "bold")   
  ) +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 26), 
               aes(group = line_type),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 6) 
plot

# Save the plot object
ggsave(filename = here::here("figures","envticks_ewda.pdf"), plot = plot, 
       width = 10.5, height = 7.5, units = "in", device = "pdf")




plot <- combined_table %>%
  filter(code_mission %in% c("Juin 2023")) %>%
  ggplot(aes(x = effectif_tick, y = factor(line_type))) + 
  facet_grid(factor(source) ~ .) +
  geom_boxplot(aes(fill = line_type, color = line_type), 
               position = position_dodge(0.75), 
               size = 2.5,       
               width = 0.6,      
               alpha = 0.8,
               outlier.size = 5) + 
  scale_fill_manual(values = c(type_palette, darker_palette)) + 
  scale_color_manual(values = darker_palette) +  
  labs(x = "Number of collected ticks", y = "Line type") +
  theme_classic(base_size = 20) + 
  theme(
    axis.title = element_text(size = 24, face = "bold"),  
    axis.text = element_text(size = 20),                  
    legend.title = element_text(size = 22),               
    legend.text = element_text(size = 20),                
    plot.title = element_text(size = 28, face = "bold")   
  ) +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 26), 
               aes(group = line_type),
               geom = "text", 
               position = position_dodge(0.8),
               color = "black", 
               size = 6) 

plot

# Save the plot object
ggsave(filename = here::here("figures","combinedticks_ewda.pdf"), plot = plot, 
       width = 11.5, height = 7.5, units = "in", device = "pdf")



