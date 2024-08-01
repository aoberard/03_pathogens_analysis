
# /!\ WARNING ----
#THIS SCRIPT IS TEMPORARY AND SHOULD NOT BE CONSERVED IN THIS R PROJECT LATER ----

#the script extracting the host data from bpm has to be run before this script
#this script should only be used to analyse 16S or rpob/glta data from spleen, in order to identify pathogenic bacteria
#this script is dependent on the generation of bpm_beprep datafram from temp_import_survey and bpm_host_macroparasite from temp_tique to join dataframes

# /!\ WARNING ----



library(ggplot2)
library(dplyr)
library(data.table)
library(here)

# GENERATE DATA ----

## Import data ----
# Import hosts and line modalities file
d_host <- readr::read_csv( here::here("data/", "raw-data/", "host_data", "20240731_bpm_modalities.csv") )

# Import 16S filtered file
file16s_run00_01 <- data.table::fread(file = here::here( "data","raw-data/","16s_run00-01","Run00-01_16S_filtered_postfrogs.txt"))

# Import rodent macroparasite
d_macroparasite <- fread(file = here::here("data/", "derived-data/", "raw-ticks", "rodents_tick", "20240731_macroparasite.csv") )

## Generate a file containing samples of interest ----

# Extract id of small mammals caught in Beprep
sm_id <- unique(d_host %>% 
                  filter(!is.na(numero_centre))%>%
                  pull(numero_centre)
)


# Identifiy taxonomy columns
taxo_name <- colnames( file16s_run00_01 [, 1:which(colnames(file16s_run00_01) == "seed_sequence") ])

# Calculate total run read number per cluster
file16s_run00_01 <- file16s_run00_01 %>%
  mutate(clean_totalrunreads = rowSums(select(., -taxo_name)))%>%
  relocate(clean_totalrunreads, .after = names(file16s_run00_01)[which(colnames(file16s_run00_01) == "seed_sequence")])

# Reattribute taxonomy columns
taxo_name <- colnames( file16s_run00_01 [, 1:which(colnames(file16s_run00_01) == "clean_totalrunreads") ])

# Take away problematic samples - if needed
#file16s_run00_01 <- file16s_run00_01 |>


# Generate new file containing only Spleen BePrep's samples
#by selecting those samples
beprep_sample_16s_sp <- file16s_run00_01 |>
  select( (contains("JPRA") | contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 
#and bind them to the taxonomy columns
beprep16s_sp <- file16s_run00_01 |>
  select( taxo_name ) |>
  cbind(beprep_sample_16s_sp)

# Delete the PCZymo positive after visual control
beprep16s_sp <- beprep16s_sp |>
  select(!contains("PCzymo"))

# Calculate total number of reads for each cluster considering chosen samples
beprep16s_sp <- beprep16s_sp %>%
  mutate(totalreads = rowSums(select(., -taxo_name)))%>%
  relocate(totalreads, .after = names(beprep16s_sp)[which(colnames(file16s_run00_01) == "clean_totalrunreads")])

# Order rows by new total read number
beprep16s_sp <- data.table::setorder(beprep16s_sp, -totalreads)

# Delete cluster with no reads for the selected individuals
beprep16s_sp <- beprep16s_sp |> 
  filter(totalreads >0)

# Write a file containing 16S Spleen OTUs containing our samples only :
data.table::fwrite(beprep16s_sp, here::here("data/", "derived-data/","16s_run00-01", "20240326_beprep16s.txt") )





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
  filter(genus %in% putative_pathos_genra | family %in% putative_pathos_family)

# Write a file containing 16S Spleen OTUs containing our samples only AND only putative pathogens :
data.table::fwrite(pathos16s, here::here("data/", "derived-data/","16s_run00-01", "20240326_beprep16s_pathos.txt") )


## Seek for species level identification for specific taxa (Mycoplasma) throught external manipulation (phylogenic tree) ----
# A new datafile is created by hand to add the checked species information and then imported
# If probable pseudo-gene affiliation is observed in data file they may also get deleted at this step
pathos16s_checked <- data.table::fread(file = here::here( "data","derived-data/","16s_run00-01","20240326_beprep16s_pathos_species-adjusted.txt"))
pathos16s_checked$species_adjusted[pathos16s_checked$species_adjusted == ""] <- NA

# Change species value when there is a species adjusted value (generated by operator manually ) and drop the species adjusted column
pathos16s_checked <- pathos16s_checked %>%
  mutate(species = if_else(condition =  !is.na(species_adjusted), species_adjusted, species )) %>%
  select(!species_adjusted)


## Create new column for best identified taxa ----

# If there are still unknown in species column, change it to NA
pathos16s_checked <- pathos16s_checked %>%
  mutate(species = if_else(species %in% c("unknown species", "Multi-affiliation"), NA , species )) 

# Creation of a new column containing best available identity (between species or genus here)
pathos16s_checked <- pathos16s_checked %>%
  mutate(identification16s = if_else(condition =  !is.na(species), species, genus) ) %>%
  relocate(identification16s, .after = "species")

# Replace space in taxa name 
pathos16s_checked <- pathos16s_checked %>%
  mutate(identification16s = gsub(" ", "_", identification16s ))


## Simplification, compress 16s information by best identified taxa ----

# Keep only best taxa and samples
taxa_pathos16s <- pathos16s_checked %>%
  select(c("identification16s", names(pathos16s_checked)[(which(names(pathos16s_checked) == "totalreads") + 1):ncol(pathos16s_checked)]))

# Sum rows by best taxa
taxa_pathos16s <- taxa_pathos16s |>
  group_by(identification16s) |>
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# Transpose the table
taxa_transposed <- taxa_pathos16s |>
  data.table::transpose(
    keep.names = "numero_centre_16s",
    make.names = 1
  )

# Finish sample name transformation to join to host data
taxa_transposed <- taxa_transposed |>
  mutate(numero_centre_16s = gsub(".SP", "", numero_centre_16s ))


## Join to macroparasite for complete pathogen data frame ----

# Select columns for macroparasite of interest (data from temp_tique script) and transform NA to 0 i effectif_tick
d_macroparasite <- d_macroparasite %>%
  select( c("effectif_tick", "numero_centre", "Ixodida", "Siphonaptera"))

# Join macroparasite on host trapping data, keep every lines from trapping data
d_macroparasite <- left_join(d_host, d_macroparasite, by = "numero_centre" )

# Transform NA to 0
d_macroparasite <- d_macroparasite %>%
  mutate( across(c("Ixodida", "Siphonaptera"), ~if_else(. == 0, NA, .) ) ) %>%
  mutate( across(c("Ixodida", "Siphonaptera"), ~if_else( is.na(.), 0, 1) ) ) %>%
  mutate( effectif_tick = if_else( is.na(effectif_tick), 0, effectif_tick) )

# Join 16s pathogens on host trapping/macroparasite data - complete join, keep only 16s screened individuals /!\
rodent_pathos <- merge(x = d_macroparasite, y = taxa_transposed, by.x = "numero_centre",  by.y = "numero_centre_16s")

# Extract Pathogens name
pathos_name <- rodent_pathos |>
  select(names(rodent_pathos)[(which(names(rodent_pathos) == "effectif_tick") + 1):ncol(rodent_pathos)]) %>%
  colnames()
pathos_name





# ANALYSIS (ONLY FOR APODEMUS FOR NOW) ----

## GLM preliminary ----

### Generate data ----
# Keep only Apodemus for analysis
d_apo_pathos_glm <- rodent_pathos %>%
  filter(taxon_mamm == "Apodemus sylvaticus")


# So filter empty pathogen for apodemus and change pathos name accordingly
d_apo_pathos_glm <- d_apo_pathos_glm %>%
  mutate(across(all_of(pathos_name), ~if(sum(.) != 0) . else NULL) )
pathos_name_apo <- d_apo_pathos_glm |>
  select(names(d_apo_pathos_glm)[(which(names(d_apo_pathos_glm) == "effectif_tick") + 1):ncol(d_apo_pathos_glm)]) %>%
  colnames()
setdiff(pathos_name, pathos_name_apo)

# Generate the dataframe for logistic analysis (0/1)
data_for_m <- d_apo_pathos_glm %>%
  mutate(across(all_of(pathos_name_apo), ~ replace(., . > 0, 1)))

# Add new variable : pathogen richness
data_for_m <- data_for_m %>%
  mutate(number_pathos = rowSums(across(all_of(pathos_name_apo))) ) 

# Identify pathos with prevalence >=10 for whole data
patho10_apo<- d_apo_pathos_glm %>%
  summarise( across(all_of(pathos_name_apo), ~ sum(. > 0) / n() ))  %>%
  unlist()
names(patho10_apo[patho10_apo >= 0.10])

# Take away rows with NA
data_for_m <- data_for_m %>%
  filter(!is.na(sexe))


### Variable exploration ----
data_for_m %>% count(Neoehrlichia_mikurensis, line_treatment)

ggplot(data_for_m, aes(x = as.factor(Haemobartonella_muris), y = poids)) +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = sexe), position=position_dodge(0.4)) + 
  coord_flip() +
  stat_summary(fun = mean, geom="point", shape=18,
               size=5, color="purple") +
  theme_minimal()

hist(data_for_m$number_pathos)


### Model : Mycoplasma haemomuris ----
rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Haemobartonella_muris ~ broadleaved_status * connectivity + code_mission + poids + sexe +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)

DHARMa::simulateResiduals(m_mycoplasma_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_mycoplasma_r_best <- lme4::glmer(
  formula = Haemobartonella_muris ~ connectivity + poids +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
DHARMa::simulateResiduals(m_mycoplasma_r_best, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_mycoplasma_r_best,.~.,test="Chisq")
summary(m_mycoplasma_r_best)

em_myco <- emmeans::emmeans(m_mycoplasma_r_best, specs = pairwise ~ connectivity, adjust = "Tukey", type = "response" )
em_myco$contrasts

plot(em_myco, comparisons = TRUE)


### Model : Mycoplasma coccoides ----

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ broadleaved_status * connectivity + code_mission + poids + sexe +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
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



### Model : Bartonella ----
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ broadleaved_status * connectivity + code_mission + poids + sexe +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_r <- lme4::glmer(
  formula = Bartonella ~ code_mission + poids  +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(m_barto_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_barto_r,.~.,test="Chisq")
summary(m_barto_r)



### Model : Neoehrlichia_mikurensis  ----
# Pas d'interaction dans le modele pour 2023 car segregation trop grande /!\
rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ broadleaved_status * connectivity + code_mission + poids + sexe +(1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ broadleaved_status  + code_mission + poids  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(m_neoeh_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_neoeh_r,.~.,test="Chisq")
summary(m_neoeh_r)

em <- emmeans::emmeans(m_neoeh_r, specs = pairwise ~ broadleaved_status, adjust = "Tukey", type = "response" )
em$contrasts

plot(em, comparisons = TRUE)


#### test graph ----

# Extract fixed effects coefficients
fixed_effects <- coef(summary(m_neoeh_r))[, "Estimate"] # ou : fixed_effects <- lme4::fixef(m_neoeh_r)

# Create data for plotting
plot_data <- expand.grid(
  broadleaved_status = unique(data_for_m$broadleaved_status),
  code_mission = unique(data_for_m$code_mission),
  poids = seq(min(data_for_m$poids), max(data_for_m$poids), length.out = 100)  
)

# Predict response variable
plot_data$predicted_prob <- predict(m_neoeh_r, newdata = plot_data, type = "response", re.form = NA)

# Plot
library(ggplot2)
ggplot(plot_data, aes(x = poids, y = predicted_prob, color = broadleaved_status)) +
  geom_line() +
  facet_wrap(~ code_mission) +
  labs(x = "Poids", y = "Predicted Probability", color = "Broadleaved Status")





### Model : Ixodida  ----
m_ixod_r <- lme4::glmer(
  formula = Ixodida ~ broadleaved_status * connectivity + code_mission + poids + sexe +(1|numero_ligne),
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
  formula = Siphonaptera ~ broadleaved_status * connectivity + code_mission + poids + sexe +(1|numero_ligne),
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


### arefaire car pas poisson + dataframe a changer car mulot ciblé |Model : Pathogen number  ----
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ broadleaved_status * connectivity + code_mission + poids  + sexe +(1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
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




### Model : Parasite Diversity for all rodents  ----

#### Generate alpha diversity index ----

alpha_div_pathos <- rodent_pathos %>%
  mutate(across(all_of(pathos_name), ~ replace(., . > 0, 1))) %>%
  group_by(numero_ligne, code_mission) %>%
  summarise(connectivity = unique(connectivity),
            broadleaved_status = unique(broadleaved_status),
            across(all_of(pathos_name), ~ sum(. > 0))
  ) %>%
  mutate( richness = rowSums(across(all_of(pathos_name)) > 0)) %>%
  mutate (shannon = vegan::diversity(across(all_of(pathos_name)) ))

hist(alpha_div_pathos$richness)
hist(alpha_div_pathos$shannon)



#### Poisson regression ----

m_rich_pathos <- stats::glm(
  formula = richness  ~ broadleaved_status * connectivity ,
  family = poisson(link="log"),
  data = alpha_div_pathos,
  na.action="na.fail"
)

DHARMa::simulateResiduals(m_rich_pathos, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

summary(m_rich_pathos)

performance::check_zeroinflation(m_rich_pathos) 


#### zero_inflated poisson regression ----
m_richness_zip <-pscl::zeroinfl(richness ~ connectivity + broadleaved_status + code_mission, 
                                data = alpha_div_pathos, dist = "poisson", na.action = "na.fail")

summary(m_richness_zip)
drop1(m_richness_zip,.~.,test="Chisq")


SelectionModels<- MuMIn::dredge(m_richness_zip, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_richness_zip <-pscl::zeroinfl(richness ~ connectivity, 
                                data = alpha_div_pathos, dist = "poisson", na.action = "na.fail")

summary(m_richness_zip)
drop1(m_richness_zip,.~.,test="Chisq")

## Graph making heatmap matrix ----

# Prevalence matrix calculation (all rodents)
#for line_treatment
graph_sequence <- c("CT_LB", "CT_HB", "NC_LB", "NC_HB", "C_LB", "C_HB" )

matrix_pathoss <- d_apo_pathos_glm %>%
  group_by(connectivity, broadleaved_status) %>%
  summarise(
    across(all_of(pathos_name_apo), 
           ~ round(sum(. > 0) / n(), digits = 2)) ) |>
  mutate(combined_name = paste( connectivity, broadleaved_status, sep = "_")) |>
  arrange(factor(combined_name, levels = graph_sequence)) |>
  tibble::column_to_rownames("combined_name") |>
  select(all_of(pathos_name_apo)) |>
  as.matrix() 

effectif_df <- d_apo_pathos_glm %>%
  mutate(combined_name = paste(connectivity, broadleaved_status, sep = "_")) %>%
  group_by(connectivity, broadleaved_status) %>%
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



# TICK EFFECTIF ANALYSIS ----

hist(rodent_pathos$effectif_tick)

m_tiqueeff <- lme4::glmer.nb(
  formula = effectif_tick ~ broadleaved_status * connectivity + (1|numero_ligne),
  data = rodent_pathos,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

DHARMa::simulateResiduals(m_tiqueeff, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

m_tiqueeff |>
  RVAideMemoire::overdisp.glmer()

SelectionModels<- MuMIn::dredge(m_tiqueeff, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels




# BARTONELLA RPOB ----
filerpobtempo <- data.table::fread(file = here::here( "data","raw-data/","rpob_run01","20240325.nobim_corrected.txt"))

#by selecting those samples
ech <- filerpobtempo |>
  select( (contains("JPRA") | contains("NCHA") & contains(".SP")) | contains("PCzymo") ) |>
  select( matches(paste(sm_id, collapse = "|")) | contains("PCzymo") ) 

# Identifiy taxonomy columns
taxo_ech <- colnames( filerpobtempo [, 1:10 ])

#and bind them to the taxonomy columns
beprep16s_sp <- file16s_run00_01 |>
  select( taxo_name ) |>
  cbind(beprep_sample_16s_sp)



