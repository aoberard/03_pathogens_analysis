#dataset used ----

# data_for_m_noforests
#beware we use only above 10% in only hedge + pine edge
# and some species were removed for richness calculation (not sure if pathogens)



#au moment du csi 2 : 
# soit faire que sur données 2023-2024
# soit faire toutes années mais que données 16S
#parce que lepto et rpob barto pas fait à ce moment là

#liste pathos associés

listwith2025 <- c(
    "Bartonella",
    "Neoehrlichia_mikurensis",
    "Streptobacillus",
    "Mycoplasma_haemomuris",
    "Mycoplasma_coccoides",
    "Mycoplasma_penetrans",
    "Ehrlichia",
    "Neisseria",
    "Spiroplasma",
    "Treponema",
    "Borrelia_miyamotoi",
    "Borreliella",
    "Mycoplasmopsis_columbina",
    "Metamycoplasma_sualvi",
    "Mycoplasma")

list2025_10percent <- intersect(listwith2025, patho10_apo)

list20232024 <-c(
  "Bartonella_taylorii",
  "Bartonella_grahamii",
  "Bartonella_doshiae",
  "Bartonella_birtlesii",
  "Bartonella_elizabethae",
  "Leptospira",
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
  "Borrelia_miyamotoi",
  "Borreliella",
  "Mycoplasmopsis_columbina",
  "Metamycoplasma_sualvi",
  "Mycoplasma")

list20232024_10percent <- intersect(list20232024, patho10_apo)

data_for_m_2023224 <- data_for_m %>%
  filter(year != 2025)

#Associations explorations ----


## ACM ----
process_data_for_m_2025 <- data_for_m %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  mutate(across(all_of(listwith2025), ~ as.factor(.)))


process_data_for_m_2023224 <- data_for_m_2023224 %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  mutate(across(all_of(list20232024), ~ as.factor(.)))

### 2025 ----

# df est ton tableau de présence/absence
acm2025<- FactoMineR::MCA(process_data_for_m_2025[, list2025_10percent], graph = FALSE)

# visualiser les agents pathogènes
factoextra::fviz_mca_var(acm2025,                     
                         axes = c(1, 2),
                         repel = F,
                         col.var = "contrib",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         pointsize = 3,
                         labelsize = 3,
                         ggtheme = theme_minimal(base_size = 14) 
)+
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "grey30")
  )

# visualiser les individus
factoextra::fviz_mca_ind(acm2025, repel = F)




### 20232024 ----
# df est ton tableau de présence/absence
acm20232024<- FactoMineR::MCA(process_data_for_m_2023224[, list20232024_10percent], graph = FALSE)

# visualiser les agents pathogènes
factoextra::fviz_mca_var(acm20232024,   
                         axes = c(1, 2),
                         repel = F,
                         col.var = "contrib",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         pointsize = 3,
                         labelsize = 3,
                         ggtheme = theme_minimal(base_size = 14) 
)+
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "grey30")
  )


# visualiser les individus
factoextra::fviz_mca_ind(acm20232024, repel = F)





# GLM especes majoritaires ----

##Model : Bartonella spp. ----

data_for_m %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_barto_r")) rm(m_barto_r)
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_r)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_barto_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification
if (exists("m_barto_r")) rm(m_barto_r)
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ season + as.factor(year) + scale(poids) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_r)

#Best model validation
DHARMa::simulateResiduals(m_barto_r) %>%
  DHARMa::testResiduals()

#Best model collinearity
performance::check_collinearity(m_barto_r)



##Model : Mycoplasma coccoides ----

data_for_m %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Mycoplasma_coccoides >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_mycoplasmacoco_r")) rm(m_mycoplasmacoco_r)
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_mycoplasmacoco_r)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification
if (exists("m_mycoplasmacoco_r")) rm(m_mycoplasmacoco_r)
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ season + scale(poids) + scale(`Apodemus sylvaticus`) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_mycoplasmacoco_r) %>%
  DHARMa::testResiduals()

#Best model collinearity
performance::check_collinearity(m_mycoplasmacoco_r)




##Model : Mycoplasma haemomuris ----
#data 2025

data_for_m %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Mycoplasma_haemomuris >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_mycoplasma_r")) rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_mycoplasma_r)


#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification
if (exists("m_mycoplasma_r")) rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ scale(poids)  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_mycoplasma_r) %>%
  DHARMa::testResiduals()

#Best model collinearity
performance::check_collinearity(m_mycoplasma_r)



## Model : Neoehrlichia mikurensis  ----

data_for_m %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Neoehrlichia_mikurensis >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_neoeh_r")) rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_neoeh_r)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_neoeh_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

# #Best model specification
# if (exists("m_neoeh_r")) rm(m_neoeh_r)
# m_neoeh_r <- lme4::glmer(
#   formula = Neoehrlichia_mikurensis ~ category * season + as.factor(year) + scale(poids) + sexe + (1|numero_ligne),
#   family = binomial(link = "logit"),
#   data = data_for_m,
#   na.action = "na.fail",                                  
#   control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
# )

#Best model validation
DHARMa::simulateResiduals(m_neoeh_r) %>%
  DHARMa::testResiduals()

#Best model collinearity
performance::check_collinearity(m_neoeh_r)



##Model : Mycoplasma sp. ----

data_for_m %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Mycoplasma >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_mycoplasmasp")) rm(m_mycoplasmasp)
m_mycoplasmasp <- lme4::glmer(
  formula = Mycoplasma ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_mycoplasmasp)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_mycoplasmasp, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification
if (exists("m_mycoplasmasp")) rm(m_mycoplasmasp)
m_mycoplasmasp <- lme4::glmer(
  formula = Mycoplasma ~ season + scale(poids) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_mycoplasmasp) %>%
  DHARMa::testResiduals()

#Best model collinearity
performance::check_collinearity(m_mycoplasmasp)



# ici données 2023 2024 ----

##Model : Bartonella taylorii  ----

data_for_m_2023224 %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella_taylorii >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_barto_t_r")) rm(m_barto_t_r)
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_2023224,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_t_r)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_barto_t_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification





##Model : Bartonella birtlesii  ----

data_for_m_2023224 %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella_birtlesii >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_barto_b_r")) rm(m_barto_b_r)
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_2023224,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_b_r)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_barto_b_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification
if (exists("m_barto_b_r")) rm(m_barto_b_r)
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ scale(poids) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_2023224,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_barto_b_r) %>%
  DHARMa::testResiduals()




##Model : Bartonella grahamii  ----

data_for_m_2023224 %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella_grahamii >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_barto_g_r")) rm(m_barto_g_r)
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category * season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_2023224,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_g_r)

#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_barto_g_r, rank = "AICc")
model_selection %>% filter(delta <2)

#Take a look at best models averaging
avg_mod <- MuMIn::model.avg(model_selection, subset = delta < 10)
summary(avg_mod)
MuMIn::sw(avg_mod)

#Best model specification
if (exists("m_barto_g_r")) rm(m_barto_g_r)
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ season + as.factor(year) + scale(poids) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_2023224,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_barto_g_r) %>%
  DHARMa::testResiduals()

