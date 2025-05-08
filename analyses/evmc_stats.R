#dataset used 

# data_for_m_noforests
#beware we use only above 10% in only hedge + pine edge
# and some species were removed for richness calculation (not sure if pathogens)

#dataset with only hedgerows
data_for_m_hedges <- data_for_m_noforests %>%
  filter(category == "hedgerows")
  
  

#Model : Bartonella_taylorii ----

## category first ----

table(data_for_m_noforests$category, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella_taylorii >0), .groups = "drop")

m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~  category + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
  
car::vif(m_barto_t_r)
SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_barto_t_r)




## broadleaved class ----

table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(broadleaved_class, season) %>%
  summarise(n_infected = sum(Bartonella_taylorii >0), .groups = "drop")

m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~  broadleaved_class * season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_barto_t_r)
SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_barto_t_r)



## connectivity first ----
#attention a prendre bon dataset


#Model : Mycoplasma haemomuris ----

## category first ----

table(data_for_m_noforests$category, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Mycoplasma_haemomuris >0), .groups = "drop")

rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ category + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_mycoplasma_r)


## broadleaved class ----

table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(broadleaved_class, season) %>%
  summarise(n_infected = sum(Mycoplasma_haemomuris >0), .groups = "drop")

rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ broadleaved_class * season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_mycoplasma_r)




#Model : Bartonella_birtlesii ----

## category first ----

table(data_for_m_noforests$category, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella_birtlesii >0), .groups = "drop")

m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ category + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_b_r)


## broadleaved class ----

table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(broadleaved_class, season) %>%
  summarise(n_infected = sum(Bartonella_birtlesii >0), .groups = "drop")

m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ broadleaved_class + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_b_r)

#Model : Mycoplasma coccoides ----

## category first ----

table(data_for_m_noforests$category, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Mycoplasma_coccoides >0), .groups = "drop")


m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ category + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_mycoplasmacoco_r)


## broadleaved class ----

table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(broadleaved_class, season) %>%
  summarise(n_infected = sum(Mycoplasma_coccoides >0), .groups = "drop")

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ broadleaved_class + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_mycoplasmacoco_r)



# Model : Neoehrlichia_mikurensis  ----

## category first ----
table(data_for_m_noforests$category, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Neoehrlichia_mikurensis >0), .groups = "drop")

rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category + season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_neoeh_r)
SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_neoeh_r)

## broadleaved class ----

table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(broadleaved_class, season) %>%
  summarise(n_infected = sum(Neoehrlichia_mikurensis >0), .groups = "drop")

rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ broadleaved_class * season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_neoeh_r)
SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_neoeh_r)

avg_model <- MuMIn::model.avg(TopModels)
summary(avg_model)


#Model : Bartonella_grahamii ----

## category first ----
table(data_for_m_noforests$category, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Bartonella_grahamii >0), .groups = "drop")


m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category + season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_barto_g_r)

## broadleaved class ----

table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$season)

data_for_m_noforests %>%
  group_by(broadleaved_class, season) %>%
  summarise(n_infected = sum(Bartonella_grahamii >0), .groups = "drop")


m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ broadleaved_class * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels



m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ broadleaved_class * season  + scale(poids)  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


DHARMa::simulateResiduals(m_barto_g_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_barto_g_r,.~.,test="Chisq")
summary(m_barto_g_r)

ggstats::ggcoef_model(m_barto_g_r)

em_bgra_brd <- emmeans::emmeans(m_barto_g_r, ~ broadleaved_class | season, type = "response")
plot(em_bgra_brd, comparisons = TRUE)
emmeans::contrast(em_bgra_brd, "pairwise", adjust = "Tukey")





# Model : Pathogen Richness ----

## category first ----
table(data_for_m_noforests$category, data_for_m_noforests$number_pathos)

m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ category * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

avg_model <- MuMIn::model.avg(TopModels)
summary(avg_model)

m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ category * season  + scale(poids)  + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

m_pathnumber_r %>%
  RVAideMemoire::overdisp.glmer()
