
####### REFAIRE MODELES AVEC ACP AXES :
- enlever interations entre axes,
- ajouter broadleavec class dans les model potetre (car paysage, pas truc local) 
- ajout interaction axe * season en priorité
- et interactions axe * annee et annee *saison dans l'optimal' 




## GLMMs ----

### Model : Neoehrlichia_mikurensis  ----

# 1. Considering treatment (C, NC, CT) - no broadleaved forests
rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ treatment + broadleaved_class + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_neoeh_r)

SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ treatment + broadleaved_class  + as.factor(year) + season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

DHARMa::simulateResiduals(m_neoeh_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_neoeh_r,.~.,test="Chisq")
summary(m_neoeh_r)

em <- emmeans::emmeans(m_neoeh_r, specs = pairwise ~ treatment, adjust = "Tukey", type = "response" )
em$contrasts
plot(em, comparisons = TRUE)

ggstats::ggcoef_model(m_neoeh_r)


# 2. Considering category (pine, hedgerows) - no broadleaved forests
rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category + broadleaved_class + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_neoeh_r)

SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category + broadleaved_class + as.factor(year) + season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
ggstats::ggcoef_model(m_neoeh_r)
  


# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled

rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_neoeh_r)


# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_neoeh_r)

SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_neoeh_r)


# 5. Considering extended broadleaved_class (LB, HB, B)
rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ broadleaved_class + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

ggstats::ggcoef_model(m_neoeh_r)


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


# 1. Considering treatment (C, NC, CT) - no broadleaved forests
rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ treatment * broadleaved_class + scale(poids) + season * as.factor(year) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_mycoplasma_r_best <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ scale(poids) + treatment + as.factor(year) + (1|numero_ligne),
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

# 2. Considering category (pine, hedgerows) - no broadleaved forests
rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ category * broadleaved_class + scale(poids) + season * as.factor(year) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


m_mycoplasma_r_best <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ scale(poids) + category + broadleaved_class + as.factor(year) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


ggstats::ggcoef_model(m_mycoplasma_r_best)


# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_mycoplasma_r)


# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_mycoplasma_r)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_mycoplasma_r)







### Model : Mycoplasma coccoides ----


# 1. Considering treatment (C, NC, CT) - no broadleaved forests
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ treatment * broadleaved_class + scale(poids) + as.factor(year)*season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ treatment + broadleaved_class + scale(poids) + as.factor(year) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


DHARMa::simulateResiduals(m_mycoplasmacoco_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

drop1(m_mycoplasmacoco_r,.~.,test="Chisq")
summary(m_mycoplasmacoco_r)


ggstats::ggcoef_model(m_mycoplasmacoco_r)

gtsummary::tbl_regression(m_mycoplasmacoco_r)


# 2. Considering category (pine, hedgerows) - no broadleaved forests
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ category * broadleaved_class + scale(poids) + as.factor(year)*season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_mycoplasmacoco_r)


# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

ggstats::ggcoef_model(m_mycoplasmacoco_r)


# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_mycoplasmacoco_r)
m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_mycoplasmacoco_r)

SelectionModels<- MuMIn::dredge(m_mycoplasmacoco_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_mycoplasmacoco_r)

### Model : Bartonella ----

# 1. Considering treatment (C, NC, CT) - no broadleaved forests
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ treatment * broadleaved_class + scale(poids) + season* as.factor(year) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_r_best <- lme4::glmer(
  formula = Bartonella ~ scale(poids) + season* as.factor(year) + broadleaved_class + (1|numero_ligne),
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


# 2. Considering category (pine, hedgerows) - no broadleaved forests
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ category + scale(poids) + season* as.factor(year) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_r)


# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_r)


# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_barto_r)
m_barto_r <- lme4::glmer(
  formula = Bartonella ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_r)

SelectionModels<- MuMIn::dredge(m_barto_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_r)


### Model : Bartonella_taylorii ----

# 1. Considering treatment (C, NC, CT) - no broadleaved forests
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ treatment * broadleaved_class + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ treatment + broadleaved_class + scale(poids) + as.factor(year) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

ggstats::ggcoef_model(m_barto_t_r)

gtsummary::tbl_regression(m_barto_t_r)

# 2. Considering category (pine, hedgerows) - no broadleaved forests
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ category * broadleaved_class + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_t_r)


# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
rm(m_barto_t_r)
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_t_r)

# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_barto_t_r)
m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_t_r)

SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_t_r)

### Model : Bartonella_grahamii ----

# 1. Considering treatment (C, NC, CT) - no broadleaved forests
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ treatment * broadleaved_class + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~  treatment + broadleaved_class + scale(poids) + as.factor(year) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

ggstats::ggcoef_model(m_barto_g_r)

gtsummary::tbl_regression(m_barto_g_r)


# 2. Considering category (pine, hedgerows) - no broadleaved forests
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category * broadleaved_class + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category * broadleaved_class + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

ggstats::ggcoef_model(m_barto_g_r)


## 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ category + scale(poids) + season  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

ggstats::ggcoef_model(m_barto_g_r)


# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_barto_g_r)
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_g_r)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_g_r)


### Model : Bartonella_birtlesii ----
# 1. Considering treatment (C, NC, CT) - no broadleaved forests
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ treatment * broadleaved_class + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
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


# 2. Considering category (pine, hedgerows) - no broadleaved forests
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ category + scale(poids) + as.factor(year) * season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_b_r)


# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ category + scale(poids) + season + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_b_r)

# 4. Considering small scale landscape characteristics PCA axis (no broadleaved forests) 
rm(m_barto_b_r)
m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ PCA_axis1 * PCA_axis2 + as.factor(year) * season  + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests_pca,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)
car::vif(m_barto_b_r)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

ggstats::ggcoef_model(m_barto_b_r)



















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


# 1. Considering treatment (C, NC, CT) - no broadleaved forests
# 2. Considering category (pine, hedgerows) - no broadleaved forests
# 3. Considering category (pine, hedgerows, broadleaved forsts) only in years with forests sampled
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ treatment * broadleaved_class + scale(poids) + season*year + sexe + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~  season + broadleaved_class + scale(poids) + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

m_pathnumber_r %>%
  RVAideMemoire::overdisp.glmer()

DHARMa::simulateResiduals(m_pathnumber_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_pathnumber_r,.~.,test="Chisq")
summary(m_pathnumber_r)

ggstats::ggcoef_model(m_pathnumber_r)




# Alternative
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ category + scale(poids) + code_mission + sexe + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_forestsyear,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


# Test for negative binomial
m_pathnumber_r<- MASS::glm.nb(formula = number_pathos ~ treatment * broadleaved_class + scale(poids) + season*year + sexe + (1|numero_ligne),
             data = data_for_m_noforests,
             na.action = "na.fail")

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_pathnumber_r<- MASS::glm.nb(formula = number_pathos ~ treatment * broadleaved_class + scale(poids) + season*year + sexe + (1|numero_ligne),
                              data = data_for_m_noforests,
                              na.action = "na.fail")

DHARMa::simulateResiduals(m_pathnumber_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

summary(m_pathnumber_r)




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

