#dataset used ----

# data_for_m_noforests
#beware we use only above 10% in only hedge + pine edge
# and some species were removed for richness calculation (not sure if pathogens)

#dataset with only hedgerows
data_for_m_hedges <- data_for_m_noforests %>%
  filter(category == "hedgerows")
  
#new palette 
brd_palette <- c("LB" = "#D2B48C",  # light brown (tan)
                 "HB" = "#8B4513")

conect_palette <- c(
  "C" = "#B07C9E",       # Muted mauve
  "NC" = "#7A91A8"    # Muted blue-gray
)

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
table(data_for_m_hedges$treatment, data_for_m_hedges$season)

data_for_m_hedges %>%
  group_by(treatment, season) %>%
  summarise(n_infected = sum(Bartonella_taylorii >0), .groups = "drop")

m_barto_t_r <- lme4::glmer(
  formula = Bartonella_taylorii ~  treatment * season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_barto_t_r)
SelectionModels<- MuMIn::dredge(m_barto_t_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_barto_t_r)




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



## connectivity first ----
table(data_for_m_hedges$treatment, data_for_m_hedges$season)

data_for_m_hedges %>%
  group_by(treatment, season) %>%
  summarise(n_infected = sum(Mycoplasma_haemomuris >0), .groups = "drop")

rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ treatment * season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_mycoplasma_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels
ggstats::ggcoef_model(m_mycoplasma_r)

avg_model <- MuMIn::model.avg(TopModels)
summary(avg_model)

rm(m_mycoplasma_r)
m_mycoplasma_r <- lme4::glmer(
  formula = Mycoplasma_haemomuris ~ treatment + scale(poids) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


DHARMa::simulateResiduals(m_mycoplasma_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_mycoplasma_r,.~.,test="Chisq")
summary(m_mycoplasma_r)





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



## connectivity first ----
table(data_for_m_hedges$treatment, data_for_m_hedges$season)

data_for_m_hedges %>%
  group_by(treatment, season) %>%
  summarise(n_infected = sum(Bartonella_birtlesii >0), .groups = "drop")

m_barto_b_r <- lme4::glmer(
  formula = Bartonella_birtlesii ~ treatment + season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_b_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels




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


## connectivity first ----
table(data_for_m_hedges$treatment, data_for_m_hedges$season)

data_for_m_hedges %>%
  group_by(treatment, season) %>%
  summarise(n_infected = sum(Mycoplasma_coccoides >0), .groups = "drop")

m_mycoplasmacoco_r <- lme4::glmer(
  formula = Mycoplasma_coccoides ~ broadleaved_class * season  + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
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
  formula = Neoehrlichia_mikurensis ~ category + season + as.factor(year) + scale(poids) + sexe  + (1|numero_ligne),
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
  formula = Neoehrlichia_mikurensis ~ broadleaved_class * season + as.factor(year) + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_neoeh_r)
SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels


avg_model <- MuMIn::model.avg(TopModels)
summary(avg_model)



## connectivity first ----
table(data_for_m_hedges$treatment, data_for_m_hedges$season)

data_for_m_hedges %>%
  group_by(treatment, season) %>%
  summarise(n_infected = sum(Neoehrlichia_mikurensis >0), .groups = "drop")


rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ treatment * season + as.factor(year) + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

car::vif(m_neoeh_r)
SelectionModels<- MuMIn::dredge(m_neoeh_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ treatment  + as.factor(year) + scale(poids) + sexe + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)


DHARMa::simulateResiduals(m_neoeh_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_neoeh_r,.~.,test="Chisq")
summary(m_neoeh_r)


em_neo_conect <- emmeans::emmeans(m_neoeh_r, ~ treatment, type = "response")
plot(em_neo_conect, comparisons = TRUE)
emmeans::contrast(em_neo_conect, "pairwise", adjust = "Tukey")

em_summary_neo_conect <- as.data.frame(em_neo_conect)

em_neo_conect_plot <- ggplot(em_summary_neo_conect, aes(y = prob, x = reorder(treatment, prob), color = treatment)) +
  geom_point(size = 4.5) +
  geom_segment(aes(y = asymp.LCL, yend = asymp.UCL, x = treatment, xend = treatment), linewidth = 1.5) +
  scale_color_manual(values = conect_palette) + 
  theme_minimal(base_size = 16) +
  labs(
    y = "Estimated something probability",
    x = "BRD",
    color = "BRD"
  ) +
  theme(
    axis.text = element_text(color = "gray20", size = 11),
    axis.title = element_text(color = "gray20", size = 12),
    strip.text = element_text(color = "gray20", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)
  ) +
  guides(
    color = "none" 
  )
em_neo_conect_plot

ggsave(filename = here::here("figures","neo_conect_evmc.pdf"),
       plot = em_neo_conect_plot, 
       width = 1100, height = 1100, units = "px", device = "pdf", dpi = 300, bg = NULL)



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

rm(m_barto_g_r)
m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ broadleaved_class * season + as.factor(year) + scale(poids) + sexe  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_noforests,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

avg_model <- MuMIn::model.avg(TopModels)
summary(avg_model)

m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ broadleaved_class * season + as.factor(year) + scale(poids)  + (1|numero_ligne),
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

em_summary_bgra_brd <- as.data.frame(em_bgra_brd)

em_bgra_brd_plot <- ggplot(em_summary_bgra_brd, aes(y = prob, x = reorder(broadleaved_class, prob), color = broadleaved_class)) +
  geom_point(size = 4.5) +
  geom_segment(aes(y = asymp.LCL, yend = asymp.UCL, x = broadleaved_class, xend = broadleaved_class), linewidth = 1.5) +
  scale_color_manual(values = brd_palette) + 
  facet_wrap(season ~ ., ncol = 2) + 
  theme_minimal(base_size = 16) +
  labs(
    y = "Estimated something probability",
    x = "BRD",
    color = "BRD"
  ) +
  theme(
    axis.text = element_text(color = "gray20", size = 11),
    axis.title = element_text(color = "gray20", size = 12),
    strip.text = element_text(color = "gray20", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)
  ) +
  guides(
    color = "none" 
  )
em_bgra_brd_plot

ggsave(filename = here::here("figures","bgrahamii_brd_evmc.pdf"),
       plot = em_bgra_brd_plot, 
       width = 1100, height = 1100, units = "px", device = "pdf", dpi = 300, bg = NULL)






## connectivity first ----
table(data_for_m_hedges$treatment, data_for_m_hedges$season)

data_for_m_hedges %>%
  group_by(treatment, season) %>%
  summarise(n_infected = sum(Bartonella_grahamii >0), .groups = "drop")


m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ treatment * season + as.factor(year) + scale(poids)  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_barto_g_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

m_barto_g_r <- lme4::glmer(
  formula = Bartonella_grahamii ~ treatment + season + as.factor(year) + scale(poids)  + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

DHARMa::simulateResiduals(m_barto_g_r, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_barto_g_r,.~.,test="Chisq")
summary(m_barto_g_r)

ggstats::ggcoef_model(m_barto_g_r)






# Model : Pathogen Richness ----

## category first ----
table(data_for_m_noforests$category, data_for_m_noforests$number_pathos)

rm(m_pathnumber_r)
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ category * season + as.factor(year) + scale(poids) + sexe  + (1|numero_ligne),
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

m_pathnumber_r %>%
  RVAideMemoire::overdisp.glmer()

ggstats::ggcoef_model(m_pathnumber_r)

#Post-hoc tests
em_pathrich <- emmeans::emmeans(m_pathnumber_r, ~ category | season, type = "response")
plot(em_pathrich, comparisons = TRUE)
emmeans::contrast(em_pathrich, "pairwise", adjust = "Tukey")





## broadleaved class ----
table(data_for_m_noforests$broadleaved_class, data_for_m_noforests$number_pathos)


rm(m_pathnumber_r)
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ broadleaved_class * season + scale(poids) + sexe  + (1|numero_ligne),
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

m_pathnumber_r %>%
  RVAideMemoire::overdisp.glmer()

ggstats::ggcoef_model(m_pathnumber_r)

#Post-hoc tests
em_pathrich <- emmeans::emmeans(m_pathnumber_r, ~ broadleaved_class | season, type = "response")
plot(em_pathrich, comparisons = TRUE)
emmeans::contrast(em_pathrich, "pairwise", adjust = "Tukey")



## connectivity first ----
table(data_for_m_hedges$treatment, data_for_m_hedges$number_pathos)

rm(m_pathnumber_r)
m_pathnumber_r <- lme4::glmer(
  formula = number_pathos ~ treatment * season + scale(poids) + sexe  + (1|numero_ligne),
  family = poisson(link = "log"),
  data = data_for_m_hedges,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

SelectionModels<- MuMIn::dredge(m_pathnumber_r, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels





#Entick test ----

# GLM making (poisson) - broadleaved_class
rm(glm_environment_tick)
glm_environment_tick <- lme4::glmer(
  formula = effectif_tick ~ broadleaved_class + season + (1|numero_ligne),
  data = dglm_environment_tick,
  family = poisson(link = "log"),  
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)


SelectionModels<- MuMIn::dredge(glm_environment_tick, rank = "AICc")              
TopModels<-subset(SelectionModels, delta < 2)
TopModels

glm_environment_tick |>
  RVAideMemoire::overdisp.glmer()


drop1(glm_environment_tick,.~.,test="Chisq")
summary(glm_environment_tick)


# GLM making (poisson) - connectivity

dglm_environment_tick_hedge <- dglm_environment_tick %>%
  filter(category == "hedgerows") %>%
  filter(season == "Spring")

rm(glm_environment_tick)
glm_environment_tick <- lme4::glmer(
  formula = effectif_tick ~ treatment + (1|numero_ligne),
  data = dglm_environment_tick_hedge,
  family = poisson(link = "log"),  
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)


SelectionModels<- MuMIn::dredge(glm_environment_tick, rank = "AICc")              
TopModels<-subset(SelectionModels, delta < 2)
TopModels

glm_environment_tick |>
  RVAideMemoire::overdisp.glmer()


drop1(glm_environment_tick,.~.,test="Chisq")
summary(glm_environment_tick)


glm_environment_tick_nb <- glmmTMB::glmmTMB(effectif_tick ~ treatment  + (1 | numero_ligne),
                   data = dglm_environment_tick_hedge,
                   family = glmmTMB::nbinom2,
                   na.action = "na.fail")

SelectionModels<- MuMIn::dredge(glm_environment_tick_nb, rank = "AICc")
TopModels<-subset(SelectionModels, delta < 2)
TopModels

summary(glm_environment_tick_nb)

