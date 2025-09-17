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


#Arrange dataframe per numero_ligne, otherwise there are problem for autocorrelation testing
data_for_m <- data_for_m %>%
  arrange(numero_ligne)

data_for_m_2023224 <- data_for_m_2023224 %>%
  arrange(numero_ligne)




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


## Screening association ----

source(here::here("R", "FctTestScreenENV.txt"))
source(here::here("R", "SCN.R"))


#Matrices des parasites

m_screen_20232024 <- data_for_m_2023224  %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(all_of(list20232024_10percent)) %>%              
  as.matrix()

m_screen_2025 <- data_for_m %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(all_of(list2025_10percent)) %>%              
  as.matrix()

#indicateur que le nombre d'hôtes est suffisant pour le nombre de combinaisons de parasites considérés
ifelse(nrow(m_screen_20232024) < (2^ncol(m_screen_20232024)), "Peu d'hôtes pour le nombre de combinaisons", "BINGO, nombre suffisants d'hôtes")
ifelse(nrow(m_screen_2025) < (2^ncol(m_screen_2025)), "Peu d'hôtes pour le nombre de combinaisons", "BINGO, nombre suffisants d'hôtes")

# Application de la fonction et struture de sauvegarde des résultats :
ResSCRENV20232024 <- FctTestScreenENV(m_screen_20232024)
ResSCRENV2025 <- FctTestScreenENV(m_screen_2025)










#Fait pour données 2023 2024

### Attribution des variables du graphique:
xNC <- length(ResSCRENV20232024[[3]]) # Nombre de combinaisons
x <- (1:xNC)
yinf <- ResSCRENV20232024[[8]] # les différentes bornes des intervalles de confiance des différentes combinaisons
ysup <- ResSCRENV20232024[[7]]
yobs <- ResSCRENV20232024[[6]]


### Graphique : 
plot(x,yobs,type="n",ylim=range(c(0.0,yinf,ysup,yobs)),
     xlab="indice de combinaison",ylab="effectif",
     main="Enveloppes de confiance 95%")                   # fond de graph
points(x,yinf,lty=1,col=4,lwd=2,pch=19)                           # borne inférieure de chaque combinaison
points(x,ysup,lty=1,col=3,lwd=2,pch=19)                           # borne supérieure de chaque combinaison
points(x,yobs,col=2,pch=19)                                       # Observation de chaque combinaison

valx <- x[(ResSCRENV20232024[[4]]!=0)]                    # combinaisons hors de l'enveloppes 

for (i in 1:length(valx) ) {                      # ligne verticale pour toutes les combinaisons hors enveloppes
  abline(v=valx[i],lty=3,col=1,lwd=1)
}



# --- Données parasites 2023-2024 ---
t_screen_20232024 <- data_for_m_2023224 %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(all_of(listwith2025))

# --- Paramètres de base ---
xNH <- nrow(t_screen_20232024)  # nombre d’hôtes
xNP <- ncol(t_screen_20232024)  # nombre de parasites
xpv <- apply(t_screen_20232024, 2, sum) / xNH # prévalence
xNC <- 2^xNP                     # nb total de combinaisons possibles

# --- Construction de la matrice des combinaisons parasites (0/1) ---
DatXMul <- array(0, c(xNC, xNP))
for (k in 1:xNP) {
  k1 <- 2^(xNP - k)
  k2 <- 2^(k - 1)
  DatXMul[, k] <- rep(c(rep(1, k1), rep(0, k1)), k2)
}

# Identifiants uniques pour chaque combinaison
DatFMul <- rep(0, xNC)
for (k in 1:xNC) {
  x <- 0
  for (j in 1:xNP) {
    if (DatXMul[k, j] == 1) {
      x <- x + (10^(xNP - j))
    }
  }
  DatFMul[k] <- x
}
rownames(DatXMul) <- DatFMul

# Nettoyage et mise en forme
DatXMul <- t(DatXMul)
rownames(DatXMul) <- names(t_screen_20232024)
DatXMul <- t(DatXMul)
DatXMul <- DatXMul[rev(1:nrow(DatXMul)), ]
DatXMul <- t(DatXMul)

# --- Heatmap des combinaisons ---
melted_cormat <- reshape2::melt(DatXMul)
melted_cormat$Var2 <- as.character(melted_cormat$Var2) # sera refactorisé plus bas
melted_cormat$Var1 <- as.factor(melted_cormat$Var1)
melted_cormat$value <- as.factor(melted_cormat$value)

# --- Données issues de ResSCRENV (déjà calculé avant) ---
combo_res <- as.character(ResSCRENV20232024[[3]])
yinf <- ResSCRENV20232024[[8]]
ysup <- ResSCRENV20232024[[7]]
yobs <- ResSCRENV20232024[[6]]

# --- Ordre commun pour les deux graphiques ---
combo_levels <- combo_res  # on choisit l’ordre de ResSCRENV comme référence

# Re-factoriser melted_cormat avec ce même ordre
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = combo_levels)

# Préparer dataframe pour le graphique IC/obs
df_scn <- data.frame(
  combo = factor(combo_res, levels = combo_levels),
  yinf  = yinf,
  ysup  = ysup,
  yobs  = yobs,
  stringsAsFactors = FALSE
)

# --- Identifier les combinaisons hors enveloppes ---
valx <- which(ResSCRENV20232024[[4]] != 0) # indices dans combo_levels

# --- Graphique IC / Observations ---
scn <- ggplot(df_scn, aes(x = combo)) +
  geom_errorbar(aes(ymin = yinf, ymax = ysup), width = 0.5) +
  geom_point(aes(y = yobs, color = "Effectif observé")) +
  geom_vline(xintercept = valx, lty = 3, col = 1, lwd = 1) +
  theme_minimal() +
  labs(y = "Nombre d'individus", x = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_color_manual(values = c("Effectif observé" =
                                  viridis::viridis_pal(option = "rocket", direction = 1, begin = 0.5)(1)))

# --- Heatmap ---
h <- ggplot(data = melted_cormat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  coord_equal() +
  theme_minimal() +
  labs(y = NULL, x = "Combinaisons exclusives de parasites") +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = viridis::viridis_pal(option = "rocket", direction = -1, begin = 0.45, end = 1)(length(unique(melted_cormat$value))))

# --- Alignement final ---
aligned_plots <- cowplot::align_plots(scn, h, align = "v", axis = "x")
cowplot::plot_grid(plotlist = aligned_plots, ncol = 1,
                   rel_heights = c(2, 1), rel_widths = c(1, 1))


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


#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_barto_r)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m[!duplicated(data_for_m$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)


#Best model collinearity
performance::check_collinearity(m_barto_r)

#Best model look
summary(m_barto_r)
drop1(m_barto_r,.~.,test="Chisq") 


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
  formula = Mycoplasma_coccoides ~ scale(poids) +as.factor(year) + (1|numero_ligne),
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",                                  
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_mycoplasmacoco_r) %>%
  DHARMa::testResiduals()

#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_mycoplasmacoco_r)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m[!duplicated(data_for_m$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)


#Best model collinearity
performance::check_collinearity(m_mycoplasmacoco_r)

#Best model look
summary(m_mycoplasmacoco_r)
drop1(m_mycoplasmacoco_r,.~.,test="Chisq") 



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

#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_mycoplasma_r)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m[!duplicated(data_for_m$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)

#Best model collinearity
performance::check_collinearity(m_mycoplasma_r)

#Best model look
summary(m_mycoplasma_r)
drop1(m_mycoplasma_r,.~.,test="Chisq") 


## Model : Neoehrlichia mikurensis  ----

data_for_m %>%
  group_by(category, season) %>%
  summarise(n_infected = sum(Neoehrlichia_mikurensis >0),
            n = n(), .groups = "drop")

#Global model specification
if (exists("m_neoeh_r")) rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~ category + season + as.factor(year) + scale(poids) + scale(`Apodemus sylvaticus`) + sexe + (1|numero_ligne),
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

#Best model specification
if (exists("m_neoeh_r")) rm(m_neoeh_r)
m_neoeh_r <- lme4::glmer(
  formula = Neoehrlichia_mikurensis ~  season + as.factor(year) + scale(poids) + sexe + (1|numero_ligne),   # Attention deux meilleurs modèles
  family = binomial(link = "logit"),
  data = data_for_m,
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))
)

#Best model validation
DHARMa::simulateResiduals(m_neoeh_r) %>%
  DHARMa::testResiduals()

#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_neoeh_r)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m[!duplicated(data_for_m$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)


#Best model collinearity
performance::check_collinearity(m_neoeh_r)

em_neo <- emmeans::emmeans(m_neoeh_r, ~ category, type = "response")
plot(em_neo, comparisons = TRUE)
emmeans::contrast(em_neo, "pairwise", adjust = "Tukey")

#Best model look
summary(m_neoeh_r)
drop1(m_neoeh_r,.~.,test="Chisq") 



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

#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_mycoplasmasp)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m[!duplicated(data_for_m$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)

#Best model collinearity
performance::check_collinearity(m_mycoplasmasp)

#Best model look
summary(m_mycoplasmasp)
drop1(m_mycoplasmasp,.~.,test="Chisq") 


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

#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_barto_b_r)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m_2023224$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m_2023224[!duplicated(data_for_m_2023224$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)

#Best model look
summary(m_barto_b_r)
drop1(m_barto_b_r,.~.,test="Chisq") 




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

#Test spatial autocorrelation
#BEWARE, as in our case we have to recalculate residuals at the numero_ligne scale
#but our lines are missing sometimes from mission (ex n°6 in 2023 spring) it can mess
#up the order in which groups are accounted for in the below script,
#so the safest option is to reorder data used for modelling first according to the spatial unit 
#used for recalculateResiduals (here done at beggining of script)
# 1. Simulated residuals
sim_res <- DHARMa::simulateResiduals(m_barto_g_r)
# 2. Aggregate residuals per line
sim_res_line <- DHARMa::recalculateResiduals(sim_res, group = data_for_m_2023224$numero_ligne)
# 3. Recover order of group used for rescaling (it is there that problem can occur as we aggregate order based on order that may not be consistent)
group_order <- unique(sim_res_line$group)
# 4. Extract coordinate of line based on order
lines_sf <- data_for_m_2023224[!duplicated(data_for_m_2023224$numero_ligne),
                      c("numero_ligne", "longitude_lambert93", "latitude_lambert93")]
lines_sf <- lines_sf[match(group_order, lines_sf$numero_ligne), ]
# 5. Test spatial correlation
DHARMa::testSpatialAutocorrelation(
  sim_res_line,
  x = lines_sf$longitude_lambert93,
  y = lines_sf$latitude_lambert93
)

#Best model look
summary(m_barto_g_r)
drop1(m_barto_g_r,.~.,test="Chisq") 


# Richness temp ----

##individual scale ----

#recalculate richness values
data_for_m <- data_for_m %>%
  mutate(number_pathos = rowSums(across(all_of(list2025_10percent), ~ (. > 0))))  
hist(data_for_m$number_pathos)

data_for_m_2023224 <- data_for_m %>%
  mutate(number_pathos = rowSums(across(all_of(list20232024_10percent), ~ (. > 0))))  
hist(data_for_m_2023224$number_pathos)


#2025
#Global model specification
if (exists("m_rich2025")) rm(m_rich2025)
m_rich2025 <- glmmTMB::glmmTMB(
  number_pathos ~ category * season + as.factor(year) +
    scale(poids) + scale(`Apodemus sylvaticus`) + sexe +
    (1 | numero_ligne),
  data = data_for_m,
  family = glmmTMB::compois(),                  # Distribution Conway-Maxwell-Poisson
  na.action = "na.fail",
)


#Model selection
if (exists("model_selection")) rm(model_selection)
model_selection <- MuMIn::dredge(m_rich2025, rank = "AICc")
model_selection %>% filter(delta <2)


#Best model specification
if (exists("m_rich2025")) rm(m_rich2025)





#Best model validation
DHARMa::simulateResiduals(m_rich2025) %>%
  DHARMa::testResiduals()




##site_mission scale ----







# beta diversity temp ----
#on 2025 data

### Generate matrices (beta diversity) ----

#Subset community matrix
m_apo_pathos <- d_apo_pathos %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(all_of(list2025_10percent)) %>%                              #Only on 10% pathos just because was not working other
  as.matrix()

#Remove rows where all species are 0 (empty rows)
m_apo_pathos_clean <- m_apo_pathos[rowSums(m_apo_pathos) > 0, ]

#Generate dataframe with keptrows
kept_rows <- which(rowSums(m_apo_pathos) > 0)
metadata_clean <- d_apo_pathos[kept_rows, ]
rm(kept_rows)

#Generat Jaccard matrix
m_apo_jacc <- vegan::vegdist(m_apo_pathos_clean, method = "jaccard", binary = TRUE)

summary(rowSums(m_apo_pathos_clean))
summary(colSums(m_apo_pathos_clean))
nrow(unique(m_apo_pathos_clean))

#PERMANOVA
vegan::adonis2(formula = m_apo_jacc ~ category, data = metadata_clean) # ACTUELLEMENT CHOIX JACCARD CAR QUE DONNEES DETECTIONS 0/1

#BETADISPER
vegan::betadisper(m_apo_jacc, metadata_clean$category) %>%   # ACTUELLEMENT CHOIX JACCARD CAR QUE DONNEES DETECTIONS 0/1
  vegan::permutest(permutations = 9999)


### nMDS ----

#Jaccard nMDS on cleaned matrix
nmds_jacc <- vegan::metaMDS(
  comm = m_apo_pathos_clean,
  distance = "jaccard",
  k = 2,
  autotransform = FALSE
)

cat(paste0(
  "Final NMDS has a stress of ", round(nmds_jacc$stress, 3), " — ",
  dplyr::case_when(
    nmds_jacc$stress < 0.05 ~ "excellent goodness of fit (inferences very reliable)",
    nmds_jacc$stress < 0.1  ~ "good goodness of fit (inferences confident)",
    nmds_jacc$stress < 0.2  ~ "fair goodness of fit (some distances misleading)",
    TRUE                    ~ "poor goodness of fit (risks in interpretation)"
  ), "\n"
))

vegan::stressplot(nmds_jacc)
# Extract scores
nmdspoint <- vegan::scores(nmds_jacc)$sites %>%
  as_tibble(rownames = "numero_centre")
nmdsvariable <- vegan::scores(nmds_jacc)$species %>%
  as_tibble(rownames = "species")

nmdspoint %>% 
  left_join(metadata_clean) %>%
  mutate(category = forcats::fct_relevel(category, category_order)) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = as.factor(category))) +
  geom_point(size = 2) +
  stat_ellipse(show.legend = FALSE, linewidth = 2) +
  geom_text(data = nmdsvariable, 
            aes(x = NMDS1, y = NMDS2, label = species), 
            colour = "grey20") +
  scale_color_manual(values = category_palette) +
  theme_minimal()


nmdspoint %>% 
  left_join(metadata_clean) %>%
  mutate(category = forcats::fct_relevel(category, category_order)) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = as.factor(category))) +
  geom_jitter(width = 0.1, height = 0.1, size = 2, alpha = 0.8) +  # jitter + transparency
  stat_ellipse(show.legend = FALSE, linewidth = 2) +
  geom_text(data = nmdsvariable, 
            aes(x = NMDS1, y = NMDS2, label = species), 
            inherit.aes = FALSE,
            colour = "grey20") +
  scale_color_manual(values = category_palette) +
  theme_minimal()


## test db-rda ----


# refait pour que ce soit que sur data sans feuillus
# et jaccard, donc sans commu 0

datatest_dbrda <- d_apo_pathos %>%
  filter(category != "broadleaved_forest") %>%
  mutate(nonzero = rowSums(across(all_of(list2025_10percent))) > 0) %>%           #Attenion ici tableau provisoire pourritos
  filter(nonzero) %>%
  select(-nonzero) %>%
  tibble::column_to_rownames("numero_centre")

if (exists("db_rda")) rm(db_rda)
db_rda <- vegan::dbrda(
  datatest_dbrda[, list2025_10percent] ~ 
    bl_buff100 +
    PCA_lidarsuperb_axis1 +
    PCA_lidarsuperb_axis2 +
    PCA_microclim_axis1 +
    PCA_microclim_axis2 +
    plant_richness,
  data = datatest_dbrda,
  distance = "jaccard"
)


datatest_dbrda$
summary(db_rda)

# Test significance of model and terms
anova(db_rda, by = "term", permutations = 999)

# Plot
plot(db_rda)
text(db_rda, display = "sites", cex = 0.7)     # site labels
text(db_rda, display = "bp", col = "red")  

#remotes::install_github("gavinsimpson/ggvegan")

ggvegan::autoplot(db_rda, scaling = 2) +
  ggplot2::theme_minimal()









# db-RDA with environmental constraints
db_rda <- vegan::dbrda(formula = m_apo_jacc ~ 
                     herb_cover +
                     shrub_cover +
                     tree_cover,
                   data = metadata_clean)

# Test significance of model and terms
anova(db_rda, by = "term", permutations = 999)

# Plot
plot(db_rda)






names(metadata_clean)

# Summarise
summary(db_rda)


# test posthoc capscale : multiconstrained