#!/usr/bin/env Rscript

library("RPostgreSQL")
library("data.table")
library("dplyr")
library("ggplot2")

# EXTRACT BPM  ----

# database and user parameters 
host <- "147.99.64.40" ; port <- 5432 ; db   <- "rongeurs" ; user <- "mus" ; pwd  <- "musculus" # read-only account (safe !)
conn <- dbConnect(PostgreSQL(), host=host, port=port, user=user, password=pwd, dbname=db)


sql_collected_ticks <- "SELECT mi.code_mission, 
       lo.localite,
       li.numero_ligne,
       MAX(co.date_collecte) AS date_collecte,
       MAX(co.observations) AS observations
FROM t_collectes AS co
LEFT JOIN t_localites AS lo ON co.id_localite = lo.id
LEFT JOIN t_missions AS mi ON co.id_mission = mi.id
LEFT JOIN t_lignes AS li ON (li.id_mission = mi.id AND li.id_localite = lo.id)
WHERE mi.code_mission ~* 'BePrep'
GROUP BY mi.code_mission, lo.localite, li.numero_ligne
ORDER BY date_collecte DESC, li.numero_ligne;"


# Get query (wait a few seconds)
bpm_site_ticks <- as.data.table(dbGetQuery(conn, sql_collected_ticks))

# Save results
xlsx::write.xlsx(bpm_site_ticks, here::here("data/", "raw-data/", "raw-ticks", "collect", "export_site_ticks_alois20240730.xlsx") )

dbDisconnect(conn)




# Data management ----

## Import data file ---- 
site_tique <- readxl::read_excel(here::here("data", "raw-data", "raw-ticks", "collect", "20240730_collect_tick.xlsx"))

## Attributing line modality ----

### Line characteristics  ----

# Glossary:
# C: Connected hedgerow / NC: Non-connected hedgerow
# CT: Control pine edge
# LB: Low Broadleaved density / HB: High Broadleaved
# B : Broadleaved forests

C <- c(2, 6, 20, 24, 30, 7, 11, 12, 28, 33, 36, 17)
CT <- c(8, 14, 18, 23, 26, 34, 5, 10, 15, 22, 32, 35)
NC <- c(3, 4, 13, 25, 27, 29, 1, 9, 16, 19, 21, 31)

LB <- c(7, 11, 12, 28, 33, 36, 5, 10, 15, 22, 32, 35, 1, 9, 16, 19, 21, 31)
HB <- c(17, 2, 6, 20, 24, 30, 8, 14, 18, 23, 26, 34, 3, 4, 13, 25, 27, 29)

B <- c(37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53)

site_tique <- site_tique %>%
  mutate(
    connectivity = case_when(
      line %in% C ~ "C",
      line %in% CT ~ "CT",
      line %in% NC ~ "NC",
      TRUE ~ NA
    ),
    broadleaved_status = case_when(
      line %in% LB ~ "LB",
      line %in% HB ~ "HB",
      TRUE ~ NA
    ),
    line_treatment = case_when(     
      line %in% C  ~ paste0(connectivity, "-", broadleaved_status),
      line %in% CT  ~ paste0(connectivity, "-", broadleaved_status),
      line %in% NC  ~ paste0(connectivity, "-", broadleaved_status),
      line %in% B  ~ "B",
      TRUE ~ NA
    )
  ) %>%
  mutate(
    line_type = case_when(
      line %in% CT ~ "pine_edge",
      line %in% B ~ "broadleaved_forest",
      line %in% c(NC, C) ~ "hedgerows",
      TRUE ~ "unknown"
    )
  )

### Mission characteristics ----

#generate year and month variable
site_tique <- site_tique %>%
  mutate(
    year = stringr::str_extract(string = mission, pattern = "\\d*$"),
    month = stringr::str_extract(string = mission, pattern = "(?<=\\s-\\sBePrep\\s-\\s)\\w+")
  )

#choose which month correspond to each season
Spring <- c("Juin")
Autumn <- c("Octobre")

#attribute season and mission code accordingly
site_tique <- site_tique %>%
  mutate(
    season = case_when(
      month %in% Spring ~ "Spring",
      month %in% Autumn ~ "Autumn",
      TRUE ~ "unknown"
    )
  ) %>%
  relocate(
    year, .after = which(names(site_tique) == "mission") 
  ) %>%
  relocate(
    season, .after = which(names(site_tique) == "mission") +1
  ) %>%
  mutate(
    mission = paste(month, year)
  ) %>%
  select(!month)

# Error check for modality attribution
unknown_rows <- site_tique %>% filter(is.na(line_treatment) | line_type == "unknown" | season == "unknown" )
if (nrow(unknown_rows) > 0) {
  warning("There are rows with unknown values. Check 'unknown_rows' for details.")
}


# Data exploration ----
site_tique %>%
  group_by(mission) %>%
  filter(!is.na(x = line)) %>%
  summarise(
    tick_number = sum(tick_number)
  )

hist(site_tique$tick_number)

graph_sequence <- c("pine_edge", "hedgerows", "broadleaved_forest")
color_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
mission_order <- c("Juin 2023", "Octobre 2023", "Juin 2024")
site_tique %>%
  ggplot(aes (x = factor(line_type, levels = graph_sequence), y = tick_number, , fill = line_type)) + 
  facet_grid(~ factor(mission, levels = mission_order) ) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_jitter(width = 0.05), dotsize = 0.5) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Type de ligne", y = "Nbr de tiques collectées au drap par ligne") +
  guides(fill = guide_legend(title = "Type de ligne")) +
  theme_minimal()


         

# Statistical analysis ----

# Making data frame for glm model /!\ spring 2023 only for poster EWDA
d_site_tique_glm <- site_tique %>%
  filter(mission == "Juin 2023")

## GLM making (poisson) ----
rm(m_site_tique)
m_site_tique <- lme4::glmer(
  formula = tick_number ~ connectivity * broadleaved_status + (1|line),
  data = d_site_tique_glm,
  family = poisson(link = "log"),  
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)

DHARMa::simulateResiduals(m_site_tique, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)

m_site_tique |>
  RVAideMemoire::overdisp.glmer()

SelectionModels<- MuMIn::dredge(m_site_tique, rank = "AICc")              
TopModels<-subset(SelectionModels, delta<2)
TopModels

rm(m_site_tique)
m_site_tique <- lme4::glmer(
  formula = tick_number ~ connectivity + (1|line),
  data = d_site_tique_glm,
  family = poisson(link = "log"), 
  na.action = "na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7))
)
DHARMa::simulateResiduals(m_site_tique, n = 250, refit = F, integerResponse = NULL, plot = T, seed = 123) |>
  plot(rank = TRUE)
drop1(m_site_tique,.~.,test="Chisq")
summary(m_site_tique)

em <- emmeans::emmeans(m_site_tique, "connectivity", type = "response")
em
plot(em, comparisons = TRUE)

emmeans::contrast(em, "pairwise", adjust = "Tukey")


