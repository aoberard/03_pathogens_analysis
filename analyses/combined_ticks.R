

# Data management  ----

## Imports ----
# Import hosts and line modalities file
hosts <- readr::read_csv( here::here("data/", "raw-data/", "host_data", "20240731_bpm_modalities.csv") )

# Import macroparasites from host file
macroparasite <- fread(file = here::here("data/", "derived-data/", "raw-ticks", "rodents_tick", "20240731_macroparasite.csv") )

# Import Collected ticks information
#


## Join ----

# Convert columns to character type (if necessary) and filter hosts data
hosts_processed <- hosts %>%
  filter(code_resultat == 1) %>%
  select(numero_centre, taxon_mamm, numero_ligne, code_mission, year, season, connectivity, broadleaved_status, line_treatment, line_type) %>%
  mutate(across(everything(), as.character))

# Convert columns to character type (if necessary) for macroparasite data
macroparasite_processed <- macroparasite %>%
  mutate(across(everything(), as.character))

# Perform the full join
host_macro <- full_join(macroparasite_processed, hosts_processed)


# Transform NA effectif_tick into 0, and make it numeric
host_macro$effectif_tick[ is.na(host_macro$effectif_tick)] <- 0
host_macro$effectif_tick <- as.numeric(host_macro$effectif_tick)


# Graph try ----




type_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
type_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
treatment_order <- c("CT-LB", "CT-HB", "NC-LB", "C-LB", "NC-HB", "C-HB", "B" )
mission_order <- c("Juin 2023", "Octobre 2023", "Juin 2024")
mission_palette <- c("Juin 2023" = "#A7DBA7", "Octobre 2023" = "#FFB07C", "Juin 2024" = "#99D8C9")

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


host_macro %>%
  filter(code_mission %in% c("Juin 2023", "Octobre 2023") )  %>%
  ggplot( aes (x = factor(line_treatment, levels = treatment_order), y = log(effectif_tick+1), fill = line_type)) +
  facet_grid(~ factor(code_mission, levels = mission_order) ) +
  geom_boxplot(position = position_dodge(1)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_jitter(width = 0.0001), dotsize = 0.5) +
  scale_fill_manual(values = type_palette) +
  labs(x = "Type de ligne", y = "Nbr de tiques collectées sur les hôtes") +
  guides(fill = guide_legend(title = "Type de ligne")) +
  theme_minimal()

# Statistical analyses




# USE INFORMATIONS from rodent tick TO COMBINED GRAPH USE ----

## Data management ----

# Integrate line modality information from site_tique
macroparasite <- left_join(x = macroparasite %>%
                             mutate(numero_ligne = as.character(numero_ligne)), 
                           y = site_tique %>%
                             select(line, connectivity, broadleaved_status, line_treatment, line_type) %>%
                             distinct(), 
                           by = c("numero_ligne" = "line") )



# GROS BLEMPRO, y'a que les rongeurs avec des parasites ici EUHHHHH les 0 y sont po là lol n'importe quoi aloies
