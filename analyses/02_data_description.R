# Script parameters ----

## Library ----
library("dplyr")
library("ggplot2")


## Graphical parameters ----

category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
category_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
type_order <- c("CT_LB", "CT_HB", "NC_LB", "NC_HB", "C_LB", "C_HB", "B" )
brd_palette <- c("LB" = "#FFB3BA", "HB" = "#B3E2B3")
mission_order <- c("Juin 2023", "Octobre 2023", "Juin 2024", "Septembre 2024")
mission_color <- c("Juin 2023" = "#66c2a5", "Octobre 2023" = "#fc8d62", "Juin 2024" = "#8da0cb", "Septembre 2024" = "#ab7a82")


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
    select(all_of(pathos_name)) %>%
    select(where(~ any(. > 0))) %>%
    colnames()
}

list_pathos_per_species





































# Exploration (for apodemus only) ----

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


