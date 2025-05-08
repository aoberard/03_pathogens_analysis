# Script parameters ----

## Library ----
library("dplyr")
library("ggplot2")


## Graphical parameters ----

category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
category_palette <- c("pine_edge" = "#FFB3BA", "hedgerows" = "#FFDFBA", "broadleaved_forest" = "#B3E2B3")
type_order <- c("CT_LB", "CT_HB", "NC_LB", "NC_HB", "C_LB", "C_HB", "B" )
type_palette <- c("CT_LB" = "#FFB3BA", "CT_HB" = "#FFB3BA" , "NC_LB" = "#FFDFBA", "NC_HB" = "#FFDFBA" , "C_LB" = "#FFDFBA" , "C_HB"= "#FFDFBA" , "B" = "#B3E2B3")
brd_palette <- c("LB" = "#FFB3BA", "HB" = "#B3E2B3")
mission_order <- c("Juin_2023", "Octobre_2023", "Juin_2024", "Septembre_2024")
mission_color <- c("Juin_2023" = "#66c2a5", "Octobre_2023" = "#fc8d62", "Juin_2024" = "#8da0cb", "Septembre_2024" = "#ab7a82")

if (!requireNamespace("cividis", quietly = TRUE)) {
  devtools::install_github("marcosci/cividis")
}
cividis_palette <- cividis::cividis(2)


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

#Number infected per taxa
rodent_pathos %>%
  mutate(across(all_of(pathos_name), ~ replace(., . > 0, 1))) %>%
  group_by(taxon_mamm) %>%
  summarise(across(all_of(pathos_name), sum))

#List of pathogen per species
list_pathos_per_species <- list()

for (i in unique(rodent_pathos$taxon_mamm)) {
  list_pathos_per_species[[i]] <- rodent_pathos %>%
    filter(taxon_mamm == i) %>%
    select(all_of(pathos_name)) %>%
    select(where(~ any(. > 0))) %>%
    colnames()
}

list_pathos_per_species



# Apodemus data exploration ----

## Descriptive numbers ----

#PREVALENCE ALL MISSION
data_for_m %>%
  summarise(
    across(all_of(pathos_name_apo), list(
      prevalence = ~ sum(. > 0) * 100 / n(),
      positives = ~ sum(. > 0)
    )),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("pathogens", ".value"),
    names_pattern = "^(.*)_(prevalence|positives)$"
  ) %>% 
  arrange(desc(positives))

#PREVALENCE PER MISSION
data_for_m %>%
  group_by(code_mission) %>%
  summarise(across(all_of(pathos_name_apo), ~ sum(.) * 100 / n()), .groups = "drop") 

#Pathogen richness 


# Calculate percentage
percent_at_least_2 <- mean(data_for_m$number_pathos >= 2) * 100
label_text <- paste0(round(percent_at_least_2, 1), "% with ≥2 pathogens")


## Graph divers anciens ----

#PREVALENCES PER SITE (beware, sometiques really few individuals per site rise prevalence + line with no individuals are not integrated)
#Init pathogen prevalence per grouping factor 
grouping_prevalence_factor <- c("numero_ligne")

# Number of positive individuals per grouping factor
pp <- data_for_m %>%
  group_by(across(all_of(grouping_prevalence_factor))) %>%
  summarise(
    effectif = n(),
    category = unique(category),
    broadleaved_class = unique(broadleaved_class),
    treatment = unique(treatment),
    type = unique(type),
    across(all_of(pathos_name_apo), ~ sum(. > 0)),  
    .groups = "drop"
  )

# Calculate the total number of individuals and prevalence for each pathogen
pp_prevalence <- pp %>%
  group_by(across(all_of(grouping_prevalence_factor))) %>%
  summarise(
    category = unique(category),
    broadleaved_class = unique(broadleaved_class),
    treatment = unique(treatment),
    type = unique(type),
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
  ggplot(aes(x = prevalence, y = pathos, fill = category)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.6) +
  geom_jitter(aes(color = category), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = category_palette) +
  scale_color_manual(values = category_palette) + 
  labs(x = "Prevalence par ligne de piegeage", y = "Pathogènes") +
  guides(fill = guide_legend(title = "category de ligne"), 
         color = guide_legend(title = "category de ligne")) +
  theme_minimal()

pp_prevalence %>%  
  ggplot(aes(x = prevalence, y = pathos, fill = type)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.6) +
  geom_jitter(aes(color = type), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = type_palette) +
  scale_color_manual(values = type_palette) + 
  labs(x = "Prevalence par ligne de piegeage", y = "Pathogènes") +
  guides(fill = guide_legend(title = "type de ligne"), 
         color = guide_legend(title = "type de ligne")) +
  theme_minimal() +
  stat_summary(fun.data = function(x) count_summary(x, y_position = 1.1), 
               aes(group = type),
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


data_for_m %>% count(Neoehrlichia_mikurensis, category)

ggplot(data_for_m, aes(x = as.factor(Mycoplasma_haemomuris), y = poids)) +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = sexe),
               dotsize = 0.7 ,
               position=position_dodge(0.4)) + 
  coord_flip() +
  stat_summary(fun = mean, geom="point", shape=18,
               size=5, color="purple") +
  theme_minimal()

hist(data_for_m$number_pathos)






## Graph heatmap matrix ----

# Prevalence matrix calculation (only apodemus)
#for line_treatment

matrix_pathos <- d_apo_pathos %>%
  group_by(code_mission) %>%
  summarise(
    across(all_of(pathos_name_apo),
           ~ round(sum(. > 0) / n(), digits = 2)) ) |>
  arrange(factor(code_mission, levels = mission_order)) |>
  tibble::column_to_rownames("code_mission") |>
  select(all_of(pathos_name_apo)) |>
  as.matrix()


effectif_df <- d_apo_pathos %>%
  group_by(code_mission) %>%
  summarise(
    effectif = n(),
  ) %>%
  arrange(factor(code_mission, levels = mission_order))

# Check the alignment by comparing combined names
if(!all(rownames(matrix_pathos) == effectif_df$code_mission)) {
  stop("Mismatch between matrix_pathos row names and effectif combined names")
}

effectif <- effectif_df %>%
  pull(effectif) 




superheat::superheat(X = matrix_pathos,
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



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")


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
ComplexHeatmap::Heatmap(matrix_pathos,
                        name = "Prévalence",
                        row_names_side = "left",
                        col = cividis_palette,
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        cluster_columns = TRUE,
                        show_column_dend = FALSE,
                        cluster_rows = FALSE,
                        row_order = mission_order,
                        right_annotation = row_ha,
) 



## Diversity analysis - to move later ----

### Graph richness ----

# Calculate percentage
percent_at_least_2 <- mean(data_for_m$number_pathos >= 2) * 100
label_text <- paste0(round(percent_at_least_2, 1), "% with ≥2 pathogens")

# Create histogram with annotation
pathos_richness_plot <- ggplot(data_for_m, aes(x = number_pathos)) +
  geom_histogram(binwidth = 1, fill = "#fddc64", color = "#7a6a45", linewidth = 1.2) +
  labs(
    title = "Distribution of Pathogen Richness per Individual",
    x = "Number of Pathogens Detected",
    y = "Number of Individuals"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.1, vjust = 2,
    size = 5,
    fontface = "italic",
    color = "#7a6a45"
  )
pathos_richness_plot

ggsave(filename = here::here("figures","pathos_richness_plot.pdf"),
       plot = pathos_richness_plot, 
       width = 1600, height = 1500, units = "px", device = "pdf", dpi = 300, bg = NULL)


### Generate matrices (beta diversity) ----

#Subset community matrix
m_apo_pathos <- d_apo_pathos %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(all_of(patho10_apo)) %>%                              #Only on 10% pathos just because was not working other
  as.matrix()

#Remove rows where all species are 0 (empty rows)
m_apo_pathos_clean <- m_apo_pathos[rowSums(m_apo_pathos) > 0, ]

#Generate dataframe with keptrows
kept_rows <- which(rowSums(m_apo_pathos) > 0)
metadata_clean <- d_apo_pathos[kept_rows, ]
rm(kept_rows)

#Generat Jaccard matrix
m_apo_jacc <- vegan::vegdist(m_apo_pathos_clean, method = "jaccard", binary = TRUE)

#Generate Bray-Curtis matrix
m_apo_bray <- vegan::vegdist(m_apo_pathos_clean, method = "bray" , binary = FALSE)

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
  stat_ellipse(show.legend = FALSE, size = 2) +
  geom_text(data = nmdsvariable, 
            aes(x = NMDS1, y = NMDS2, label = species), 
            colour = "grey20") +
  scale_color_manual(values = category_palette) +
  theme_minimal()








































# -------------BELOW is TEST GRAPH FOR EWDA POSTER -------

category <- c("pine_edge" = "#8FC08C", "hedgerows" = "#A95738", "broadleaved_forest" = "#B3E2B3")


# Generate a darker version of the fill colors
darker_palette <- sapply(category, darken_color, amount = 0.25)

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
  scale_fill_manual(values = c(category, darker_palette)) + 
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
  scale_fill_manual(values = c(category, darker_palette)) + 
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
  scale_fill_manual(values = c(category, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall prevalence per line", y = "Line category") +
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
  scale_fill_manual(values = c(category, darker_palette)) + 
  scale_color_manual(values = darker_palette) + 
  labs(x = "Overall prevalence per line", y = "Line category") +
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
  scale_fill_manual(values = c(category, darker_palette)) + 
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


