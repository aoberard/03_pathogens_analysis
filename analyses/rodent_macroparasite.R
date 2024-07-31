#!/usr/bin/env Rscript

library("RPostgreSQL")
library("data.table")
library("dplyr")


# EXTRACT BPM  ----

# database and user parameters 
host <- "147.99.64.40" ; port <- 5432 ; db   <- "rongeurs" ; user <- "mus" ; pwd  <- "musculus" # read-only account (safe !)
conn <- dbConnect(PostgreSQL(), host=host, port=port, user=user, password=pwd, dbname=db)


sql_rodent_fast_parasite <- "SELECT mi.code_mission, 
lo.localite, 
li.date_2 AS date_pose, li.numero_ligne,
re.date_releve_2 AS date_releve,
pi.code_resultat, dis.numero_centre,
dis.sexe, dis.poids, par.effectif AS effectif_para, txp.taxon_name AS parasite
FROM t_missions AS mi
LEFT JOIN t_lignes AS li ON li.id_mission = mi.id
LEFT JOIN t_localites AS lo ON li.id_localite = lo.id
LEFT JOIN t_releves AS re ON re.id_ligne = li.id
LEFT JOIN t_pieges AS pg ON pg.id_ligne = li.id
LEFT JOIN t_piegeages AS pi ON (pi.id_piege = pg.id AND pi.id_releve = re.id)
LEFT JOIN t_types_resultats AS ttr ON pi.code_resultat = ttr.code_resultat
LEFT JOIN t_captures AS ca ON ca.id_piegeage = pi.id
LEFT JOIN t_taxons_nomenclature AS txnc ON ca.id_nomenclature = txnc.id
LEFT JOIN t_taxons AS txc ON txnc.id_taxon = txc.id
LEFT JOIN t_dissections AS dis ON dis.id_capture = ca.id
LEFT JOIN t_parasites AS par ON par.id_dissection = dis.id
LEFT JOIN t_taxons_nomenclature AS txnp ON par.id_nomenclature = txnp.id
LEFT JOIN t_taxons AS txp ON txnp.id_taxon = txp.id
WHERE mi.code_mission ~* 'BePrep'
AND pi.code_resultat = 1
ORDER BY mi.date_debut_2, li.numero_ligne, re.date_releve_2, pg.numero, dis.numero_centre
;"

# Get query (wait a few seconds)
bpm_fast_parasite <- as.data.table(dbGetQuery(conn, sql_rodent_fast_parasite))
# Save results
fwrite(bpm_fast_parasite, here::here("data/", "raw-data/", "raw-ticks", "rodents_tick", "export_fast_parasite_alois20240328.csv") )

dbDisconnect(conn)



# Data management ----

## Import data file ---- 
bpm_host_macroparasite <- data.table::fread(file = here::here( "data/", "raw-data/", "raw-ticks", "rodents_tick", "export_fast_parasite_alois20240328.csv"))










# DELETE DUPLICATE (WOULD NOT HAVE TO BE DONE IF QUERY WORKED CORRECTLY)
bpm_host_macroparasite <- bpm_host_macroparasite %>%
  distinct()

# Replace blank space in data
bpm_host_macroparasite$parasite[bpm_host_macroparasite$parasite == ""] <- NA

# Rename columns
bpm_host_macroparasite <- bpm_host_macroparasite %>%
  mutate( parasite = gsub(" ", "_", parasite))

## Generate tidy dataframe with macroparasite (0/1) ----
bpm_host_macroparasite_tidy <- bpm_host_macroparasite %>%
  filter( !is.na(parasite)) %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    names_from = parasite,
    values_from = value,
    values_fill = list(value = 0),
    id_cols = numero_centre
  )

# Add tick number when available
bpm_host_macroparasite_tidy <- left_join( x = bpm_host_macroparasite_tidy,
                                          y = bpm_host_macroparasite %>%
                                            filter(parasite == "Ixodida") %>%
                                            select( c("numero_centre", "effectif_para")),
                                          by = "numero_centre")

## Explore Macroparasite data from BPM ----

# How many ticks per season for rodent in BPM 
bpm_host_macroparasite %>%
  group_by(mission) %>%
  summarise(total_effectif_para = sum(effectif_para, na.rm = TRUE))

