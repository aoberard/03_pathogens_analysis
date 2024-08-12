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



# Data control ----

# /!\ The final data file is currently made outside of R :

## Import data file ---- 
site_tique <- readxl::read_excel(here::here("data", "raw-data", "raw-ticks", "collect", "20240730_collect_tick.xlsx"))





