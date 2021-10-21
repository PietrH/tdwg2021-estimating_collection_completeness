# Databricks notebook source
# MAGIC %md 
# MAGIC 
# MAGIC # Poster Figures using rarefaction
# MAGIC 
# MAGIC Making two figures:
# MAGIC 
# MAGIC 1. comparing a number of families held by meise against the rest of the world
# MAGIC 1. comparing the african taxa held by meise with the african taxa held by the world

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC setup env

# COMMAND ----------

# load libraries

library(sparklyr)
if(!require("countrycode")){install.packages("countrycode")}
if(!require("iNEXT")){install.packages("iNEXT")}
library(iNEXT)
if(!require("furrr")){install.packages("furrr")}

# COMMAND ----------

# declare spark cluster
sc <- sparklyr::spark_connect(method = "databricks")

# COMMAND ----------

# get a country to continent lookup table to spark
sparklyr::copy_to(sc,countrycode::codelist,"codelist")

# COMMAND ----------

# create spark dataframe with only columns we need
codelist_spark <- sparklyr::sdf_sql(sc,'SELECT iso2c,continent FROM codelist')

# COMMAND ----------

# lets put the occurrences in cache for speed later on (takes a while), this might not be worth the wait. It also takes a while for the cluster to scale up.
tbl_cache(sc,"pieter.occurrence")

# COMMAND ----------

# reading in the taxon table (gbif taxonomic backbone)

taxon_table_spark <- sparklyr::spark_read_csv(sc,"dbfs:/FileStore/tables/Taxon/Taxon_tsv.gz",delimiter="\t")

dplyr::glimpse(head(taxon_table_spark,50))

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC ### Decide if we want to calculate curves individually
# MAGIC 
# MAGIC Depending on the size of the abundance datasets (amount of taxa in a family), and the amount of families, you might want to calculate every cuve individually instead of all at once. 

# COMMAND ----------

calculate_curves_individually = FALSE

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC ## Comparing a number of families

# COMMAND ----------

# setting the families of interest

interesting_families <- c(7681,
                          8798,
                          4170,
                          6669) # all families are interesting really
family_info <- filter(taxon_table_spark,taxonID %in% interesting_families)
family_names <- dplyr::pull(family_info,canonicalName)


# lets make sure we have the right families selected
family_info %>%
  SparkR::as.data.frame() %>%
  display

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC Let's get some numbers for these families just to check our filters, we expect around 250 species for the Connaraceae (APG3: 180, CoL:240, GBIF: 289)

# COMMAND ----------


sdf_sql(sc,
          'SELECT acceptedtaxonkey,familykey FROM pieter.occurrence WHERE familykey = 8798 OR familykey = 4170 OR familykey = 6669 OR familykey = 7681') %>%
  sparklyr::left_join(select(taxon_table_spark,taxonID, taxonRank,canonicalName,kingdom,scientificName,taxonomicStatus), by = c("acceptedtaxonkey" = "taxonID")) %>%
    filter(taxonomicStatus == "accepted") %>%
    #filter(kingdom == "Plantae") %>%
    filter(taxonRank == "species") %>%
    filter(familykey == 6669) %>%
    dplyr::group_by(acceptedtaxonkey) %>%
    dplyr::tally() %>%
    sparklyr::sdf_nrow()


# COMMAND ----------

# MAGIC %md 
# MAGIC 
# MAGIC ### For Meise
# MAGIC 
# MAGIC Let's do the taxa within these families held by Meise Botanic Garden Herbarium first, since that's a smaller set to work on

# COMMAND ----------

# only accepted taxa, species level, from the herbarium (including fungi because we'll need them!)

# query occurrence table via function so we don't repeat ourselves when we do this again for the whole of gbif
get_occ <- function(only_meise=FALSE){
  sdf_sql(sc,
          sprintf("SELECT acceptedtaxonkey,familykey FROM pieter.occurrence%s",
                  ifelse(only_meise,
                         ' WHERE datasetkey = "b740eaa0-0679-41dc-acb7-990d562dfa37"',
                         ' WHERE datasetkey != "b740eaa0-0679-41dc-acb7-990d562dfa37"')
                 )
         ) %>%
    sparklyr::left_join(select(taxon_table_spark,taxonID, taxonRank,canonicalName,kingdom,scientificName,taxonomicStatus), by = c("acceptedtaxonkey" = "taxonID")) %>%
    filter(taxonomicStatus == "accepted") %>%
    #filter(kingdom == "Plantae") %>% # include Fungi and others! 
    filter(taxonRank == "species") %>%
    filter(familykey %in% interesting_families)
}


meise_occ <- get_occ(only_meise = TRUE)


# COMMAND ----------

# function to get abundance counts for a certain family key

get_abundance <- function(iter_family_key, occ_table){

filter(occ_table,familykey == iter_family_key) %>%
  dplyr::group_by(acceptedtaxonkey) %>%
  dplyr::tally() %>%
  dplyr::pull(n)
}

# COMMAND ----------

inext_conn <- iNEXT(get_abundance(6669,meise_occ),q = 0, datatype = "abundance") #calc rarefaction curve for Connaraceae

# COMMAND ----------

# plot for the Connaraceae
ggiNEXT(inext_conn, type = 1)

# COMMAND ----------

# create a named list of abundance counts for our interesting families

#meise_abundance <- purrr::map(interesting_families,get_abundance,meise_occ)
#names(meise_abundance) <- family_names

meise_abundance <- purrr::map(dplyr::pull(family_info,taxonID),get_abundance,meise_occ)
names(meise_abundance) <- dplyr::pull(family_info,canonicalName)


# COMMAND ----------

if(calculate_curves_individually){

library(furrr)

plan("multisession", workers = 8) # double the amount of threads of a single worker, doesn't propagate over the cluster

meise_inext <- furrr::future_map(meise_abundance, function(x){
  inext_out <-
    iNEXT(x,
       q = 0,
       datatype = "abundance" #,endpoint = calc_endpoint
       )
  out_list <-
    list(DataInfo = inext_out$DataInfo,
       iNextEst = inext_out$iNextEst,
       AsyEst = inext_out$AsyEst)
  return(out_list)
},
                    .options = furrr_options(seed = 123)
          )
names(meise_inext) <- names(meise_abundance)}

# COMMAND ----------

## plotting one of the calculated curves

if(calculate_curves_individually){

inext <- meise_inext[[3]]
class(inext) <- "iNEXT"

ggiNEXT(inext, type = 1)}



# COMMAND ----------

# plotting multiple curves
inext_meise <-
  iNEXT(meise_abundance,
        q = 0,
        endpoint = 2 * max(unlist(purrr::map(meise_abundance,sum))), # stop extrapolating at double the number of occurrences of the largest family
        conf = 0.99
       )
ggiNEXT(inext_meise, type = 1)

# COMMAND ----------

#comparing 4 assemblages

inext_meise

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC ### For the whole of GBIF

# COMMAND ----------

gbif_occ <- get_occ(only_meise = FALSE)

# COMMAND ----------

# create a named list of abundance counts for our interesting families
gbif_abundance <- purrr::map(dplyr::pull(family_info,taxonID),get_abundance,gbif_occ)
names(gbif_abundance) <- dplyr::pull(family_info,canonicalName)

# COMMAND ----------

# calculating rarefaction curves individually, we don't need to do this if the cluster can handle them all at once. 

if(calculate_curves_individually){

# we set our future plan earlier
gbif_inext <- furrr::future_map(gbif_abundance, function(x){
  inext_out <-
    iNEXT(x,
       q = 0,
       datatype = "abundance" #,endpoint = calc_endpoint
       )
  out_list <-
    list(DataInfo = inext_out$DataInfo,
       iNextEst = inext_out$iNextEst,
       AsyEst = inext_out$AsyEst)
  return(out_list)
},
                    .options = furrr_options(seed = 123)
          )
names(gbif_inext) <- names(gbif_abundance)}

# COMMAND ----------

# plotting multiple curves
inext_gbif <- 
  iNEXT(gbif_abundance,
        q = 0,
        endpoint = 2 * max(unlist(purrr::map(meise_abundance,sum))), # stop extrapolating at double the number of occurrences of the largest family
        conf = 0.99
       )
ggiNEXT(inext_gbif,type = 1)

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC ### Plotting interesting families in Meise vs GBIF on the same graph

# COMMAND ----------

# combining the abundance lists, setting names to family names underscore gbif or meise

both_abundance <- append(purrr::map(dplyr::pull(family_info,taxonID),get_abundance,meise_occ),
                         purrr::map(dplyr::pull(family_info,taxonID),get_abundance,gbif_occ))
#both_abundance <- append(meise_abundance,gbif_abundance)

names(both_abundance) <- unlist(purrr::map(c("meise","gbif"),function(x,y){paste(y,x,sep = "_")},dplyr::pull(family_info,canonicalName)))

# COMMAND ----------

both_inext <- 
  iNEXT(both_abundance,
        q = 0,
        endpoint = 2 * max(unlist(purrr::map(meise_abundance,sum))), # stop extrapolating at double the number of occurrences of meise
        conf = 0.99
       )

both_inext_no_endpoint <- 
    iNEXT(both_abundance,
        q = 0,
        conf = 0.99
       )

# COMMAND ----------

# inspect abundance table
dplyr::glimpse(both_abundance)

# COMMAND ----------

ggiNEXT(both_inext,type = 1)

# COMMAND ----------

ggiNEXT(both_inext,type = 2)

# COMMAND ----------

# MAGIC %md 
# MAGIC 
# MAGIC ## All African Plant Taxa: Meise vs GBIF except Meise. For equal collection size, who's more rich?

# COMMAND ----------

# Get all species level records, for plants, for both meise, and all of gbif-except-meise

get_all_occ <- function(only_meise=FALSE){
  sdf_sql(sc,
          sprintf("SELECT acceptedtaxonkey,familykey,countrycode FROM pieter.occurrence%s",
                  ifelse(only_meise,
                         ' WHERE datasetkey = "b740eaa0-0679-41dc-acb7-990d562dfa37"',
                         ' WHERE datasetkey != "b740eaa0-0679-41dc-acb7-990d562dfa37"')
                 )
         ) %>%
    sparklyr::left_join(select(taxon_table_spark,taxonID, taxonRank,canonicalName,kingdom,scientificName,taxonomicStatus), by = c("acceptedtaxonkey" = "taxonID")) %>%
    filter(taxonomicStatus == "accepted") %>%
    filter(kingdom == "Plantae") %>%
    filter(taxonRank == "species") %>%
    sparklyr::left_join(codelist_spark, by = c("countrycode" = "iso2c")) %>%
    filter(continent == "Africa")
  }

gbif_all_occ <- get_all_occ(only_meise = FALSE)
meise_all_occ <- get_all_occ(only_meise = TRUE)



# COMMAND ----------

# convert to abundance table

get_abundance_scientific_name <- function(occurrence_table){
  occurrence_table %>%
    dplyr::group_by(acceptedtaxonkey) %>%
    dplyr::tally() %>%
    dplyr::pull(n)
}

gbif_all_abundance <- get_abundance_scientific_name(gbif_all_occ)
meise_all_abundance <- get_abundance_scientific_name(meise_all_occ)

# COMMAND ----------

both_all_inext <- iNEXT(
  list(gbif = gbif_all_abundance, meise = meise_all_abundance),
  q = 0,
  endpoint =  1.5e8,
  nboot = 10, #we should increase this for the final figure
  conf = 0.99
)

# COMMAND ----------

both_all_inext_plot <- 
  ggiNEXT(both_all_inext, type = 1)

# COMMAND ----------

both_all_inext_plot

# COMMAND ----------

# function to extract the highest extrapolated richness value from a facet from an inext object
get_max_extrapolated <- function(inext_result, facet){
  inext_result$iNextEst[facet] %>% 
  data.table::rbindlist() %>% 
  dplyr::pull(qD) %>% 
  max()
}


