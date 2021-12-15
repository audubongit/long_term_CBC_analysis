
# Long-term Christmas Bird Count Analysis -----------------------------------------------------------
# Code developed by T. Meehan and S. Saunders 2020 - 2021
# Analysis associated with Global Change Biology paper entitled "Unraveling a century of global change
# impacts on winter bird distributions in the eastern United States"
# Time period: 1930s to 2010s
# Analysis: Generalized additive mixed models (GAMMs) and analysis of deviance (ANODEV)
# Analysis: 9 bird groups comprising 89 species
#  ------------------------------------------------------------------------------------------------


# set up -----------------------------------------------------------------------
# libraries
library(mapview)
library(rgeos)
library(RODBC)
library(psych)
library(sf)
library(mgcViz)
library(cowplot)
library(mgcv)
library(dplyr)
library(tidyr)

# options
options(scipen = 9999999)

# plot theme
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# model selction
slcpo <- function(m, na.rm = TRUE) {
  -2 * sum(log(m$cpo$cpo), na.rm = na.rm)
}

# helpers
diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in%
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in%
                                                          colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp))
    stop('ramp should be either the name of a RColorBrewer palette, ',
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p[[grep('^legend', names(p))]][[1]]$args$key$col <-
    ramp(1000)[zlim[-length(zlim)]]
  p$panel.args.common$col.regions <- ramp(1000)[zlim[-length(zlim)]]
  p
}

# set wd
root_dir <- "~/GitHub/QuantitativeMetrics/vlt_cbc"
setwd(root_dir)
# ------------------------------------------------------------------------------


# get site by decade covariate data --------------------------------------------
# site data
focal_sites <- read_sf("deep_cbc_sites.shp" ) %>%
  rename(full_circle=site)
focal_area <- read_sf("study_area.shp" ) %>%
  dplyr::select(country=ISO)
focal_states <- read_sf("study_states.shp" ) %>%
  dplyr::select(country=ISO)
mapview(focal_area) + mapview(focal_states) + mapview(focal_sites)

# site and effort data
effort_data <- read.csv("deep_cbc_sites_effort_by_year.csv")

# tmin data
tmin_per_site <- read.csv("VLT_CBC_climate_summary.csv") %>%
  dplyr::select(-POINT_FID, -lon, -lat, -geometry) %>% 
  pivot_longer(-site, names_to="clim_par", values_to="value") %>%
  separate(clim_par, into = c("clim_par", "year"), sep="_") %>%
  mutate(year=as.numeric(gsub(pattern = "s", "", year)) + 9,
         value=as.numeric(value)) %>%
  filter(clim_par=="tmin") %>%
  dplyr::select(full_circle=site, year, tmin=value) %>%
  arrange(full_circle, year)

# precip data
precip_per_site <- read.csv("VLT_CBC_climate_summary.csv") %>%
  dplyr::select(-POINT_FID, -lon, -lat, -geometry) %>% 
  pivot_longer(-site, names_to="clim_par", values_to="value") %>%
  separate(clim_par, into = c("clim_par", "year"), sep="_") %>%
  mutate(year=as.numeric(gsub(pattern = "s", "", year)) + 9,
         value=as.numeric(value)) %>%
  filter(clim_par=="ppt") %>%
  dplyr::select(full_circle=site, year, precip=value) %>%
  arrange(full_circle, year)

# lulc data
prop_crop <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="crop") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, POINT_FID, prop_crop=MEAN, year=YEAR)
prop_decid <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="forest-decid") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_decid=MEAN, -YEAR) %>%
  cbind(prop_crop)
prop_evrgn <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="forest-evrgn") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_evrgn=MEAN, -YEAR) %>%
  cbind(prop_decid)
prop_mixed <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="forest-mixed") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_mixed=MEAN, -YEAR) %>%
  cbind(prop_evrgn)
prop_grass <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="grass") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_grass=MEAN, -YEAR) %>%
  cbind(prop_mixed)
prop_pasture <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="pasture") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_pasture=MEAN, -YEAR) %>%
  cbind(prop_grass)
prop_shrub <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="shrub") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_shrub=MEAN, -YEAR) %>%
  cbind(prop_pasture)
prop_urban <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="urban") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_urban=MEAN, -YEAR) %>%
  cbind(prop_shrub)
prop_water <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="water") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_water=MEAN, -YEAR) %>%
  cbind(prop_urban)
prop_herbs <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="wetland-herbs") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_herbs=MEAN, -YEAR) %>%
  cbind(prop_water)
prop_woody <- read.csv("VLT_CBC_LULC_summary.csv") %>% 
  filter(CLASS=="wetland-woody") %>%
  dplyr::select(-COUNT, -AREA, -CLASS, -POINT_FID, prop_woody=MEAN, -YEAR) %>%
  cbind(prop_herbs) %>% dplyr::select(POINT_FID, year, everything())
lulc_props <- prop_woody %>% 
  dplyr::select(POINT_FID, year, everything()) %>% 
  mutate(sum_props=rowSums(.[3:13])) %>% 
  arrange(POINT_FID, year) %>%
  mutate(full_circle=precip_per_site$full_circle,
         year2=precip_per_site$year) %>%
  dplyr::select(full_circle, year,  everything(), -POINT_FID, -year2)

# combine like types
lulc_props <- lulc_props %>% 
  mutate(prop_decid_evrgn_mixed=prop_decid+prop_evrgn+prop_mixed, 
         prop_decid_evrgn_mixed_shrub=prop_decid+prop_evrgn+prop_mixed+prop_shrub, 
         prop_decid_evrgn_mixed_shrub_pasture=prop_decid+prop_evrgn+prop_mixed+prop_shrub+prop_pasture, 
         prop_crop_pasture=prop_crop+prop_pasture, 
         prop_woody_herb=prop_woody+prop_herbs, 
         prop_crop_pasture_urban=prop_crop+prop_urban+prop_pasture) %>%
  dplyr::select(full_circle, year, prop_decid_evrgn_mixed, prop_decid_evrgn_mixed_shrub, prop_decid_evrgn_mixed_shrub_pasture, 
                prop_woody_herb, prop_grass, prop_crop_pasture, prop_crop_pasture_urban, prop_urban)

# join covars
site_year_covs <- covs <- left_join(tmin_per_site, precip_per_site) %>% 
  left_join(lulc_props) %>% rename(period=year) %>%
  arrange(full_circle, period)

# clean
rm(list=ls(pattern="prop_"))

# get norms for anomalies
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_tmin=mean(tmin, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_precip=mean(precip, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_forest=mean(prop_decid_evrgn_mixed, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>%  
  summarise(site_mean_forest_shrub=mean(prop_decid_evrgn_mixed_shrub, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_forest_shrub_pasture=mean(prop_decid_evrgn_mixed_shrub_pasture, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_wetland=mean(prop_woody_herb, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_grass=mean(prop_grass, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_ag=mean(prop_crop_pasture, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_anthro=mean(prop_crop_pasture_urban, na.rm=T)) %>%
  right_join(covs)
covs <- site_year_covs %>% group_by(full_circle) %>% 
  summarise(site_mean_urban=mean(prop_urban, na.rm=T)) %>%
  right_join(covs)

# create anomalies
names(covs)
site_year_covs <- covs %>% 
  mutate(site_urban_anom=prop_urban-site_mean_urban, 
         site_anthro_anom=prop_crop_pasture_urban-site_mean_anthro,
         site_ag_anom=prop_crop_pasture-site_mean_ag, 
         site_grass_anom=prop_grass-site_mean_grass, 
         site_wetland_anom=prop_woody_herb-site_mean_wetland, 
         site_forest_anom=prop_decid_evrgn_mixed-site_mean_forest, 
         site_forest_shrub_anom=prop_decid_evrgn_mixed_shrub-site_mean_forest_shrub, 
         site_forest_shrub_pasture_anom=prop_decid_evrgn_mixed_shrub_pasture-site_mean_forest_shrub_pasture, 
         site_tmin_anom=tmin-site_mean_tmin,
         site_precip_anom=precip-site_mean_precip) %>%
  dplyr::select(full_circle, period, 
                site_mean_tmin, site_tmin_anom,
                site_mean_precip, site_precip_anom,
                site_mean_urban, site_urban_anom,
                site_mean_anthro, site_anthro_anom,
                site_mean_ag, site_ag_anom,
                site_mean_grass, site_grass_anom,
                site_mean_wetland, site_wetland_anom,
                site_mean_forest_shrub, site_forest_shrub_anom,
                site_mean_forest_shrub_pasture, site_forest_shrub_pasture_anom,
                site_mean_forest, site_forest_anom)

# prep for modeling
names(site_year_covs)
site_year_covs <- site_year_covs %>% 
  mutate(decade=as.numeric(factor(period))-
           min(as.numeric(factor(period))) + 1,
         site_int=as.numeric(factor(full_circle)),
         tmin_scl=as.numeric(scale(site_mean_tmin)),
         tmin_anom_scl=as.numeric(scale(site_tmin_anom)),
         precip_scl=as.numeric(scale(site_mean_precip)),
         precip_anom_scl=as.numeric(scale(site_precip_anom)),
         urban_scl=as.numeric(scale(site_mean_urban)),
         urban_anom_scl=as.numeric(scale(site_urban_anom)),
         anthro_scl=as.numeric(scale(site_mean_anthro)),
         anthro_anom_scl=as.numeric(scale(site_anthro_anom)),
         ag_scl=as.numeric(scale(site_mean_ag)),
         ag_anom_scl=as.numeric(scale(site_ag_anom)),
         forest_shrub_scl=as.numeric(scale(site_mean_forest_shrub)),
         forest_shrub_anom_scl=as.numeric(scale(site_forest_shrub_anom)),
         forest_shrub_pasture_scl=as.numeric(scale(site_mean_forest_shrub_pasture)),
         forest_shrub_pasture_anom_scl=as.numeric(scale(site_forest_shrub_pasture_anom)),
         grass_scl=as.numeric(scale(site_mean_grass)),
         grass_anom_scl=as.numeric(scale(site_grass_anom)),
         wetland_scl=as.numeric(scale(site_mean_wetland)),
         wetland_anom_scl=as.numeric(scale(site_wetland_anom)),
         forest_scl=as.numeric(scale(site_mean_forest)),
         forest_anom_scl=as.numeric(scale(site_forest_anom)))
summary(site_year_covs)
# ------------------------------------------------------------------------------

# choose species ---------------------------------------------------------------
#grassland birds
species_list <- c("Savannah Sparrow", "Vesper Sparrow", 
                  "Eastern Meadowlark", "American Pipit") 

#shrubland birds
species_list <- c("Eastern Towhee", "Brown Thrasher", "Field Sparrow", "Fox Sparrow",
                  "Northern Bobwhite", "Northern Cardinal", "Northern Mockingbird", 
                  "American Tree Sparrow", "Loggerhead Shrike", "Harris's Sparrow")

#mixed habitat birds 
species_list <- c("Common Grackle", "Brown-headed Cowbird",
                  "Red-winged Blackbird", "American Crow",
                  "American Goldfinch", "White-crowned Sparrow",
                  "Killdeer", "Mourning Dove",
                  "Song Sparrow", "House Finch", "Orange-crowned Warbler")

#wetland birds: waterfowl
species_list <- c("American Black Duck", "Wood Duck", "Hooded Merganser", 
                  "Canvasback", "Ring-necked Duck", "Snow Goose",
                  "Cackling/Canada Goose", "Lesser Scaup", "Redhead", "Ruddy Duck",
                  "American Wigeon", "Greater Scaup", "Bufflehead")

#wetland birds: shorebirds gulls and waders
species_list <- c("American Woodcock", "Ring-billed Gull", "Bonaparte's Gull", "Common/Wilson's Snipe", "American Coot",
                  "Sandhill Crane", "Double-crested Cormorant", "Fish Crow",
                  "Great Blue Heron", "Bald Eagle", "Great Black-backed Gull", "Herring Gull")

#wetland birds: passerines
species_list <- c("Belted Kingfisher", "Rusty Blackbird","Swamp Sparrow","Marsh Wren",
                  "Nelson's/Saltmarsh Sparrow (Sharp-tailed Sparrow)","LeConte's Sparrow")

# large forest birds
species_list <- c("Red-shouldered Hawk","Wild Turkey", "Red-tailed Hawk","Cooper's Hawk")

# forest birds: woodpeckers
species_list <- c("Red-cockaded Woodpecker", "Red-headed Woodpecker",
                  "Red-bellied Woodpecker","Yellow-bellied Sapsucker", 
                  "Pileated Woodpecker","Northern Flicker",
                  "Downy Woodpecker", "Hairy Woodpecker")

# forest birds: passerines
species_list <- c("Tufted Titmouse","Pine Warbler","Carolina Chickadee",
                  "Carolina Wren", "Winter Wren","Eastern Phoebe","Eastern Bluebird","Blue Jay", "Hermit Thrush",
                  "Yellow-rumped Warbler","Ruby-crowned Kinglet","Golden-crowned Kinglet","Brown Creeper",
                  "Brown-headed Nuthatch", "Purple Finch","White-throated Sparrow",
                  "Dark-eyed Junco","American Robin", "White-breasted Nuthatch","Pine Siskin","Red-breasted Nuthatch")

# ------------------------------------------------------------------------------

# alter query ---------------------------------------------------------------------------

spnames <- read.csv("Z:/GitHub/QuantitativeMetrics/cbc_analysis_library/code_2019/cluster/data/cbc_trends_spp_list_66_to_119.csv")

#head(spnames)
querylist <- spnames %>% filter(ebird_com_name %in% species_list) %>% dplyr::select(query, species_code)

# get count data and merge with cov data -------------------------------------
count_data <- c()
for(i in 1:nrow(querylist)){
  # define query
  species <- querylist[i,2]
  focal_sp <- querylist[i,1] #paste0("ref_species.com_name = '", species, "'")
  first_count_num <- 30
  last_count_num <- 119
  query_text = paste("SELECT loc_circle.subnational_code, loc_circle.name,",
                     "cnt_submission.count_yr, ref_species.sci_name,",
                     "cnt_observation.how_many FROM cnt_submission",  
                     "FULL JOIN loc_circle ON loc_circle.circle_id =",
                     "cnt_submission.circle_id FULL JOIN cnt_observation",
                     "ON cnt_submission.submission_id = cnt_observation.submission_id",
                     "FULL JOIN ref_species ON cnt_observation.species_id = ",
                     "ref_species.species_id WHERE", focal_sp, "AND",
                     "cnt_submission.count_yr BETWEEN", first_count_num, 
                     "AND", last_count_num)
  # establish a connection
  channel <- odbcConnect(dsn = "CBC database", uid = "xx", 
                         pwd = "xx")
  # query database
  cbc_dat <- sqlQuery(channel, query_text)
  # close the connection
  odbcClose(channel)
  # clean and merge count data
  d1 <- cbc_dat %>% dplyr::select(circle_name=name, subnational_code, 
                                  species=sci_name,
                                  count_year=count_yr, count=how_many) %>%
    mutate(full_circle=paste(circle_name, subnational_code)) %>%
    dplyr::select(-c(circle_name, subnational_code)) %>%
    as.data.frame()
  d2 <- merge(effort_data, d1, by=c("full_circle", "count_year"), all.x=T, all.y=F)
  d2$species <- species
  d2$count[is.na(d2$count)] <- 0
  count_dat <- d2 %>%
    dplyr::select(full_circle, lon, lat, count_year, period, distance_km, species,
                  count) %>%
    group_by(full_circle, lon, lat, period, species) %>%
    summarise(count=sum(count), distance_km=sum(distance_km)) %>%
    mutate(present=if_else(count>0, 1, 0)) %>%
    arrange(desc(lat), period)
  
  # get proportion 0s
  detection_info <- count_dat %>% group_by(full_circle) %>% 
    summarise(site_number_present=sum(present),
              site_possible_present=n()) %>% 
    mutate(fraction_present=site_number_present/site_possible_present)
  
  # merge detection info
  count_dat <- count_dat %>% left_join(detection_info)
  
  # stack data
  count_data <- rbind(count_data, count_dat)
  print(paste("completed", species))
}

# merge counts and covs
model_data <- as.data.frame(arrange(count_data, full_circle, period)) %>%
  left_join(site_year_covs) %>%
  mutate(period=period-4, 
         log_distance_1k_km=log(distance_km/1000),
         species_int=as.integer(factor(species)),
         intercept=1)
#-------------------------------------------------------------------------------

# filter data to presence absence -----------------------------------------------------------
filtered_data <- model_data %>% filter(fraction_present >= 0 & 
                                         fraction_present <= 1) %>% 
  as.data.frame()
# ------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# run full model, model fit stats, then model subsets to determine relative importance of variables
# mgcv spatial gamm ------------------------------------------------------------
gam1 <- gam(present ~ s(log_distance_1k_km) +
              s(tmin_scl, tmin_anom_scl) +
              s(precip_scl, precip_anom_scl) +
              s(forest_scl, forest_anom_scl) + 
              s(urban_scl, urban_anom_scl) +
              s(site_int, bs="re") + 
              s(species_int, bs="re") + 
              s(decade, bs="re") +
              s(lon, lat, k=6),
            data=filtered_data, 
            family="binomial")

summary(gam1)

#check model fit
par(mfrow = c(2,2))
gam.check(gam1)

#check plot of effects
b <- getViz(gam1)

gridPrint(plot(sm(b, 2)) + l_fitRaster() + l_fitContour(size=0.6, colour="lightgray") + l_rug() + labs(title=NULL, y="Min temp anomaly",x="Min temp avg") + scale_fill_viridis(name="s(x)"),
          plot(sm(b, 3)) + l_fitRaster() + l_fitContour(size=0.6, colour="lightgray") + l_rug() + labs(title=NULL, y="Precip anomaly",x="Precip avg") + scale_fill_viridis(name="s(x)"),
          plot(sm(b, 4)) + l_fitRaster() + l_fitContour(size=0.6, colour="lightgray") + l_rug() + labs(title=NULL, y="Forest anomaly",x="Forest avg") + scale_fill_viridis(name="s(x)"),
          plot(sm(b, 5)) + l_fitRaster() + l_fitContour(size=0.6, colour="lightgray") + l_rug() + labs(title=NULL, y="Urban anomaly",x="Urban avg") + scale_fill_viridis(name="s(x)") ,ncol = 2)

# Two dimensional residual checks
b <- getViz(b, nsim = 50)
ck1 <- check2D(b, x1 = "lon", x2 = "lat")
res1 <- ck1 + l_gridCheck2D(gridFun = mean)

ck1 <- check2D(b, x1 = "tmin_scl", x2 = "tmin_anom_scl")
res2 <- ck1 + l_gridCheck2D(gridFun = mean)

ck1 <- check2D(b, x1 = "precip_scl", x2 = "precip_anom_scl")
res3 <- ck1 + l_gridCheck2D(gridFun = mean)

ck1 <- check2D(b, x1 = "forest_scl", x2 = "forest_anom_scl")
res4 <- ck1 + l_gridCheck2D(gridFun = mean)

ck1 <- check2D(b, x1 = "urban_scl", x2 = "urban_anom_scl")
res5 <- ck1 + l_gridCheck2D(gridFun = mean)

#-------------------------------------------------------------------------------
## Determining relative importance of LULC vs. climate/spatial vs. 
# temporal/spatial components/temporal components/interactions
full.glm <- gam(present ~ s(log_distance_1k_km) +
                  s(tmin_scl, tmin_anom_scl) +
                  s(precip_scl, precip_anom_scl) +
                  s(forest_scl, forest_anom_scl) + 
                  s(urban_scl, urban_anom_scl) +
                  s(site_int, bs="re") + 
                  s(species_int, bs="re") + 
                  s(decade, bs="re") +
                  s(lon, lat, k=6),
                data=filtered_data, 
                family="binomial")

clim.glm <- gam(present ~ s(log_distance_1k_km) +
                  s(tmin_scl, tmin_anom_scl) +
                  s(precip_scl, precip_anom_scl) +
                  s(site_int, bs="re") + 
                  s(species_int, bs="re") + 
                  s(decade, bs="re") +
                  s(lon, lat, k=6),
                data=filtered_data, 
                family="binomial")

lulc.glm <- gam(present ~ s(log_distance_1k_km) +
                  s(forest_scl, forest_anom_scl) +
                  s(urban_scl, urban_anom_scl) +
                  s(site_int, bs="re") + 
                  s(species_int, bs="re") + 
                  s(decade, bs="re") +
                  s(lon, lat, k=6),
                data=filtered_data, 
                family="binomial")

spatialfe.glm <- gam(present ~ s(log_distance_1k_km) +
                       s(tmin_scl) +
                       s(precip_scl) +
                       s(forest_scl) + 
                       s(urban_scl) +
                       s(site_int, bs="re") + 
                       s(species_int, bs="re") + 
                       s(decade, bs="re") +
                       s(lon, lat, k=6),
                     data=filtered_data, 
                     family="binomial")

spatialclim.glm <- gam(present ~ s(log_distance_1k_km) + 
                         s(tmin_scl) +
                         s(precip_scl) +
                         s(site_int, bs="re") + 
                         s(species_int, bs="re") + 
                         s(decade, bs="re") +
                         s(lon, lat, k=6),
                       data=filtered_data, 
                       family="binomial")

spatiallulc.glm <- gam(present ~ s(log_distance_1k_km) + 
                         s(forest_scl) + 
                         s(urban_scl) +
                         s(site_int, bs="re") + 
                         s(species_int, bs="re") + 
                         s(decade, bs="re") +
                         s(lon, lat, k=6),
                       data=filtered_data, 
                       family="binomial")

temporal.glm <- gam(present ~ s(log_distance_1k_km) +
                      s(tmin_anom_scl) +
                      s(precip_anom_scl) +
                      s(forest_anom_scl) + 
                      s(urban_anom_scl) +
                      s(site_int, bs="re") + 
                      s(species_int, bs="re") + 
                      s(decade, bs="re") +
                      s(lon, lat, k=6),
                    data=filtered_data, 
                    family="binomial")

temporalclim.glm <- gam(present ~ s(log_distance_1k_km) + 
                          s(tmin_anom_scl) +
                          s(precip_anom_scl) +
                          s(site_int, bs="re") + 
                          s(species_int, bs="re") + 
                          s(decade, bs="re") +
                          s(lon, lat, k=6),
                        data=filtered_data, 
                        family="binomial")

temporallulc.glm <- gam(present ~ s(log_distance_1k_km) + 
                          s(forest_anom_scl) +
                          s(urban_anom_scl) +
                          s(site_int, bs="re") + 
                          s(species_int, bs="re") + 
                          s(decade, bs="re") +
                          s(lon, lat, k=6),
                        data=filtered_data, 
                        family="binomial")

spteadd.glm <- gam(present ~ s(log_distance_1k_km) +
                     s(tmin_scl) +
                     s(tmin_anom_scl) +
                     s(precip_scl) +
                     s(precip_anom_scl) +
                     s(forest_scl) +
                     s(forest_anom_scl) +
                     s(urban_scl) +
                     s(urban_anom_scl) +
                     s(site_int, bs="re") + 
                     s(species_int, bs="re") + 
                     s(decade, bs="re") +
                     s(lon, lat, k=6),
                   data=filtered_data, 
                   family="binomial")

climinter.glm <- gam(present ~ s(log_distance_1k_km) +
                       s(tmin_scl,tmin_anom_scl) +
                       s(precip_scl,precip_anom_scl) +
                       s(forest_scl) +
                       s(forest_anom_scl) +
                       s(urban_scl) +
                       s(urban_anom_scl) +
                       s(site_int, bs="re") + 
                       s(species_int, bs="re") + 
                       s(decade, bs="re") +
                       s(lon, lat, k=6),
                     data=filtered_data, 
                     family="binomial")

lulcinter.glm <- gam(present ~ s(log_distance_1k_km) +
                       s(tmin_scl) +
                       s(tmin_anom_scl) +
                       s(precip_scl) +
                       s(precip_anom_scl) +
                       s(forest_scl,forest_anom_scl) +
                       s(urban_scl,urban_anom_scl) +
                       s(site_int, bs="re") + 
                       s(species_int, bs="re") + 
                       s(decade, bs="re") +
                       s(lon, lat, k=6),
                     data=filtered_data, 
                     family="binomial")

null.glm <- gam(present ~ s(log_distance_1k_km) +
                  s(site_int, bs="re") + 
                  s(species_int, bs="re") + 
                  s(decade, bs="re") +
                  s(lon, lat, k=6),
                data=filtered_data, 
                family="binomial")

#-----------------------------------------------------------------------------------------------------------------------
# version of null model with just decade and lon/lat splines to use
# for visualizing presence per decade empirically

null3.glm <- gam(present ~ s(log_distance_1k_km) +
                   s(site_int, bs="re") + 
                   s(species_int, bs="re") + 
                   te(lon, lat, decade, k=c(30,5),d=c(2,1), bs=c("tp","cr")), 
                 data=filtered_data, 
                 family="binomial")

# create new dataset
# clip grid to map
namap1 <- as(focal_area, "Spatial")
grid1 <- as(namap1, "sf") %>%
  st_make_grid(cellsize = 0.5, what = "centers") %>% as("Spatial")

predgrid <- coordinates(grid1)
newdata3 <- data.frame(lon=rep(predgrid[,1],3),lat=rep(predgrid[,2],3),
                       decade=rep(c(1,2,3,4,5,6,7,8,9),each=nrow(predgrid)),
                       log_distance_1k_km=0,
                       site_int=0,
                       species_int=0)

newdata3$Pred <- plogis(predict.gam(null3.glm,newdata3))

ggplot(newdata3) +
  geom_tile(aes(lon,lat, fill=Pred)) +
  facet_wrap(~decade) + scale_fill_distiller(palette = "Spectral") + 
  geom_sf(data=focal_area, fill=NA)

#fitted values per decade
fit.val <- c(null3.glm$fitted.values)
dec.data <- c(filtered_data$decade)
df <- as.data.frame(cbind(fit.val,dec.data))

# Save deviance for each model subset ------------------------------------------------------------------------
fdev <- full.glm$deviance
spdev <- spatialfe.glm$deviance
tedev <- temporal.glm$deviance
spclidev <- spatialclim.glm$deviance
splulcdev <- spatiallulc.glm$deviance
teclidev <- temporalclim.glm$deviance
telulcdev <- temporallulc.glm$deviance
nulldev <- null.glm$deviance
clidev <- clim.glm$deviance
lulcdev <- lulc.glm$deviance

mod_devs <- c(fdev,spdev,tedev,spclidev,splulcdev,teclidev,telulcdev,nulldev,clidev,lulcdev)
names(mod_devs) <- c('Full Dev', 'Sp Dev', "Tem Dev", 'Sp Cli Dev', 'Sp LULC Dev', 'Tem Cli Dev', 'Tem LULC Dev', 'Null Dev', 'Cli Dev', 'LULC Dev')
mdevs <- sort(mod_devs)

#ANODEV for each model subset --------------------------------------------------------------------------------------------------------
#deviance explained by climate
clim.devex <- (null.glm$deviance-clim.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance explained by lulc
lulc.devex <- (null.glm$deviance-lulc.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance explained by spatial component of climate & LULC 
spatial.devex <- (null.glm$deviance-spatialfe.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance expl by spatial component of climate 
spatialclim.devex <- (null.glm$deviance-spatialclim.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance expl by spatial component of lulc 
spatiallulc.devex <- (null.glm$deviance-spatiallulc.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance explained by temporal component of climate & lulc
temporal.devex <- (null.glm$deviance-temporal.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance explained by temporal component of climate 
temporalclim.devex <- (null.glm$deviance-temporalclim.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance explained by temporal component of lulc 
temporallulc.devex <- (null.glm$deviance-temporallulc.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#devaince explained by spatial + temporal additive (no interactions)
sptemadd.devex <- (null.glm$deviance-spteadd.glm$deviance)/(null.glm$deviance-full.glm$deviance)
#deviance explained by climate interaction
climinter.devex <- (null.glm$deviance-climinter.glm$deviance)/(null.glm$deviance-full.glm$deviance) - sptemadd.devex
#deviance explained by LULC interaction 
lulcinter.devex <- (null.glm$deviance-lulcinter.glm$deviance)/(null.glm$deviance-full.glm$deviance) - sptemadd.devex

devs <- c(clim.devex, lulc.devex, spatial.devex, spatialclim.devex, spatiallulc.devex,
          temporal.devex, temporalclim.devex, temporallulc.devex, sptemadd.devex, climinter.devex, lulcinter.devex)

names(devs) <- c("Climate Dev", "LULC Dev", "Spatial Dev", "Spatial Clim Dev", "Spatial LULC Dev",
                 "Temporal Dev", "Temporal Clim Dev", "Temporal LULC Dev", "Additive Dev","Clim Inter Dev", "LULC Inter Dev")

var.impt <- sort(devs) 

# Log likelihoods -----------------------------------------------------------------------------
# save log likelihood comparisons with full model
t1 <- lrtest(full.glm,spatialfe.glm)
sp.test <- t1$LogLik[2]-t1$LogLik[1]
t2 <- lrtest(full.glm,temporal.glm)
temp.test <- t2$LogLik[2]-t2$LogLik[1]
t3 <- lrtest(full.glm,clim.glm)
clim.test <- t3$LogLik[2]-t3$LogLik[1]
t4 <- lrtest(full.glm,lulc.glm)
lulc.test <- t4$LogLik[2]-t4$LogLik[1]
t5 <- lrtest(full.glm,spatialclim.glm)
spclim.test <- t5$LogLik[2]-t5$LogLik[1]
t6 <- lrtest(full.glm,spatiallulc.glm)
splulc.test <- t6$LogLik[2]-t6$LogLik[1]
t7 <- lrtest(full.glm,temporalclim.glm)
tempclim.test <- t7$LogLik[2]-t7$LogLik[1]
t8 <- lrtest(full.glm,temporallulc.glm)
templulc.test <- t8$LogLik[2]-t8$LogLik[1]
t9 <- lrtest(full.glm,null.glm)
null.test <- t9$LogLik[2]-t9$LogLik[1]

lrtests <- c(sp.test,temp.test,clim.test,lulc.test,spclim.test,splulc.test,tempclim.test,templulc.test,null.test)
names(lrtests) <- c("Sp test","Temp test","Climate test","LULC test","Sp Clim test","Sp LULC test","Temp Clim test","Temp LULC test","Null test")

lratios <- sort(lrtests)

#Likelihood ratio tests  -------------------------------------------------------------------------
library(lmtest)
options(scipen = -2)
# save LRTs comparing each submodel to the null model 
lr1 <- lrtest(temporallulc.glm,null.glm)
tlulc_null <- c(lr1$Chisq[2],lr1$Df[2],lr1$`Pr(>Chisq)`[2])
lr2 <- lrtest(temporalclim.glm,null.glm)
tcli_null <- c(lr2$Chisq[2],lr2$Df[2],lr2$`Pr(>Chisq)`[2])
lr3 <- lrtest(temporal.glm,null.glm)
temp_null <- c(lr3$Chisq[2],lr3$Df[2],lr3$`Pr(>Chisq)`[2])
lr4 <- lrtest(spatialclim.glm,null.glm)
sclim_null <- c(lr4$Chisq[2],lr4$Df[2],lr4$`Pr(>Chisq)`[2])
lr5 <- lrtest(spatiallulc.glm, null.glm)
slulc_null <- c(lr5$Chisq[2],lr5$Df[2],lr5$`Pr(>Chisq)`[2])
lr6 <- lrtest(spatialfe.glm,null.glm)
spatial_null <- c(lr6$Chisq[2],lr6$Df[2],lr6$`Pr(>Chisq)`[2])
lr7 <- lrtest(clim.glm,null.glm)
clim_null <- c(lr7$Chisq[2],lr7$Df[2],lr7$`Pr(>Chisq)`[2])
lr8 <- lrtest(lulc.glm,null.glm)
lulc_null <- c(lr8$Chisq[2],lr8$Df[2],lr8$`Pr(>Chisq)`[2])
lr9 <- lrtest(full.glm,null.glm)
full_null <- c(lr9$Chisq[2],lr9$Df[2],lr9$`Pr(>Chisq)`[2])

lrts_null <- matrix(c(tlulc_null,tcli_null,temp_null,sclim_null,slulc_null,spatial_null,clim_null,lulc_null,full_null),nrow=9,ncol=3,byrow = TRUE)
rownames(lrts_null) <- c('TempLULC to Null','TempCli to Null','Temp to Null','SpCli to Null',
                         'SpLULC to Null','Spat to Null','Clim to Null','LULC to Null','Full to Null')
colnames(lrts_null) <- c('ChiSq','Delta DF','P value')
