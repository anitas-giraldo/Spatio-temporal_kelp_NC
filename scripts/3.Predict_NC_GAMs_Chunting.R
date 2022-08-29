##

## Script by Anita Giraldo, 4 May 2022
## Last modified by Anita Giraldo, 4 May 2022

##  This script predicts and projects spatially the result from the models selected in script 2

# libraries ----
library(Amelia) # to map missing data
library(here)
library(tidyverse)
library(dplyr)
library(ggplot2) 
library(caret) ## For model fitting and evaluation
library(RANN)
library(corrplot)
library(rsample)
library(yardstick) ## to evaluate a model
library(visreg)  ## For visualizing regression models
library(plotROC) ## For constructing ROC curves
library(mgcv)    ## For fitting GAM models
library(kernlab) ## Contains an example dataset
library(glmnet)  ## For fitting regularized models
library(pROC)
library(ROCR)
library(car)
library(OddsPlotty) # to plot from the caret package
library(tidymodels) # modern version of caret
library(cutpointr) # Find optimal cutoff point for binary classification
library(bestglm) # Best subset GLM : The function bestglm selects the best subset of inputs for the glm family
library(MASS) # Stepwise regression
library(reshape2)
library(sp)
library(sf)
library(rgdal)
library(raster) 
library(terra)
library(RColorBrewer)
#library(geobr)
library(ggspatial)
library(ggrepel)
library(forcats)
library(beepr)
library(gridExtra)
library(FSSgam)
library(viridis)
library(hrbrthemes)
library(gplots)
library(beepr)
library(stringr)



# Clear environment ----
rm(list=ls())



# directories ----
m.dir <- here()
d.dir <- here('data')
k.dir <- here('outputs_nc_rcca')
o.dir <- paste(k.dir, "gam_V5", sep ='/') # kelp model results
u.dir <- paste(k.dir, "gam_urchins3", sep ='/') # urchin model results
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



## Load info on years RCCA ----
years <- read.csv(paste(d.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


# get the sites from with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 64


# 1. Load RCCA data ----

df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  mutate(zone_new = case_when(
    transect == '1' ~ 'OUTER',
    transect == '2' ~ 'OUTER', 
    transect == '3' ~ 'OUTER', 
    transect == '4' ~ 'INNER',
    transect == '5' ~ 'INNER',
    transect == '6' ~ 'INNER')) %>%
  dplyr::select(-zone) %>%
  rename(zone = zone_new) %>%
  mutate_at(vars(zone), list(as.factor)) %>%
  relocate(zone, .after = transect) %>%
  glimpse() # Rows: 1,154



## get the sites for North Coast model ----
df.nc <- df %>%
  dplyr::select(-c(latitude, longitude)) %>%
  right_join(ncsites, by = c('site_name')) %>%
  droplevels() %>% #glimpse()
  #dplyr::select(-c(total.years, pre.mhw.years, during.mhw.years, post.mhw.years)) %>%
  relocate(c(latitude, longitude), .after = zone) %>%
  glimpse() # Rows: 708

length(levels(df.nc$site_name)) # 10
levels(df.nc$site_name)
any(is.na(df.nc$Max_Monthly_Anomaly_Temp))





# 2. Choose variables and transform needed ----

names(df.nc)

dat1 <- df.nc %>%
  dplyr::select(
    # Factors 
    latitude, longitude,
    site_name, year, transect, zone,
    # Bio vars
    den_NERLUE , den_MESFRAAD , den_STRPURAD , den_PYCHEL, den_HALRUF,
    # Nitrate vars 
    Days_10N, 
    Min_Monthly_Nitrate, 
    Max_Monthly_Nitrate,
    Mean_Monthly_Nitrate,
    Mean_Monthly_Upwelling_Nitrate,
    Max_Monthly_Anomaly_Nitrate,
    Mean_Monthly_Summer_Nitrate,
    # Temperature vars
    Days_16C ,
    Mean_Monthly_Temp ,
    Mean_Monthly_Summer_Temp,
    MHW_Upwelling_Days  , 
    Min_Monthly_Anomaly_Temp,
    Max_Monthly_Anomaly_Upwelling_Temp,
    Min_Monthly_Temp, 
    Mean_Monthly_Upwelling_Temp,
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc,
    # Orb vel
    UBR_Mean, UBR_Max,
    # NPP
    Mean_Monthly_NPP, Max_Monthly_NPP_Upwelling, Mean_Monthly_NPP_Upwelling, Min_Monthly_NPP,
  ) %>%
  # Bio transformations
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_den_PYCHEL = log(den_PYCHEL + 1),
         log_den_HALRUF = log(den_HALRUF + 1),
         log_mean_vrm = log(mean_vrm + 1)) %>%
  dplyr::select(-c(den_NERLUE,
                   den_MESFRAAD,
                   den_STRPURAD,
                   den_PYCHEL,
                   den_HALRUF,
                   mean_vrm)) %>%
  # Temperature transformations
  mutate(log_Days_16C = log(Days_16C + 1)) %>%
  dplyr::select(-c(Days_16C)) %>%
  # Orb vel transformations
  mutate(log_UBR_Mean = log(UBR_Mean + 1),
         log_UBR_Max = log(UBR_Max +1)) %>%
  dplyr::select(-c(UBR_Mean,
                   UBR_Max)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)



# 3. Divide data into train and test ----

inTraining <- createDataPartition(dat2$log_den_NERLUE, p = 0.75, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]




# 4. Run GAM  ----

gam1 <- gam(formula = log_den_NERLUE ~ 
                s(log_den_STRPURAD, k = 5, bs = "cr") +
                s(Max_Monthly_Nitrate, k = 5, bs = "cr") + 
                s(wh_max, k = 5, bs = "cr") +
                s(log_UBR_Max, k = 4, bs = "cr") +
                s(site_name, zone, bs = "re") + 
                s(year, bs = "re"), 
              family = tw(), data = dat2, method = "REML")


# 5. Check GAM ----
gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)

# visualize responses
par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()




# 6. Predict to compare to observed ----

testdata <- dat2 %>%
  dplyr::select("log_den_STRPURAD", 
                "log_mean_vrm",
                "log_UBR_Max",
                "Max_Monthly_Nitrate",
                "wh_max",
                "log_den_NERLUE",
                "site_name", "zone", "year")

head(testdata)


fits <- predict.gam(gam1, newdata=testdata, type='response', se.fit=T)


## predict average kelp per year --
predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit, na.rm = T),se.fit=mean(se.fit, na.rm = T))%>%
  ungroup()


ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  scale_fill_viridis(discrete = T) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

ggmod.year


###


# 7. Plot observed vs. predicted ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()



# Plot observed vs. predicted --
library(ggpmisc)

my.formula <- y ~ x

p <- ggplot(predicts.all, aes(x = fit, y = log_den_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  #scale_color_viridis(discrete = T) +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p






###


###


###


# *** PREDICT BEST MODEL ACROSS ALL YEARS AND SITE THAT I HAVE DATA FOR ----

#gam1
#Max_Monthly_Nitrate
#log_UBR_Max
#wh_max
#log_den_STRPURAD

### Get depth ----

depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

names(dat2)

depth <- rast(paste(depth.dir, "depth_mean_nc_all_wInterp_300m_30m.tif", sep ='/'))
plot(depth)

n.extent <- ext(depth)

crs1 <- "epsg:4326"
d2 <- project(depth, crs1)

n.extent <- ext(d2)


## Get rock ----
sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

dir(sub.dir)

rock <- rast(paste(sub.dir, "prob_rock_nc_MASKED_LatLong_all_300m_wInterp.tif", sep ='/'))
rock 

# crop to NC --
rock2 <- crop(rock, ext(d2))
plot(rock2)

rock3 <- resample(rock2, d2)
plot(rock3)


### Get Env predictors  ----

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data"



### Get nitrate predictors  ----
max_nit <- rast(paste(re.dir, "Nitrate", "Max_Monthly_Nitrate.tif", sep ='/'))
max_nit

# # crop to NC --
max_nit2 <- crop(max_nit, n.extent)
plot(max_nit2[[1]])

# resample predictors to bathy ----
max_nit3 <- resample(max_nit2, d2)

# mask predictors to bathy ----
max_nit4 <- mask(max_nit3, d2)
plot(max_nit4[[1]])


### Get Wave predictors  ----

w.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Wave_data"

## Max Wave height --

# load raster data --

wave.dir <- paste(w.dir, "Wave_extended2", "wh_max", sep ='/')

# load raster data --

n.files <- dir(wave.dir)
# list files in source --
n.files <- list.files(wave.dir, pattern = '.tif$', full.names = TRUE)
n.files
length(n.files)
# list names to load onto the Environment --
names.list <- list.files(wave.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
whmax.stack <- c()

# use do call to create a raster otherwise it creates a list
whmax.stack <- do.call("c", tfiles)
plot(whmax.stack[[1]])


# # crop to NC --
whmax.stack2 <- crop(whmax.stack, n.extent)
plot(whmax.stack2[[1]])

# resample predictors to bathy ----
whmax.stack3 <- resample(whmax.stack2, d2)

# mask predictors to bathy ----
whmax.stack4 <- mask(whmax.stack3, d2)
plot(whmax.stack4[[1]])


## Mean UBR MAX ----

w2.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Waves_UBR_30m"

# load raster data --

ubr <- rast(paste(w2.dir, "UBR_Max_30m_NC.tif", sep ='/'))
plot(ubr[[1]])

ubr <- project(ubr, rock)

plot(ubr[[1]])

ubr1 <- classify(ubr, cbind(0, NA))
plot(ubr1[[1]])


# # crop to NC --
ubr2 <- crop(ubr1, n.extent)
plot(ubr2[[1]])

# resample predictors to bathy ----
ubr3 <- resample(ubr2, d2)

# mask predictors to bathy ----
ubr4 <- mask(ubr3, d2)
plot(ubr4[[1]])

ubr5 <- log(ubr4)

### Get urchins ----

# load purple urchin predictions ----
urch.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Spatio_temporal_GAMs/outputs_nc_rcca/gam_urchins3/log_sp_predictions"

# load raster data --

u.files <- dir(urch.dir)
u.files <- list.files(urch.dir, pattern = '.tif')
u.files
length(u.files)


##



# stack rasters --

preds1 <- c(max_nit4[[1]], whmax.stack4[[1]], ubr5[[1]])
names(preds1)


# get year raster ----

year1998 <- classify(max_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("Max_Monthly_Nitrate" , 
                   "wh_max", 
                   "log_UBR_Max",                         
                   "year") 


# sites ----

rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)

site.raster <- rast(rdf, type = 'xyz', crs="EPSG:4326", extent = ext(year1998))
site.raster

plot(site.raster)

ext(year1998)
ext(site.raster)

site.raster2 <- extend(site.raster, year1998)

preds3 <- c(preds2, site.raster2)
names(preds3) <- c("Max_Monthly_Nitrate"  , 
                   "wh_max", 
                   "log_UBR_Max",                         
                   "year" ,   
                   "site_name") 


##

# zone ----

zone.raster <- d2
names(zone.raster) <- 'zone'
plot(zone.raster)
levels(dat2$zone)

rec.m <-  c(-Inf, -10, 2,
            -10, 0.1, 1)

rclmat <- matrix(rec.m, ncol=3, byrow=TRUE)

zone.raster2 <- classify(zone.raster, rclmat, right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("Max_Monthly_Nitrate" , 
                   "wh_max", 
                   "log_UBR_Max",                         
                   "year",    
                   "site_name",
                   "zone")




## LOOP to predict each year using data frame ----

nereo.mod <- gam1
summary(nereo.mod)

# make list of years --
year.list <- paste(2004:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(d2, cbind(-Inf, 0.1, 2004), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 



# outputs dir ----
# * use an output directory of yours
preds.dir <- paste(o2.dir, "preds", sep ='/')
preds.dir

# output for rasters scaled by rock 
# * use an output directory of yours
rock.preds.dir <- paste(o2.dir, "rock_preds", sep ='/')
rock.preds.dir




for (i in 1:length(year.list)) {
  
  # 1. get urchins
  urchin.rast <- rast(paste(urch.dir, u.files[i], sep ='/'))
  urchin.rast2 <- resample(urchin.rast, d2)
  
  # 2. stack with predictors for that year
  #env.raster <- c(d2, max_nit4[[i+6]], whmax.stack4[[i]], wymean.stack4[[i]])
  
  #env.raster <- c(d2, mean_up_T4[[i+6]], max_nit4[[i+6]], whmax.stack4[[i]], wymean.stack4[[i]])
  
  # V3
  env.raster <- c(max_nit4[[i+6]], whmax.stack4[[i]], ubr5[[i]])
  
  preds1 <- c(urchin.rast2, env.raster)
  
  # 3. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(-Inf, 0, year.no), right=FALSE)
  
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster2)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("log_den_STRPURAD",
                     "Max_Monthly_Nitrate"  , 
                     "wh_max", 
                     "log_UBR_Max",                         
                     "year", 
                     "zone",
                     "site_name")
  
  df4 <- as.data.frame(preds4, xy = T) %>%
    mutate_at(vars(year, zone, site_name), list(as.factor)) %>%
    mutate(zone = recode_factor(zone, '1' = 'INNER', '2' = 'OUTER')) %>%
    glimpse()
  
  # 5. predict
  year.pred.df <- predict.gam(nereo.mod, newdata=df4, type='response', se.fit=T)
  head(year.pred.df)
  
  # join with df for lats and lons
  preds.all <-  df4 %>% 
    data.frame(year.pred.df) %>%
    dplyr::select(x, y, fit) %>%
    glimpse()
  
  # 6. Rasterize
  crs.p <- "epsg:4326"
  year.prediction <- rast(preds.all, type = 'xyz', crs = crs.p, digits = 6)
  plot(year.prediction)
  
  # 7. save raw raster
  name.raster <- paste(year.no, "Log_Nereo_GIVE_IT_A_NAME.tif", sep = '_')
  writeRaster(year.prediction, paste(preds.dir, name.raster, sep = '/'))
  
  # 8. scale by rock
  rock4 <- resample(rock3, year.prediction)
  year.prediction2 <- rock4*year.prediction
  
  # 9. save rastr scaled by rock
  name.raster.rock <- paste(year.no, "Log_Nereo_rock__GIVE_IT_A_NAME.tif", sep = '_')
  writeRaster(year.prediction2, paste(rock.preds.dir, name.raster.rock, sep = '/'))
}




###

###

