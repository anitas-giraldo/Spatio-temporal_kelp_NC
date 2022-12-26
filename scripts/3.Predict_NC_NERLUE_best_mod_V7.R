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
o.dir <- paste(k.dir, "gam_V7", sep ='/') # kelp model results
u.dir <- paste('outputs_nc_rcca_urchins', "gam_urchins4_120m", sep ='/') # urchin model results
# rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
# dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



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


# read best modesl ----
bm <- read.csv(paste(o.dir, "best_models.csv", sep ='/'))




# 4. Run GAM  ----

gam1 <- gam(formula = log_den_NERLUE ~ 
                s(log_den_STRPURAD, k = 6, bs = "cr") +
                #s(log_Days_16C, k = 3, bs = "cr") + 
                s(Mean_Monthly_Upwelling_Temp, k = 3, bs = "cr") + 
                s(Max_Monthly_Nitrate, k = 6, bs = "cr") +
                s(log_UBR_Max, k = 10, bs = "cr") +
                s(wh_mean, k = 10, bs = "cr") +
                s(mean_depth, k = 4, bs = "cr") +
                s(Mean_Monthly_NPP, k = 6, bs = "cr") +
                s(site_name, zone, bs = "re") + 
                s(year, bs = "re"), 
              family = tw(), data = dat2, method = "REML")


# 5. Check GAM ----
gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)

# visualize rensponse 
par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()



# save plot --
pdf(file=paste(o.dir, paste(name,'logden_NERLUE',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,3),mar=c(2,4,3,1))
#plot(gam5.1,all.terms=T,pages=1,residuals=T,pch=16)
visreg(gam5.1)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


# 10. Predict to compare to observed ----

testdata <- dat2 %>%
  dplyr::select("log_den_STRPURAD", 
                #"log_Days_16C",
                'Mean_Monthly_Upwelling_Temp',
                "log_UBR_Max",
                "Max_Monthly_Nitrate",    
                "wh_mean",
                'mean_depth',
                'Mean_Monthly_NPP',
                log_den_NERLUE,
                site_name, zone, year)

head(testdata)


fits <- predict.gam(gam1, newdata=testdata, type='response', se.fit=T)


### predict average kelp per year ----
predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit, na.rm = T),se.fit=mean(se.fit, na.rm = T))%>%
  ungroup()


ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_viridis(discrete = T) +
  #scale_fill_manual(values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year


###


### Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

#write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = log_den_NERLUE, y = fit, color = zone)) +
  geom_point()


# Plot observed vs. predicted ----
library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = log_den_NERLUE, col = zone)) +        
  geom_point() +
  scale_color_viridis(discrete =T)+
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

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


# add lat lons
pred.obs.all <- predicts.all %>% #pred.obs.all %>%
  left_join(ncsites, by = 'site_name') %>%
  pivot_longer(cols = c(fit, log_den_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)

# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all, aes(x = year, y = latitude, color = values)) +
  geom_jitter(size = 4, pch = 16) +
  labs(x = 'Year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


###

###



###


###


###


# *** PREDICT BEST MODEL ACROSS ALL YEARS AND SITE THAT I HAVE DATA FOR ----

#gam5.1
# log_den_STRPURAD
# log_Days_16C
# Max_Monthly_Nitrate
# log_UBR_Max
# wh_mean


### Get depth ----

depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

names(dat2)

depth <- rast(paste(depth.dir, "depth_mean_nc_120res_30depth_latlon.tif", sep ='/'))
plot(depth)

n.extent <- ext(depth)

# crs1 <- "epsg:4326"
# d2 <- project(depth, crs1)
# 
# n.extent <- ext(d2)


## Get rock ----
sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

dir(sub.dir)

rock <- rast(paste(sub.dir, "prob_rock_mean_nc_120res_30depth_latlon.tif", sep ='/'))
rock 

# crop to NC --
test <- c(depth, rock)
# rock2 <- crop(rock, ext(depth))
# plot(rock2)
# 
# rock3 <- resample(rock2, d2)
# plot(rock3)


### Get Env predictors  ----

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"



### Get nitrate predictors  ----
max_nit <- rast(paste(re.dir, "Nitrate", "Max_Monthly_Nitrate", "Max_Monthly_Nitrate_extended_stack.tif", sep ='/'))
max_nit

# # crop to NC --
max_nit2 <- crop(max_nit, n.extent)
plot(max_nit2[[1]])

# resample predictors to bathy ----
max_nit3 <- resample(max_nit2, depth)
plot(max_nit3[[1]])

# mask predictors to bathy ----
max_nit4 <- mask(max_nit3, depth)
plot(max_nit4[[1]])


### Get NPP predictors  ----
npp <- rast(paste(re.dir, "NPP", "Mean_Monthly_NPP", "Mean_Monthly_NPP_extended_stack.tif", sep ='/'))
npp

# # crop to NC --
npp2 <- crop(npp, n.extent)
plot(npp2[[1]])

# resample predictors to bathy ----
npp3 <- resample(npp2, depth)
plot(npp3[[1]])

# mask predictors to bathy ----
npp4 <- mask(npp3, depth)
plot(npp4[[1]])


### Get Wave predictors  ----

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"

## Max Wave height --

# load raster data --

wave.dir <- paste(re.dir, "Waves", "wh_mean", sep ='/')


### Get wh_mean predictors  ----
wh <- rast(paste(wave.dir, "wh_mean_extended_stack.tif", sep ='/'))
wh

# # crop to NC --
wh2 <- crop(wh, n.extent)
plot(wh2[[1]])

# resample predictors to bathy ----
wh3 <- resample(wh2, depth)
plot(wh3[[1]])

# mask predictors to bathy ----
wh4 <- mask(wh3, depth)
plot(wh4[[1]])



## Mean UBR MAX ----

w2.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Waves_UBR_30m"

# load raster data --

ubr <- rast(paste(w2.dir, "UBR_Max_30m_NC.tif", sep ='/'))
plot(ubr[[1]])
ubr

ubr <- project(ubr, rock)

plot(ubr[[1]])

ubr1 <- classify(ubr, cbind(0, NA))
plot(ubr1[[1]])


# # crop to NC --
ubr2 <- crop(ubr1, n.extent)
plot(ubr2[[1]])

# resample predictors to bathy ----
ubr3 <- resample(ubr2, depth)
plot(ubr3[[1]])

# mask predictors to bathy ----
ubr4 <- mask(ubr3, depth)
plot(ubr4[[1]])

ubr5 <- log(ubr4)
plot(ubr5[[2]])


### Get TEMP predictors  ----

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"

## log_Days_16C --

### Get wh_mean predictors  ----
temp <- rast(paste(re.dir, "Temperature", "Mean_Monthly_Upwelling_Temp", "Mean_Monthly_Upwelling_Temp_extended_stack.tif", sep ='/'))
temp

# # crop to NC --
temp2 <- crop(temp, n.extent)
plot(temp2[[1]])

# resample predictors to bathy ----
temp3 <- resample(temp2, depth)
plot(temp3[[1]])

# mask predictors to bathy ----
temp4 <- mask(temp3, depth)
plot(temp4[[1]])

# log
# temp5 <- log(temp4 + 1)
# plot(temp5[[11]])


### Get urchins ----

# load purple urchin predictions ----
urch.dir <- paste(u.dir, "predictions", sep ='/')
urch.dir

urchins <- rast(paste(urch.dir, "log_STRPURAD_preds_NC_V4_120m_stack.tif", sep ='/'))
urchins

plot(urchins[[1]])
max.u <- max(dat2$log_den_STRPURAD) # 8.936298

# reclassify urchins to max limit
urchins2 <- classify(urchins, cbind(max.u, Inf, max.u))
plot(urchins2[[1]])

##


# YEAR ----

year1998 <- classify(max_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(2004:2021)
length(year.list)

 


# SITES ----

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




##

# ZONE ----

zone.raster <- depth
names(zone.raster) <- 'zone'
plot(zone.raster)
levels(dat2$zone)

rec.m <-  c(-Inf, -10, 2,
            -10, 0.1, 1)

rclmat <- matrix(rec.m, ncol=3, byrow=TRUE)

zone.raster2 <- classify(zone.raster, rclmat, right=FALSE)
plot(zone.raster2, col = c('red', 'blue'))




## LOOP - data frame ----



## loop to predict several years ----

nereo.mod <- gam1
summary(nereo.mod)

# make list of years --
#year.list <- paste(1998:2021)
year.list <- paste(2004:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(depth, cbind(-Inf, 0.1, 2004), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 



# outputs dir ----

#o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Spatio_temporal_GAMs/outputs_nc_rcca/gam_V4/gam_5.1"
preds.dir <- paste(o.dir, "preds_120m", sep ='/')
preds.dir

# output for rasters scaled by rock 
rock.preds.dir <- paste(o.dir, "rock_preds_120m", sep ='/')
rock.preds.dir

# Version 3: V4_5.1.1_v3 : Using urchins3, which don't have VRM
# Version sp_predictions_v3 : using gam 5.1, with urchins 3


for (i in 1:length(year.list)) {
  
  # 1. get urchins
  urchin.rast <- urchins2[[i]]
  urchin.rast2 <- resample(urchin.rast, depth)
  
  # 2. stack with predictors for that year
  #env.raster <- c(d2, max_nit4[[i+6]], whmax.stack4[[i]], wymean.stack4[[i]])
  
  #env.raster <- c(d2, mean_up_T4[[i+6]], max_nit4[[i+6]], whmax.stack4[[i]], wymean.stack4[[i]])
  
  # V3
  env.raster <- c(depth, max_nit4[[i+6]], temp4[[i+6]], wh4[[i]], ubr5[[i]], npp4[[i+1]])
  
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
                     "mean_depth",
                     "Max_Monthly_Nitrate"  , 
                     "Mean_Monthly_Upwelling_Temp", 
                     "wh_mean",  
                     "log_UBR_Max",
                     "Mean_Monthly_NPP",
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
  name.raster <- paste(year.no, "log_Nereo_preds_NC_V7_120m.tif", sep = '_')
  writeRaster(year.prediction, paste(preds.dir, name.raster, sep = '/'))
  
  # 8. scale by rock
  rock4 <- resample(rock, year.prediction)
  year.prediction2 <- rock4*year.prediction
  
  # 9. save rastr scaled by rock
  name.raster.rock <- paste(year.no, "log_Nereo_preds_rock_NC_V7_120m.tif", sep = '_')
  #writeRaster(year.prediction, paste(preds.dir, "sp_predictions", paste(year.no, "logNereo_preds_NC.tif", sep = '_'), sep ='/'))
  writeRaster(year.prediction2, paste(rock.preds.dir, name.raster.rock, sep = '/'))
}




###

###

### GET STABILITY ----



# load raster data --
# preds.dir <- paste(o.dir, "sp_predictions_5.1.1_V3_rock", sep ='/')
# preds.dir

preds.dir <- paste(o.dir, "rock_preds_120m", sep ='/')
preds.dir

#preds.dir <- rock.preds.dir

n.files <- dir(preds.dir)
# list files in source --
n.files <- list.files(preds.dir, pattern = '.tif', full.names = TRUE)
n.files
length(n.files)
#n.files <- n.files[-c(19, 20)]
# list names to load onto the Environment --
names.list <- list.files(preds.dir, pattern = '.tif')
names.list <- str_replace(names.list, ".tif", "")
length(names.list)
#names.list <- names.list[-c(19,20)]

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]]
length(tfiles)

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)
plot(preds.stack[[18]])

names(preds.stack) <- paste("Nereo", paste(2004:2021), sep ='_')
preds.stack

# save
writeRaster(preds.stack, paste(preds.dir, "rock_stack_log_Nereo_NC_V7_120m.tif", sep ='/'), overwrite = T)

# get csv --
dfn <- as.data.frame(preds.stack, xy = T) %>% glimpse()

write.csv(dfn, paste(preds.dir, "rock_stack_log_Nereo_NC_V7_120m.csv", sep ='/'))


## Cap at max log_nerlue ----
max.n <- max(dat2$log_den_NERLUE) # 6.216606
preds.stack2 <- classify(preds.stack, cbind(max.n, Inf, max.n))
plot(preds.stack2[[10]])

# save
writeRaster(preds.stack2, paste(preds.dir, "rock_stack_log_Nereo_NC_V7_120m_capped.tif", sep ='/'))

# get csv --
dfn <- as.data.frame(preds.stack2, xy = T) %>% glimpse()

write.csv(dfn, paste(preds.dir, "rock_stack_log_Nereo_NC_V7_120m_capped.csv", sep ='/'))


# get only  PRE MHW STABILITY ----
preds.stack3 <- preds.stack2[[1:10]]
preds.stack3 

### Calculate mean across years ----

stability.stack <- c()

# calculate mean across years
mean.nereo <- mean(preds.stack3, na.rm = T)
plot(mean.nereo)


# calculate sd across years
sd.nereo <- app(preds.stack3, fun = sd, na.rm = T)
plot(sd.nereo)


# calculate stability
s.nereo <- mean.nereo/sd.nereo
plot(s.nereo)

# calculate stability index
s_index.nereo <- (s.nereo*mean.nereo)
plot(s_index.nereo)

# stack 

stability.stack <- c(mean.nereo, sd.nereo, s.nereo, s_index.nereo)
names(stability.stack) <- c("mean_nereo", "sd_nereo", "s_nereo", "s_index_nereo")
plot(stability.stack)

# save 
#writeRaster(stability.stack, paste(preds.dir, "rock_Stabilty_PreMHW_calculations_Nereo_NC_V7_120m_capped.tif", sep ='/'), overwrite = T)

## get csv --

st.df <- as.data.frame(stability.stack, xy = T) %>% 
  mutate(s_nereo = replace(s_nereo, is.na(s_nereo), 0),
         s_index_nereo = replace(s_index_nereo, is.na(s_index_nereo), 0)) %>%
  glimpse()

write.csv(st.df, paste(preds.dir, "rock_Stabilty_PreMHW_calculations_Nereo_NC_V7_120m_capped.csv", sep ='/'))

# n.df <- as.data.frame(preds.stack, xy = T) %>% glimpse()
# write.csv(n.df, paste(preds.dir, "rock_log_Nereo_year_preds_V7.csv", sep ='/'))





###

###

# FOR 29 June Meeting ----

d.nc1 <- df.nc %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(den_NERLUE, na.rm = T),
            sd.nereo = sd(den_NERLUE, na.rm = T),
            n_nereo = length(den_NERLUE),
            se_nereo = sd.nereo/sqrt(n_nereo)) %>%
  glimpse()


ggplot(d.nc1, aes(x = site_name, y = mean_nereo)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin = mean_nereo-se_nereo, ymax = mean_nereo+se_nereo)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))



ncs <- vect(ncsites, geom = c("longitude", "latitude"), crs = "epsg:4326")
plot(stability.stack[[1]])
plot(ncs, add=T)

ncs2 <- extract(preds.stack, ncs)
head(ncs2)
str(ncs2)
#plot(preds.stack)

years.n <- paste('years', paste(2004:2021), sep = '_')
years.n <- paste(2004:2021)
names(ncs2) <- c("id", years.n)
names(ncs2)

ncsites
ncs2$site_name <- ncsites$site_name
head(ncs2)

ncs3 <- ncs2 %>%
  dplyr::select(-id) %>%
  pivot_longer(cols = '2004':'2021', names_to = "year", values_to = "nereo") %>% #glimpse()
  group_by(site_name) %>%
  summarise(mean_nereo2 = mean(nereo, na.rm = T), 
            sd_nereo2 = sd(nereo, na.rm = T), 
            n_nereo2 = length(nereo),
            se_nereo2 = sd_nereo2/sqrt(n_nereo2)) %>%
  glimpse()


ncs4 <- d.nc1 %>%
  left_join(ncs3, by = "site_name") %>%
  glimpse()


ggplot(ncs4) +
  geom_bar(aes(x = site_name, y = mean_nereo), stat = 'identity', fill = "skyblue") +
  geom_errorbar(aes(x = site_name, y = mean_nereo, ymin = mean_nereo-se_nereo, ymax = mean_nereo+se_nereo),
                width = 0.3, size = 0.8) +
  geom_point(aes(x = site_name, y = mean_nereo2), stat="identity", col = "red", size = 3, group = 1) +
  geom_errorbar(aes(x = site_name, y = mean_nereo2, ymin = mean_nereo2-se_nereo2, ymax = mean_nereo2+se_nereo2),
                col = 'red', width = 0.2, size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


#####


d.nc1 <- df %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(den_NERLUE, na.rm = T),
            sd.nereo = sd(den_NERLUE, na.rm = T),
            n_nereo = length(den_NERLUE),
            se_nereo = sd.nereo/sqrt(n_nereo)) %>%
  glimpse()


ggplot(d.nc1, aes(x = site_name, y = mean_nereo)) +
  geom_bar( stat = 'identity', fill = "skyblue") +
  geom_errorbar(aes(x = site_name, y = mean_nereo, ymin = mean_nereo-se_nereo, ymax = mean_nereo+se_nereo),
                width = 0.3, size = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))



ncs <- vect(years, geom = c("longitude", "latitude"), crs = "epsg:4326")
plot(stability.stack[[1]])
plot(ncs, add=T)

ncs2 <- extract(preds.stack, ncs)
head(ncs2)
str(ncs2)
#plot(preds.stack)

years.n <- paste('years', paste(2004:2021), sep = '_')
years.n <- paste(2004:2021)
names(ncs2) <- c("id", years.n)
names(ncs2)

years
ncs2$site_name <- years$site_name
head(ncs2)

ncs3 <- ncs2 %>%
  dplyr::select(-id) %>%
  pivot_longer(cols = '2004':'2021', names_to = "year", values_to = "nereo") %>% #glimpse()
  group_by(site_name) %>%
  summarise(mean_nereo2 = mean(nereo, na.rm = T), 
            sd_nereo2 = sd(nereo, na.rm = T), 
            n_nereo2 = length(nereo),
            se_nereo2 = sd_nereo2/sqrt(n_nereo2)) %>%
  glimpse()


ncs4 <- d.nc1 %>%
  left_join(ncs3, by = "site_name") %>%
  glimpse()


ggplot(ncs4) +
  geom_bar(aes(x = site_name, y = mean_nereo), stat = 'identity', fill = "skyblue") +
  geom_errorbar(aes(x = site_name, y = mean_nereo, ymin = mean_nereo-se_nereo, ymax = mean_nereo+se_nereo),
                width = 0.3, size = 0.8) +
  geom_point(aes(x = site_name, y = mean_nereo2), stat="identity", col = "red", size = 3, group = 1) +
  geom_errorbar(aes(x = site_name, y = mean_nereo2, ymin = mean_nereo2-se_nereo2, ymax = mean_nereo2+se_nereo2),
                col = 'red', width = 0.2, size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))




##### Save points -----
stability.stack
s.df <- as.data.frame(stability.stack[[1]], xy = T)
head(s.df)
write.csv(s.df, paste(preds.dir, "points_mean_nereo_preds_V2.csv", sep='/'))


# load raster data --
preds.dir <- paste(o.dir, "sp_predictions_5.1.1_V2_rock", sep ='/')
preds.dir

stability.stack <- rast(paste(preds.dir, "Stabilty_calculations_Nereo_NC_gam_V4_5.1.1_V2_rock.tif", sep ='/'))
plot(stability.stack[[1]])
