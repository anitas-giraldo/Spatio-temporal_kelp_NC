###
### Script created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 27 March 2022

## This script prepares the data for the density models of kelp in the north coast --

## resources --
# Box-Cox transformation
# https://www.statology.org/box-cox-transformation-in-r/
# https://robjhyndman.com/hyndsight/transformations/

### Load libraries ----

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

# Clear environment ----
rm(list=ls())

### Set directories ----
#w.dir<- dirname(rstudioapi::getActiveDocumentContext()$path)
m.dir <- here()
mlpa.dir <- "G:/Shared drives/MLPA L-T Kelp Forest Monitoring Project - Admin/Monitoring Data/MLPA merged datasets"
merged.data <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Merge_bio_env_vars_transect/data_tidy"
d.dir <- here("data")
raw.dir <- here("data_raw")
plots.dir <- here("plots")
o.dir <- here("outputs_gam")

## Load info on transect no. years ----
years <- read.csv(paste(merged.data, "No_survey_years_per_site.csv", sep ='/')) %>%
  glimpse()


# get the sites from Nc with preMHW data ----
scsites <- years %>%
  dplyr::select(-c(X, latitude, longitude)) %>%
  mutate_at(vars(site_campus_unique_ID), list(as.factor)) %>%
  # get north coast 
  dplyr::filter(region == 'sc') %>%
  # get only sites with PRE MHW data 
  dplyr::filter(preMHW > 0) %>%
  # get sites with 3 or more total survey years
  dplyr::filter(no.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 64

# gets lats and lons ----
sclatlons <- years %>%
  dplyr::select(-c(X)) %>%
  mutate_at(vars(site_campus_unique_ID), list(as.factor)) %>%
  # get north coast 
  dplyr::filter(region == 'sc') %>%
  # get only sites with PRE MHW data 
  dplyr::filter(preMHW > 0) %>%
  # get sites with 3 or more total survey years
  dplyr::filter(no.years > 2) %>%
  dplyr::select(site_campus_unique_ID, latitude, longitude) %>%
  droplevels() %>%
  glimpse() # Rows: 16


## Load transect data ----

df <- read.csv(paste(merged.data, "swath_upc_fish_terr_clim_env_depth_seaotter_temp_nitrate_transect.csv", sep ='/')) %>%
  mutate_at(vars(site_campus_unique_ID, survey_year), list(as.factor)) %>%
  glimpse()


## get the sites for North Coast model ----
df.sc <- df %>%
  dplyr::select(-c(X.x, Y.x)) %>%
  right_join(scsites, by = c('site_campus_unique_ID', 'mlpa_region', 'region')) %>%
  rename(X = X.y, Y = Y.y) %>%
  droplevels() %>%
  relocate(c(no.years, preMHW, duringMHW, postMHW), .after = longitude) %>%
  relocate(c(X, Y), .after = longitude) %>%
  glimpse() # Rows: 4,088

length(levels(df.sc$site_campus_unique_ID)) # 64

## Choose pre MHW years ----
levels(df.sc$survey_year)

df.sc <- df.sc %>%
  dplyr::filter(survey_year == "1999" |
                  survey_year == "2000" |
                  survey_year == "2001" |
                  survey_year == "2002" |
                  survey_year == "2003" |
                  survey_year == "2004" |
                  survey_year == "2005" |
                  survey_year == "2006" |
                  survey_year == "2007" |
                  survey_year == "2008" |
                  survey_year == "2009" |
                  survey_year == "2010" |
                  survey_year == "2011" |
                  survey_year == "2012" |
                  survey_year == "2013" ) %>%
  droplevels() %>%
  glimpse() # Rows: 2,579



### GAM V1 ----

#### Remove variables that are not necessary ----
names(df.sc)

dat <- df.sc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_NERLUE,
    #den_PANINT,
    count_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    #seaotter_dens_sm
    seaotter_pupratio:seaotter_trend5yr,
    Days_15C:Days_19C, 
    Days_23C:Degree_Days_19C, 
    Degree_Days_23C,
    Days_10N:Days_15N,
    Days_5N:Days_9N:Degree_Days_15N, 
    Degree_Days_5N:Degree_Days_9N,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    npp.mean,
    ID.x, ID.y
  )) %>%
  glimpse()

## UP TO HERE ----

#### remove correlated substrate terrain variables ----

dat2 <- dat %>%
  # transform the variable
  #mutate(mean_logvrm = log(mean_vrm)) %>%
  # remove untransformed 
  #dplyr::select(-mean_vrm) %>%
  # remove unselected variables from corr analysis
  dplyr::select(-c(
    # proportions mapped
    prop_map_prob_of_rock, prop_map_vrm, prop_map_slope
  )) %>% # 
  glimpse()


## 

#### Select predictors for this GAMs ----

names(dat2)

dat.gam1 <- dat2 %>%
  dplyr::select(# factors
    survey_year, zone, site_campus_unique_ID,
    # Substrate vars
    sd_prob_of_rock, min_prob_of_rock, mean_prob_of_rock, 
    max_depth, min_depth, depth_mean, 
    mean_slope, min_vrm, mean_vrm,
    # Nitrate vars 
    Max_Monthly_Nitrate, Mean_Monthly_Summer_Nitrate, Days_3N,
    Min_Monthly_Anomaly_Upwelling_Nitrate, Min_Monthly_Nitrate,
    # Temperature vars
    Max_Monthly_Anomaly_Summer_Temp,
    Mean_Monthly_Temp,
    Days_21C, 
    MHW_Upwelling_Days,
    # Bio vars
    den_MACPYRAD, den_MESFRAAD, den_STRPURAD, 
    #seaotter_dens_sm, # remove seaotters from correlations
    # Other vars
    npgo_mean, mei_mean, pdo_mean,
    #orb.max, orb.95,
    wh.max, wh.95) %>%
  mutate(log_min_vrm = log(min_vrm + 1),
         log_mean_vrm = log(mean_vrm + 1),
         log_Days_21C = log(Days_21C + 1),
         log_MHW_Upwelling_Days = log(MHW_Upwelling_Days + 1),
         log_den_MACPYRAD = log(den_MACPYRAD + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(min_vrm, mean_vrm, Days_21C, MHW_Upwelling_Days,
                   den_MACPYRAD, den_MESFRAAD, den_STRPURAD)) %>%
  mutate_at(vars(survey_year, zone), list(as.factor)) %>%
  glimpse() # 2,579 - 29 columns


#### Transform vars - Remove NAs ----

dat.gam1 <- dat.gam1 %>%
  # tranform biological variables
  #mutate(logden_NERLUE = log(den_NERLUE + 1)) %>% #glimpse()
  #dplyr::select(-c(den_NERLUE)) %>%
  drop_na() %>%
  droplevels() %>%
  glimpse() # Rows: 3,971

any(is.na(dat.gam1))


## Divide into train and test ----
inTraining <- createDataPartition(dat.gam1$log_den_MACPYRAD, p = 0.7, list = FALSE)
train.gam1 <- dat.gam1[ inTraining,]
test.gam1  <- dat.gam1[-inTraining,]

glimpse(train.gam1) # Rows: 1,747
glimpse(test.gam1) # Rows: 746




#### V2. Run the full subset model selection ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

names(dat.gam1)

# get train and test data sets ----
inTraining <- createDataPartition(dat.gam1$log_den_MACPYRAD, p = .5, list = FALSE)
train.gam1 <- dat.gam1[ inTraining,]
test.gam1  <- dat.gam1[-inTraining,]

name <- 'V2'
o.dir <- paste(m.dir, "outputs_sc", paste("gam", name, sep = '_'), sep = '/')

pred.vars <- c(# choose substrate vars --
  "sd_prob_of_rock" , 
  "min_prob_of_rock" , 
  "mean_prob_of_rock" ,                           
  "max_depth", 
  "min_depth" , 
  "depth_mean", 
  "mean_slope", 
  "log_min_vrm", 
  "log_mean_vrm"
  # choose other env variables --
  #"npgo_mean", "wh.95", "wh.max", "mei_mean", "pdo_mean",
  # choose temperature vars --
  #"Max_Monthly_Anomaly_Summer_Temp", 
  #"Mean_Monthly_Temp",
  #"log_Days_21C", 
  #"log_MHW_Upwelling_Days",
  #"Min_Monthly_Anomaly_Summer_Temp", 
  #"Max_Monthly_Anomaly_Temp",
  #"Mean_Monthly_Temp", 
  #"MHW_Days", "MHW_Upwelling_Days",
  #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
  #"Mean_Monthly_Upwelling_Temp",
  # choose nitrate vars --
  #"Max_Monthly_Nitrate", 
  #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
  #"Mean_Monthly_Summer_Nitrate", 
  #"Days_3N",
  #"Min_Monthly_Anomaly_Upwelling_Nitrate",
  # choose bio vars --
  #"log_den_MACPYRAD" 
  #"log_den_MESFRAAD",
  #"log_den_STRPURAD"
  )

length(pred.vars) # 10

#fact.vars <- c("survey_year")

model.v1 <- gam(log_den_MACPYRAD ~ 
                  #s(site_campus_unique_ID, bs = 're') +
                  #s(transect_unique, bs = 're') +
                  s(site_campus_unique_ID, zone, bs = 're') +
                  s(survey_year, bs = 're') ,
                  #te(depth_mean, wh.max,  bs = 'cr') +
                  #te(depth_mean, wh.95,  bs = 'cr'),
                #s(mean_slope) + 
                #s(mean_prob_of_rock) + s(depth_mean) + s(mean_logvrm) +
                #s(min_depth) + s(min_prob_of_rock) +
                #s(max_depth) 
                #s(npgo_mean) + s(npp.mean) + 
                #s(orb.95) + s(orb.max) + 
                #s(wh.95) + s(wh.max) +
                #s(Max_Monthly_Anomaly_Upwelling_Temp) + 
                #s(Max_Monthly_Temp) + 
                #s(Mean_Monthly_Summer_Temp) + 
                #s(Mean_Monthly_Temp) + 
                #s(Mean_Monthly_Upwelling_Temp) + 
                #s(Min_Monthly_Anomaly_Summer_Temp) +
                #s(Min_Monthly_Anomaly_Temp) + 
                #s(Min_Monthly_Temp) + 
                #s(logMHW_Days) + 
                #s(logMHW_Summer_Days) + 
                #s(logMHW_Upwelling_Days) +
                #s(logden_HALRUF) + 
                #s(logden_PYCHEL) + 
                #s(transect_unique, bs = 're'),
                data = train.gam1, 
                family = tw(),
                method = "REML") 

model.set <- generate.model.set(use.dat = train.gam1,
                                test.fit = model.v1,
                                pred.vars.cont = pred.vars,
                                #smooth.smooth.interactions = c(c("depth_mean", "wh.max"), c("depth_mean", "wh.95")),
                                max.predictors = 4,
                                cov.cutoff = 0.65, # cut off for correlations
                                #pred.vars.fact = fact.vars,
                                #linear.vars="Distance",
                                k=3,
                                null.terms = "s(site_campus_unique_ID, zone, bs = 're') +
                                s(survey_year, bs = 're')")

length(model.set$mod.formula) # 178
model.set$mod.formula[3]
out.list <- fit.model.set(model.set,
                          max.models= 500,
                          parallel=T)

beep()

# check results ----
out.all=list()
var.imp=list()

names(out.list)
length(out.list$success.models)
out.list$failed.models
out.list$mod.data.out[,1]
length(out.list$mod.data.out)
mod.table <- out.list$mod.data.out
test <- which(str_detect(mod.table$modname, "log_den_MACPYRAD"))
length(test)

#mod.table <- mod.table[-test,]
#length(mod.table2$modname)
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=4),]
out.i <- mod.table[c(1:5),]
nrow(out.i)
length(out.i)
names(out.i)
out.i[1]
out.i[2]
class(out.i)

out.all <- c(out.all,list(out.i))
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

out.list$failed.models
out.list$success.models[7]

length(out.list$failed.models)
length(out.list$success.models)


# plot the best models ----
# for(m in 1:nrow(out.i)){
#   best.model.name=as.character(out.i$modname[1])
#   
#   png(file=paste(o.dir, paste(name,m,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
#   if(best.model.name!="null"){
#     par(mfrow=c(3,1),mar=c(9,4,3,1))
#     best.model <- out.list$success.models[[best.model.name]]
#     plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
#     mtext(side=2,text='logden_NERLUE',outer=F)}  
#   dev.off()
# }


# Model fits and importance ----
names(out.all) <- 'logden_MACPYRAD'
names(var.imp) <- 'logden_MACPYRAD'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name,"best_models.csv",sep="_"), sep ='/'))
write.csv(all.var.imp, file=paste(o.dir, paste(name,"all.var.imp.csv",sep="_"), sep ='/'))

# Generic importance plots-
# heatmap.2(all.var.imp, notecex=0.4,  
#           #dendrogram ="none",
#           #col = colorRampPalette(c("white","yellow","red"))(10),
#           trace = "none", key.title = "", keysize=2,
#           notecol="black",key=T,
#           sepcolor = "black", margins=c(12,8), lhei=c(4,15),
#           Rowv=FALSE,Colv=FALSE)


# Part 2 - custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep = '_'), sep ='/')) %>% #from local copy
  rename(resp.var=X)%>%
  dplyr::select(-resp.var) %>% #glimpse()
  gather(key = predictor, value = importance) %>%
  arrange(importance) %>%
  mutate_at(vars(predictor), list(as.factor)) %>%
  glimpse()

## set the levels in order we want

dat.taxa %>%
  #count(predictor) %>% glimpse()
  mutate(predictor = fct_reorder(predictor, importance, .desc = TRUE)) %>%
  ggplot(aes(x = predictor, y = importance)) + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


# Manually make the most parsimonious GAM models ----

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = log_den_MACPYRAD ~ s(depth_mean, k = 3, bs = "cr") + 
              s(log_min_vrm, k = 3, bs = "cr") + 
              s(mean_slope, k = 3, bs = "cr") + 
              s(min_prob_of_rock, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + 
              s(survey_year, bs = "re"),
            family = tw(), data = train.gam1, method = "REML")


beep()
gam1$aic
gam1$deviance
summary(gam1)

29.4##

# Plot model results ----

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,1),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()
  

#### gam2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp, k = 3, bs = "cr") + 
              s(wh.max, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + 
              s(survey_year, bs = "re"),
            family = gaussian(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,1),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()
  

###

###


### GAM V2 ----

#### Remove variables that are not necessary ----

dat <- df.nc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_MACPYRAD,
    den_PANINT,
    pct_cov_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    seaotter_dens_sm:seaotter_trend5yr,
    Days_16C, Days_18C:Degree_Days_16C, Degree_Days_18C:Degree_Days_23C,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    npp.mean
  )) %>%
  glimpse()


#### remove correlated substrate terrain variables ----

dat2 <- dat %>%
  # transform the variable
  mutate(mean_logvrm = log(mean_vrm)) %>%
  # remove untransformed 
  dplyr::select(-mean_vrm) %>%
  # remove unselected variables from corr analysis
  dplyr::select(-c(
    # all medians
    median_slope, median_depth, median_vrm, median_prob_of_rock,
    # all ranges
    range_slope, range_depth, range_vrm, range_prob_of_rock,
    # mins
    min_slope, min_vrm, 
    # sd
    sd_slope, sd_vrm,
    # max
    max_slope, max_vrm,
    # proportions mapped
    prop_map_prob_of_rock, prop_map_vrm, prop_map_slope,
    # max prob rock is always 1
    max_prob_of_rock,
    # days above 17C always 0
    #Days_17C, Degree_Days_17C,
    # the same as Max_Monthly_Anomaly_Temp
    #Max_Monthly_Anomaly_Summer_Temp,
    # the same as Min_Monthly_Anomaly_Temp
    #Min_Monthly_Anomaly_Upwelling_Temp,
    # MHW Days are the same as Intensity
    #MHW_Intensity, MHW_Summer_Intensity, MHW_Upwelling_Intensity,
    # remove other correlated vars
    mei_mean, pdo_mean, 
    #Max_Monthly_Anomaly_Temp,
    sd_depth, sd_prob_of_rock
  )) %>% # 
  glimpse()


## 

#### Select predictors for this GAMs ----

names(dat2)

dat.gam1 <- dat2 %>%
  dplyr::select(# get factors 
    "survey_year", "site_campus_unique_ID", "zone", 
    # get dependent var
    "den_NERLUE",
    # choose substrate vars
    "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
    "min_prob_of_rock" , "max_depth", "depth_mean",
    # choose other env variables
    "npgo_mean", "wh.95", "wh.max", 
    # choose temperature vars
    "Max_Monthly_Anomaly_Summer_Temp", 
    "Max_Monthly_Temp",
    "Days_15C", 
    "MHW_Summer_Days",
    "Min_Monthly_Anomaly_Summer_Temp", 
    #"Max_Monthly_Anomaly_Temp",
    "Mean_Monthly_Temp", 
    #"MHW_Days", "MHW_Upwelling_Days",
    #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
    #"Mean_Monthly_Upwelling_Temp",
    # choose nitrate vars
    "Max_Monthly_Nitrate", 
    "Min_Monthly_Anomaly_Summer_Nitrate",
    "Max_Monthly_Anomaly_Summer_Nitrate",
    #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
    "Mean_Monthly_Nitrate" 
    #"Max_Monthly_Anomaly_Nitrate"
  ) %>%
  mutate_at(vars(survey_year, zone), list(as.factor)) %>%
  glimpse() # Rows: 478 - 24 columns


#### Transform vars - Remove NAs ----

dat.gam1 <- dat.gam1 %>%
  # tranform biological variables
  mutate(logden_NERLUE = log(den_NERLUE + 1)) %>% #glimpse()
  dplyr::select(-c(den_NERLUE)) %>%
  drop_na() %>%
  droplevels() %>%
  glimpse() # Rows: 460

any(is.na(dat.gam1))



#### V1. Run the full subset model selection ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

name <- 'V2'
o.dir <- paste("outputs_gam", name, sep = '_')

pred.vars <- c(# choose substrate vars
  "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
  "min_prob_of_rock" , "max_depth", "depth_mean",
  # choose other env variables
  "npgo_mean", "wh.95", "wh.max", 
  # choose temperature vars
  "Max_Monthly_Anomaly_Summer_Temp", 
  "Max_Monthly_Temp",
  "Days_15C", 
  "MHW_Summer_Days",
  "Min_Monthly_Anomaly_Summer_Temp", 
  #"Max_Monthly_Anomaly_Temp",
  "Mean_Monthly_Temp", 
  #"MHW_Days", "MHW_Upwelling_Days",
  #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
  #"Mean_Monthly_Upwelling_Temp",
  # choose nitrate vars
  "Max_Monthly_Nitrate", 
  #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
  #"Mean_Monthly_Nitrate", "Max_Monthly_Anomaly_Nitrate",
  "Min_Monthly_Anomaly_Summer_Nitrate", "Max_Monthly_Anomaly_Summer_Nitrate")

length(pred.vars) # 19

#fact.vars <- c("survey_year")

model.v1 <- gam(logden_NERLUE ~ 
                  #s(site_campus_unique_ID, bs = 're') +
                  #s(transect_unique, bs = 're') +
                  s(site_campus_unique_ID, zone, bs = 're') +
                  s(survey_year, bs = 're') ,
                #te(depth_mean, wh.max,  bs = 'cr') +
                #te(depth_mean, wh.95,  bs = 'cr'),
                #s(mean_slope) + 
                #s(mean_prob_of_rock) + s(depth_mean) + s(mean_logvrm) +
                #s(min_depth) + s(min_prob_of_rock) +
                #s(max_depth) 
                #s(npgo_mean) + s(npp.mean) + 
                #s(orb.95) + s(orb.max) + 
                #s(wh.95) + s(wh.max) +
                #s(Max_Monthly_Anomaly_Upwelling_Temp) + 
                #s(Max_Monthly_Temp) + 
                #s(Mean_Monthly_Summer_Temp) + 
                #s(Mean_Monthly_Temp) + 
                #s(Mean_Monthly_Upwelling_Temp) + 
                #s(Min_Monthly_Anomaly_Summer_Temp) +
                #s(Min_Monthly_Anomaly_Temp) + 
                #s(Min_Monthly_Temp) + 
                #s(logMHW_Days) + 
                #s(logMHW_Summer_Days) + 
                #s(logMHW_Upwelling_Days) +
                #s(logden_HALRUF) + 
                #s(logden_PYCHEL) + 
                #s(transect_unique, bs = 're'),
                data = dat.gam1, 
                family = tw(),
                method = "REML") 

model.set <- generate.model.set(use.dat = dat.gam1,
                                test.fit = model.v1,
                                pred.vars.cont = pred.vars,
                                #smooth.smooth.interactions = c(c("depth_mean", "wh.max"), c("depth_mean", "wh.95")),
                                max.predictors = 4,
                                cov.cutoff = 0.6, # cut off for correlations
                                #pred.vars.fact = fact.vars,
                                #linear.vars="Distance",
                                k=3,
                                null.terms = "s(site_campus_unique_ID, zone, bs = 're') +
                                s(survey_year, bs = 're')")


out.list <- fit.model.set(model.set,
                          max.models= 2000,
                          parallel=T)

beep()

# check results ----
out.all=list()
var.imp=list()

names(out.list)
length(out.list$success.models)
out.list$failed.models
out.list$mod.data.out[,1]
length(out.list$mod.data.out)
mod.table <- out.list$mod.data.out
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=3),]
nrow(out.i)
length(out.i)
names(out.i)
out.i[1]
out.i[2]
class(out.i)

out.all <- c(out.all,list(out.i))
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

out.list$failed.models
out.list$success.models[7]

length(out.list$failed.models)
length(out.list$success.models)


# plot the best models ----
# for(m in 1:nrow(out.i)){
#   best.model.name=as.character(out.i$modname[1])
#   
#   png(file=paste(o.dir, paste(name,m,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
#   if(best.model.name!="null"){
#     par(mfrow=c(3,1),mar=c(9,4,3,1))
#     best.model <- out.list$success.models[[best.model.name]]
#     plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
#     mtext(side=2,text='logden_NERLUE',outer=F)}  
#   dev.off()
# }


# Model fits and importance ----
names(out.all) <- 'logden_NERLUE'
names(var.imp) <- 'logden_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name,"best_models.csv",sep="_"), sep ='/'))
write.csv(all.var.imp, file=paste(o.dir, paste(name,"all.var.imp.csv",sep="_"), sep ='/'))

# Generic importance plots-
# heatmap.2(all.var.imp, notecex=0.4,  
#           #dendrogram ="none",
#           #col = colorRampPalette(c("white","yellow","red"))(10),
#           trace = "none", key.title = "", keysize=2,
#           notecol="black",key=T,
#           sepcolor = "black", margins=c(12,8), lhei=c(4,15),
#           Rowv=FALSE,Colv=FALSE)


# Part 2 - custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep = '_'), sep ='/')) %>% #from local copy
  rename(resp.var=X)%>%
  dplyr::select(-resp.var) %>% #glimpse()
  gather(key = predictor, value = importance) %>%
  arrange(importance) %>%
  mutate_at(vars(predictor), list(as.factor)) %>%
  glimpse()

## set the levels in order we want

dat.taxa %>%
  #count(predictor) %>% glimpse()
  mutate(predictor = fct_reorder(predictor, importance, .desc = TRUE)) %>%
  ggplot(aes(x = predictor, y = importance)) + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


# Manually make the most parsimonious GAM models ----

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(MHW_Summer_Days, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Nitrate,  k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + s(survey_year, bs = "re"), 
            family = tw(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,1),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


#### gam2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(mean_slope, k = 3, bs = "cr") + 
              s(MHW_Summer_Days, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + 
              s(survey_year, bs = "re"),
            family = tw(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,1),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###

###

### GAM V3 ----

#### Remove variables that are not necessary ----

dat <- df.nc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_MACPYRAD,
    den_PANINT,
    pct_cov_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    seaotter_dens_sm:seaotter_trend5yr,
    Days_16C, Days_18C:Degree_Days_16C, Degree_Days_18C:Degree_Days_23C,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    npp.mean
  )) %>%
  glimpse()


#### remove correlated substrate terrain variables ----

dat2 <- dat %>%
  # transform the variable
  mutate(mean_logvrm = log(mean_vrm)) %>%
  # remove untransformed 
  dplyr::select(-mean_vrm) %>%
  # remove unselected variables from corr analysis
  dplyr::select(-c(
    # all medians
    median_slope, median_depth, median_vrm, median_prob_of_rock,
    # all ranges
    range_slope, range_depth, range_vrm, range_prob_of_rock,
    # mins
    min_slope, min_vrm, 
    # sd
    sd_slope, sd_vrm,
    # max
    max_slope, max_vrm,
    # proportions mapped
    prop_map_prob_of_rock, prop_map_vrm, prop_map_slope,
    # max prob rock is always 1
    max_prob_of_rock,
    # days above 17C always 0
    #Days_17C, Degree_Days_17C,
    # the same as Max_Monthly_Anomaly_Temp
    #Max_Monthly_Anomaly_Summer_Temp,
    # the same as Min_Monthly_Anomaly_Temp
    #Min_Monthly_Anomaly_Upwelling_Temp,
    # MHW Days are the same as Intensity
    #MHW_Intensity, MHW_Summer_Intensity, MHW_Upwelling_Intensity,
    # remove other correlated vars
    mei_mean, pdo_mean, 
    #Max_Monthly_Anomaly_Temp,
    sd_depth, sd_prob_of_rock
  )) %>% # 
  glimpse()


## 

#### Select predictors for this GAMs ----

names(dat2)

dat.gam1 <- dat2 %>%
  dplyr::select(# get factors 
    "survey_year", "site_campus_unique_ID", "zone", 
    # get dependent var
    "den_NERLUE",
    # choose substrate vars
    "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
    "min_prob_of_rock" , "max_depth", "depth_mean",
    # choose other env variables
    "npgo_mean", "wh.95", "wh.max", 
    # choose temperature vars
    "Max_Monthly_Anomaly_Summer_Temp", 
    "Max_Monthly_Temp",
    "Days_15C", 
    "MHW_Summer_Days",
    "Min_Monthly_Anomaly_Summer_Temp", 
    #"Max_Monthly_Anomaly_Temp",
    "Mean_Monthly_Temp", 
    #"MHW_Days", "MHW_Upwelling_Days",
    #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
    #"Mean_Monthly_Upwelling_Temp",
    # choose nitrate vars
    "Max_Monthly_Nitrate", 
    "Min_Monthly_Anomaly_Summer_Nitrate",
    "Max_Monthly_Anomaly_Summer_Nitrate",
    #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
    "Mean_Monthly_Nitrate" 
    #"Max_Monthly_Anomaly_Nitrate"
  ) %>%
  mutate_at(vars(survey_year, zone), list(as.factor)) %>%
  glimpse() # Rows: 478 - 24 columns


#### Transform vars - Remove NAs ----

dat.gam1 <- dat.gam1 %>%
  # tranform biological variables
  mutate(logden_NERLUE = log(den_NERLUE + 1)) %>% #glimpse()
  dplyr::select(-c(den_NERLUE)) %>%
  drop_na() %>%
  droplevels() %>%
  glimpse() # Rows: 460

any(is.na(dat.gam1))



#### V1. Run the full subset model selection ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

name <- 'V3'
o.dir <- paste("outputs_gam", name, sep = '_')

pred.vars <- c(# choose substrate vars
  "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
  "min_prob_of_rock" , "max_depth", "depth_mean",
  # choose other env variables
  "npgo_mean", "wh.95", "wh.max", 
  # choose temperature vars
  "Max_Monthly_Anomaly_Summer_Temp", 
  "Max_Monthly_Temp",
  "Days_15C", 
  "MHW_Summer_Days",
  "Min_Monthly_Anomaly_Summer_Temp", 
  #"Max_Monthly_Anomaly_Temp",
  "Mean_Monthly_Temp", 
  #"MHW_Days", "MHW_Upwelling_Days",
  #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
  #"Mean_Monthly_Upwelling_Temp",
  # choose nitrate vars
  "Max_Monthly_Nitrate", 
  #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
  #"Mean_Monthly_Nitrate", "Max_Monthly_Anomaly_Nitrate",
  "Min_Monthly_Anomaly_Summer_Nitrate", "Max_Monthly_Anomaly_Summer_Nitrate")

length(pred.vars) # 19

#fact.vars <- c("survey_year")

model.v1 <- gam(logden_NERLUE ~ 
                  #s(site_campus_unique_ID, bs = 're') +
                  #s(transect_unique, bs = 're') +
                  s(site_campus_unique_ID, zone, bs = 're') +
                  s(survey_year, bs = 're') ,
                #te(depth_mean, wh.max,  bs = 'cr') +
                #te(depth_mean, wh.95,  bs = 'cr'),
                #s(mean_slope) + 
                #s(mean_prob_of_rock) + s(depth_mean) + s(mean_logvrm) +
                #s(min_depth) + s(min_prob_of_rock) +
                #s(max_depth) 
                #s(npgo_mean) + s(npp.mean) + 
                #s(orb.95) + s(orb.max) + 
                #s(wh.95) + s(wh.max) +
                #s(Max_Monthly_Anomaly_Upwelling_Temp) + 
                #s(Max_Monthly_Temp) + 
                #s(Mean_Monthly_Summer_Temp) + 
                #s(Mean_Monthly_Temp) + 
                #s(Mean_Monthly_Upwelling_Temp) + 
                #s(Min_Monthly_Anomaly_Summer_Temp) +
                #s(Min_Monthly_Anomaly_Temp) + 
                #s(Min_Monthly_Temp) + 
                #s(logMHW_Days) + 
                #s(logMHW_Summer_Days) + 
                #s(logMHW_Upwelling_Days) +
                #s(logden_HALRUF) + 
                #s(logden_PYCHEL) + 
                #s(transect_unique, bs = 're'),
                data = dat.gam1, 
                family = tw(),
                method = "REML") 

model.set <- generate.model.set(use.dat = dat.gam1,
                                test.fit = model.v1,
                                pred.vars.cont = pred.vars,
                                #smooth.smooth.interactions = c(c("depth_mean", "wh.max"), c("depth_mean", "wh.95")),
                                max.predictors = 5,
                                cov.cutoff = 0.6, # cut off for correlations
                                #pred.vars.fact = fact.vars,
                                #linear.vars="Distance",
                                k=3,
                                null.terms = "s(site_campus_unique_ID, zone, bs = 're') +
                                s(survey_year, bs = 're')")


out.list <- fit.model.set(model.set,
                          max.models= 2000,
                          parallel=T)

beep()

# check results ----
out.all=list()
var.imp=list()

names(out.list)
length(out.list$success.models)
out.list$failed.models
out.list$mod.data.out[,1]
length(out.list$mod.data.out)
mod.table <- out.list$mod.data.out
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=3),]
nrow(out.i)
length(out.i)
names(out.i)
out.i[1]
out.i[2]
class(out.i)

out.all <- c(out.all,list(out.i))
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

out.list$failed.models
out.list$success.models[7]

length(out.list$failed.models)
length(out.list$success.models)


# plot the best models ----
# for(m in 1:nrow(out.i)){
#   best.model.name=as.character(out.i$modname[1])
#   
#   png(file=paste(o.dir, paste(name,m,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
#   if(best.model.name!="null"){
#     par(mfrow=c(3,1),mar=c(9,4,3,1))
#     best.model <- out.list$success.models[[best.model.name]]
#     plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
#     mtext(side=2,text='logden_NERLUE',outer=F)}  
#   dev.off()
# }


# Model fits and importance ----
names(out.all) <- 'logden_NERLUE'
names(var.imp) <- 'logden_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name,"best_models.csv",sep="_"), sep ='/'))
write.csv(all.var.imp, file=paste(o.dir, paste(name,"all.var.imp.csv",sep="_"), sep ='/'))

# Generic importance plots-
# heatmap.2(all.var.imp, notecex=0.4,  
#           #dendrogram ="none",
#           #col = colorRampPalette(c("white","yellow","red"))(10),
#           trace = "none", key.title = "", keysize=2,
#           notecol="black",key=T,
#           sepcolor = "black", margins=c(12,8), lhei=c(4,15),
#           Rowv=FALSE,Colv=FALSE)


# Part 2 - custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep = '_'), sep ='/')) %>% #from local copy
  rename(resp.var=X)%>%
  dplyr::select(-resp.var) %>% #glimpse()
  gather(key = predictor, value = importance) %>%
  arrange(importance) %>%
  mutate_at(vars(predictor), list(as.factor)) %>%
  glimpse()

## set the levels in order we want

dat.taxa %>%
  #count(predictor) %>% glimpse()
  mutate(predictor = fct_reorder(predictor, importance, .desc = TRUE)) %>%
  ggplot(aes(x = predictor, y = importance)) + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


# Manually make the most parsimonious GAM models ----

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(Max_Monthly_Temp, k = 3, bs = "cr") + s(MHW_Summer_Days, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Nitrate, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + s(survey_year,  bs = "re"), 
            family = tw(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###

###


### GAM V4 ----

#### Remove variables that are not necessary ----

dat <- df.nc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_MACPYRAD,
    den_PANINT,
    pct_cov_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    seaotter_dens_sm:seaotter_trend5yr,
    Days_16C, Days_18C:Degree_Days_16C, Degree_Days_18C:Degree_Days_23C,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    npp.mean
  )) %>%
  glimpse()


#### remove correlated substrate terrain variables ----

dat2 <- dat %>%
  # transform the variable
  mutate(mean_logvrm = log(mean_vrm)) %>%
  # remove untransformed 
  dplyr::select(-mean_vrm) %>%
  # remove unselected variables from corr analysis
  dplyr::select(-c(
    # all medians
    median_slope, median_depth, median_vrm, median_prob_of_rock,
    # all ranges
    range_slope, range_depth, range_vrm, range_prob_of_rock,
    # mins
    min_slope, min_vrm, 
    # sd
    sd_slope, sd_vrm,
    # max
    max_slope, max_vrm,
    # proportions mapped
    prop_map_prob_of_rock, prop_map_vrm, prop_map_slope,
    # max prob rock is always 1
    max_prob_of_rock,
    # days above 17C always 0
    #Days_17C, Degree_Days_17C,
    # the same as Max_Monthly_Anomaly_Temp
    #Max_Monthly_Anomaly_Summer_Temp,
    # the same as Min_Monthly_Anomaly_Temp
    #Min_Monthly_Anomaly_Upwelling_Temp,
    # MHW Days are the same as Intensity
    #MHW_Intensity, MHW_Summer_Intensity, MHW_Upwelling_Intensity,
    # remove other correlated vars
    mei_mean, pdo_mean, 
    #Max_Monthly_Anomaly_Temp,
    sd_depth, sd_prob_of_rock
  )) %>% # 
  glimpse()


## 

#### Select predictors for this GAMs ----

names(dat2)

dat.gam1 <- dat2 %>%
  dplyr::select(# get factors 
    "survey_year", "site_campus_unique_ID", "zone", 
    # get dependent var
    "den_NERLUE",
    # choose substrate vars
    "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
    "min_prob_of_rock" , "max_depth", "depth_mean",
    # choose other env variables
    "npgo_mean", "wh.95", "wh.max", 
    # choose temperature vars
    "Max_Monthly_Anomaly_Summer_Temp", 
    "Max_Monthly_Temp",
    "Days_15C", 
    "MHW_Summer_Days",
    "Min_Monthly_Anomaly_Summer_Temp", 
    #"Max_Monthly_Anomaly_Temp",
    "Mean_Monthly_Temp", 
    #"MHW_Days", "MHW_Upwelling_Days",
    #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
    #"Mean_Monthly_Upwelling_Temp",
    # choose nitrate vars
    "Max_Monthly_Nitrate", 
    "Min_Monthly_Anomaly_Summer_Nitrate",
    "Max_Monthly_Anomaly_Summer_Nitrate",
    #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
    "Mean_Monthly_Nitrate" 
    #"Max_Monthly_Anomaly_Nitrate"
  ) %>%
  mutate_at(vars(survey_year, zone), list(as.factor)) %>%
  glimpse() # Rows: 478 - 24 columns


#### Transform vars - Remove NAs ----

dat.gam1 <- dat.gam1 %>%
  # tranform biological variables
  mutate(logden_NERLUE = log(den_NERLUE + 1)) %>% #glimpse()
  dplyr::select(-c(den_NERLUE)) %>%
  drop_na() %>%
  droplevels() %>%
  glimpse() # Rows: 460

any(is.na(dat.gam1))



#### V1. Run the full subset model selection ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

name <- 'V4'
o.dir <- paste("outputs_gam", name, sep = '_')

pred.vars <- c(# choose substrate vars
  "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
  "min_prob_of_rock" , "max_depth", "depth_mean",
  # choose other env variables
  "npgo_mean", "wh.95", "wh.max", 
  # choose temperature vars
  "Max_Monthly_Anomaly_Summer_Temp", 
  "Max_Monthly_Temp",
  "Days_15C", 
  "MHW_Summer_Days",
  "Min_Monthly_Anomaly_Summer_Temp", 
  #"Max_Monthly_Anomaly_Temp",
  "Mean_Monthly_Temp", 
  #"MHW_Days", "MHW_Upwelling_Days",
  #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
  #"Mean_Monthly_Upwelling_Temp",
  # choose nitrate vars
  "Max_Monthly_Nitrate", 
  #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
  #"Mean_Monthly_Nitrate", "Max_Monthly_Anomaly_Nitrate",
  "Min_Monthly_Anomaly_Summer_Nitrate", "Max_Monthly_Anomaly_Summer_Nitrate")

length(pred.vars) # 19

#fact.vars <- c("survey_year")

model.v1 <- gam(logden_NERLUE ~ 
                  #s(site_campus_unique_ID, bs = 're') +
                  #s(transect_unique, bs = 're') +
                  s(site_campus_unique_ID, zone, bs = 're') +
                  s(survey_year, bs = 're') ,
                #te(depth_mean, wh.max,  bs = 'cr') +
                #te(depth_mean, wh.95,  bs = 'cr'),
                #s(mean_slope) + 
                #s(mean_prob_of_rock) + s(depth_mean) + s(mean_logvrm) +
                #s(min_depth) + s(min_prob_of_rock) +
                #s(max_depth) 
                #s(npgo_mean) + s(npp.mean) + 
                #s(orb.95) + s(orb.max) + 
                #s(wh.95) + s(wh.max) +
                #s(Max_Monthly_Anomaly_Upwelling_Temp) + 
                #s(Max_Monthly_Temp) + 
                #s(Mean_Monthly_Summer_Temp) + 
                #s(Mean_Monthly_Temp) + 
                #s(Mean_Monthly_Upwelling_Temp) + 
                #s(Min_Monthly_Anomaly_Summer_Temp) +
                #s(Min_Monthly_Anomaly_Temp) + 
                #s(Min_Monthly_Temp) + 
                #s(logMHW_Days) + 
                #s(logMHW_Summer_Days) + 
                #s(logMHW_Upwelling_Days) +
                #s(logden_HALRUF) + 
                #s(logden_PYCHEL) + 
                #s(transect_unique, bs = 're'),
                data = dat.gam1, 
                family = tw(),
                method = "REML") 

model.set <- generate.model.set(use.dat = dat.gam1,
                                test.fit = model.v1,
                                pred.vars.cont = pred.vars,
                                #smooth.smooth.interactions = c(c("depth_mean", "wh.max"), c("depth_mean", "wh.95")),
                                max.predictors = 6,
                                cov.cutoff = 0.6, # cut off for correlations
                                #pred.vars.fact = fact.vars,
                                #linear.vars="Distance",
                                k=3,
                                null.terms = "s(site_campus_unique_ID, zone, bs = 're') +
                                s(survey_year, bs = 're')")


out.list <- fit.model.set(model.set,
                          max.models= 2000,
                          parallel=T)

beep()

# check results ----
out.all=list()
var.imp=list()

names(out.list)
length(out.list$success.models)
out.list$failed.models
out.list$mod.data.out[,1]
length(out.list$mod.data.out)
mod.table <- out.list$mod.data.out
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=3),]
nrow(out.i)
length(out.i)
names(out.i)
out.i[1]
out.i[2]
class(out.i)

out.all <- c(out.all,list(out.i))
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

out.list$failed.models
out.list$success.models[7]

length(out.list$failed.models)
length(out.list$success.models)


# plot the best models ----
# for(m in 1:nrow(out.i)){
#   best.model.name=as.character(out.i$modname[1])
#   
#   png(file=paste(o.dir, paste(name,m,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
#   if(best.model.name!="null"){
#     par(mfrow=c(3,1),mar=c(9,4,3,1))
#     best.model <- out.list$success.models[[best.model.name]]
#     plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
#     mtext(side=2,text='logden_NERLUE',outer=F)}  
#   dev.off()
# }


# Model fits and importance ----
names(out.all) <- 'logden_NERLUE'
names(var.imp) <- 'logden_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name,"best_models.csv",sep="_"), sep ='/'))
write.csv(all.var.imp, file=paste(o.dir, paste(name,"all.var.imp.csv",sep="_"), sep ='/'))

# Generic importance plots-
# heatmap.2(all.var.imp, notecex=0.4,  
#           #dendrogram ="none",
#           #col = colorRampPalette(c("white","yellow","red"))(10),
#           trace = "none", key.title = "", keysize=2,
#           notecol="black",key=T,
#           sepcolor = "black", margins=c(12,8), lhei=c(4,15),
#           Rowv=FALSE,Colv=FALSE)


# Part 2 - custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep = '_'), sep ='/')) %>% #from local copy
  rename(resp.var=X)%>%
  dplyr::select(-resp.var) %>% #glimpse()
  gather(key = predictor, value = importance) %>%
  arrange(importance) %>%
  mutate_at(vars(predictor), list(as.factor)) %>%
  glimpse()

## set the levels in order we want

dat.taxa %>%
  #count(predictor) %>% glimpse()
  mutate(predictor = fct_reorder(predictor, importance, .desc = TRUE)) %>%
  ggplot(aes(x = predictor, y = importance)) + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


# Manually make the most parsimonious GAM models ----

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(Max_Monthly_Temp, k = 3, bs = "cr") + 
              s(mean_slope, k = 3, bs = "cr") + 
              s(MHW_Summer_Days, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Nitrate, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone,  bs = "re") + s(survey_year, bs = "re"), 
            family = tw(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()

###


#### gam2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(Max_Monthly_Temp, k = 3, bs = "cr") + 
              s(MHW_Summer_Days, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Nitrate, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp, k = 3, bs = "cr") + 
              s(wh.max, k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + s(survey_year, bs = "re"), 
            family = tw(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###


#### gam3 ----

subname <- "3"

best.model.name=as.character(out.i$modname[3])
best.model <- out.list$success.models[[best.model.name]]
best.model

gam1 <- gam(formula = logden_NERLUE ~ s(depth_mean, k = 3, bs = "cr") + 
              s(Max_Monthly_Temp, k = 3, bs = "cr") + 
              s(mean_prob_of_rock, k = 3, bs = "cr") + 
              s(MHW_Summer_Days, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Nitrate, k = 3, bs = "cr") + 
              s(Min_Monthly_Anomaly_Summer_Temp,  k = 3, bs = "cr") + 
              s(site_campus_unique_ID, zone, bs = "re") + s(survey_year, bs = "re"), 
            family = tw(), data = dat.gam1, method = "REML")


gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###

#### gam4 ----

subname <- "4"

best.model.name=as.character(out.i$modname[4])
best.model <- out.list$success.models[[best.model.name]]
best.model


# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


# predict - gam5 ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(depth_mean=mean(mod$model$depth_mean),
                        Max_Monthly_Temp=mean(mod$model$Max_Monthly_Temp),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        MHW_Summer_Days=mean(mod$model$MHW_Summer_Days),
                        Min_Monthly_Anomaly_Summer_Nitrate=mean(mod$model$Min_Monthly_Anomaly_Summer_Nitrate),
                        Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
                        site_campus_unique_ID=(mod$model$site_campus_unique_ID),
                        zone=(mod$model$zone),
                        survey_year=(mod$model$survey_year))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.year = testdata%>%data.frame(fits)%>%
  group_by(survey_year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=survey_year,y=response,fill=survey_year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
                    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year 


###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- dat.gam1 %>%
  dplyr::select(survey_year, site_campus_unique_ID, zone, logden_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  left_join(ner.obs, by = c('survey_year', 'site_campus_unique_ID', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = logden_NERLUE, y = fit, color = zone)) +
  geom_point() 



# add lat lons
pred.obs.all <- pred.obs.all %>%
  left_join(nclatlons, by = 'site_campus_unique_ID') %>%
  pivot_longer(cols = c(fit, logden_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)

# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all, aes(x = zone, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

ggplot(pred.obs.all, aes(x = survey_year, y = latitude, color = values, shape = type)) +
  geom_jitter(size = 2.5) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~zone) +
  theme_bw()

head(preds)



##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###


#### gam5 ----

subname <- "5"

best.model.name=as.character(out.i$modname[5])
best.model <- out.list$success.models[[best.model.name]]
best.model
class(best.model)


# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


# predict - gam5 ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(depth_mean=mean(mod$model$depth_mean),
                        Max_Monthly_Temp=mean(mod$model$Max_Monthly_Temp),
                        mean_logvrm=mean(mod$model$mean_logvrm),
                        MHW_Summer_Days=mean(mod$model$MHW_Summer_Days),
                        Min_Monthly_Anomaly_Summer_Nitrate=mean(mod$model$Min_Monthly_Anomaly_Summer_Nitrate),
                        Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
                        site_campus_unique_ID=(mod$model$site_campus_unique_ID),
                        zone=(mod$model$zone),
                        survey_year=(mod$model$survey_year))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.year = testdata%>%data.frame(fits)%>%
  group_by(survey_year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=survey_year,y=response,fill=survey_year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
                    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year 


###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- dat.gam1 %>%
  dplyr::select(survey_year, site_campus_unique_ID, zone, logden_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  left_join(ner.obs, by = c('survey_year', 'site_campus_unique_ID', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = logden_NERLUE, y = fit, color = zone)) +
  geom_point() 



# add lat lons
pred.obs.all <- pred.obs.all %>%
  left_join(nclatlons, by = 'site_campus_unique_ID') %>%
  pivot_longer(cols = c(fit, logden_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)

# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all, aes(x = zone, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

ggplot(pred.obs.all, aes(x = survey_year, y = latitude, color = values, shape = type)) +
  geom_jitter(size = 2.5) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~zone) +
  theme_bw()

head(preds)


###


###     ###     ###     ###


###


### GAM V5 ----

#### Remove variables that are not necessary ----

dat <- df.nc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_MACPYRAD,
    den_PANINT,
    pct_cov_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    seaotter_dens_sm:seaotter_trend5yr,
    Days_16C, Days_18C:Degree_Days_16C, Degree_Days_18C:Degree_Days_23C,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    npp.mean
  )) %>%
  glimpse()


#### remove correlated substrate terrain variables ----

dat2 <- dat %>%
  # transform the variable
  mutate(mean_logvrm = log(mean_vrm)) %>%
  # remove untransformed 
  dplyr::select(-mean_vrm) %>%
  # remove unselected variables from corr analysis
  dplyr::select(-c(
    # all medians
    median_slope, median_depth, median_vrm, median_prob_of_rock,
    # all ranges
    range_slope, range_depth, range_vrm, range_prob_of_rock,
    # mins
    min_slope, min_vrm, 
    # sd
    sd_slope, sd_vrm,
    # max
    max_slope, max_vrm,
    # proportions mapped
    prop_map_prob_of_rock, prop_map_vrm, prop_map_slope,
    # max prob rock is always 1
    max_prob_of_rock,
    # days above 17C always 0
    #Days_17C, Degree_Days_17C,
    # the same as Max_Monthly_Anomaly_Temp
    #Max_Monthly_Anomaly_Summer_Temp,
    # the same as Min_Monthly_Anomaly_Temp
    #Min_Monthly_Anomaly_Upwelling_Temp,
    # MHW Days are the same as Intensity
    #MHW_Intensity, MHW_Summer_Intensity, MHW_Upwelling_Intensity,
    # remove other correlated vars
    mei_mean, pdo_mean, 
    #Max_Monthly_Anomaly_Temp,
    sd_depth, sd_prob_of_rock
  )) %>% # 
  glimpse()


## 

#### Select predictors for this GAMs ----

names(dat2)

dat.gam1 <- dat2 %>%
  dplyr::select(# get factors 
    "survey_year", "site_campus_unique_ID", "zone", 
    # get dependent var
    "den_NERLUE",
    # choose bio variables
    "den_MESFRAAD", "den_STRPURAD", "den_PYCHEL",
    # choose substrate vars
    "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
    "min_prob_of_rock" , "max_depth", "depth_mean",
    # choose other env variables
    "npgo_mean", "wh.95", "wh.max", 
    # choose temperature vars
    "Max_Monthly_Anomaly_Summer_Temp", 
    "Max_Monthly_Temp",
    "Days_15C", 
    "MHW_Summer_Days",
    "Min_Monthly_Anomaly_Summer_Temp", 
    #"Max_Monthly_Anomaly_Temp",
    "Mean_Monthly_Temp", 
    #"MHW_Days", "MHW_Upwelling_Days",
    #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
    #"Mean_Monthly_Upwelling_Temp",
    # choose nitrate vars
    "Days_8N", 
    "Max_Monthly_Nitrate", 
    "Min_Monthly_Anomaly_Summer_Nitrate",
    "Max_Monthly_Anomaly_Summer_Nitrate",
    #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
    "Mean_Monthly_Nitrate" 
    #"Max_Monthly_Anomaly_Nitrate"
  ) %>%
  mutate_at(vars(survey_year, zone), list(as.factor)) %>%
  glimpse() # Rows: 478 - 27 columns


#### Transform vars - Remove NAs ----

dat.gam1 <- dat.gam1 %>%
  # tranform biological variables
  mutate(logden_NERLUE = log(den_NERLUE + 1),
         logden_MESFRAAD = log(den_MESFRAAD + 1),
         logden_STRPURAD = log(den_STRPURAD + 1),
         logden_PYCHEL = log(den_PYCHEL)) %>% #glimpse()
  dplyr::select(-c(den_NERLUE, den_MESFRAAD, den_STRPURAD, den_PYCHEL)) %>%
  drop_na() %>%
  droplevels() %>%
  glimpse() # Rows: 460

any(is.na(dat.gam1))



#### V1. Run the full subset model selection ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

name <- 'V5'
o.dir <- paste("outputs_gam", name, sep = '_')

pred.vars <- c(# choose substrate vars
  "mean_slope" , "mean_prob_of_rock" , "min_depth" ,  "mean_logvrm",                         
  "min_prob_of_rock" , "max_depth", "depth_mean",
  # choose other env variables
  "npgo_mean", "wh.95", "wh.max", 
  # choose temperature vars
  "Max_Monthly_Anomaly_Summer_Temp", 
  "Max_Monthly_Temp",
  "Days_15C", 
  "MHW_Summer_Days",
  "Min_Monthly_Anomaly_Summer_Temp", 
  #"Max_Monthly_Anomaly_Temp",
  "Mean_Monthly_Temp", 
  #"MHW_Days", "MHW_Upwelling_Days",
  #"Min_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp", 
  #"Mean_Monthly_Upwelling_Temp",
  # choose nitrate vars
  "Days_8N",
  "Max_Monthly_Nitrate", 
  #"Max_Monthly_Anomaly_Upwelling_Nitrate", 
  #"Mean_Monthly_Nitrate", "Max_Monthly_Anomaly_Nitrate",
  "Min_Monthly_Anomaly_Summer_Nitrate", "Max_Monthly_Anomaly_Summer_Nitrate",
  # choose bio variables
  "logden_MESFRAAD", "logden_STRPURAD", "logden_PYCHEL")

length(pred.vars) # 23

#fact.vars <- c("survey_year")

model.v1 <- gam(logden_NERLUE ~ 
                  #s(site_campus_unique_ID, bs = 're') +
                  #s(transect_unique, bs = 're') +
                  s(site_campus_unique_ID, zone, bs = 're') +
                  s(survey_year, bs = 're') ,
                #te(depth_mean, wh.max,  bs = 'cr') +
                #te(depth_mean, wh.95,  bs = 'cr'),
                #s(mean_slope) + 
                #s(mean_prob_of_rock) + s(depth_mean) + s(mean_logvrm) +
                #s(min_depth) + s(min_prob_of_rock) +
                #s(max_depth) 
                #s(npgo_mean) + s(npp.mean) + 
                #s(orb.95) + s(orb.max) + 
                #s(wh.95) + s(wh.max) +
                #s(Max_Monthly_Anomaly_Upwelling_Temp) + 
                #s(Max_Monthly_Temp) + 
                #s(Mean_Monthly_Summer_Temp) + 
                #s(Mean_Monthly_Temp) + 
                #s(Mean_Monthly_Upwelling_Temp) + 
                #s(Min_Monthly_Anomaly_Summer_Temp) +
                #s(Min_Monthly_Anomaly_Temp) + 
                #s(Min_Monthly_Temp) + 
                #s(logMHW_Days) + 
                #s(logMHW_Summer_Days) + 
                #s(logMHW_Upwelling_Days) +
                #s(logden_HALRUF) + 
                #s(logden_PYCHEL) + 
                #s(transect_unique, bs = 're'),
                data = dat.gam1, 
                family = tw(),
                method = "REML") 

model.set <- generate.model.set(use.dat = dat.gam1,
                                test.fit = model.v1,
                                pred.vars.cont = pred.vars,
                                #smooth.smooth.interactions = c(c("depth_mean", "wh.max"), c("depth_mean", "wh.95")),
                                max.predictors = 6,
                                cov.cutoff = 0.6, # cut off for correlations
                                #pred.vars.fact = fact.vars,
                                #linear.vars="Distance",
                                k=3,
                                null.terms = "s(site_campus_unique_ID, zone, bs = 're') +
                                s(survey_year, bs = 're')")


out.list <- fit.model.set(model.set,
                          max.models= 2000,
                          parallel=T)

beep()

# check results ----
out.all=list()
var.imp=list()

names(out.list)
length(out.list$success.models)
out.list$failed.models
out.list$mod.data.out[,1]
length(out.list$mod.data.out)
mod.table <- out.list$mod.data.out
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=3),]
nrow(out.i)
length(out.i)
names(out.i)
out.i[1]
out.i[2]
class(out.i)

out.all <- c(out.all,list(out.i))
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

out.list$failed.models
out.list$success.models[7]

length(out.list$failed.models)
length(out.list$success.models)


# plot the best models ----
# for(m in 1:nrow(out.i)){
#   best.model.name=as.character(out.i$modname[1])
#   
#   png(file=paste(o.dir, paste(name,m,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
#   if(best.model.name!="null"){
#     par(mfrow=c(3,1),mar=c(9,4,3,1))
#     best.model <- out.list$success.models[[best.model.name]]
#     plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
#     mtext(side=2,text='logden_NERLUE',outer=F)}  
#   dev.off()
# }


# Model fits and importance ----
names(out.all) <- 'logden_NERLUE'
names(var.imp) <- 'logden_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name,"best_models.csv",sep="_"), sep ='/'))
write.csv(all.var.imp, file=paste(o.dir, paste(name,"all.var.imp.csv",sep="_"), sep ='/'))

# Generic importance plots-
# heatmap.2(all.var.imp, notecex=0.4,  
#           #dendrogram ="none",
#           #col = colorRampPalette(c("white","yellow","red"))(10),
#           trace = "none", key.title = "", keysize=2,
#           notecol="black",key=T,
#           sepcolor = "black", margins=c(12,8), lhei=c(4,15),
#           Rowv=FALSE,Colv=FALSE)


# Part 2 - custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep = '_'), sep ='/')) %>% #from local copy
  rename(resp.var=X)%>%
  dplyr::select(-resp.var) %>% #glimpse()
  gather(key = predictor, value = importance) %>%
  arrange(importance) %>%
  mutate_at(vars(predictor), list(as.factor)) %>%
  glimpse()

## set the levels in order we want

dat.taxa %>%
  #count(predictor) %>% glimpse()
  mutate(predictor = fct_reorder(predictor, importance, .desc = TRUE)) %>%
  ggplot(aes(x = predictor, y = importance)) + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


# Manually make the most parsimonious GAM models ----

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()

###


#### gam2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###


#### gam3 ----

subname <- "3"

best.model.name=as.character(out.i$modname[3])
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###

#### gam4 ----

subname <- "4"

best.model.name=as.character(out.i$modname[4])
best.model <- out.list$success.models[[best.model.name]]
best.model


# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


# predict - gam5 ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(depth_mean=mean(mod$model$depth_mean),
                        Max_Monthly_Temp=mean(mod$model$Max_Monthly_Temp),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        MHW_Summer_Days=mean(mod$model$MHW_Summer_Days),
                        Min_Monthly_Anomaly_Summer_Nitrate=mean(mod$model$Min_Monthly_Anomaly_Summer_Nitrate),
                        Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
                        site_campus_unique_ID=(mod$model$site_campus_unique_ID),
                        zone=(mod$model$zone),
                        survey_year=(mod$model$survey_year))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.year = testdata%>%data.frame(fits)%>%
  group_by(survey_year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=survey_year,y=response,fill=survey_year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
                    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year 


###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- dat.gam1 %>%
  dplyr::select(survey_year, site_campus_unique_ID, zone, logden_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  left_join(ner.obs, by = c('survey_year', 'site_campus_unique_ID', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = logden_NERLUE, y = fit, color = zone)) +
  geom_point() 



# add lat lons
pred.obs.all <- pred.obs.all %>%
  left_join(nclatlons, by = 'site_campus_unique_ID') %>%
  pivot_longer(cols = c(fit, logden_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)

# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all, aes(x = zone, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

ggplot(pred.obs.all, aes(x = survey_year, y = latitude, color = values, shape = type)) +
  geom_jitter(size = 2.5) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~zone) +
  theme_bw()

head(preds)



##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


###


#### gam5 ----

subname <- "5"

best.model.name=as.character(out.i$modname[5])
best.model <- out.list$success.models[[best.model.name]]
best.model
class(best.model)


# get formula of best model --
bm <- print(best.model)
bm.form <- as.formula(paste("logden_NERLUE ", paste(bm, collapse= " ")))
bm.form
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = dat.gam1, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

##

# Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
par(mfrow=c(3,3),mar=c(9,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()


# predict - gam5 ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(depth_mean=mean(mod$model$depth_mean),
                        Max_Monthly_Temp=mean(mod$model$Max_Monthly_Temp),
                        mean_logvrm=mean(mod$model$mean_logvrm),
                        MHW_Summer_Days=mean(mod$model$MHW_Summer_Days),
                        Min_Monthly_Anomaly_Summer_Nitrate=mean(mod$model$Min_Monthly_Anomaly_Summer_Nitrate),
                        Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
                        site_campus_unique_ID=(mod$model$site_campus_unique_ID),
                        zone=(mod$model$zone),
                        survey_year=(mod$model$survey_year))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.year = testdata%>%data.frame(fits)%>%
  group_by(survey_year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=survey_year,y=response,fill=survey_year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
                    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year 


###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- dat.gam1 %>%
  dplyr::select(survey_year, site_campus_unique_ID, zone, logden_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  left_join(ner.obs, by = c('survey_year', 'site_campus_unique_ID', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = logden_NERLUE, y = fit, color = zone)) +
  geom_point() 



# add lat lons
pred.obs.all <- pred.obs.all %>%
  left_join(nclatlons, by = 'site_campus_unique_ID') %>%
  pivot_longer(cols = c(fit, logden_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)

# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all, aes(x = zone, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

ggplot(pred.obs.all, aes(x = survey_year, y = latitude, color = values, shape = type)) +
  geom_jitter(size = 2.5) +
  labs(x = 'Transect', y = 'Latitude') +
  facet_wrap(~zone) +
  theme_bw()

head(preds)

