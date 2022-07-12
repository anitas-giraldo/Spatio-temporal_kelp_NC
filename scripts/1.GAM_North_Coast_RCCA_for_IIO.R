###
### Script based on one created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 19 April 2022

## This script prepares the data for the density models of kelp in the north coast with RCCA data --

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
library(beepr)
library(stringr)

# Clear environment ----
rm(list=ls())

### Set directories ----
#w.dir<- dirname(rstudioapi::getActiveDocumentContext()$path)
m.dir <- here()
mlpa.dir <- "G:/Shared drives/MLPA L-T Kelp Forest Monitoring Project - Admin/Monitoring Data/MLPA merged datasets"
merged.data <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Merge_bio_env_vars_transect/data_tidy"
d.dir <- here("data")
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"
raw.dir <- here("data_raw")
plots.dir <- here("plots")
o.dir <- here("For_IIO")
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"

## Load info on years RCCA ----
years <- read.csv(paste(rcca.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


# get the sites from with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 64


## Load RCCA data ----

df <- read.csv(paste(dd.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse()


## get the sites for North Coast model ----
df.nc <- df %>%
  dplyr::select(-c(latitude, longitude)) %>%
  right_join(ncsites, by = c('site_name')) %>%
  droplevels() %>% #glimpse()
  dplyr::select(-c(total.years, pre.mhw.years, during.mhw.years, post.mhw.years)) %>%
  relocate(c(latitude, longitude), .after = zone) %>%
  glimpse() # Rows: 708

length(levels(df.nc$site_name)) # 10
levels(df.nc$site_name)
any(is.na(df.nc$Max_Monthly_Anomaly_Temp))

## Choose variables and transform needed ----
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
    #Mean_Monthly_Upwelling_Nitrate,
    Max_Monthly_Anomaly_Nitrate,
    #Mean_Monthly_Summer_Nitrate,
    # Temperature vars
    Days_16C ,
    Mean_Monthly_Temp ,
    #Mean_Monthly_Summer_Temp,
    MHW_Upwelling_Days  , 
    #Min_Monthly_Anomaly_Temp,
    Max_Monthly_Anomaly_Upwelling_Temp,
    #Min_Monthly_Temp, 
    Mean_Monthly_Upwelling_Temp,
    # Other vars
    npp.mean,                             
    #wh.95 ,   wh.max,
    wh_max, mean_waveyear,
    #npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
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
  # Other Env transformations
  mutate(log_npp.mean = log(npp.mean + 1)) %>%
  dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708


# Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 507

## Divide Pre and Post ----
levels(dat2$year)

# Pre-MHW --
PreMHW <- dat2 %>%
  dplyr::filter(year == "2006" |
                  year == "2007" |
                  year == "2008" | 
                  year == "2009" |
                  year == "2010" |
                  year == "2011" |
                  year == "2012" |
                  year == "2013" ) %>%
  dplyr::filter(year != "2006",
                year != "2007") %>%
  droplevels() %>%
  glimpse() # Rows: 282

levels(PreMHW$year)
length(levels(PreMHW$year)) # 8
nrow(PreMHW)

# Pre-MHW --
PostMHW <- dat2 %>%
  dplyr::filter(year == "2014" |
                  year == "2015" |
                  year == "2016" | 
                  year == "2017" |
                  year == "2018" |
                  year == "2019" |
                  year == "2020" |
                  year == "2021" ) %>%
  droplevels() %>%
  glimpse() # Rows: 225

levels(PostMHW$year)
length(levels(PostMHW$year)) # 7


###

###


### GAM V1 ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

# Use Pre MHW data to fit the model --

# 1. Select predictors for this GAM ----
names(PreMHW)


# 2. Divide data into train and test ----

inTraining <- createDataPartition(PreMHW$log_den_NERLUE, p = 0.8, list = FALSE)
train.gam <- PreMHW[ inTraining,]
test.gam  <- PreMHW[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'V2_premhw2'
o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c("Days_10N",
               "Min_Monthly_Nitrate",
               "Max_Monthly_Nitrate",
               "Mean_Monthly_Nitrate",
               #"Mean_Monthly_Upwelling_Nitrate",
               "Max_Monthly_Anomaly_Nitrate" ,
               #"Mean_Monthly_Summer_Nitrate"  ,      
               "Mean_Monthly_Temp"  ,
               #"Mean_Monthly_Summer_Temp" ,
               "MHW_Upwelling_Days"  ,               
               #"Min_Monthly_Anomaly_Temp"   ,        
               "Max_Monthly_Anomaly_Upwelling_Temp",
               #"Min_Monthly_Temp"  ,                 
               "Mean_Monthly_Upwelling_Temp",        
               #"wh.95",                             
               #"wh.max",                             
               #"npgo_mean",
               #"mei_mean",                    
               "log_den_MESFRAAD" ,
               "log_den_STRPURAD" ,                 
               "log_den_PYCHEL" ,                   
               #"log_den_HALRUF",
               "log_Days_16C"                      
               #"log_npp.mean" 
               )

length(pred.vars) # 24



# 5. Define Null model ----

#fact.vars <- c("survey_year")

model.v1 <- gam(log_den_NERLUE ~ 
                  s(site_name, zone, bs = 're') +
                  s(year, bs = 're') ,
                data = train.gam, 
                family = tw(),
                method = "REML") 

# 6. Define model set up ----

model.set <- generate.model.set(use.dat = train.gam,
                                test.fit = model.v1,
                                pred.vars.cont = pred.vars,
                                #smooth.smooth.interactions = c("depth_mean", "wh.max", "wh.95"),
                                max.predictors = 6,
                                cov.cutoff = 0.7, # cut off for correlations
                                #pred.vars.fact = fact.vars,
                                #linear.vars="Distance",
                                k=4,
                                null.terms = "s(site_name, zone, bs = 're') +
                                s(year, bs = 're')")


# 7. Run the full subset model selection ---- 

out.list <- fit.model.set(model.set,
                          max.models= 500,
                          parallel=T)

beep()

# 8. Model fits and importance ----
out.all=list()
var.imp=list()

#names(out.list)
# length(out.list$success.models)
# out.list$failed.models
# out.list$mod.data.out[,1]
# length(out.list$mod.data.out)

# Put results in a table --
mod.table <- out.list$mod.data.out
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=3),]
nrow(out.i)

out.all <- c(out.all,list(out.i)) 
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

# 9. Save model fits and importance ----
names(out.all) <- 'Nereo.density'
names(var.imp) <- 'Nereo.density'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "all.var.imp.csv", sep="/"))


# 10. plot the best models ----
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[i])
  
  png(file=paste(o2.dir, paste(name,m,'Nereo.density',"mod_fits.png",sep="_"), sep ='/'))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model <- out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text='log_den_NERLUE',outer=F)}
  dev.off()
}



# 11. custom plot of importance scores----
dat.taxa <- read.csv(paste(o2.dir, "all.var.imp.csv", sep ='/')) %>% #from local copy
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

namep <- paste(name,"var_imp.png", sep ='_')

ggsave(namep, device = 'png', path = o2.dir)


 # Manually make the most parsimonious GAM models ----

out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
  glimpse()

#### gam1 ----

subname <- "8"

best.model.name=as.character(out.i$modname[8])
best.model.name
best.model <- out.list$success.models[[best.model.name]]
best.model <- print(out.i$formula[1])

gam1 <- gam(formula = log_den_NERLUE ~ #s(Days_10N, k = 6, bs = "cr") + 
              s(Max_Monthly_Anomaly_Upwelling_Temp, k = 3, bs = "cr") + 
              s(Max_Monthly_Anomaly_Nitrate, k = 6, bs = "cr") + 
              s(Min_Monthly_Nitrate, k = 6, bs = "cr") +
              #s(Min_Monthly_Nitrate, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

# png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
# par(mfrow=c(3,3),mar=c(2,4,3,1))
# visreg(gam1)
# dev.off()

##

# Plot model results ----

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o2.dir, paste(name,subname,'logden_NERLUE',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()

##

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(Max_Monthly_Anomaly_Upwelling_Temp=mean(mod$model$Max_Monthly_Anomaly_Upwelling_Temp),
                        Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        #Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE,
    Max_Monthly_Anomaly_Upwelling_Temp,
                Max_Monthly_Anomaly_Nitrate, Min_Monthly_Nitrate,
                #Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)


predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
                    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  scale_fill_viridis(direction = -1, discrete = T) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = log_den_NERLUE, y = fit, color = zone)) +
  geom_point() 


# Plot observed vs. predicted ----
library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = log_den_NERLUE, col = zone)) +        
  geom_point() +
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
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o2.dir)


# ner.obs <- dat.gam1 %>%
#   dplyr::select(survey_year, site_campus_unique_ID, zone, log_den_MACPYRAD) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('survey_year', 'site_campus_unique_ID', 'zone')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))
# 
# # plot pred vs. observed ----
# pred.obs.all %>%
#   ggplot(aes(x = log_den_MACPYRAD, y = fit, color = zone)) +
#   geom_point() 
# 
# pred.obs.all %>%
#   ggplot(aes(x = log_den_MACPYRAD, y = fit, color = zone)) +
#   geom_point()



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
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)



###

## PREDICT ON POSTHEATWAVE ----

glimpse(PostMHW)


## PREDict ----
head(mod$model)
testdata <- PostMHW %>%
  dplyr::select(log_den_NERLUE,
                Max_Monthly_Anomaly_Upwelling_Temp, 
                Max_Monthly_Anomaly_Nitrate, 
                Min_Monthly_Nitrate, 
                #Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  scale_fill_viridis(direction = -1, discrete = T) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- PostMHW %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = log_den_NERLUE, y = fit, color = zone)) +
  geom_point() 


# Plot observed vs. predicted ----
library(ggpmisc)

glimpse(predicts.all)

predicts.all2 <- predicts.all %>%
  #dplyr::filter(zone != 'OUTER') %>%
  dplyr::filter(fit < 10) %>%
  droplevels() %>%
  glimpse()

p <- ggplot(predicts.all2, aes(x = fit, y = log_den_NERLUE, col = zone)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all2, aes(x = fit, y = log_den_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_postMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)






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
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

namep <- paste(name, subname, "pred-obs_latitude_postMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)








####

####

## PREDICT NEREOCYSTIS ----

## get model again ----

gam1 <- gam(formula = log_den_NERLUE ~ 
              s(Max_Monthly_Anomaly_Upwelling_Temp, k = 3, bs = "cr") + 
              s(Max_Monthly_Anomaly_Nitrate, k = 6, bs = "cr") + 
              s(Min_Monthly_Nitrate, k = 6, bs = "cr") + 
              s(mean_depth, k = 3, bs = "cr") +
              #s(MHW_Upwelling_Days, k = 7, bs = "cr") +
              #s(wh.95, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)

visreg(gam1)


nereo.mod <- gam1
nereo.mod$formula



# 1. Get rasters of predictors for every year ----

# set directory of rasters --

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data"

# check variables ---
nereo.mod$formula

# log_den_STRPURAD 
# Min_Monthly_Temp
# Mean_Monthly_Summer_Nitrate 
# site_name + zone + year


# load susbtrate predictors ----

# set extent
e.nc <- extent(-124.60, -120.6956, 37.62184, 42.00996)
# bathy
sp.dir <- here('spatial')
nc.mask <- rast(paste(sp.dir, "mask_NC.tif", sep ='/'))

sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

depth <- rast(paste(sub.dir, "depth_mean_nc_all_wInterp_300m_30m.tif", sep ='/'))
depth

crs1 <- "epsg:4326"
depth2 <- project(depth, crs1)
plot(depth2)




# load temperature predictors ----


Max_An_UP <- rast(paste(re.dir, 'Temperature', "Max_Monthly_Anomaly_Upwelling_Temp.tif", sep ='/'))
Max_An_UP

Max_An_UP2 <- crop(Max_An_UP, e.nc)


# resample predictors to bathy ----
Max_An_UP3 <- resample(Max_An_UP, depth2)

# mask predictors to bathy ----
Max_An_UP4 <- mask(Max_An_UP3, depth2)
plot(Max_An_UP4[[1]])



###


# load nitrate predictors ----

## Max_Monthly_Anomaly_Nitrate --

max_nit <- rast(paste(re.dir, "Nitrate", "Max_Monthly_Anomaly_Nitrate.tif", sep ='/'))
max_nit

# # crop to NC --
# mean_sum_nit2 <- crop(mean_sum_nit, extent(e.nc))
# plot(mean_sum_nit2[[1]])

# resample predictors to bathy ----
max_nit3 <- resample(max_nit, depth2)

# mask predictors to bathy ----
max_nit4 <- mask(max_nit3, depth2)
plot(max_nit4)


##

# load nitrate predictors ----

min_nit <- rast(paste(re.dir, "Nitrate", "Min_Monthly_Nitrate.tif", sep ='/'))
min_nit

# # crop to NC --
# mean_sum_nit2 <- crop(mean_sum_nit, extent(e.nc))
# plot(mean_sum_nit2[[1]])

# resample predictors to bathy ----
min_nit3 <- resample(min_nit, depth2)

# mask predictors to bathy ----
min_nit4 <- mask(min_nit3, depth2)
plot(min_nit4[[4]])



##



# stack and try to predict one year ----

preds1 <- c(depth2, Max_An_UP4[[1]], max_nit4[[1]], min_nit4[[1]])
names(preds1)

year1998 <- classify(depth2[[1]], cbind(-Inf, 0.1, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("mean_depth",  "Max_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Nitrate"  ,
                   "Min_Monthly_Nitrate",                          
                   'year') 


##


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
names(preds3) <- c("mean_depth",  "Max_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Nitrate"  ,
                   "Min_Monthly_Nitrate",                          
                   'year'     ,
                   "site_name") 

##

# zone ----

zone.raster <- depth2
names(zone.raster) <- 'zone'

plot(zone.raster)

m <- c(-Inf, -10, 2,
       -10, 0, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

zone.raster2 <- classify(zone.raster, rclmat, include.lowest=TRUE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("mean_depth",  "Max_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Nitrate"  ,
                   "Min_Monthly_Nitrate",                          
                   'year'     ,
                   "site_name",
                   "zone")


##



### LOOP ----

## loop to predict several years ----

# make list of years --
year.list <- paste(1998:2021)
length(year.list)

# make template raster of year ----

year.raster <- classify(depth2[[1]], cbind(-Inf, 0.1, 1998), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 

plot(zone.raster2)

# make sites raster ----

plot(site.raster2)

# outputs dir ----

preds.dir <- paste(o2.dir,  "sp_preds", sep ='/')


for (i in 1:length(year.list)) {
  
  
  # 2. stack with predictors for that year
  env.raster <- c(depth2, Max_An_UP4[[i]], max_nit4[[i]], min_nit4[[i]])
  
  preds1 <- env.raster
  
  # 3. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster2)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("mean_depth",
                     "Max_Monthly_Anomaly_Upwelling_Temp",  
                     "Max_Monthly_Anomaly_Nitrate", 
                     "Min_Monthly_Nitrate"  , 
                     #"Max_Monthly_Nitrate",                          
                     "year"     ,
                     "zone",
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, nereo.mod)
  plot(year.prediction)
  
  # 6. Change the log
  year.prediction2 <- exp(year.prediction)
  plot(year.prediction2)
  
  # save
  name.raster <- paste(year.no, "logNereo_preds_NC.tif", sep = '_')
  #writeRaster(year.prediction, paste(preds.dir, "sp_predictions", paste(year.no, "logNereo_preds_NC.tif", sep = '_'), sep ='/'))
  writeRaster(year.prediction2, paste(preds.dir, name.raster, sep = '/'))
}




###

###

### GET STABILITY ----

## UP TO HERE ----

# load raster data --

n.files <- dir(preds.dir)
# list files in source --
n.files <- list.files(preds.dir, pattern = '.tif', full.names = TRUE)
n.files
length(u.files)
# list names to load onto the Environment --
names.list <- list.files(preds.dir, pattern = '.tif')
names.list <- str_replace(names.list, ".tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)

### Calculate mean across years ----

stability.stack <- c()

# calculate mean across years
mean.nereo <- mean(preds.stack, na.rm = T)
plot(mean.nereo)


# calculate sd across years
sd.nereo <- app(preds.stack, fun = sd, na.rm = T)
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
writeRaster(stability.stack, paste(preds.dir, "Stabilty_calculations_Nereo_NC_gam_V3.tif", sep ='/'), overwrite = T)


###

###

# save points ----

s_index.n.points <- as.data.frame(s_index.nereo, xy = T)
head(s_index.n.points)

write.csv(s_index.n.points, paste(preds.dir, "S_index_points_Nereo_NC_gam_pre.csv", sep ='/'))



####

####

## TEST 2 ----

#### gam1 ----

# try same vars as GAM v4
gam2 <- gam(formula = log_den_NERLUE ~ #s(Days_10N, k = 6, bs = "cr") + 
              s(log_den_STRPURAD, k = 4, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 6, bs = "cr") + 
              s(mean_waveyear, k = 6, bs = "cr") +
              s(wh_max, k = 6, bs = "cr") +
              s(mean_depth, k = 4, bs = "cr") +
              #s(Min_Monthly_Nitrate, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = PreMHW, method = "REML")



gam2$aic
gam2$deviance
summary(gam2)
gam.check(gam2)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2)
dev.off()

# png(file=paste(o.dir, paste(name,subname,'logden_NERLUE',"mod_fits.png",sep="_"), sep ='/'))
# par(mfrow=c(3,3),mar=c(2,4,3,1))
# visreg(gam1)
# dev.off()

##

# Plot model results ----

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o2.dir, paste(name,subname,'logden_NERLUE',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='logden_NERLUE',outer=F)
dev.off()

##

## PREDICT  gam ----
mod<-gam2
head(mod$model)
testdata <- expand.grid(Max_Monthly_Anomaly_Upwelling_Temp=mean(mod$model$Max_Monthly_Anomaly_Upwelling_Temp),
                        Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        #Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE,
                log_den_STRPURAD,
                Max_Monthly_Nitrate,
                mean_waveyear, wh_max,
                mean_depth,
                #Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)


predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  scale_fill_viridis(direction = -1, discrete = T) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = log_den_NERLUE, y = fit, color = zone)) +
  geom_point() 


# Plot observed vs. predicted ----
library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = log_den_NERLUE, col = zone)) +        
  geom_point() +
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
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o2.dir)


# ner.obs <- dat.gam1 %>%
#   dplyr::select(survey_year, site_campus_unique_ID, zone, log_den_MACPYRAD) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('survey_year', 'site_campus_unique_ID', 'zone')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))
# 
# # plot pred vs. observed ----
# pred.obs.all %>%
#   ggplot(aes(x = log_den_MACPYRAD, y = fit, color = zone)) +
#   geom_point() 
# 
# pred.obs.all %>%
#   ggplot(aes(x = log_den_MACPYRAD, y = fit, color = zone)) +
#   geom_point()



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
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)



###

## PREDICT ON POSTHEATWAVE ----

glimpse(PostMHW)


## PREDict ----
head(mod$model)
testdata <- PostMHW %>%
  dplyr::select(log_den_NERLUE,
                log_den_STRPURAD,
                Max_Monthly_Nitrate,
                mean_waveyear, wh_max,
                mean_depth,
                #Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()


# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
    values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  scale_fill_viridis(direction = -1, discrete = T) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- PostMHW %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = log_den_NERLUE, y = fit, color = zone)) +
  geom_point() 


# Plot observed vs. predicted ----
library(ggpmisc)

glimpse(predicts.all)

predicts.all2 <- predicts.all %>%
  #dplyr::filter(zone != 'OUTER') %>%
  dplyr::filter(fit < 10) %>%
  droplevels() %>%
  glimpse()

p <- ggplot(predicts.all2, aes(x = fit, y = log_den_NERLUE, col = zone)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all2, aes(x = fit, y = log_den_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_postMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)






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
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  facet_wrap(~type) +
  theme_bw()

namep <- paste(name, subname, "pred-obs_latitude_postMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)








####

####

## PREDICT NEREOCYSTIS ----

## get model again ----

gam1 <- gam2



gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)

visreg(gam1)


nereo.mod <- gam1
nereo.mod$formula



### Get depth ----

depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

names(dat3)

depth <- rast(paste(depth.dir, "depth_mean_nc_all_wInterp_300m_30m.tif", sep ='/'))
plot(depth)

n.extent <- ext(depth)

crs1 <- "epsg:4326"

d2 <- project(depth, crs1)

n.extent <- ext(d2)


# env directory
re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data"



###


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

## WH max ----

# load raster data --

w.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Wave_data"

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



## Mean Wave year --

wave.dir <- paste(w.dir, "Wave_extended2", "wy_mean", sep ='/')

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
wymean.stack <- c()

# use do call to create a raster otherwise it creates a list
wymean.stack <- do.call("c", tfiles)
plot(wymean.stack[[1]])


# # crop to NC --
wymean.stack2 <- crop(wymean.stack, n.extent)
plot(wymean.stack2[[1]])

# resample predictors to bathy ----
wymean.stack3 <- resample(wymean.stack2, d2)

# mask predictors to bathy ----
wymean.stack4 <- mask(wymean.stack3, d2)
plot(wymean.stack4[[1]])





### Get urchins ----

# load purple urchin predictions ----
urch.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Spatio_temporal_GAMs/outputs_nc_rcca/gam_urchins2/sp_predictions"

# load raster data --

u.files <- dir(urch.dir)
u.files <- list.files(urch.dir, pattern = '.tif')
u.files
length(u.files)


##



# stack and try to predict one year ----

preds1 <- c(d2,  max_nit4[[1]], whmax.stack4[[1]], wymean.stack4[[1]])
names(preds1)

year1998 <- classify(max_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("mean_depth",  "Max_Monthly_Nitrate"  , "wh_max", "mean_waveyear",                         
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
names(preds3) <- c("mean_depth", "Max_Monthly_Nitrate"  , "wh_max", "mean_waveyear",                         
                   "year"     ,
                   "site_name") 


##

# zone ----

zone.raster <- d2
names(zone.raster) <- 'zone'
plot(zone.raster)

rec.m <-  c(-Inf, -10, 2,
            -10, 0.1, 1)

rclmat <- matrix(rec.m, ncol=3, byrow=TRUE)

zone.raster2 <- classify(zone.raster, rclmat, right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("mean_depth", "Max_Monthly_Nitrate"  , "wh_max", "mean_waveyear",                         
                   "year"      ,
                   "site_name",
                   "zone")


##



### LOOP ----

## loop to predict several years ----

nereo.mod <- gam_pre_V2


# make list of years --
#year.list <- paste(1998:2021)
year.list <- paste(2004:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(d2, cbind(-Inf, 0, 2004), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 

zone.raster <- d2
names(zone.raster) <- 'zone'
plot(zone.raster)

rec.m <-  c(-Inf, -10, 2,
            -10, 0.1, 1)

rclmat <- matrix(rec.m, ncol=3, byrow=TRUE)

zone.raster2 <- classify(zone.raster, rclmat, right=FALSE)
plot(zone.raster2)
names(zone.raster2) <- 'zone'

# make sites raster ----

# each site is latitude --
rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)
# make raster
site.raster <- rast(rdf, type = 'xyz', crs="EPSG:4326", extent = ext(year1998))
site.raster
# extend to fit other rasters
site.raster2 <- extend(site.raster, year1998)

# outputs dir ----

preds.dir <- paste(o2.dir, "sp_preds2", sep ='/')


for (i in 1:length(year.list)) {
  
  # 1. get urchins
  urchin.rast <- rast(paste(urch.dir, u.files[6+i], sep ='/'))
  urchin.rast2 <- resample(urchin.rast, d2)
  
  # 2. stack with predictors for that year
  env.raster <- c(d2,  max_nit4[[6+i]], whmax.stack4[[i]], wymean.stack4[[i]])
  
  preds1 <- c(urchin.rast2, env.raster)
  
  # 3. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(-Inf, 0, year.no), right=FALSE)
  
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("log_den_STRPURAD",
                     "mean_depth",  
                     #"Mean_Monthly_Upwelling_Temp", 
                     "Max_Monthly_Nitrate"  , 
                     "wh_max", 
                     "mean_waveyear",
                     "year"     ,
                     "zone",
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, nereo.mod)
  plot(year.prediction)
  
  # 6. Change the log
  year.prediction2 <- exp(year.prediction)
  plot(year.prediction2)
  
  # save
  name.raster <- paste(year.no, "Nereo_preds_NC_PreV2.tif", sep = '_')
  #writeRaster(year.prediction, paste(preds.dir, "sp_predictions", paste(year.no, "logNereo_preds_NC.tif", sep = '_'), sep ='/'))
  writeRaster(year.prediction2, paste(preds.dir, name.raster, sep = '/'))
}


###

###

### GET STABILITY ----


# load raster data --

n.files <- dir(preds.dir)
# list files in source --
n.files <- list.files(preds.dir, pattern = '.tif', full.names = TRUE)
n.files
length(u.files)
# list names to load onto the Environment --
names.list <- list.files(preds.dir, pattern = '.tif')
names.list <- str_replace(names.list, ".tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)

### Calculate mean across years ----

stability.stack <- c()

# calculate mean across years
mean.nereo <- mean(preds.stack, na.rm = T)
plot(mean.nereo)


# calculate sd across years
sd.nereo <- app(preds.stack, fun = sd, na.rm = T)
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
#writeRaster(stability.stack, paste(preds.dir, "Stabilty_calculations_Nereo_NC_pre_V2.tif", sep ='/'), overwrite = T)


# save points ----

s_index.n.points <- as.data.frame(s_index.nereo, xy = T)
head(s_index.n.points)

write.csv(s_index.n.points, paste(preds.dir, "S_index_points_Nereo_NC_gam_pre_V2.csv", sep ='/'))

#

s.n.points <- as.data.frame(s.nereo, xy = T)
head(s.n.points)

write.csv(s.n.points, paste(preds.dir, "S_points_Nereo_NC_gam_pre_V2.csv", sep ='/'))


#

mean.n.points <- as.data.frame(mean.nereo, xy = T)
head(mean.n.points)

write.csv(mean.n.points, paste(preds.dir, "mean_points_Nereo_NC_gam_pre_V2.csv", sep ='/'))
