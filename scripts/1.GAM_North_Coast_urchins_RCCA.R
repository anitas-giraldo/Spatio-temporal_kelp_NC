###
### Script based on one created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 11 July 2022 - with orb vel and Npp from Tom

## This script runs models

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
o.dir <- here("outputs_nc_rcca_urchins")
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"

#o2.dir <- paste(o.dir, "gam_urchins3", sep ='/')

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


## Load RCCA data ----

df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
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
    # Other vars
    #npp.mean,                             
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc,
    # Orb vel
    UBR_Mean,
    # NPP
    Mean_Monthly_NPP, Max_Monthly_NPP_Upwelling, Mean_Monthly_NPP_Upwelling, Min_Monthly_NPP) %>%
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
  mutate(log_UBR_Mean = log(UBR_Mean + 1)) %>%
  dplyr::select(-c(UBR_Mean)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
  glimpse() # Rows: 708


# Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)

# 1. Select predictors for this GAM ----
names(dat2)



# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat2$log_den_STRPURAD, p = 0.75, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

# URCHINS 4 ----

name <- 'urchins4'

o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c("Days_10N",
               "Min_Monthly_Nitrate",
               "Max_Monthly_Nitrate",
               "Mean_Monthly_Nitrate",
               "Mean_Monthly_Upwelling_Nitrate",
               "Max_Monthly_Anomaly_Nitrate" ,
               "Mean_Monthly_Summer_Nitrate"  ,      
               "Mean_Monthly_Temp"  ,
               "Mean_Monthly_Summer_Temp" ,
               "MHW_Upwelling_Days"  ,               
               "Min_Monthly_Anomaly_Temp"   ,        
               "Max_Monthly_Anomaly_Upwelling_Temp",
               "Min_Monthly_Temp"  ,                 
               "Mean_Monthly_Upwelling_Temp",        
               #"wh.95",                             
               #"wh.max",                             
               #"npgo_mean",
               #"mei_mean",                    
               #"log_den_MESFRAAD" ,
               #"log_den_STRPURAD" ,                 
               "log_den_PYCHEL" ,                   
               #"log_den_HALRUF",
               "log_Days_16C",                      
               #"log_npp.mean",
               "log_mean_vrm",
               "mean_slope",
               "mean_depth", 
               "mean_prob_of_rock",
               "wh_max",
               "wh_mean",
               "mean_waveyear",
               "wh_95prc",
               "log_UBR_Mean",                      
               "Mean_Monthly_NPP",
               "Max_Monthly_NPP_Upwelling",
               "log_Mean_Monthly_NPP_Upwelling",    
               "log_Min_Monthly_NPP")

length(pred.vars) # 29

# 5. Define Null model ----

#fact.vars <- c("survey_year")

model.v1 <- gam(log_den_STRPURAD ~ 
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
names(out.all) <- 'log_den_STRPURAD'
names(var.imp) <- 'log_den_STRPURAD'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "STRPURAD_all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "STRPURAD_best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "STRPURAD_all.var.imp.csv", sep="/"))


# 10. plot the best models ----
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(o2.dir, paste(name,m,'Nereo.density',"mod_fits.png",sep="_"), sep ='/'))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model <- out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text='log_den_NERLUE',outer=F)}
  dev.off()
}



# 11. custom plot of importance scores----
dat.taxa <- read.csv(paste(o2.dir, "STRPURAD_all.var.imp.csv", sep ='/')) %>% #from local copy
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



### Check BEST MODELS ----

# Manually make the most parsimonious GAM models ----

out.i <- read.csv(paste(o2.dir, "STRPURAD_best_models.csv", sep ='/')) %>%
  glimpse()

#### gam3 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model.name
#best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
#best.model


gam1 <- gam(formula = log_den_STRPURAD ~ 
              s(Days_10N, k = 6, bs = "cr") +
              #s(log_Days_16C, k = 6, bs = "cr") +
              #log_mean_vrm +
              #s(log_Mean_Monthly_NPP_Upwelling, k = 3, bs = "cr") + 
              s(log_Min_Monthly_NPP, k = 6, bs = "cr") + 
              mean_depth +
              s(Max_Monthly_Anomaly_Upwelling_Temp, k = 6, bs = "cr") +
              mean_prob_of_rock +
              #s(wh_95prc, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = dat2, method = "GCV")



gam1$aic
gam1$gcv.ubre
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
pdf(file=paste(o2.dir, paste(name,subname,'log_den_STRPURAD',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='log_den_STRPURAD',outer=F)
dev.off()

##

## PREDICT  gam ----
mod<-gam1

glimpse(test.gam)


testdata <- dat2 %>%
  dplyr::select(log_den_STRPURAD, 
                log_Min_Monthly_NPP,
                Days_10N, 
                #wh_95prc,
                mean_depth, 
                Max_Monthly_Anomaly_Upwelling_Temp,
                mean_prob_of_rock, 
                #log_mean_vrm,
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
  ylab("Log S. purpuratus")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"),limits = rev(levels(predicts.year$survey_year)))+
  #scale_x_discrete(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020")) +
  geom_bar(stat = "identity")+
  scale_fill_viridis(discrete=T, option = 'E', direction = -1) +
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year 


### Plot mean yearly kelp predicted vs observed ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

predicts.year3 <- predicts.all %>%
  group_by(year)%>% #only change here
  summarise(response = mean(fit, na.rm = T),
            se.fit=mean(se.fit, na.rm = T),
            observed = mean(log_den_STRPURAD, na.rm = T),
            sd.obs = sd(log_den_STRPURAD, na.rm = T),
            n.obs = length(log_den_STRPURAD),
            se.obs=sd.obs/(sqrt(n.obs)))%>%
  dplyr::filter(year != "2006") %>% # remove 2006 because only 1 observation
  ungroup() %>%
  glimpse()


ggmod.year <- ggplot( data=predicts.year3) +
  geom_bar(aes(x=year,y=response,fill=year), stat = "identity")+
  scale_fill_viridis(discrete = T, direction = 1, option = "G") +
  geom_errorbar(aes(x=year,y=response, ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  geom_line(aes(x=year, y= observed), group = 1, color = 'red', size = 1.5) +
  ylab("Log density S.purpuratus")+
  xlab('survey_year')+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year



predicts.all2 <- predicts.all %>%
  mutate(fit.exp = exp(fit)) %>%
  glimpse()

predicts.year2 <- predicts.all2 %>%
  group_by(year)%>% #only change here
  summarise(response = mean(fit.exp, na.rm = T),
            sd.fit = sd(fit.exp, na.rm = T),
            n.fit = length(fit.exp),
            se.fit=sd.fit/(sqrt(n.fit)))%>%
  ungroup() %>%
  glimpse()


ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year2) +
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


#write.csv(predicts.all,paste(o2.dir, paste(name, subname, "urchins_preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_STRPURAD) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_STRPURAD) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

#write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "urchins_preds-obs.csv", sep = '-'), sep ='/'))

# plot pred vs. observed ----
pred.obs.all %>%
  ggplot(aes(x = log_den_STRPURAD, y = fit, color = zone)) +
  geom_point() 


# Plot observed vs. predicted ----
library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = log_den_STRPURAD, col = zone)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = log_den_STRPURAD)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               #aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               aes(label = paste(..rr.label..)),
               parse = TRUE, size = 4) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='S. purpuratus') +
  theme_bw()
p

namep <- paste(name, subname, "urchins_pred-obs.png", sep ='_')
#ggsave(namep, plot = last_plot(), device = 'png', path = o2.dir)


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
  pivot_longer(cols = c(fit, log_den_STRPURAD), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)

# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all, aes(x = year, y = latitude, color = values)) +
  geom_jitter(size = 3.5, pch = 16) +
  labs(x = 'Year', y = 'Latitude', color = "Log S.purpuratus") +
  facet_wrap(~type) +
  scale_color_viridis(option = 'E', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "urchins_pred-obs_latitude.png", sep ='_')
#ggsave(namep, plot = last_plot(), device = 'png', path = o2.dir)


###

### PREDICT URCHINS SPATIALLY FOR EVERY YEAR ----


# 1. Get rasters of predictors for every year ----

# set directory of rasters --

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data"

# check variables ---
mod$pred.formula

# Days_10N + Max_Monthly_Anomaly_Upwelling_Temp + Mean_Monthly_Summer_Temp + 
# Min_Monthly_Nitrate + site_name + zone + year




# set extent
e.nc <- extent(-124.60, -120.6956, 37.62184, 42.00996)
# bathy
sp.dir <- here('spatial')
nc.mask <- rast(paste(sp.dir, "mask_NC.tif", sep ='/'))

# Get substrate ----
sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

# get Depth ----
dir(sub.dir)

d <- rast(paste(sub.dir, "depth_mean_nc.all_wInterp_300m_30m.tif", sep ='/'))
d

# project depth --
crs1 <- "epsg:4326"

# transform  --
d1 <- project(d, crs1)


# crop to NC --
d2 <- crop(d1, extent(e.nc))
plot(d2)


##

# Get rock ----

dir(sub.dir)

rock <- rast(paste(sub.dir, "prob_rock_nc_MASKED_LatLong_all_300m_wInterp.tif", sep ='/'))
rock 

# project depth --
#crs1 <- "epsg:4326"

# transform  --
#d1 <- project(d, crs1)


# crop to NC --
rock2 <- crop(rock , extent(e.nc))
plot(rock2)

rock3 <- resample(rock2, d2)
plot(rock3)



# load temperature predictors ----

max_an_up_temp <- rast(paste(re.dir, 'Temperature', "Max_Monthly_Anomaly_Upwelling_Temp.tif", sep ='/'))
max_an_up_temp

# crop to NC --
max_an_up_temp2 <- crop(max_an_up_temp, extent(e.nc))
plot(max_an_up_temp2[[1]])

# resample predictors to bathy ----
#max_an_up_temp3 <- resample(max_an_up_temp2, nc.mask)
max_an_up_temp3 <- resample(max_an_up_temp2, d2)

# mask predictors to bathy ----
max_an_up_temp4 <- mask(max_an_up_temp3, d2)
plot(max_an_up_temp4[[1]])

###



# load nitrate predictors ----

days10n <- rast(paste(re.dir, "Nitrate", "Days_10N.tif", sep ='/'))
days10n

# crop to NC --
days10n2 <- crop(days10n, extent(e.nc))
plot(days10n2[[1]])

# resample predictors to bathy ----
days10n3 <- resample(days10n2, d2)

# mask predictors to bathy ----
days10n4 <- mask(days10n3, d2)
plot(days10n4[[1]])


###

## load WH_95prc ----

# load raster data --

w.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Wave_data"

wave.dir <- paste(w.dir, "Wave_extended2", "wave_95prc", sep ='/')

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
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
wave.stack <- do.call("c", tfiles)
plot(wave.stack[[1]])


# # crop to NC --
wave.stack2 <- crop(wave.stack, e.nc)
plot(wave.stack2[[1]])

# resample predictors to bathy ----
wave.stack3 <- resample(wave.stack2, d2)

# mask predictors to bathy ----
wave.stack4 <- mask(wave.stack3, d2)
plot(wave.stack4[[1]])



###



##



# stack and try to predict one year ----

preds1 <- c(d2, rock3, max_an_up_temp4[[1]], days10n4[[1]],wave.stack4[[1]])
names(preds1)

year1998 <- classify(rock3[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(2004:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("mean_depth", 
                   "mean_prob_of_rock",
                   "Max_Monthly_Anomaly_Upwelling_Temp",  
                   "Days_10N",  
                   "wh_95prc",
                   'year') 

year.raster <- year1998

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
names(preds3) <- c("mean_depth", 
                   "mean_prob_of_rock",
                   "Max_Monthly_Anomaly_Upwelling_Temp",  
                   "Days_10N",  
                   "wh_95prc",
                   'year',
                   "site_name") 

##

# zone ----

zone.raster <- d2
names(zone.raster) <- 'zone'


rec.m <-  c(-31, -10, 2,
            -10, 0.1, 1)

rclmat <- matrix(rec.m, ncol=3, byrow=TRUE)

zone.raster2 <- classify(zone.raster, rclmat, right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("mean_depth", 
                   "mean_prob_of_rock",
                   "Max_Monthly_Anomaly_Upwelling_Temp",  
                   "Days_10N",  
                   "wh_95prc",
                   'year',
                   "site_name",
                   "zone") 


# LOOP PREDICT raster ----

# set predictios directory --
#preds.dir <- paste(o2.dir, "log_sp_predictions", sep ='/')
preds.dir <- paste(o.dir, "gam_urchins3", "log_sp_predictions", sep ='/')

# loop --

for (i in 1:length(year.list)) {
  
  # 1. stack predictors for that year
  preds1 <- c(d2, rock3, max_an_up_temp4[[i+6]], days10n4[[i+6]], wave.stack4[[i]])
  
  # 2. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  preds2 <- c(preds1, year.r)

  # 3. stack zone
  preds3 <- c(preds2, zone.raster2)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("mean_depth", 
                     "mean_prob_of_rock",
                     "Max_Monthly_Anomaly_Upwelling_Temp",  
                     "Days_10N",  
                     "wh_95prc",
                     'year',
                     "zone",
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, mod)
  plot(year.prediction)
  
  # 6. transform 
  #year.prediction2 <- exp(year.prediction)
  
  # save
  writeRaster(year.prediction, paste(preds.dir, paste(year.no, "Log_purple_urchins_NC_V3.tif", sep = '_'), sep ='/'))
  
}

plot(year.prediction2)


###

###



###

###


## LOOP - data frame ----



## loop to predict several years ----

urch.mod <- gam1


# make list of years --
#year.list <- paste(1998:2021)
year.list <- paste(2004:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(d2, cbind(-Inf, 0.1, 2004), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 




#outputs dir
preds.dir <- paste(o.dir, "gam_urchins3", "log_sp_predictions", sep ='/')
preds.dir <- paste(o.dir, "gam_urchins3", "sp_predictions", sep ='/')


# Version 3: V4_5.1.1_v3 : Using urchins3, which don't have VRM
# Version sp_predictions_v3 : using gam 5.1, with urchins 3


for (i in 1:length(year.list)) {
  
  # 1. stack predictors for that year
  preds1 <- c(d2, rock3, max_an_up_temp4[[i+6]], days10n4[[i+6]], wave.stack4[[i]])
  
  # 2. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster2)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("mean_depth", 
                     "mean_prob_of_rock",
                     "Max_Monthly_Anomaly_Upwelling_Temp",  
                     "Days_10N",  
                     "wh_95prc",
                     'year',
                     "zone",
                     "site_name")
  
  # get data frame
  df4 <- as.data.frame(preds4, xy = T) %>%
    mutate_at(vars(year, zone, site_name), list(as.factor)) %>%
    mutate(zone = recode_factor(zone, '1' = 'INNER', '2' = 'OUTER')) %>%
    glimpse()
  
  # 5. predict
  year.pred.df <- predict.gam(urch.mod, newdata=df4, type='response', se.fit=T)
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
  
  yp <- year.prediction-1
  plot(yp)
  
  # 7. Change the log
  year.prediction2 <- exp(yp)
  year.prediction2 <- exp(year.prediction)
  plot(year.prediction2)
  
  plot(year.prediction2, breaks = c(0,1,10,50,100,200,300, 2000))
  
  # save
  name.raster <- paste(year.no, "log_purple_urchins_NC_V3.tif", sep = '_')
  #writeRaster(year.prediction, paste(preds.dir, "sp_predictions", paste(year.no, "logNereo_preds_NC.tif", sep = '_'), sep ='/'))
  writeRaster(year.prediction, paste(preds.dir, name.raster, sep = '/'))
}






###

### GET STABILITY ----

#preds.dir <- paste(o2.dir, "sp_predictions", sep ='/')

# load raster data --

n.files <- dir(preds.dir)
# list files in source --
n.files <- list.files(preds.dir, pattern = '.tif', full.names = TRUE)
n.files
length(n.files)
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

# Calculate mean and Stability for urchins ----

u.dir <- paste(o.dir, "gam_urchins2", "sp_predictions", sep ='/')

# load raster data --

n.files <- dir(u.dir)
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
plot(preds.stack)

# transform the log 
preds.stack2 <- exp(preds.stack)
plot(preds.stack2)

### Calculate mean across years ----

stability.stack <- c()

# calculate mean across years
mean.urch <- mean(preds.stack, na.rm = T)
plot(mean.urch)


# calculate sd across years
sd.urch <- app(preds.stack, fun = sd, na.rm = T)
plot(sd.urch)


# calculate stability
s.urch <- mean.urch/sd.urch
plot(s.urch)

# calculate stability index
s_index.urch <- (s.urch*mean.urch)
plot(s_index.urch)

# stack 

stability.stack <- c(mean.urch, sd.urch, s.urch, s_index.urch)
names(stability.stack) <- c("mean_urch", "sd_urch", "s_urch", "s_index_urch")
plot(stability.stack)

# save 
writeRaster(stability.stack, paste(u.dir, "Stabilty_calculations_Urchins_NC_gam_V2.tif", sep ='/'), overwrite = T)



### Save as points ----

preds.stack[[10]]

# 2012
u2012 <- as.data.frame(preds.stack[[9]], xy = T)
head(u2012)
names(u2012)  <- c("x"  ,  "y"  ,  "u.2012")

write.csv(u2012, paste(preds.dir, "Sp_preds_urchins_2012_NC_gam_V3.csv", sep ='/'))

# 2013
u2013 <- as.data.frame(preds.stack[[10]], xy = T)
head(u2013)
names(u2013)  <- c("x"  ,  "y"  ,  "u.2013")

write.csv(u2013, paste(preds.dir, "Sp_preds_urchins_2013_NC_gam_V3.csv", sep ='/'))


# 2014
u2014 <- as.data.frame(preds.stack[[11]], xy = T)
head(u2014)
names(u2014)  <- c("x"  ,  "y"  ,  "u.2014")

write.csv(u2014, paste(preds.dir, "Sp_preds_urchins_2014_NC_gam_V3.csv", sep ='/'))


# 2015
u2015 <- as.data.frame(preds.stack[[12]], xy = T)
head(u2015)
names(u2015)  <- c("x"  ,  "y"  ,  "u.2015")

write.csv(u2015, paste(preds.dir, "Sp_preds_urchins_2015_NC_gam_V3.csv", sep ='/'))


# 2016
u2016 <- as.data.frame(preds.stack[[13]], xy = T)
head(u2016)
names(u2016)  <- c("x"  ,  "y"  ,  "u.2016")

write.csv(u2016, paste(preds.dir, "Sp_preds_urchins_2016_NC_gam_V3.csv", sep ='/'))


