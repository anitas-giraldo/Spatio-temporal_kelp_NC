###
### Script based on one created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 11 July 2022 with added orb velocity and NPP from Tom

## This script run the GAM of north coast with RCCA data --



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
o.dir <- here("outputs_nc_rcca")
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"

## Load info on years RCCA ----
years <- read.csv(paste(d.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


# get the sites from with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(total.years > 2) %>%
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


## Load RCCA data ----

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
    #wh.95 ,   wh.max,
    #npgo_mean , mei_mean,
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
         log_UBR_Max = log(UBR_Max + 1)) %>%
  dplyr::select(-c(UBR_Mean, UBR_Max)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
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
  droplevels() %>%
  glimpse() # Rows: 282

levels(PreMHW$year)
length(levels(PreMHW$year)) # 8

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

inTraining <- createDataPartition(dat2$log_den_NERLUE, p = 0.75, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'V7'
o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')
o2.dir

names(train.gam)
length(names(train.gam))

# 4. Define predictor variables ----

pred.vars <- c(#"Days_10N" ,
               #"Min_Monthly_Nitrate" ,
               "Max_Monthly_Nitrate",
               #"Mean_Monthly_Nitrate",
               #"Mean_Monthly_Upwelling_Nitrate" ,
               #"Max_Monthly_Anomaly_Nitrate"  ,     
               #"Mean_Monthly_Summer_Nitrate" ,
               #"Mean_Monthly_Temp" ,
               #"Mean_Monthly_Summer_Temp" ,
               #"MHW_Upwelling_Days" ,               
               "Min_Monthly_Anomaly_Temp" ,
               #"Max_Monthly_Anomaly_Upwelling_Temp",
               #"Min_Monthly_Temp" ,
               #"Mean_Monthly_Upwelling_Temp" ,      
               "mean_depth" ,
               #"mean_prob_of_rock" ,
               #"mean_slope"  ,
               #"wh_max" ,
               "wh_mean"  ,
               #"mean_waveyear"  ,
               #"wh_95prc" ,                       
               "Mean_Monthly_NPP" ,
               #"Max_Monthly_NPP_Upwelling" ,
               #"log_den_NERLUE" ,
               #"log_den_MESFRAAD" ,
               "log_den_STRPURAD" ,
               #"log_den_PYCHEL"  ,
               #"log_den_HALRUF" ,
               #"log_mean_vrm" ,                     
               #"log_Days_16C"  ,
               #"log_UBR_Mean"  ,
               #"log_Mean_Monthly_NPP_Upwelling" ,   
               #"log_Min_Monthly_NPP",
               "log_UBR_Max")

length(pred.vars) # 29



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
                                max.predictors = 7,
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
out.i <- mod.table
nrow(out.i)

out.all <- c(out.all,list(out.i)) 
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

out.list$variable.importance

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
  best.model.name=as.character(out.i$modname[2])
  
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

nrow(out.i)

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[4])
best.model <- out.list$success.models[[best.model.name]]
best.model <- print(out.i$formula[1])

gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_UBR_Max, k = 4, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 6, bs = "cr") + 
              s(mean_depth, k = 6, bs = "cr") + 
              s(MHW_Upwelling_Days, k = 6, bs = "cr") +
              s(mean_prob_of_rock, k = 6, bs = "cr") +
              s(wh_max, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")


gam1$var.summary$log_Days_16C
gam1$null.deviance
gam1$aic
gam1$coefficients
gam1$deviance
summary(gam1$var.summary)
summary(gam1$aic)
summary(gam1$coefficients)
summary(gam1$deviance)
s <- summary(gam1)
s$cov.unscaled
s$cov.scaled
s$s.table
gam.check(gam1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()


# GAM 2 ----
best.model <- print(out.i$formula[2])

gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_UBR_Max, k = 4, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 6, bs = "cr") + 
              s(mean_depth, k = 6, bs = "cr") + 
              s(Mean_Monthly_NPP, k = 6, bs = "cr") +
              s(mean_prob_of_rock, k = 6, bs = "cr") +
              s(wh_max, k = 6, bs = "cr") +
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


# GAM 3 ----
best.model <- print(out.i$formula[3])

gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_UBR_Max, k = 4, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 6, bs = "cr") + 
              s(mean_depth, k = 6, bs = "cr") + 
              s(Mean_Monthly_NPP, k = 6, bs = "cr") +
              s(mean_prob_of_rock, k = 6, bs = "cr") +
              s(wh_max, k = 6, bs = "cr") +
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
testdata <- expand.grid(log_den_PYCHEL=mean(mod$model$log_den_PYCHEL),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Mean_Monthly_Summer_Temp=mean(mod$model$Mean_Monthly_Summer_Temp),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, 
                log_Days_16C, 
                log_den_STRPURAD,
                Max_Monthly_Nitrate, 
                Mean_Monthly_NPP, 
                mean_waveyear,
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
  scale_fill_viridis(option = "A", discrete = T) +
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green", "purple")) +
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
  dplyr::select(log_den_NERLUE, log_den_PYCHEL, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)



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





#######

###

#### gam2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$success.models[[best.model.name]]
best.model <- print(out.i$formula[2])

gam1 <- gam(formula = log_den_NERLUE ~ s(log_den_PYCHEL, k = 5, bs = "cr") + 
              s(log_den_STRPURAD, k = 8, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 5, bs = "cr") + 
              s(Mean_Monthly_Summer_Temp, k = 5, bs = "cr") + 
              s(Mean_Monthly_Temp, k = 5, bs = "cr") +
              #s(Min_Monthly_Nitrate, k = 3, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
par(mfrow=c(3,3),mar=c(2,4,3,1))
gam.check(gam1)
dev.off()

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
pdf(file=paste(o2.dir, paste(name,subname,'log_den_NERLUE',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='log_den_NERLUE',outer=F)
dev.off()

##

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_PYCHEL=mean(mod$model$log_den_PYCHEL),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        Mean_Monthly_Summer_Temp=mean(mod$model$Mean_Monthly_Summer_Temp),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        #Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, log_den_PYCHEL, log_den_STRPURAD, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, 
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
  dplyr::select(log_den_NERLUE, log_den_PYCHEL, log_den_STRPURAD, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, 
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)



###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds_postMHW.csv", sep = '-'), sep ='/'))

ner.obs <- PostMHW %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs_postMHW.csv", sep = '-'), sep ='/'))

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

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

###

#######

###

#### gam3 ----

subname <- "3"

best.model.name=as.character(out.i$modname[3])
best.model <- out.list$success.models[[best.model.name]]
best.model <- print(out.i$formula[3])

gam1 <- gam(formula = log_den_NERLUE ~ s(log_den_PYCHEL, k = 6, bs = "cr") + 
              s(Days_10N, k = 10, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 5, bs = "cr") + 
              s(Mean_Monthly_Summer_Temp, k = 5, bs = "cr") + 
              s(Mean_Monthly_Temp, k = 6, bs = "cr") +
              s(Min_Monthly_Nitrate, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
par(mfrow=c(3,3),mar=c(2,4,3,1))
gam.check(gam1)
dev.off()

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
pdf(file=paste(o2.dir, paste(name,subname,'log_den_NERLUE',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='log_den_NERLUE',outer=F)
dev.off()

##

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_PYCHEL=mean(mod$model$log_den_PYCHEL),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        Days_10N=mean(mod$model$Days_10N),
                        Mean_Monthly_Summer_Temp=mean(mod$model$Mean_Monthly_Summer_Temp),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, log_den_PYCHEL, Days_10N, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, Min_Monthly_Nitrate,
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
  dplyr::select(log_den_NERLUE, log_den_PYCHEL, Days_10N, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)



###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds_postMHW.csv", sep = '-'), sep ='/'))

ner.obs <- PostMHW %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs_postMHW.csv", sep = '-'), sep ='/'))

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
  scale_color_viridis() +
  theme_bw()

namep <- paste(name, subname, "pred-obs_latitude_postMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)



#######

###

#### gam4 ----

# testing one of the best models without PYCHEL

subname <- "4"

# best.model.name=as.character(out.i$modname[3])
# best.model <- out.list$success.models[[best.model.name]]
# best.model <- print(out.i$formula[3])

gam1 <- gam(formula = log_den_NERLUE ~ s(log_den_STRPURAD, k = 5, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 5, bs = "cr") + 
              s(Mean_Monthly_Summer_Temp, k = 5, bs = "cr") + 
              s(Mean_Monthly_Temp, k = 5, bs = "cr") + 
              s(Min_Monthly_Nitrate, k = 5, bs = "cr") + 
              s(site_name, zone, bs = "re") + s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
par(mfrow=c(3,3),mar=c(2,4,3,1))
gam.check(gam1)
dev.off()

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
pdf(file=paste(o2.dir, paste(name,subname,'log_den_NERLUE',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='log_den_NERLUE',outer=F)
dev.off()

##

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        #Days_10N=mean(mod$model$Days_10N),
                        Mean_Monthly_Summer_Temp=mean(mod$model$Mean_Monthly_Summer_Temp),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, Min_Monthly_Nitrate,
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
  dplyr::select(log_den_NERLUE, log_den_STRPURAD, Max_Monthly_Nitrate,  Mean_Monthly_Summer_Temp,
                Mean_Monthly_Temp, Min_Monthly_Nitrate,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)



###

# Plot latitudinally ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds_postMHW.csv", sep = '-'), sep ='/'))

ner.obs <- PostMHW %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs_postMHW.csv", sep = '-'), sep ='/'))

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
  scale_color_viridis() +
  theme_bw()

namep <- paste(name, subname, "pred-obs_latitude_postMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###



###

### GAM V2 ----

# Use PRE and Post MHW data to fit the model

glimpse(dat2)
levels(dat2$year)

# 1. Select predictors for this GAM ----
names(dat2)

dat3 <- dat2 %>%
  dplyr::select(site_name, year, transect, zone, 
                log_den_NERLUE, log_den_STRPURAD, log_den_PYCHEL,
                Days_10N, Max_Monthly_Nitrate, Mean_Monthly_Summer_Nitrate,
                log_Days_16C, Mean_Monthly_Temp, MHW_Upwelling_Days, Min_Monthly_Temp, 
                Mean_Monthly_Upwelling_Temp, 
                wh.max, wh.95, npgo_mean, mei_mean)


# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat3$log_den_NERLUE, p = 0.7, list = FALSE)
train.gam <- dat3[ inTraining,]
test.gam  <- dat3[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'V2'
o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c("log_den_STRPURAD", "log_den_PYCHEL",
               "Days_10N", "Max_Monthly_Nitrate", "Mean_Monthly_Summer_Nitrate",
               "log_Days_16C", "Mean_Monthly_Temp", "MHW_Upwelling_Days","Min_Monthly_Temp", 
               "Mean_Monthly_Upwelling_Temp", 
               "wh.max", "wh.95", "npgo_mean", "mei_mean" )

length(pred.vars) # 14



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
                                k=3,
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
names(out.all) <- 'log_den_NERLUE'
names(var.imp) <- 'log_den_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "all.var.imp.csv", sep="/"))


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

# out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
#   glimpse()

#### gam1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
best.model


gam1 <- gam(formula = log_den_NERLUE ~ s(log_den_STRPURAD, k = 6, bs = "cr") + 
              s(Mean_Monthly_Summer_Nitrate, k = 5, bs = "cr") + 
              #s(Mean_Monthly_Summer_Temp, k = 6, bs = "cr") + 
              #s(Mean_Monthly_Temp, k = 6, bs = "cr") +
              s(Min_Monthly_Temp, k = 6, bs = "cr") +
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
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Min_Monthly_Temp=mean(mod$model$Min_Monthly_Temp),
                        #Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        #Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD,
                Mean_Monthly_Summer_Nitrate, Min_Monthly_Temp,
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
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'S of N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 2.5) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


###

###

###

## PREDICT ON ENTIRE DATA SET ----

glimpse(df)
names(df)

df.pred <- df %>%
  dplyr::select(year, site_name, transect, zone,
                longitude, latitude,
                den_NERLUE, den_STRPURAD, 
                Mean_Monthly_Summer_Nitrate,
                Min_Monthly_Temp) %>%
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(den_NERLUE, den_STRPURAD)) %>%
  glimpse()

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Min_Monthly_Temp=mean(mod$model$Min_Monthly_Temp),
                        #Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        #Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- df.pred %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD,
                Mean_Monthly_Summer_Nitrate, Min_Monthly_Temp,
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

glimpse(predicts.year)

# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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

glimpse(years)

years2 <- years %>%
  dplyr::select(site_name, longitude, latitude) %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  glimpse()

# add lat lons
pred.obs.all <- predicts.all %>% #pred.obs.all %>%
  left_join(years2, by = 'site_name') %>%
  pivot_longer(cols = c(fit, log_den_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)
glimpse(pred.obs.all)
# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all %>% dplyr::filter(type =='fit'), aes(x = year, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  #facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 4) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


####


####


### GAM V3 ####


# Use PRE and Post MHW data to fit the model

glimpse(dat2)
levels(dat2$year)

# 1. Select predictors for this GAM ----
names(dat2)

dat3 <- dat2 %>%
  dplyr::select(site_name, year, transect, zone, 
                log_den_NERLUE, log_den_STRPURAD, log_den_PYCHEL,
                Days_10N, Max_Monthly_Nitrate, Mean_Monthly_Summer_Nitrate,
                log_Days_16C, Mean_Monthly_Temp, MHW_Upwelling_Days, Min_Monthly_Temp, 
                Mean_Monthly_Upwelling_Temp, 
                wh.max, wh.95, npgo_mean, mei_mean,
                mean_depth, mean_prob_of_rock, log_mean_vrm, mean_slope)


# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat3$log_den_NERLUE, p = 0.7, list = FALSE)
train.gam <- dat3[ inTraining,]
test.gam  <- dat3[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'V3'
o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c("log_den_STRPURAD", "log_den_PYCHEL",
               "Days_10N", "Max_Monthly_Nitrate", "Mean_Monthly_Summer_Nitrate",
               "log_Days_16C", "Mean_Monthly_Temp", "MHW_Upwelling_Days","Min_Monthly_Temp", 
               "Mean_Monthly_Upwelling_Temp", 
               "wh.max", "wh.95", "npgo_mean", "mei_mean",
               "mean_depth", "mean_prob_of_rock", "log_mean_vrm", "mean_slope")

length(pred.vars) # 18



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
names(out.all) <- 'log_den_NERLUE'
names(var.imp) <- 'log_den_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "all.var.imp.csv", sep="/"))


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

# out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
#   glimpse()



#### gam1 ----

# load models if needed ----

out.i <- read.csv(paste(file=paste(o2.dir, "best_models.csv", sep="/"))) %>% glimpse()


subname <- "1"

best.model.name=as.character(out.i$modname[3])
best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
best.model


gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 5, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 7, bs = "cr") + 
              s(mean_prob_of_rock, k = 8, bs = "cr") +
              s(MHW_Upwelling_Days, k = 7, bs = "cr") +
              #s(wh.95, k = 6, bs = "cr") +
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
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        log_Days_16C=mean(mod$model$log_Days_16C),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        mean_prob_of_rock=mean(mod$model$mean_prob_of_rock),
                        MHW_Upwelling_Days=mean(mod$model$MHW_Upwelling_Days),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD, log_Days_16C, 
                Max_Monthly_Nitrate, mean_prob_of_rock, MHW_Upwelling_Days,
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
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'S of N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 2.5) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


###

###

###

## PREDICT ON ENTIRE DATA SET ----

glimpse(df)
names(df)

df.pred <- df %>%
  dplyr::select(year, site_name, transect, zone,
                longitude, latitude,
                den_NERLUE, den_STRPURAD, 
                Days_16C, 
                Max_Monthly_Nitrate,
                mean_prob_of_rock,
                MHW_Upwelling_Days) %>%
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_Days_16C = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(den_NERLUE, den_STRPURAD, Days_16C)) %>%
  glimpse()

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        log_Days_16C=mean(mod$model$log_Days_16C),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        mean_prob_of_rock=mean(mod$model$mean_prob_of_rock),
                        MHW_Upwelling_Days=mean(mod$model$MHW_Upwelling_Days),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- df.pred %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD,
                log_Days_16C, Max_Monthly_Nitrate,
                mean_prob_of_rock, MHW_Upwelling_Days,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)


predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit, na.rm = T),se.fit=mean(se.fit, na.rm = T))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()

glimpse(predicts.year)

# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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

glimpse(years)

years2 <- years %>%
  dplyr::select(site_name, longitude, latitude) %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  glimpse()

# add lat lons
pred.obs.all <- predicts.all %>% #pred.obs.all %>%
  left_join(years2, by = 'site_name') %>%
  pivot_longer(cols = c(fit, log_den_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)
glimpse(pred.obs.all)
# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all %>% dplyr::filter(type =='fit'), aes(x = year, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  #facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 4) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


####

####


#### 

### GAM 4 ####

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
o.dir <- here("outputs_nc_rcca")
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
    wh_max, wh_mean, mean_waveyear, wh_95prc
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
  #mutate(log_npp.mean = log(npp.mean + 1)) %>%
  #dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708


# Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)

# 1. Select predictors for this GAM ----
names(dat2)

dat3 <- dat2 %>%
  dplyr::select(site_name, year, transect, zone, 
                # bio
                log_den_NERLUE, log_den_STRPURAD, log_den_PYCHEL,
                # nitrate
                Days_10N, 
                Max_Monthly_Nitrate, Min_Monthly_Nitrate,
                Mean_Monthly_Nitrate, Max_Monthly_Anomaly_Nitrate, 
                Mean_Monthly_Summer_Nitrate,
                # temperature
                log_Days_16C, Mean_Monthly_Temp, MHW_Upwelling_Days, 
                Min_Monthly_Temp, Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_Anomaly_Upwelling_Temp,
                # climatic indices
                npgo_mean, mei_mean,
                # substrate
                mean_depth, mean_prob_of_rock, log_mean_vrm, mean_slope,
                # waves
                wh_max, wh_mean, mean_waveyear, wh_95prc)

length(dat3) # 29

# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat3$log_den_NERLUE, p = 0.7, list = FALSE)
train.gam <- dat3[ inTraining,]
test.gam  <- dat3[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'V4'
o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c(# bio
  'log_den_STRPURAD', 'log_den_PYCHEL',
  # nitrate
  'Days_10N', 
  'Max_Monthly_Nitrate', 'Min_Monthly_Nitrate',
  'Mean_Monthly_Nitrate', 'Max_Monthly_Anomaly_Nitrate', 
  'Mean_Monthly_Summer_Nitrate',
  # temperature
  'log_Days_16C', 'Mean_Monthly_Temp', 'MHW_Upwelling_Days', 
  'Min_Monthly_Temp', 'Mean_Monthly_Upwelling_Temp', 
  'Max_Monthly_Anomaly_Upwelling_Temp',
  # climatic indices
  'npgo_mean', 'mei_mean',
  # substrate
  'mean_depth', 'mean_prob_of_rock', 'log_mean_vrm', 'mean_slope',
  # waves
  'wh_max', 'wh_mean', 'mean_waveyear', 'wh_95prc')

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
names(out.all) <- 'log_den_NERLUE'
names(var.imp) <- 'log_den_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "all.var.imp.csv", sep="/"))


#out.i <- read.csv(paste(o2.dir, "best_models.csv", sep="/"))

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

# out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
#   glimpse()



#### gam1 ----

# load models if needed ----

out.i <- read.csv(paste(file=paste(o2.dir, "best_models.csv", sep="/"))) %>% glimpse()


subname <- "1"

best.model.name=as.character(out.i$modname[3])
best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
best.model.name


gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 5, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 4, bs = "cr") + 
              s(mean_waveyear, k = 8, bs = "cr") +
              #s(mean_prob_of_rock, k = 4, bs = "cr") +
              #s(mean_depth, k = 4, bs = "cr") +
              s(wh_max, k = 4, bs = "cr") +
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
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        log_Days_16C=mean(mod$model$log_Days_16C),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        mean_waveyear=mean(mod$model$mean_waveyear),
                        wh_max=mean(mod$model$wh_max),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD, log_Days_16C, 
                Max_Monthly_Nitrate, mean_waveyear, wh_max,
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
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'S of N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 2.5) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


###

###

###

## PREDICT ON ENTIRE DATA SET ----

glimpse(df)
names(df)

df.pred <- df %>%
  dplyr::select(year, site_name, transect, zone,
                longitude, latitude,
                den_NERLUE, den_STRPURAD, 
                Days_16C, 
                Max_Monthly_Nitrate,
                mean_waveyear,
                wh_max) %>%
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_Days_16C = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(den_NERLUE, den_STRPURAD, Days_16C)) %>%
  glimpse()

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        log_Days_16C=mean(mod$model$log_Days_16C),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        mean_waveyear=mean(mod$model$mean_waveyear),
                        wh_max=mean(mod$model$wh_max),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- df.pred %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD,
                log_Days_16C, Max_Monthly_Nitrate,
                mean_waveyear, wh_max,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)


predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit, na.rm = T),se.fit=mean(se.fit, na.rm = T))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()

glimpse(predicts.year)

# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs_all.csv", sep = '-'), sep ='/'))

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

namep <- paste(name, subname, "pred-obs_all.png", sep ='_')
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

glimpse(years)

years2 <- years %>%
  dplyr::select(site_name, longitude, latitude) %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  glimpse()

# add lat lons
pred.obs.all <- predicts.all %>% #pred.obs.all %>%
  left_join(years2, by = 'site_name') %>%
  pivot_longer(cols = c(fit, log_den_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)
glimpse(pred.obs.all)
# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all %>% dplyr::filter(type =='fit'), aes(x = year, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  #facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_all.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 4) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


# save --
namep <- paste(name, subname, "s_latitude_all.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o2.dir)



####

#### MODEL URCHINS ####

# Use PRE and Post MHW data to fit the model

glimpse(dat3)
levels(dat3$year)

# 1. Select predictors for this GAM ----
names(dat3)

# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat3$log_den_STRPURAD, p = 0.8, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'urchins3'

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
               "wh_95prc")

length(pred.vars) # 24

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
                                k=3,
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


# Manually make the most parsimonious GAM models ----

# out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
#   glimpse()

#### gam1 ----

subname <- "3"

best.model.name=as.character(out.i$modname[22])
best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
best.model


gam1 <- gam(formula = log_den_STRPURAD ~ s(Days_10N, k = 30, bs = "cr") +
              #s(log_mean_vrm, k = 10, bs = "cr") +
              log_mean_vrm +
              s(Max_Monthly_Anomaly_Upwelling_Temp, k = 30, bs = "cr") + 
              #s(mean_depth, k = 10, bs = "cr") + 
              mean_depth +
              #s(mean_prob_of_rock, k = 10, bs = "cr") +
              mean_prob_of_rock +
              s(wh_95prc, k = 30, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "GCV")



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
head(mod$model)
testdata <- expand.grid(Days_10N=mean(mod$model$Days_10N),
                        wh_95prc=mean(mod$model$wh_95prc),
                        mean_depth=mean(mod$model$mean_depth),
                        Max_Monthly_Anomaly_Upwelling_Temp=mean(mod$model$Max_Monthly_Anomaly_Upwelling_Temp),
                        mean_prob_of_rock=mean(mod$model$mean_prob_of_rock),
                        log_mean_vrm=mean(mod$model$log_mean_vrm),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)
glimpse(testdata)


testdata <- test.gam %>%
  dplyr::select(log_den_STRPURAD, Days_10N, wh_95prc,
                mean_depth, Max_Monthly_Anomaly_Upwelling_Temp,
                mean_prob_of_rock, log_mean_vrm,
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
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "urchins_preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_STRPURAD) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_STRPURAD) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "urchins_preds-obs.csv", sep = '-'), sep ='/'))

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
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "urchins_pred-obs.png", sep ='_')
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
  labs(x = 'Zone', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "urchins_pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


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

# bathy
sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

vrm <- rast(paste(sub.dir, "vrm_nc_all-mapped_300m_wInterp.tif", sep ='/'))

# transform to log vrm --
vrm1 <- vrm + 1
plot(vrm1)

log_vrm <- log(vrm1)
plot(log_vrm)

# transform vrm --
log_vrm1 <- project(log_vrm, "EPSG:4326")


# crop to NC --
log_vrm2 <- crop(log_vrm1, extent(e.nc))
plot(log_vrm2)

log_vrm3 <- project(log_vrm2, "EPSG:26910")




# load temperature predictors ----

max_an_up_temp <- rast(paste(re.dir, 'Temperature', "Max_Monthly_Anomaly_Upwelling_Temp.tif", sep ='/'))
max_an_up_temp

# crop to NC --
max_an_up_temp2 <- crop(max_an_up_temp, extent(e.nc))
plot(max_an_up_temp2[[1]])

# resample predictors to bathy ----
#max_an_up_temp3 <- resample(max_an_up_temp2, nc.mask)
max_an_up_temp3 <- resample(max_an_up_temp2, log_vrm3)

# mask predictors to bathy ----
max_an_up_temp4 <- mask(max_an_up_temp3, log_vrm3)
plot(max_an_up_temp4)

###

Days_16C <- rast(paste(re.dir, 'Temperature', "Days_16C.tif", sep ='/'))
Days_16C
plot(Days_16C)

# log transform it if needed --
Days_16_1 <- Days_16C + 1
plot(Days_16_1)

log_days16C <- log(Days_16_1)
plot(log_days16C)

# crop to NC --
log_days16C2 <- crop(log_days16C, extent(e.nc))
plot(log_days16C2[[1]])

# resample predictors to bathy ----
log_days16C3 <- resample(log_days16C2, log_vrm3)

# mask predictors to bathy ----
log_days16C4 <- mask(log_days16C3, log_vrm3)
plot(log_days16C4)


# load nitrate predictors ----

days10n <- rast(paste(re.dir, "Nitrate", "Days_10N.tif", sep ='/'))
days10n

# crop to NC --
days10n2 <- crop(days10n, extent(e.nc))
plot(days10n2[[1]])

# resample predictors to bathy ----
days10n3 <- resample(days10n2, log_vrm3)

# mask predictors to bathy ----
days10n4 <- mask(days10n3, log_vrm3)
plot(days10n4)


###

min_nit <- rast(paste(re.dir, "Nitrate", "Min_Monthly_Nitrate.tif", sep ='/'))
min_nit

# crop to NC --
min_nit2 <- crop(min_nit, extent(e.nc))
plot(min_nit2[[2]])

# resample predictors to bathy ----
min_nit3 <- resample(min_nit2, log_vrm3)

# mask predictors to bathy ----
min_nit4 <- mask(min_nit3, log_vrm3)
plot(min_nit4)


###



##



# stack and try to predict one year ----

preds1 <- c(log_vrm3, max_an_up_temp4[[1]], log_days16C4[[1]], days10n4[[1]], min_nit4[[1]])
names(preds1)

year1998 <- classify(min_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("log_mean_vrm", "Max_Monthly_Anomaly_Upwelling_Temp", "log_Days_16C"  , "Days_10N",                          
"Min_Monthly_Nitrate", 'year') 


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
names(preds3) <- c("log_mean_vrm", "Max_Monthly_Anomaly_Upwelling_Temp", "log_Days_16C"  , "Days_10N",                          
                   "Min_Monthly_Nitrate", 
                   "year",
                   "site_name") 

##

# zone ----

zone.raster <- year1998
names(zone.raster) <- 'zone'

zone.raster2 <- classify(zone.raster, cbind(0, Inf, 1), right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("log_mean_vrm", "Max_Monthly_Anomaly_Upwelling_Temp", "log_Days_16C"  , "Days_10N",                          
                   "Min_Monthly_Nitrate", 
                   "year",
                   "site_name",
                   "zone") 


##

# test predict ----

testp <- predict(preds4, mod)
plot(testp)

#testp2 <- approximate(testp, method = 'linear')

## aggragate with threshold ----
# mean.thresh <- function(x) {ifelse(length(which(is.na(x) == T)) >= 70, NA, mean(x, na.rm = T))} 
# 
# sd.thresh <- function(x) {ifelse(length(which(is.na(x) == T)) >= 70, NA, sd(x, na.rm = T))}


### LOOP ----

## loop to predict several years ----

# make list of years --
year.list <- paste(1998:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(min_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 

zone.raster <- year.raster
year.raster <- classify(year.raster, cbind(0, Inf, 1), right=FALSE)
names(zone.raster) <- 'zone'

# make sites raster ----

# each site is latitude --
rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)
# make raster
site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year1998))
site.raster
plot(site.raster)
# extend to fit other rasters
site.raster2 <- extend(site.raster, year1998)
plot(site.raster2)


for (i in 1:length(year.list)) {
  
  # 1. stack predictors for that year
  preds1 <- c(log_vrm3, max_an_up_temp4[[i]], log_days16C4[[i]], days10n4[[i]], min_nit4[[i]])
  
  # 2. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  preds2 <- c(preds1, year.r)

  # 3. stack zone
  preds3 <- c(preds2, zone.raster)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("log_mean_vrm",
                     "Max_Monthly_Anomaly_Upwelling_Temp",
                     "log_Days_16C",
                     "Days_10N"  ,
                     "Min_Monthly_Nitrate",
                     "year"     ,
                     "zone" ,
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, mod)
  plot(year.prediction)
  
  # save
  writeRaster(year.prediction, paste(o2.dir, "sp_predictions", paste(year.no, "purple_urchins_NC.tif", sep = '_'), sep ='/'))
  
}

plot(year.prediction)

r1 <- rast(paste(o2.dir, "sp_predictions", "2007_purple_urchins_NC.tif", sep ='/'))
plot(r1)

plot(year.r)


####

####

## PREDICT NEREOCYSTIS ----

## get model again ----

gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 5, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 7, bs = "cr") + 
              s(mean_prob_of_rock, k = 8, bs = "cr") +
              s(MHW_Upwelling_Days, k = 7, bs = "cr") +
              #s(wh.95, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)


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

prob_rock <- rast(paste(sub.dir, "prob_rock_nc.all_300m_wInterp.tif", sep ='/'))
prob_rock



# load temperature predictors ----


Days_16C <- rast(paste(re.dir, 'Temperature', "Days_16C.tif", sep ='/'))
Days_16C

# transform to log --
days_16C1 <- Days_16C + 1

log_days_16C <- log(days_16C1)

# crop to NC --
# Days_16C2 <- crop(Days_16C, extent(prob_rock))
# plot(min_temp2[[1]])

# resample predictors to bathy ----
Days_16C3 <- resample(log_days_16C, prob_rock)

# mask predictors to bathy ----
Days_16C4 <- mask(Days_16C3, prob_rock)
plot(Days_16C4)

###

# MHW 

mhw_up <- rast(paste(re.dir, 'Temperature', "MHW_Upwelling_Days.tif", sep ='/'))
mhw_up

# resample predictors to bathy ----
mhw_up3 <- resample(mhw_up, prob_rock)

# mask predictors to bathy ----
mhw_up4 <- mask(mhw_up3, prob_rock)
plot(mhw_up4)





# load nitrate predictors ----

max_nit <- rast(paste(re.dir, "Nitrate", "Max_Monthly_Nitrate.tif", sep ='/'))
max_nit

# # crop to NC --
# mean_sum_nit2 <- crop(mean_sum_nit, extent(e.nc))
# plot(mean_sum_nit2[[1]])

# resample predictors to bathy ----
max_nit3 <- resample(max_nit, prob_rock)

# mask predictors to bathy ----
max_nit4 <- mask(max_nit3, prob_rock)
plot(max_nit4)


# load purple urchin predictions ----
urch.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Spatio_temporal_GAMs/outputs_nc_rcca/gam_urchins2/sp_predictions"

# load raster data --

u.files <- dir(urch.dir)
u.files <- list.files(urch.dir, pattern = '.tif')
u.files
length(u.files)




##



# stack and try to predict one year ----

preds1 <- c(prob_rock, Days_16C4[[1]], mhw_up4[[1]], max_nit4[[1]])
names(preds1)

year1998 <- classify(Days_16C4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("mean_prob_of_rock",  "Log_Days_16C", "MHW_Upwelling_Days"  , "Max_Monthly_Nitrate",                          
                   'year') 


##


# sites ----

rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)

site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year1998))
site.raster

plot(site.raster)

ext(year1998)
ext(site.raster)

site.raster2 <- extend(site.raster, year1998)

preds3 <- c(preds2, site.raster2)
names(preds3) <- c("mean_prob_of_rock",  "Log_Days_16C", "MHW_Upwelling_Days"  , "Max_Monthly_Nitrate",                          
                   "year"     ,
                   "site_name") 

##

# zone ----

zone.raster <- year1998
names(zone.raster) <- 'zone'

zone.raster2 <- classify(zone.raster, cbind(0, Inf, 1), right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("mean_prob_of_rock",  "Log_Days_16C", "MHW_Upwelling_Days"  , "Max_Monthly_Nitrate",                          
                   "year"     ,
                   "site_name",
                   "zone")


##



### LOOP ----

## loop to predict several years ----

# make list of years --
year.list <- paste(1998:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(prob_rock, cbind(0, Inf, 1998), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 

zone.raster <- year.raster
year.raster <- classify(year.raster, cbind(0, Inf, 1), right=FALSE)
names(zone.raster) <- 'zone'

# make sites raster ----

# each site is latitude --
rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)
# make raster
site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year1998))
site.raster
# extend to fit other rasters
site.raster2 <- extend(site.raster, year1998)

# outputs dir ----

preds.dir <- paste(o.dir, "gam_V3", "sp_predictions", sep ='/')


for (i in 1:length(year.list)) {
  
  # 1. get urchins
  urchin.rast <- rast(paste(urch.dir, u.files[i], sep ='/'))
  urchin.rast2 <- resample(urchin.rast, prob_rock)
  
  # 2. stack with predictors for that year
  env.raster <- c(prob_rock, Days_16C4[[i]], mhw_up4[[i]], max_nit4[[i]])
  
  preds1 <- c(urchin.rast2, env.raster)
  
  # 3. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("log_den_STRPURAD",
                     "mean_prob_of_rock",  
                     "log_Days_16C", 
                     "MHW_Upwelling_Days"  , 
                     "Max_Monthly_Nitrate",                          
                     "year"     ,
                     "zone",
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, nereo.mod)
  plot(year.prediction)
  
  # 6. Change the log
  year.prediction2 <- exp(year.prediction)
  
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



###


###


## GAM V5 ----

## ADDED ORB VEL AND NPP FROM TOM --

# 3. Set parameters to save outputs ----

name <- 'V5'
o2.dir <- paste(o.dir, paste("gam", name, sep = '_'), sep ='/')

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c(# bio
  'log_den_STRPURAD', 'log_den_PYCHEL',
  # nitrate
  'Days_10N', 
  'Max_Monthly_Nitrate', 'Min_Monthly_Nitrate',
  'Mean_Monthly_Nitrate', 'Max_Monthly_Anomaly_Nitrate', 
  'Mean_Monthly_Summer_Nitrate',
  # temperature
  'log_Days_16C', 'Mean_Monthly_Temp', 'MHW_Upwelling_Days', 
  'Min_Monthly_Temp', 'Mean_Monthly_Upwelling_Temp', 
  'Max_Monthly_Anomaly_Upwelling_Temp',
  # climatic indices
  #'npgo_mean', 'mei_mean',
  # substrate
  'mean_depth', 'mean_prob_of_rock', 'log_mean_vrm', 'mean_slope',
  # waves
  'wh_max', 'wh_mean', 'mean_waveyear', 'wh_95prc', 
  # Orb vel
  'log_UBR_Mean',
  # NPP
  'Mean_Monthly_NPP', 'Max_Monthly_NPP_Upwelling', 'log_Mean_Monthly_NPP_Upwelling', 'log_Min_Monthly_NPP')

length(pred.vars) # 27



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

#beep()

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
names(out.all) <- 'log_den_NERLUE'
names(var.imp) <- 'log_den_NERLUE'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "all.var.imp.csv", sep="/"))


#out.i <- read.csv(paste(o2.dir, "best_models.csv", sep="/"))

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

# out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
#   glimpse()



#### gam1 ----

# load models if needed ----

out.i <- read.csv(paste(file=paste(o2.dir, "best_models.csv", sep="/"))) %>% glimpse()


subname <- "1"

best.model.name=as.character(out.i$modname[3])
best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
best.model.name


gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 5, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 4, bs = "cr") + 
              s(mean_waveyear, k = 8, bs = "cr") +
              s(mean_prob_of_rock, k = 4, bs = "cr") +
              #s(mean_depth, k = 4, bs = "cr") +
              s(Mean_Monthly_NPP, k = 4, bs = "cr") +
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
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        log_Days_16C=mean(mod$model$log_Days_16C),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        mean_waveyear=mean(mod$model$mean_waveyear),
                        wh_max=mean(mod$model$wh_max),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, 
                log_den_STRPURAD, 
                log_Days_16C, 
                Max_Monthly_Nitrate, 
                mean_waveyear,
                mean_prob_of_rock,
                Mean_Monthly_NPP,
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
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'S of N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 2.5) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


###

###

###

## PREDICT ON ENTIRE DATA SET ----

glimpse(df)
names(df)

df.pred <- df %>%
  dplyr::select(year, site_name, transect, zone,
                longitude, latitude,
                den_NERLUE, den_STRPURAD, 
                Days_16C, 
                Max_Monthly_Nitrate,
                mean_waveyear,
                wh_max) %>%
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_Days_16C = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(den_NERLUE, den_STRPURAD, Days_16C)) %>%
  glimpse()

## PREDICT  gam ----
mod<-gam1
head(mod$model)
testdata <- expand.grid(log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
                        log_Days_16C=mean(mod$model$log_Days_16C),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        mean_waveyear=mean(mod$model$mean_waveyear),
                        wh_max=mean(mod$model$wh_max),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- df.pred %>%
  dplyr::select(log_den_NERLUE, log_den_STRPURAD,
                log_Days_16C, Max_Monthly_Nitrate,
                mean_waveyear, wh_max,
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)


predicts.year = testdata%>%data.frame(fits)%>%
  group_by(year)%>% #only change here
  summarise(response=mean(fit, na.rm = T),se.fit=mean(se.fit, na.rm = T))%>%
  ungroup()
# write.csv(predicts.year,paste(o.dir, "predicts_mod1.csv", sep ='/')) #there is some BUG in dplyr - that this fixes
# predicts.year<-read.csv(paste(o.dir, "predicts_mod1.csv", sep ='/')) %>%
#   mutate_at(vars(survey_year), list(as.factor)) %>%
#   glimpse()

glimpse(predicts.year)

# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_NERLUE) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_NERLUE) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "preds-obs_all.csv", sep = '-'), sep ='/'))

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

namep <- paste(name, subname, "pred-obs_all.png", sep ='_')
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

glimpse(years)

years2 <- years %>%
  dplyr::select(site_name, longitude, latitude) %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  glimpse()

# add lat lons
pred.obs.all <- predicts.all %>% #pred.obs.all %>%
  left_join(years2, by = 'site_name') %>%
  pivot_longer(cols = c(fit, log_den_NERLUE), names_to = 'type', values_to = 'values') %>%
  mutate_at(vars(zone), list(as.character)) %>%
  #mutate(zone = strsplit(transect_unique, "[_]")[[1]]) %>%
  mutate(zone = sub("\\_.*", "", zone)) %>%
  mutate_at(vars(type, zone), list(as.factor)) %>%
  glimpse()



levels(pred.obs.all$zone)
glimpse(pred.obs.all)
# PLOT LATITUDINALLY
## not a map
ggplot(pred.obs.all %>% dplyr::filter(type =='fit'), aes(x = year, y = latitude, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'Zone', y = 'Latitude') +
  #facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_all.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# Calculate Stability ----
glimpse(pred.obs.all)

stability <- pred.obs.all %>%
  dplyr::filter(type == 'fit') %>%
  group_by(site_name) %>%
  summarise(mean_nereo = mean(values, na.rm = T),
            sd_nereo = sd(values, na.rm = T),
            stability = mean_nereo/sd_nereo,
            s.index = stability*mean_nereo, 
            latitude = max(latitude),
            longitude = max(longitude)) %>%
  glimpse()


## PLOT on map ----

shp.dir <- here('spatial')

# load California shapefile --
# https://cengel.github.io/R-spatial/mapping.html
ca <- st_read(paste(shp.dir, "North_CA.shp", sep = '/'))
plot(ca)
ca

cax <- st_transform(ca, CRS("+proj=longlat +datum=NAD83 +no_defs"))
cax
plot(cax, col = 'blue')



dfsp <- st_as_sf(x = stability,
                 coords = c('longitude', 'latitude'), crs = CRS("+proj=longlat +datum=NAD83 +no_defs"))



# plot with ggplot --

# plot train dataset --
p <- ggplot() +
  geom_sf(data = cax, aes(fill='Id')) +
  labs(title = 'N. luetkeana') +
  geom_sf(data = dfsp, aes(color = s.index), size = 4) +
  scale_fill_manual(name = NULL, values = c('Id' = 'light yellow'), labels = 'California') +
  #scale_color_manual(name = 'CV density of N. luetkeana') +
  scale_color_viridis(direction = -1) +
  #facet_wrap(~zone)+
  theme_bw() +
  theme(legend.title = element_text(size= 12, face = 'bold'),
        legend.text = element_text(size= 12))

p


# save --
namep <- paste(name, subname, "s_latitude_all.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o2.dir)



####

#### MODEL URCHINS ####

# Use PRE and Post MHW data to fit the model

glimpse(dat3)
levels(dat3$year)

# 1. Select predictors for this GAM ----
names(dat3)

# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat3$log_den_STRPURAD, p = 0.8, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'urchins3'

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
               "wh_95prc")

length(pred.vars) # 24

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
                                k=3,
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


# Manually make the most parsimonious GAM models ----

# out.i <- read.csv(paste(o2.dir, "best_models.csv", sep ='/')) %>%
#   glimpse()

#### gam1 ----

subname <- "3"

best.model.name=as.character(out.i$modname[22])
best.model <- out.list$success.models[[best.model.name]]
#best.model <- out.list$formula[[best.model.name]]
best.model


gam1 <- gam(formula = log_den_STRPURAD ~ s(Days_10N, k = 30, bs = "cr") +
              #s(log_mean_vrm, k = 10, bs = "cr") +
              log_mean_vrm +
              s(Max_Monthly_Anomaly_Upwelling_Temp, k = 30, bs = "cr") + 
              #s(mean_depth, k = 10, bs = "cr") + 
              mean_depth +
              #s(mean_prob_of_rock, k = 10, bs = "cr") +
              mean_prob_of_rock +
              s(wh_95prc, k = 30, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "GCV")



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
head(mod$model)
testdata <- expand.grid(Days_10N=mean(mod$model$Days_10N),
                        wh_95prc=mean(mod$model$wh_95prc),
                        mean_depth=mean(mod$model$mean_depth),
                        Max_Monthly_Anomaly_Upwelling_Temp=mean(mod$model$Max_Monthly_Anomaly_Upwelling_Temp),
                        mean_prob_of_rock=mean(mod$model$mean_prob_of_rock),
                        log_mean_vrm=mean(mod$model$log_mean_vrm),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)
glimpse(testdata)


testdata <- test.gam %>%
  dplyr::select(log_den_STRPURAD, Days_10N, wh_95prc,
                mean_depth, Max_Monthly_Anomaly_Upwelling_Temp,
                mean_prob_of_rock, log_mean_vrm,
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
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green")) +
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

write.csv(predicts.all,paste(o2.dir, paste(name, subname, "urchins_preds.csv", sep = '-'), sep ='/'))

ner.obs <- test.gam %>%
  dplyr::select(year, site_name, zone, log_den_STRPURAD) %>%
  glimpse()

pred.obs.all <- predicts.all %>%
  dplyr::select(-log_den_STRPURAD) %>%
  left_join(ner.obs, by = c('year', 'site_name', 'zone')) %>%
  glimpse()

write.csv(pred.obs.all ,paste(o2.dir, paste(name, subname, "urchins_preds-obs.csv", sep = '-'), sep ='/'))

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
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "urchins_pred-obs.png", sep ='_')
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
  labs(x = 'Zone', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "urchins_pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


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

# bathy
sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

vrm <- rast(paste(sub.dir, "vrm_nc_all-mapped_300m_wInterp.tif", sep ='/'))

# transform to log vrm --
vrm1 <- vrm + 1
plot(vrm1)

log_vrm <- log(vrm1)
plot(log_vrm)

# transform vrm --
log_vrm1 <- project(log_vrm, "EPSG:4326")


# crop to NC --
log_vrm2 <- crop(log_vrm1, extent(e.nc))
plot(log_vrm2)

log_vrm3 <- project(log_vrm2, "EPSG:26910")




# load temperature predictors ----

max_an_up_temp <- rast(paste(re.dir, 'Temperature', "Max_Monthly_Anomaly_Upwelling_Temp.tif", sep ='/'))
max_an_up_temp

# crop to NC --
max_an_up_temp2 <- crop(max_an_up_temp, extent(e.nc))
plot(max_an_up_temp2[[1]])

# resample predictors to bathy ----
#max_an_up_temp3 <- resample(max_an_up_temp2, nc.mask)
max_an_up_temp3 <- resample(max_an_up_temp2, log_vrm3)

# mask predictors to bathy ----
max_an_up_temp4 <- mask(max_an_up_temp3, log_vrm3)
plot(max_an_up_temp4)

###

Days_16C <- rast(paste(re.dir, 'Temperature', "Days_16C.tif", sep ='/'))
Days_16C
plot(Days_16C)

# log transform it if needed --
Days_16_1 <- Days_16C + 1
plot(Days_16_1)

log_days16C <- log(Days_16_1)
plot(log_days16C)

# crop to NC --
log_days16C2 <- crop(log_days16C, extent(e.nc))
plot(log_days16C2[[1]])

# resample predictors to bathy ----
log_days16C3 <- resample(log_days16C2, log_vrm3)

# mask predictors to bathy ----
log_days16C4 <- mask(log_days16C3, log_vrm3)
plot(log_days16C4)


# load nitrate predictors ----

days10n <- rast(paste(re.dir, "Nitrate", "Days_10N.tif", sep ='/'))
days10n

# crop to NC --
days10n2 <- crop(days10n, extent(e.nc))
plot(days10n2[[1]])

# resample predictors to bathy ----
days10n3 <- resample(days10n2, log_vrm3)

# mask predictors to bathy ----
days10n4 <- mask(days10n3, log_vrm3)
plot(days10n4)


###

min_nit <- rast(paste(re.dir, "Nitrate", "Min_Monthly_Nitrate.tif", sep ='/'))
min_nit

# crop to NC --
min_nit2 <- crop(min_nit, extent(e.nc))
plot(min_nit2[[2]])

# resample predictors to bathy ----
min_nit3 <- resample(min_nit2, log_vrm3)

# mask predictors to bathy ----
min_nit4 <- mask(min_nit3, log_vrm3)
plot(min_nit4)


###



##



# stack and try to predict one year ----

preds1 <- c(log_vrm3, max_an_up_temp4[[1]], log_days16C4[[1]], days10n4[[1]], min_nit4[[1]])
names(preds1)

year1998 <- classify(min_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("log_mean_vrm", "Max_Monthly_Anomaly_Upwelling_Temp", "log_Days_16C"  , "Days_10N",                          
                   "Min_Monthly_Nitrate", 'year') 


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
names(preds3) <- c("log_mean_vrm", "Max_Monthly_Anomaly_Upwelling_Temp", "log_Days_16C"  , "Days_10N",                          
                   "Min_Monthly_Nitrate", 
                   "year",
                   "site_name") 

##

# zone ----

zone.raster <- year1998
names(zone.raster) <- 'zone'

zone.raster2 <- classify(zone.raster, cbind(0, Inf, 1), right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("log_mean_vrm", "Max_Monthly_Anomaly_Upwelling_Temp", "log_Days_16C"  , "Days_10N",                          
                   "Min_Monthly_Nitrate", 
                   "year",
                   "site_name",
                   "zone") 


##

# test predict ----

testp <- predict(preds4, mod)
plot(testp)

#testp2 <- approximate(testp, method = 'linear')

## aggragate with threshold ----
# mean.thresh <- function(x) {ifelse(length(which(is.na(x) == T)) >= 70, NA, mean(x, na.rm = T))} 
# 
# sd.thresh <- function(x) {ifelse(length(which(is.na(x) == T)) >= 70, NA, sd(x, na.rm = T))}


### LOOP ----

## loop to predict several years ----

# make list of years --
year.list <- paste(1998:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(min_nit4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 

zone.raster <- year.raster
year.raster <- classify(year.raster, cbind(0, Inf, 1), right=FALSE)
names(zone.raster) <- 'zone'

# make sites raster ----

# each site is latitude --
rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)
# make raster
site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year1998))
site.raster
plot(site.raster)
# extend to fit other rasters
site.raster2 <- extend(site.raster, year1998)
plot(site.raster2)


for (i in 1:length(year.list)) {
  
  # 1. stack predictors for that year
  preds1 <- c(log_vrm3, max_an_up_temp4[[i]], log_days16C4[[i]], days10n4[[i]], min_nit4[[i]])
  
  # 2. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("log_mean_vrm",
                     "Max_Monthly_Anomaly_Upwelling_Temp",
                     "log_Days_16C",
                     "Days_10N"  ,
                     "Min_Monthly_Nitrate",
                     "year"     ,
                     "zone" ,
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, mod)
  plot(year.prediction)
  
  # save
  writeRaster(year.prediction, paste(o2.dir, "sp_predictions", paste(year.no, "purple_urchins_NC.tif", sep = '_'), sep ='/'))
  
}

plot(year.prediction)

r1 <- rast(paste(o2.dir, "sp_predictions", "2007_purple_urchins_NC.tif", sep ='/'))
plot(r1)

plot(year.r)


####

####

## PREDICT NEREOCYSTIS ----

## get model again ----

gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 5, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 7, bs = "cr") + 
              s(mean_prob_of_rock, k = 8, bs = "cr") +
              s(MHW_Upwelling_Days, k = 7, bs = "cr") +
              #s(wh.95, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)


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

prob_rock <- rast(paste(sub.dir, "prob_rock_nc.all_300m_wInterp.tif", sep ='/'))
prob_rock



# load temperature predictors ----


Days_16C <- rast(paste(re.dir, 'Temperature', "Days_16C.tif", sep ='/'))
Days_16C

# transform to log --
days_16C1 <- Days_16C + 1

log_days_16C <- log(days_16C1)

# crop to NC --
# Days_16C2 <- crop(Days_16C, extent(prob_rock))
# plot(min_temp2[[1]])

# resample predictors to bathy ----
Days_16C3 <- resample(log_days_16C, prob_rock)

# mask predictors to bathy ----
Days_16C4 <- mask(Days_16C3, prob_rock)
plot(Days_16C4)

###

# MHW 

mhw_up <- rast(paste(re.dir, 'Temperature', "MHW_Upwelling_Days.tif", sep ='/'))
mhw_up

# resample predictors to bathy ----
mhw_up3 <- resample(mhw_up, prob_rock)

# mask predictors to bathy ----
mhw_up4 <- mask(mhw_up3, prob_rock)
plot(mhw_up4)





# load nitrate predictors ----

max_nit <- rast(paste(re.dir, "Nitrate", "Max_Monthly_Nitrate.tif", sep ='/'))
max_nit

# # crop to NC --
# mean_sum_nit2 <- crop(mean_sum_nit, extent(e.nc))
# plot(mean_sum_nit2[[1]])

# resample predictors to bathy ----
max_nit3 <- resample(max_nit, prob_rock)

# mask predictors to bathy ----
max_nit4 <- mask(max_nit3, prob_rock)
plot(max_nit4)


# load purple urchin predictions ----
urch.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Spatio_temporal_GAMs/outputs_nc_rcca/gam_urchins2/sp_predictions"

# load raster data --

u.files <- dir(urch.dir)
u.files <- list.files(urch.dir, pattern = '.tif')
u.files
length(u.files)




##



# stack and try to predict one year ----

preds1 <- c(prob_rock, Days_16C4[[1]], mhw_up4[[1]], max_nit4[[1]])
names(preds1)

year1998 <- classify(Days_16C4[[1]], cbind(0, Inf, 1998), right=FALSE)
plot(year1998)
names(year1998) <- 'year'

year.list <- paste(1998:2021)
length(year.list)

preds2 <- c(preds1, year1998)
names(preds2)

names(preds2) <- c("mean_prob_of_rock",  "Log_Days_16C", "MHW_Upwelling_Days"  , "Max_Monthly_Nitrate",                          
                   'year') 


##


# sites ----

rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)

site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year1998))
site.raster

plot(site.raster)

ext(year1998)
ext(site.raster)

site.raster2 <- extend(site.raster, year1998)

preds3 <- c(preds2, site.raster2)
names(preds3) <- c("mean_prob_of_rock",  "Log_Days_16C", "MHW_Upwelling_Days"  , "Max_Monthly_Nitrate",                          
                   "year"     ,
                   "site_name") 

##

# zone ----

zone.raster <- year1998
names(zone.raster) <- 'zone'

zone.raster2 <- classify(zone.raster, cbind(0, Inf, 1), right=FALSE)
plot(zone.raster2)

preds4 <- c(preds3, zone.raster2)
names(preds4)

names(preds4) <- c("mean_prob_of_rock",  "Log_Days_16C", "MHW_Upwelling_Days"  , "Max_Monthly_Nitrate",                          
                   "year"     ,
                   "site_name",
                   "zone")


##



### LOOP ----

## loop to predict several years ----

# make list of years --
year.list <- paste(1998:2021)
length(year.list)

# make template raster of year ----
year.raster <- classify(prob_rock, cbind(0, Inf, 1998), right=FALSE)
plot(year.raster)
names(year.raster) <- 'year'

# make zone raster ---- 

zone.raster <- year.raster
year.raster <- classify(year.raster, cbind(0, Inf, 1), right=FALSE)
names(zone.raster) <- 'zone'

# make sites raster ----

# each site is latitude --
rdf <- as.data.frame(year1998, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)
# make raster
site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year1998))
site.raster
# extend to fit other rasters
site.raster2 <- extend(site.raster, year1998)

# outputs dir ----

preds.dir <- paste(o.dir, "gam_V3", "sp_predictions", sep ='/')


for (i in 1:length(year.list)) {
  
  # 1. get urchins
  urchin.rast <- rast(paste(urch.dir, u.files[i], sep ='/'))
  urchin.rast2 <- resample(urchin.rast, prob_rock)
  
  # 2. stack with predictors for that year
  env.raster <- c(prob_rock, Days_16C4[[i]], mhw_up4[[i]], max_nit4[[i]])
  
  preds1 <- c(urchin.rast2, env.raster)
  
  # 3. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(0, Inf, year.no), right=FALSE)
  
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # name predictors 
  names(preds4) <- c("log_den_STRPURAD",
                     "mean_prob_of_rock",  
                     "log_Days_16C", 
                     "MHW_Upwelling_Days"  , 
                     "Max_Monthly_Nitrate",                          
                     "year"     ,
                     "zone",
                     "site_name")
  
  # 5. predict
  year.prediction <- predict(preds4, nereo.mod)
  plot(year.prediction)
  
  # 6. Change the log
  year.prediction2 <- exp(year.prediction)
  
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

###

###

###

## For July 13th


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
    UBR_Mean, UBRYear_Mean, UBRYear_Max, UBR_Max,
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
         log_UBR_Max = log(UBR_Max + 1),
         log_UBRYear_Mean = log(UBRYear_Mean + 1),
         log_UBRYear_Max = log(UBRYear_Max + 1)) %>%
  #dplyr::select(-c(UBR_Mean)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
  glimpse() # Rows: 708


# Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 507

## Divide Pre and Post ----
levels(dat2$year)
names(dat2)


# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat2$log_den_NERLUE, p = 0.8, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]


# run gam ----

gam1 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 4, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 6, bs = "cr") + 
              s(Mean_Monthly_NPP, k = 8, bs = "cr") + 
              #s(mean_prob_of_rock, k = 6, bs = "cr") +
              s(log_UBR_Max, k = 6, bs = "cr") +
              #s(mean_waveyear, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



summary(gam1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

gam2 <- gam(formula = log_den_NERLUE ~ s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 4, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 6, bs = "cr") + 
              #s(Mean_Monthly_NPP, k = 8, bs = "cr") + 
              #s(mean_prob_of_rock, k = 6, bs = "cr") +
              s(log_UBR_Max, k = 6, bs = "cr") +
              s(mean_waveyear, k = 6, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = train.gam, method = "REML")



summary(gam2)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2)
dev.off()



gam1$aic
gam1$deviance
summary(gam1)
gam.check(gam1)



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
testdata <- expand.grid(log_den_PYCHEL=mean(mod$model$log_den_PYCHEL),
                        Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
                        #mean_logvrm=mean(mod$model$mean_logvrm),
                        Mean_Monthly_Summer_Temp=mean(mod$model$Mean_Monthly_Summer_Temp),
                        Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
                        Min_Monthly_Nitrate=mean(mod$model$Min_Monthly_Nitrate),
                        site_name=(mod$model$site_name),
                        zone=(mod$model$zone),
                        year=(mod$model$year))%>%
  distinct()%>%
  glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(log_den_NERLUE, 
                log_Days_16C, 
                log_den_STRPURAD,
                Max_Monthly_Nitrate, 
                Mean_Monthly_NPP, 
                log_UBR_Max,
                mean_waveyear,
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
  scale_fill_viridis(option = "A", discrete = T) +
  #   ggtitle(substitute(italic(name)))+
  # scale_fill_manual(#labels = c("2010", "2011", "2013", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange", "green", "purple")) +
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

###

###


###

###

## For 20 July meeting ----

gam1 <- gam(formula = log_den_NERLUE ~ 
              #s(log_Days_16C, k = 3, bs = "cr") + 
              s(log_den_STRPURAD, k = 8, bs = "cr") + 
              s(Max_Monthly_Nitrate, k = 4, bs = "cr") + 
              s(Mean_Monthly_Upwelling_Temp, k = 3, bs = "cr") +
              #Mean_Monthly_Upwelling_Temp +
              #s(mean_prob_of_rock, k = 4, bs = "cr") +
              #s(mean_depth, k = 4, bs = "cr") +
              #s(mean_waveyear, k = 6, bs = "cr") +
              s(log_UBR_Max, k = 3, bs = "cr") +
              s(Mean_Monthly_NPP, k = 6, bs = "cr") +
              #s(Mean_Monthly_Summer_Nitrate, k = 6, bs = "cr") +
              s(wh_max, k = 4, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = tw(), data = dat2, method = "REML")




summary(gam1)
gam.check(gam1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

## Test ----
library(ggpmisc)

fits <- predict.gam(gam1, newdata=test.gam, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  test.gam %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = log_den_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1, color = 'red'), show.legend = F) +
  labs(x='Predicted', y='Observed', title='Log density N. luetkeana') +
  theme_bw()
p

## bar plot ----

predicts.year = test.gam%>%data.frame(fits)%>% #glimpse()
  group_by(year)%>% #only change here
  summarise(response=mean(fit),
            se.fit=mean(se.fit),
            observed = mean(log_den_NERLUE, na.rm = T),
            sd.observed = sd(log_den_NERLUE, na.rm = T),
            n.observed = length(log_den_NERLUE),
            se.observed = sd.observed/(sqrt(n.observed)))%>%
  dplyr::filter(year != "2006") %>%
  ungroup() %>%
  glimpse()

ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  scale_fill_viridis(discrete =T)+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  geom_line(aes(x=year, y= observed), group = 1, color = 'red', size = 1.5) +
  geom_errorbar(aes(x=year,y=observed, ymin = observed-se.observed,ymax = observed+se.observed),width = 0.5, col = 'red') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

ggmod.year 

