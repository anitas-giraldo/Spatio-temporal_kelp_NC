###
### Script created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 13 April 2022

## This script runs GAMs for North Coast using Landsat data

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

# Clear environment ----
rm(list=ls())

### Set directories ----
#w.dir<- dirname(rstudioapi::getActiveDocumentContext()$path)
m.dir <- here()
d.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.ls.outputs"
out.dir <- here("outputs_nc_ls")
p.dir <- paste(o.dir, 'plots', sep ='/')


## 1. Load data ----

df <- read.csv(paste(d.dir, "NC_subsamp_Landsat_99-21_merged_vars_noNAs.csv", sep = '/')) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  mutate(depth = depth*(-1)) %>% # make depth positive
  rename(kelp.area = max.area) %>%
  glimpse()


## 2. Divide data in PRE and POST MHW ----
levels(df$year)

## Choose pre MHW years --

preMHW <- df %>%
  dplyr::filter(year == "1999" |
                  year == "2000" |
                  year == "2001" |
                  year == "2002" |
                  year == "2003" |
                  year == "2004" |
                  year == "2005" |
                  year == "2006" |
                  year == "2007" |
                  year == "2008" |
                  year == "2009" |
                  year == "2010" |
                  year == "2011" |
                  year == "2012" |
                  year == "2013" ) %>%
  droplevels() %>%
  glimpse() # Rows: 2,232

## Choose post MHW years --

postMHW <- df %>%
  dplyr::filter(year == "2014" |
                  year == "2015" |
                  year == "2016" |
                  year == "2017" |
                  year == "2018" |
                  year == "2019" |
                  year == "2020" |
                  year == "2021" ) %>%
  droplevels() %>%
  glimpse() # Rows: 919

###

###

# Use PRE-MHW data to fit the model

# 3. Select all possible variables ----
glimpse(preMHW)

preMHW2 <- preMHW %>%
  dplyr::select(c( # factors and coords
    year,
    x,y,
    # dependent
    kelp.area,
    # substrate
    depth, slope, vrm, prob_of_rock,
    # other env
    npgo_mean, mei_mean,
    # nitrate
    Max_Monthly_Anomaly_Summer_Nitrate,
    Mean_Monthly_Summer_Nitrate,
    Min_Monthly_Anomaly_Nitrate,
    Min_Monthly_Anomaly_Upwelling_Nitrate,
    Max_Monthly_Anomaly_Nitrate,
    Days_3N,
    Max_Monthly_Anomaly_Upwelling_Nitrate, 
    Max_Monthly_Nitrate,
    Mean_Monthly_Upwelling_Nitrate,
    # temperature
    MHW_Intensity,
    Min_Monthly_Anomaly_Summer_Temp,
    Min_Monthly_Anomaly_Temp,
    Mean_Monthly_Temp,
    Days_16C,
    Mean_Monthly_Summer_Temp
  )) %>%
  glimpse() # Rows: 219


# 4. Transform necessary variables ----
glimpse(preMHW2)

preMHW3 <- preMHW2 %>%
  mutate(kelp.area_sqrt = sqrt(kelp.area),
         kelp.area_log = log(kelp.area + 1),
         depth_sqrt = sqrt(depth),
         slope_log = log(slope + 1), 
         vrm_sqrt = sqrt(vrm),
         prob_of_rock_log = log(prob_of_rock + 1),
         MHW_Intensity_log = log(MHW_Intensity + 1),
         Days_16C_log = log(Days_16C + 1)) %>%
  dplyr::select(-c(depth, slope, vrm, prob_of_rock,
                   MHW_Intensity, Days_16C)) %>%
  relocate(kelp.area_sqrt, .after = y) %>%
  glimpse()
         


## 

### GAM V1 ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

# 1. Select predictors for this GAM ----

names(preMHW3)

# get data frame --

dat.gam <- preMHW3 %>%
  dplyr::select(
    # factors
    year,
    # coords
    x,y,
    # dependent 
    #kelp.area_log,
    kelp.area,
    # environmental
    npgo_mean, mei_mean,
    # nitrate
    Mean_Monthly_Summer_Nitrate,
    Max_Monthly_Anomaly_Nitrate,
    # temperature
    Days_16C_log,
    MHW_Intensity_log,
    Min_Monthly_Anomaly_Summer_Temp,
    ) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse() # Rows: 478 - 24 columns

##

# 2. Divide data into train and test ----
inTraining <- createDataPartition(dat.gam$kelp.area, p = 0.75, list = FALSE)
train.gam <- dat.gam[ inTraining,]
test.gam  <- dat.gam[-inTraining,]

glimpse(train.gam) # Rows: 1,675
glimpse(test.gam) # Rows: 557

##

# 3. Set parameters to save outputs ----

name <- 'V1'
o.dir <- paste(out.dir, paste("gam", name, sep = '_'), sep ='/')

names(dat.gam)
glimpse(dat.gam)

##

# 4. Define predictor variables ---- 

pred.vars <- c(# environmental
  'npgo_mean', 'mei_mean',
  # nitrate
  'Mean_Monthly_Summer_Nitrate',
  'Max_Monthly_Anomaly_Nitrate',
  # temperature
  'Days_16C_log',
  'MHW_Intensity_log',
  'Min_Monthly_Anomaly_Summer_Temp'
  )

length(pred.vars) # 7

##

# 5. Define Null model ----

#fact.vars <- c("survey_year")

model.v1 <- gam(kelp.area ~ 
                  s(year, bs = 're'),
                data = train.gam, 
                family = tw(),
                method = "REML") 
##

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
                                null.terms = "s(year, bs = 're')")



# # To remove unwanted interactions : --

# length(model.set$mod.formula) # 12539
# test <- model.set$mod.formula
# remove.these <- which(str_detect(test, "depth_mean, wh.max, wh.95"))
# length(remove.these)
# # remove those
# model.set$mod.formula[remove.these] <- NULL
# #
# remove.these <- which(str_detect(test, "wh.max, wh.95"))
# length(remove.these)
# # remove those
# model.set$mod.formula[remove.these] <- NULL
# length(model.set$mod.formula) # 10710

## 

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
names(out.all) <- 'kelp.area'
names(var.imp) <- 'kelp.area'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(all.mod.fits, file = paste(o.dir, "all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o.dir, "best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o.dir, "all.var.imp.csv", sep="/"))


# 10. plot the best models ----
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[1])

  png(file=paste(o.dir, paste(name,m,'Kelp_area',"mod_fits.png",sep="_"), sep ='/'))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model <- out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text='logden_NERLUE',outer=F)}
  dev.off()
}



# 11. custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, "all.var.imp.csv", sep ='/')) %>% #from local copy
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

ggsave(namep, device = 'png', path = o.dir)

###

###

#### CHECK BEST MODELS ----

# Manually make the most parsimonious GAM models --

# 1.1 Gam 1.1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$formula[[best.model.name]]
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
# bm <- print(best.model$formula)
# bm.form <- as.formula(paste("kelp.area_sqrt ", paste(bm, collapse= " ")))
bm.form <- print(best.model$formula)
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = train.gam, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

# 1.2. Check gam 
gam.check(gam1)

# 1.3. Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o.dir, paste(name,subname,'kelp.area',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='kelp.area',outer=F)
dev.off()

##

# 1.4. PREDICT gam ----

# mod<-gam1
# head(mod$model)
# testdata <- expand.grid(Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
#                         Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
#                         MHW_Intensity_log=mean(mod$model$MHW_Intensity_log),
#                         Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
#                         npgo_mean=mean(mod$model$npgo_mean),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()
# 

glimpse(test.gam)

# 1.4.1. Get test data ----

testdata <- test.gam %>%
  dplyr::select(x, y, kelp.area,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity_log,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  glimpse()

glimpse(testdata)

length(testdata)

# 1.4.2. Predict on test data ----

fits <- predict.gam(mod, newdata=testdata[,4:9], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

predicts.all <- cbind(testdata, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,10), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# 1.7. PREDICT ON POST-MHW ----

glimpse(postMHW)

# 1.7.1  select and transform variables ----
testdata2 <- postMHW %>%
  dplyr::select(x, y, kelp.area,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  # transform needed
  mutate(MHW_Intensity_log = log(MHW_Intensity + 1)) %>%
  dplyr::select(- MHW_Intensity) %>%
  glimpse()


# 1.7.2. Predict on POST-MHW data ----
glimpse(testdata2)

fits <- predict.gam(mod, newdata=testdata2[,4:9], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata2, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

# save
namep <- paste(name, subname, "year_preds_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

# Plot latitudinal ----

predicts.all <- cbind(testdata2, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs_PostMHW.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,10), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

###

# 2.1. Gam 1.2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$formula[[best.model.name]]
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
# bm <- print(best.model$formula)
# bm.form <- as.formula(paste("kelp.area_sqrt ", paste(bm, collapse= " ")))
bm.form <- print(best.model$formula)
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = train.gam, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

# 1.2. Check gam 
gam.check(gam1)

# 1.3. Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o.dir, paste(name,subname,'kelp.area',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='kelp.area',outer=F)
dev.off()

##

# 1.4. PREDICT gam ----

# mod<-gam1
# head(mod$model)
# testdata <- expand.grid(Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
#                         Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
#                         MHW_Intensity_log=mean(mod$model$MHW_Intensity_log),
#                         Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
#                         npgo_mean=mean(mod$model$npgo_mean),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()
# 

glimpse(test.gam)

# 1.4.1. Get test data ----

testdata <- test.gam %>%
  dplyr::select(x, y, kelp.area,
                mei_mean,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                Days_16C_log,
                MHW_Intensity_log,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  glimpse()

glimpse(testdata)

length(testdata)

# 1.4.2. Predict on test data ----

fits <- predict.gam(mod, newdata=testdata[,4:11], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

predicts.all <- cbind(testdata, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,12), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# 1.7. PREDICT ON POST-MHW ----

glimpse(postMHW)

# 1.7.1  select and transform variables ----
testdata2 <- postMHW %>%
  dplyr::select(x, y, kelp.area,
                mei_mean,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                Days_16C,
                MHW_Intensity,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  # transform needed
  mutate(MHW_Intensity_log = log(MHW_Intensity + 1),
         Days_16C_log = log(Days_16C + 1)) %>%
  dplyr::select(- c(MHW_Intensity, Days_16C)) %>%
  glimpse()


# 1.7.2. Predict on POST-MHW data ----
glimpse(testdata2)
names(testdata2)

fits <- predict.gam(mod, newdata=testdata2[,4:11], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata2, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

# save
namep <- paste(name, subname, "year_preds_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

# Plot latitudinal ----

predicts.all <- cbind(testdata2, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs_PostMHW.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,12), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

###

## 

### GAM V2 ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

# 1. Select predictors for this GAM ----

names(preMHW3)

# get data frame --

dat.gam <- preMHW3 %>%
  dplyr::select(
    # factors
    year,
    # coords
    x,y,
    # dependent 
    #kelp.area_log,
    kelp.area,
    # environmental
    npgo_mean, mei_mean,
    # nitrate
    Mean_Monthly_Summer_Nitrate,
    Min_Monthly_Anomaly_Nitrate,
    # temperature
    Min_Monthly_Anomaly_Temp,
    Mean_Monthly_Temp,
    Days_16C_log,
    MHW_Intensity_log,
    Min_Monthly_Anomaly_Summer_Temp,
  ) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse() # Rows: 478 - 24 columns

##

# 2. Divide data into train and test ----
inTraining <- createDataPartition(dat.gam$kelp.area, p = 0.75, list = FALSE)
train.gam <- dat.gam[ inTraining,]
test.gam  <- dat.gam[-inTraining,]

glimpse(train.gam) # Rows: 1,675
glimpse(test.gam) # Rows: 557

##

# 3. Set parameters to save outputs ----

name <- 'V2'
o.dir <- paste(out.dir, paste("gam", name, sep = '_'), sep ='/')

names(dat.gam)
glimpse(dat.gam)

##

# 4. Define predictor variables ---- 

pred.vars <- c(# environmental
  'npgo_mean', 'mei_mean',
  # nitrate
  'Mean_Monthly_Summer_Nitrate',
  'Min_Monthly_Anomaly_Nitrate',
  # temperature
  'Days_16C_log',
  'MHW_Intensity_log',
  'Min_Monthly_Anomaly_Summer_Temp',
  'Min_Monthly_Anomaly_Temp',
  'Mean_Monthly_Temp'
)

length(pred.vars) # 9

##

# 5. Define Null model ----

#fact.vars <- c("survey_year")

model.v1 <- gam(kelp.area ~ 
                  s(year, bs = 're'),
                data = train.gam, 
                family = tw(),
                method = "REML") 
##

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
                                null.terms = "s(year, bs = 're')")



# # To remove unwanted interactions : --

# length(model.set$mod.formula) # 12539
# test <- model.set$mod.formula
# remove.these <- which(str_detect(test, "depth_mean, wh.max, wh.95"))
# length(remove.these)
# # remove those
# model.set$mod.formula[remove.these] <- NULL
# #
# remove.these <- which(str_detect(test, "wh.max, wh.95"))
# length(remove.these)
# # remove those
# model.set$mod.formula[remove.these] <- NULL
# length(model.set$mod.formula) # 10710

## 

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
names(out.all) <- 'kelp.area'
names(var.imp) <- 'kelp.area'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(all.mod.fits, file = paste(o.dir, paste(name, "all.mod.fits.csv", sep ='_'), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name, "best_models.csv", sep ='_'), sep="/"))
write.csv(all.var.imp, file=paste(o.dir, paste(name, "all.var.imp.csv", sep ='_'), sep="/"))


# 10. plot the best models ----
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[1])
  
  png(file=paste(o.dir, paste(name,m,'Kelp_area',"mod_fits.png",sep="_"), sep ='/'))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model <- out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text='logden_NERLUE',outer=F)}
  dev.off()
}



# 11. custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep ='_'), sep ='/')) %>% #from local copy
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

ggsave(namep, device = 'png', path = o.dir)

###

###

#### CHECK BEST MODELS ----

# Manually make the most parsimonious GAM models --

# 1.1 Gam 1.1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$formula[[best.model.name]]
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
# bm <- print(best.model$formula)
# bm.form <- as.formula(paste("kelp.area_sqrt ", paste(bm, collapse= " ")))
bm.form <- print(best.model$formula)
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = train.gam, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

# 1.2. Check gam 
par(mfrow=c(3,3),mar=c(2,4,3,1))
gam.check(gam1)

# 1.3. Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o.dir, paste(name,subname,'kelp.area',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='kelp.area',outer=F)
dev.off()

##

# 1.4. PREDICT gam ----

# mod<-gam1
# head(mod$model)
# testdata <- expand.grid(Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
#                         Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
#                         MHW_Intensity_log=mean(mod$model$MHW_Intensity_log),
#                         Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
#                         npgo_mean=mean(mod$model$npgo_mean),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()
# 

glimpse(test.gam)
gam1

# 1.4.1. Get test data ----

testdata <- test.gam %>%
  dplyr::select(x, y, kelp.area,
                mei_mean,
                #Min_Monthly_Anomaly_Nitrate,
                #Mean_Monthly_Summer_Nitrate,
                Min_Monthly_Anomaly_Temp,
                Mean_Monthly_Temp,
                Days_16C_log,
                MHW_Intensity_log,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  glimpse()

glimpse(testdata)

length(testdata)

# 1.4.2. Predict on test data ----

mod <- gam1

fits <- predict.gam(mod, newdata=testdata[,4:11], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

predicts.all <- cbind(testdata, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,12), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# 1.7. PREDICT ON POST-MHW ----

glimpse(postMHW)

# 1.7.1  select and transform variables ----
testdata2 <- postMHW %>%
  dplyr::select(x, y, kelp.area,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  # transform needed
  mutate(MHW_Intensity_log = log(MHW_Intensity + 1)) %>%
  dplyr::select(- MHW_Intensity) %>%
  glimpse()


# 1.7.2. Predict on POST-MHW data ----
glimpse(testdata2)

fits <- predict.gam(mod, newdata=testdata2[,4:9], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata2, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

# save
namep <- paste(name, subname, "year_preds_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

# Plot latitudinal ----

predicts.all <- cbind(testdata2, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs_PostMHW.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,10), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

###

# 2.1. Gam 1.2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$formula[[best.model.name]]
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
# bm <- print(best.model$formula)
# bm.form <- as.formula(paste("kelp.area_sqrt ", paste(bm, collapse= " ")))
bm.form <- print(best.model$formula)
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = train.gam, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

# 1.2. Check gam 
gam.check(gam1)

# 1.3. Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o.dir, paste(name,subname,'kelp.area',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='kelp.area',outer=F)
dev.off()

##

# 1.4. PREDICT gam ----

# mod<-gam1
# head(mod$model)
# testdata <- expand.grid(Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
#                         Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
#                         MHW_Intensity_log=mean(mod$model$MHW_Intensity_log),
#                         Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
#                         npgo_mean=mean(mod$model$npgo_mean),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()
# 

glimpse(test.gam)

# 1.4.1. Get test data ----

testdata <- test.gam %>%
  dplyr::select(x, y, kelp.area,
                mei_mean,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                Days_16C_log,
                MHW_Intensity_log,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  glimpse()

glimpse(testdata)

length(testdata)

# 1.4.2. Predict on test data ----

fits <- predict.gam(mod, newdata=testdata[,4:11], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

predicts.all <- cbind(testdata, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,12), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# 1.7. PREDICT ON POST-MHW ----

glimpse(postMHW)

# 1.7.1  select and transform variables ----
testdata2 <- postMHW %>%
  dplyr::select(x, y, kelp.area,
                mei_mean,
                #Min_Monthly_Anomaly_Nitrate,
                #Mean_Monthly_Summer_Nitrate,
                Min_Monthly_Anomaly_Temp,
                Mean_Monthly_Temp,
                Days_16C,
                MHW_Intensity,
                Min_Monthly_Anomaly_Summer_Temp,
                npgo_mean,
                year) %>%
  # transform needed
  mutate(MHW_Intensity_log = log(MHW_Intensity + 1),
         Days_16C_log = log(Days_16C + 1)) %>%
  dplyr::select(- c(MHW_Intensity, Days_16C)) %>%
  glimpse()


# 1.7.2. Predict on POST-MHW data ----
glimpse(testdata2) # 919 obs
names(testdata2) # Rows: 918

# remove outlier 
testdata2 <- testdata2 %>%
  dplyr::filter(kelp.area < 899) %>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata2[,4:11], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata2, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

# save
namep <- paste(name, subname, "year_preds_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

# Plot latitudinal ----

predicts.all <- cbind(testdata2, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs_PostMHW.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,12), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

###

## 

### GAM V3 ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R

# 1. Select predictors for this GAM ----

names(preMHW3)

# get data frame --

dat.gam <- preMHW3 %>%
  dplyr::select(
    # factors
    year,
    # coords
    x,y,
    # dependent 
    #kelp.area_log,
    kelp.area,
    # environmental
    #npgo_mean, mei_mean,
    # nitrate
    Mean_Monthly_Summer_Nitrate,
    Min_Monthly_Anomaly_Nitrate,
    # temperature
    Min_Monthly_Anomaly_Temp,
    Mean_Monthly_Temp,
    Days_16C_log,
    MHW_Intensity_log,
    Min_Monthly_Anomaly_Summer_Temp,
  ) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse() # Rows: 478 - 24 columns

##

# 2. Divide data into train and test ----
inTraining <- createDataPartition(dat.gam$kelp.area, p = 0.75, list = FALSE)
train.gam <- dat.gam[ inTraining,]
test.gam  <- dat.gam[-inTraining,]

glimpse(train.gam) # Rows: 1,675
glimpse(test.gam) # Rows: 557

##

# 3. Set parameters to save outputs ----

name <- 'V3'
o.dir <- paste(out.dir, paste("gam", name, sep = '_'), sep ='/')

names(dat.gam)
glimpse(dat.gam)

##

# 4. Define predictor variables ---- 

pred.vars <- c(# environmental
  #'npgo_mean', 'mei_mean',
  # nitrate
  'Mean_Monthly_Summer_Nitrate',
  'Min_Monthly_Anomaly_Nitrate',
  # temperature
  'Days_16C_log',
  'MHW_Intensity_log',
  'Min_Monthly_Anomaly_Summer_Temp',
  'Min_Monthly_Anomaly_Temp',
  'Mean_Monthly_Temp'
)

length(pred.vars) # 7

##

# 5. Define Null model ----

#fact.vars <- c("survey_year")

model.v1 <- gam(kelp.area ~ 
                  s(year, bs = 're'),
                data = train.gam, 
                family = tw(),
                method = "REML") 
##

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
                                null.terms = "s(year, bs = 're')")



# # To remove unwanted interactions : --

# length(model.set$mod.formula) # 12539
# test <- model.set$mod.formula
# remove.these <- which(str_detect(test, "depth_mean, wh.max, wh.95"))
# length(remove.these)
# # remove those
# model.set$mod.formula[remove.these] <- NULL
# #
# remove.these <- which(str_detect(test, "wh.max, wh.95"))
# length(remove.these)
# # remove those
# model.set$mod.formula[remove.these] <- NULL
# length(model.set$mod.formula) # 10710

## 

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
names(out.all) <- 'kelp.area'
names(var.imp) <- 'kelp.area'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(all.mod.fits, file = paste(o.dir, paste(name, "all.mod.fits.csv", sep ='_'), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(name, "best_models.csv", sep ='_'), sep="/"))
write.csv(all.var.imp, file=paste(o.dir, paste(name, "all.var.imp.csv", sep ='_'), sep="/"))


# 10. plot the best models ----
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[1])
  
  png(file=paste(o.dir, paste(name,m,'Kelp_area',"mod_fits.png",sep="_"), sep ='/'))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model <- out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text='logden_NERLUE',outer=F)}
  dev.off()
}



# 11. custom plot of importance scores----
dat.taxa <- read.csv(paste(o.dir, paste(name, "all.var.imp.csv", sep ='_'), sep ='/')) %>% #from local copy
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

ggsave(namep, device = 'png', path = o.dir)

###

###

#### CHECK BEST MODELS ----

# Manually make the most parsimonious GAM models --

# 1.1 Gam 1.1 ----

subname <- "1"

best.model.name=as.character(out.i$modname[1])
best.model <- out.list$formula[[best.model.name]]
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
# bm <- print(best.model$formula)
# bm.form <- as.formula(paste("kelp.area_sqrt ", paste(bm, collapse= " ")))
bm.form <- print(best.model$formula)
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = train.gam, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

# 1.2. Check gam 
par(mfrow=c(3,3),mar=c(2,4,3,1))
gam.check(gam1)

# 1.3. Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o.dir, paste(name,subname,'kelp.area',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='kelp.area',outer=F)
dev.off()

##

# 1.4. PREDICT gam ----

# mod<-gam1
# head(mod$model)
# testdata <- expand.grid(Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
#                         Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
#                         MHW_Intensity_log=mean(mod$model$MHW_Intensity_log),
#                         Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
#                         npgo_mean=mean(mod$model$npgo_mean),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()
# 

glimpse(test.gam)
gam1

# 1.4.1. Get test data ----

testdata <- test.gam %>%
  dplyr::select(x, y, kelp.area,
                #mei_mean,
                #Min_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity_log,
                Min_Monthly_Anomaly_Summer_Temp,
                Min_Monthly_Anomaly_Temp,
                #npgo_mean,
                year) %>%
  glimpse()

glimpse(testdata)

length(testdata)

# 1.4.2. Predict on test data ----

mod <- gam1

fits <- predict.gam(mod, newdata=testdata[,4:8], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

predicts.all <- cbind(testdata, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,9), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# 1.7. PREDICT ON POST-MHW ----

glimpse(postMHW)

# 1.7.1  select and transform variables ----
testdata2 <- postMHW %>%
  dplyr::select(x, y, kelp.area,
                #mei_mean,
                #Min_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity,
                Min_Monthly_Anomaly_Summer_Temp,
                Min_Monthly_Anomaly_Temp,
                #npgo_mean,
                year) %>%
  # transform needed
  mutate(MHW_Intensity_log = log(MHW_Intensity + 1)) %>%
  dplyr::select(- MHW_Intensity) %>%
  glimpse()


# 1.7.2. Predict on POST-MHW data ----
glimpse(testdata2)

fits <- predict.gam(mod, newdata=testdata2[,4:8], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata2, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

# save
namep <- paste(name, subname, "year_preds_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

# Plot latitudinal ----

predicts.all <- cbind(testdata2, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs_PostMHW.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,9), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

###

# 2.1. Gam 1.2 ----

subname <- "2"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$formula[[best.model.name]]
best.model <- out.list$success.models[[best.model.name]]
best.model

# get formula of best model --
# bm <- print(best.model$formula)
# bm.form <- as.formula(paste("kelp.area_sqrt ", paste(bm, collapse= " ")))
bm.form <- print(best.model$formula)
class(bm.form)

# run gam of best model --
gam1 <- gam(formula = bm.form, family = tw(), data = train.gam, method = "REML")

gam1$aic
gam1$deviance
summary(gam1)

# 1.2. Check gam 
gam.check(gam1)

# 1.3. Plot model results ----

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam1)
dev.off()

##

par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1)
dev.off()

# save plot --
pdf(file=paste(o.dir, paste(name,subname,'kelp.area',"fits.pdf",sep="_"), sep ='/'))
#par(mfrow=c(3,2),mar=c(9,4,3,1))
par(mfrow=c(3,2),mar=c(2,4,3,1))
plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
mtext(side=2,text='kelp.area',outer=F)
dev.off()

##

# 1.4. PREDICT gam ----

# mod<-gam1
# head(mod$model)
# testdata <- expand.grid(Max_Monthly_Anomaly_Nitrate=mean(mod$model$Max_Monthly_Anomaly_Nitrate),
#                         Mean_Monthly_Summer_Nitrate=mean(mod$model$Mean_Monthly_Summer_Nitrate),
#                         MHW_Intensity_log=mean(mod$model$MHW_Intensity_log),
#                         Min_Monthly_Anomaly_Summer_Temp=mean(mod$model$Min_Monthly_Anomaly_Summer_Temp),
#                         npgo_mean=mean(mod$model$npgo_mean),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()
# 

glimpse(test.gam)

# 1.4.1. Get test data ----
mod <- gam1
mod

testdata <- test.gam %>%
  dplyr::select(x, y, kelp.area,
                Days_16C_log,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity_log,
                Min_Monthly_Anomaly_Summer_Temp,
                Min_Monthly_Anomaly_Temp,
                year) %>%
  glimpse()

glimpse(testdata)

length(testdata)

# 1.4.2. Predict on test data ----

fits <- predict.gam(mod, newdata=testdata[,4:9], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

predicts.all <- cbind(testdata, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,10), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


###

# 1.7. PREDICT ON POST-MHW ----

glimpse(postMHW)

# 1.7.1  select and transform variables ----
testdata2 <- postMHW %>%
  dplyr::select(x, y, kelp.area,
                Days_16C,
                Mean_Monthly_Summer_Nitrate,
                MHW_Intensity,
                Min_Monthly_Anomaly_Summer_Temp,
                Min_Monthly_Anomaly_Temp,
                year) %>%
  # transform needed
  mutate(MHW_Intensity_log = log(MHW_Intensity + 1),
         Days_16C_log = log(Days_16C + 1)) %>%
  dplyr::select(- c(MHW_Intensity, Days_16C)) %>%
  glimpse()


# 1.7.2. Predict on POST-MHW data ----
glimpse(testdata2) # 919 obs
names(testdata2) # Rows: 918

# remove outlier 
testdata2 <- testdata2 %>%
  dplyr::filter(kelp.area < 899) %>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata2[,4:9], type='response', se.fit=T)
fits <- as.data.frame(fits)
glimpse(fits)

predicts.year <- cbind(testdata2, fits) %>%
  group_by(year) %>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit)) %>%
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
  # scale_fill_manual(labels = c("2010", "2011", "2016", "2017", "2018", "2019", "2020"), 
  #                   values=c("light blue", "blue", "dark blue", "black", "dark red", "red", "orange")) +
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

# save
namep <- paste(name, subname, "year_preds_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)

###

# Plot latitudinal ----

predicts.all <- cbind(testdata2, fits) %>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

write.csv(predicts.all,paste(o.dir, paste(name, subname, "preds-obs_PostMHW.csv", sep = '-'), sep ='/'))

# ner.obs <- dat.gam %>%
#   dplyr::select(year, kelp.area) %>%
#   glimpse()
# 
# pred.obs.all <- predicts.all %>%
#   left_join(ner.obs, by = c('year')) %>%
#   glimpse()
# 
# write.csv(pred.obs.all ,paste(o.dir, paste(name, subname, "preds-obs.csv", sep = '-'), sep ='/'))

# 1.5. plot predicted vs. observed ----
pred.obs.all <- predicts.all

library(ggpmisc)

p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +        
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = kelp.area)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Area of N. luetkeana') +
  theme_bw()
p

namep <- paste(name, subname, "pred-obs_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)


# 1.6 Plot latitudinally ----
glimpse(pred.obs.all)
names(pred.obs.all)

pred.obs.all2 <- pred.obs.all %>%
  pivot_longer(cols = c(3,10), names_to = "type", values_to = "values") %>%
  mutate_at(vars(type), list(as.factor)) %>%
  glimpse()

## not a map
ggplot(pred.obs.all2, aes(x = year, y = y, color = values)) +
  geom_jitter(size = 2.5, pch = 16) +
  labs(x = 'year', y = 'Latitude') +
  facet_wrap(~type) +
  scale_color_viridis(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

namep <- paste(name, subname, "pred-obs_latitude_PostMHW.png", sep ='_')
ggsave(namep, plot = last_plot(), device = 'png', path = o.dir)
