###
### Script based on one created by : Anita Giraldo on 18 July 2022 2022
### Script last updated by : Anita Giraldo on 18 July 2022 with added orb velocity and NPP from Tom

## This script calculates the percentage maximum for each site and runs gams accordingly --



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
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 64


## Load RCCA data ----

df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse() # Rows: 1,154 rows

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

# get just kelp and factors ----

just.kelp <- df.nc %>%
  dplyr::select(site_name, year, transect, zone, den_NERLUE) %>%
  glimpse() # 708 rows


# Calculate proportion of maximum ----

# get maximum per site and zone ----
mean.site <- just.kelp %>%
  group_by(year, site_name, zone) %>%
  summarise(mean.site.zone = mean(den_NERLUE, na.rm = T)) %>%
  ungroup() %>%
  glimpse() # Rows: 236

max.site <- mean.site %>%
  group_by(site_name, zone) %>%
  summarise(max.site.zone = max(mean.site.zone, na.rm = T)) %>%
  glimpse() # Rows: 20


prop.max <- mean.site %>%
  left_join(max.site, by = c("site_name", "zone")) %>%
  mutate(prop.max_NERLUE = mean.site.zone/max.site.zone) %>%
  #dplyr::select(-den_NERLUE) %>%
  glimpse() # Rows: 236

prop.max2 <- just.kelp %>%
  left_join(prop.max, by = c("site_name", "year", "zone")) %>%
  dplyr::select(-den_NERLUE) %>%
  glimpse()


prop.max3 <- prop.max2

prop.max3$prop.max_NERLUE[prop.max3$prop.max_NERLUE == 1] <- 0.999
prop.max3$prop.max_NERLUE[prop.max3$prop.max_NERLUE == 0] <- 0.001

glimpse(prop.max3)


## joing back with df.nc ----

df.nc2 <- prop.max3 %>%
  left_join(df.nc, by = c("site_name", "year", "transect", "zone")) %>%
  glimpse()



## Choose variables and transform needed ----
names(df.nc2)

dat1 <- df.nc2 %>%
  dplyr::select(
    # Factors 
    latitude, longitude,
    site_name, year, transect, zone,
    # Bio vars
    prop.max_NERLUE, den_NERLUE , den_MESFRAAD , den_STRPURAD , den_PYCHEL, den_HALRUF,
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
         log_prop.max_NERLUE = log(prop.max_NERLUE + 0.001),
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
  glimpse() # Rows: 507



###

###


### GAM V1 ----

# https://stats.stackexchange.com/questions/391912/how-to-fit-a-longitudinal-gam-mixed-model-gamm
# https://github.com/beckyfisher/FSSgam/blob/master/case_study2_soft_sediment.R


# 1. Select predictors for this GAM ----
names(dat2)


# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat2$log_prop.max_NERLUE, p = 0.75, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]

#### Select predictors for this GAMs ----

names(train.gam)



# 3. Set parameters to save outputs ----

name <- 'V1'
o2.dir <- paste(o.dir, paste("gam_prop_max", name, sep = '_'), sep ='/')
o2.dir

names(train.gam)

# 4. Define predictor variables ----

pred.vars <- c("Min_Monthly_Nitrate" ,
               "Max_Monthly_Nitrate"  ,             
               "Mean_Monthly_Nitrate" ,
               "Mean_Monthly_Upwelling_Nitrate",
               "Mean_Monthly_Summer_Nitrate" ,     
               "Mean_Monthly_Temp"    ,
               "Mean_Monthly_Summer_Temp"  ,
               "MHW_Upwelling_Days"   ,
               "Max_Monthly_Anomaly_Upwelling_Temp",
               "Min_Monthly_Temp"  ,                
               "Mean_Monthly_Upwelling_Temp",
               "mean_depth" ,
               "mean_prob_of_rock" ,
               "wh_max" ,
               "wh_mean" ,
               "mean_waveyear"   ,
               "wh_95prc"   ,
               "Mean_Monthly_NPP" ,
               "Max_Monthly_NPP_Upwelling"  ,
               "log_den_MESFRAAD"  ,
               "log_den_STRPURAD" ,                    
               "log_UBR_Mean"  ,
               "log_Mean_Monthly_NPP_Upwelling" ,
               "log_Min_Monthly_NPP"
               )

length(pred.vars) # 24



# 5. Define Null model ----

#fact.vars <- c("survey_year")

model.v1 <- gam(prop.max_NERLUE ~ 
                  s(site_name, zone, bs = 're') +
                  s(year, bs = 're') ,
                data = train.gam, 
                family = betar(link = "logit", eps=.Machine$double.eps*100),
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
names(out.all) <- 'prop_max_Nereo'
names(var.imp) <- 'prop_max_Nereo'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)
#write.csv(all.mod.fits[,-2], file=paste(o.dir, paste(name,"all.mod.fits.csv",sep="_"), sep ='/'))
write.csv(mod.table, file = paste(o2.dir, "prop_max_all.mod.fits.csv", sep ='/'))
write.csv(out.i, file=paste(o2.dir, "prop_max_best_models.csv", sep="/"))
write.csv(all.var.imp, file=paste(o2.dir, "prop_max_all.var.imp.csv", sep="/"))


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

#### gam1 ----

subname <- "4"

best.model.name=as.character(out.i$modname[2])
best.model <- out.list$success.models[[best.model.name]]
best.model <- print(out.i$formula[2])

gam1 <- gam(formula = prop.max_NERLUE ~ 
              #s(log_den_PYCHEL, k = 5, bs = "cr") +
              s(log_den_STRPURAD, k = 5, bs = "cr") + 
              #s(Min_Monthly_Nitrate, k = 5, bs = "cr") +
              s(Max_Monthly_Nitrate, k = 5, bs = "cr") +
              #s(log_Days_16C, k = 5, bs = "cr") +
              #s(Mean_Monthly_Summer_Temp, k = 5, bs = "cr") + 
              #s(Mean_Monthly_Upwelling_Temp, k = 4, bs = "cr") +
              s(Mean_Monthly_NPP, k = 5, bs = "cr") +
              s(mean_waveyear, k = 5, bs = "cr") +
              s(log_UBR_Max, k = 5, bs = "cr") +
              #s(Min_Monthly_Nitrate, k = 3, bs = "cr") +
              s(site_name, zone, bs = "re") + 
              s(year, bs = "re"), 
            family = betar(), data = train.gam, method = "REML")



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
# pdf(file=paste(o2.dir, paste(name,subname,'log_den_NERLUE',"fits.pdf",sep="_"), sep ='/'))
# #par(mfrow=c(3,2),mar=c(9,4,3,1))
# par(mfrow=c(3,2),mar=c(2,4,3,1))
# plot(gam1,all.terms=T,pages=1,residuals=T,pch=16)
# mtext(side=2,text='log_den_NERLUE',outer=F)
# dev.off()

##

## PREDICT  gam ----
mod<-gam1
head(mod$model)
# testdata <- expand.grid(log_prop.max_NERLUE = mean(mod$model$prop.max_NERLUE),
#                         mean_waveyear=mean(mod$model$mean_waveyear),
#                         Max_Monthly_Nitrate=mean(mod$model$Max_Monthly_Nitrate),
#                         log_den_STRPURAD=mean(mod$model$log_den_STRPURAD),
#                         Mean_Monthly_NPP=mean(mod$model$Mean_Monthly_NPP),
#                         Mean_Monthly_Temp=mean(mod$model$Mean_Monthly_Temp),
#                         log_UBR_Max=mean(mod$model$log_UBR_Max),
#                         site_name=(mod$model$site_name),
#                         zone=(mod$model$zone),
#                         year=(mod$model$year))%>%
#   distinct()%>%
#   glimpse()

glimpse(test.gam)



testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                Mean_Monthly_Upwelling_Temp, 
                mean_waveyear,
                log_UBR_Max, 
                site_name, zone, year)

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)


predicts.year = testdata%>%data.frame(fits)%>% #glimpse()
  group_by(year)%>% #only change here
  summarise(response=mean(fit),
            se.fit=mean(se.fit),
            observed = mean(prop.max_NERLUE, na.rm = T),
            sd.observed = sd(prop.max_NERLUE, na.rm = T),
            n.observed = length(prop.max_NERLUE),
            se.observed = sd.observed/(sqrt(n.observed)))%>%
  ungroup() %>%
  glimpse()



# PLOTS for model for year  ----
ggmod.year <- ggplot(aes(x=year,y=response,fill=year), data=predicts.year) +
  ylab(" ")+
  xlab('survey_year')+
  scale_fill_viridis(discrete =T)+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  geom_line(aes(x=year, y= observed), group = 1, color = 'red', size = 1.5) +
  geom_errorbar(aes(x=year,y=observed, ymin = observed-se.observed,ymax = observed+se.observed),width = 0.5, col = 'red') +
  theme_classic()
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year 


## plot pred.vs obs ----

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = prop.max_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  labs(x='Predicted', y='Observed', title='Proportion of maximum density N. luetkeana') +
  theme_bw()
p

###


## Plot mean yearly kelp predicted vs observed ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()



predicts.year2 <- predicts.all %>%
  rename(log_fit = fit,
         log_se.fit = se.fit) %>%
  mutate(prop.max_NERLUE = (exp(log_prop.max_NERLUE) - 0.001),
         fit = (exp(log_fit) - 0.001),
         se.fit = (exp(log_se.fit) - 0.001)) %>%
  # mutate(den_NERLUE = (10^(log_den_NERLUE) - 1),
  #        fit = (10^(log_fit) - 1),
  #        se.fit = (10^(log_se.fit) - 1)) %>%
  glimpse()

predicts.year3 <- predicts.year2 %>%
  group_by(year) %>% # only change here
  summarise(response = mean(fit, na.rm = T),
            se.fit2 = mean(se.fit, na.rm = T),
            sd.fit = sd(fit, na.rm = T),
            n.fit = length(fit),
            se.fit = sd.fit/(sqrt(n.fit)),
            log_response = mean(log_fit, na.rm = T),
            log_se.fit = mean(log_se.fit, na.rm = T),
            observed = mean(prop.max_NERLUE, na.rm = T),
            sd.obs = sd(prop.max_NERLUE, na.rm = T),
            n.obs = length(prop.max_NERLUE),
            se.obs = sd.obs/(sqrt(n.obs)),
            log_observed = mean(log_prop.max_NERLUE, na.rm = T),
            log_sd.obs = sd(log_prop.max_NERLUE, na.rm = T),
            log_n.obs = length(log_prop.max_NERLUE),
            log_se.obs = log_sd.obs/(sqrt(log_n.obs))) %>%
  dplyr::filter(year != "2006") %>% # remove 2006 because only 1 observation
  ungroup() %>%
  glimpse()

# plot test in log scale
ggmod.year <- ggplot( data=predicts.year3) +
  geom_bar(aes(x=year,y=log_response,fill=year), stat = "identity")+
  scale_fill_viridis(discrete = T, direction = 1, option = "G") +
  geom_errorbar(aes(x=year,y=log_response, ymin = log_response-log_se.fit,ymax = log_response+log_se.fit),width = 0.5) +
  geom_line(aes(x=year, y= log_observed), group = 1, color = 'red', size = 1.5) +
  geom_errorbar(aes(x=year,y=log_observed, ymin = log_observed-log_se.obs,ymax = log_observed+log_se.obs),width = 0.5, col = 'red') +
  ylab("Log proportion of maximum density N. luetkeana")+
  xlab('survey_year')+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year

# plot test actual density
ggmod.year <- ggplot(data=predicts.year3) +
  geom_bar(aes(x=year,y=response,fill=year), stat = "identity")+
  scale_fill_viridis(discrete = T, direction = 1, option = "G") +
  geom_errorbar(aes(x=year,y=response, ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  geom_line(aes(x=year, y= observed), group = 1, color = 'red', size = 1.5) +
  geom_errorbar(aes(x=year,y=observed, ymin = observed-se.obs,ymax = observed+se.obs),width = 0.5, col = 'red') +
  ylab("Proportion of maximum density of N. luetkeana")+
  xlab('survey_year')+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
#Theme1+
#annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
#annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.year


# prop pred vs. observed ----
glimpse(predicts.all)

my.formula <- y ~ x


p <- ggplot(predicts.year2, aes(x = log_fit, y = log_prop.max_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               #aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               aes(label = paste(..rr.label..)),
               parse = TRUE, size = 4) +         
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1, color = 'red'), show.legend = F) +
  #scale_color_viridis(discrete = T) +
  labs(x='Predicted', y='Observed', title='Log proportion of maximum density N. luetkeana') +
  theme_classic()
p

p <- ggplot(predicts.year2, aes(x = fit, y = prop.max_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               #aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               aes(label = paste(..rr.label..)),
               parse = TRUE, size = 4) +         
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1, color = 'red'), show.legend = F) +
  #scale_color_viridis(discrete = T) +
  labs(x='Predicted', y='Observed', title='Proportion of maximum density N. luetkeana') +
  theme_classic()
p







