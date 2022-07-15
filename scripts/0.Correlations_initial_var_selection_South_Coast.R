### test 

### Script created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 23 March 2022

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
  glimpse() # Rows: 64


## Load transect data ----

df <- read.csv(paste(merged.data, "swath_upc_fish_terr_clim_env_depth_seaotter_temp_nitrate_transect.csv", sep ='/')) %>%
  mutate_at(vars(site_campus_unique_ID), list(as.factor)) %>%
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



## Remove variables that are not necessary ----
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


## remove correlated substrate terrain variables ----

dat2 <- dat %>%
  # transform the variable
  mutate(mean_logvrm = log(mean_vrm)) %>%
  # remove untransformed 
  dplyr::select(-mean_vrm) %>%
  # remove unselected variables from corr analysis
  dplyr::select(-c(
    # all medians
    #median_slope, median_depth, median_vrm, median_prob_of_rock,
    # all ranges
    #range_slope, range_depth, range_vrm, range_prob_of_rock,
    # mins
    #min_slope, min_vrm, 
    # sd
    #sd_slope, sd_vrm,
    # max
    #max_slope, max_vrm,
    # proportions mapped
    #prop_map_prob_of_rock, prop_map_vrm, prop_map_slope,
    # max prob rock is always 1
    #max_prob_of_rock,
    # days above 17C always 0
    #Days_22C, Degree_Days_22C,
    # the same as Max_Monthly_Anomaly_Temp
    #Max_Monthly_Anomaly_Summer_Temp,
    # the same as Min_Monthly_Anomaly_Temp
    #Min_Monthly_Anomaly_Upwelling_Temp,
    # MHW Days are the same as Intensity
    #MHW_Intensity, MHW_Summer_Intensity, MHW_Upwelling_Intensity,
    # remove other correlated vars
    #mei_mean, pdo_mean, 
    #Max_Monthly_Anomaly_Temp,
    #sd_depth, sd_prob_of_rock
  )) %>% # 
  glimpse()


###

###


# Substrate distribution of variables  ----

names(dat2)

predictors <- dat2[,c(6:85)]

# nitrate --
sub.preds <- c("depth_mean", "mean_slope"  , "mean_prob_of_rock",                    
               "sd_slope", "sd_depth", "sd_prob_of_rock",                      
               "sd_vrm" , "min_slope", "min_depth",                            
               "min_prob_of_rock" , "min_vrm" ,  "max_slope" ,                          
               "max_depth" , "max_prob_of_rock" , "max_vrm",                              
               "range_slope", "range_depth" , "range_prob_of_rock" ,                  
               "range_vrm", "median_slope",  "median_depth",                         
               "median_prob_of_rock" , "median_vrm" , "prop_map_slope",                       
               "prop_map_prob_of_rock" , "prop_map_vrm")

pred.sub <- predictors %>%
  dplyr::select(sub.preds) %>% #glimpse()
  pivot_longer(cols = depth_mean:prop_map_vrm, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()


psub <- pred.sub %>%
  #mutate(text = fct_reorder(Variable, Value)) %>%
  ggplot( aes(x=Value, color=Variable, fill=Variable)) +
  #geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
  geom_histogram(position="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~Variable, scales = 'free') +
  #theme_ipsum() +
  # theme(
  #   legend.position="none",
  #   panel.spacing = unit(0.1, "lines"),
  #   strip.text.x = element_text(size = 8)
  # ) +
  xlab("") +
  ylab("Frequency") + 
  theme(legend.position="bottom") 

psub

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in nit.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}




## Check correlations ----

## Substrate Only 1 ----

names(dat2)

predictors <- dat2 %>%
  dplyr::select("depth_mean", "mean_slope"  , "mean_prob_of_rock",                    
                "sd_slope", "sd_depth", "sd_prob_of_rock",                      
                "sd_vrm" , "min_slope", "min_depth",                            
                "min_prob_of_rock" , "min_vrm" ,  "max_slope" ,                          
                "max_depth" , "max_prob_of_rock" , "max_vrm",                              
                "range_slope", "range_depth" , "range_prob_of_rock" ,                  
                "range_vrm", "median_slope",  "median_depth",                         
                "median_prob_of_rock" , "median_vrm" , "prop_map_slope",                       
                "prop_map_prob_of_rock" , "prop_map_vrm") %>%
  glimpse()

names(predictors)

length(names(predictors)) # 23

length(which(is.na(predictors))) # 58
predictors <- na.omit(predictors)
nrow(predictors) #  4030
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C

# compute the p-value of correlations --
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


## Substrate Only 2 ----

# remove: some vars 

names(predictors)

predictors <- dat2 %>%
  dplyr::select("depth_mean", "mean_slope"  , "mean_prob_of_rock",                    
                #"sd_slope", 
                #"sd_depth", 
                "sd_prob_of_rock",                      
                #"sd_vrm" , 
                #"min_slope", 
                "min_depth",                            
                "min_prob_of_rock" , 
                "min_vrm" ,  
                #"max_slope" ,                          
                "max_depth" , "max_prob_of_rock" , 
                #"max_vrm",                              
                #"range_slope", 
                #"range_depth" , 
                #"range_prob_of_rock" ,                  
                #"range_vrm", 
                #"median_slope", 
                "median_depth"                         
                #"median_prob_of_rock" , 
                #"median_vrm" , "prop_map_slope",                       
                #"prop_map_prob_of_rock" , 
                #"prop_map_vrm"
                ) %>%
  glimpse()

names(predictors)

length(names(predictors)) # 10

length(which(is.na(predictors))) # 58
predictors <- na.omit(predictors)
nrow(predictors) #  4030
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.7, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


###


## Nitrate - distribution of variables  ----


names(dat2)

#predictors <- dat2[,c(35:77)]

# nitrate --
nit.preds <- c("Days_2N", "Days_3N", "Days_3N",
               "Max_Monthly_Anomaly_Nitrate" , 
               "Max_Monthly_Anomaly_Summer_Nitrate",
               "Max_Monthly_Anomaly_Upwelling_Nitrate", 
               "Max_Monthly_Nitrate"  , 
               "Mean_Monthly_Nitrate" ,
               "Mean_Monthly_Summer_Nitrate" ,          
               "Mean_Monthly_Upwelling_Nitrate" ,      
               "Min_Monthly_Anomaly_Nitrate" ,         
               "Min_Monthly_Anomaly_Summer_Nitrate" ,   
               "Min_Monthly_Anomaly_Upwelling_Nitrate", 
               "Min_Monthly_Nitrate"  ,
               "Min_Monthly_Temp_Nitrate" )

predictors <- pred.nit <- dat2 %>%
  dplyr::select(nit.preds) %>% 
  glimpse()

pred.nit <- dat2 %>%
  dplyr::select(nit.preds) %>% #glimpse()
  pivot_longer(cols = Days_2N:Min_Monthly_Temp_Nitrate, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()


pnit <- pred.nit %>%
  #mutate(text = fct_reorder(Variable, Value)) %>%
  ggplot( aes(x=Value, color=Variable, fill=Variable)) +
  #geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
  geom_histogram(position="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~Variable, scales = 'free') +
  #theme_ipsum() +
  # theme(
  #   legend.position="none",
  #   panel.spacing = unit(0.1, "lines"),
  #   strip.text.x = element_text(size = 8)
  # ) +
  xlab("") +
  ylab("Frequency") + 
  theme(legend.position="bottom") 

pnit

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,3))
for (i in nit.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}






###




## Check correlations ----
 
## Nitrate Only 1 ----

names(dat2)

predictors <- dat2 %>%
  dplyr::select(Days_1N, Days_2N, Days_3N,                        
                Days_4N , Degree_Days_1N , Degree_Days_2N ,Degree_Days_3N,                       
                Degree_Days_4N,  Max_Monthly_Anomaly_Nitrate , Max_Monthly_Anomaly_Summer_Nitrate,
                Max_Monthly_Anomaly_Upwelling_Nitrate, Max_Monthly_Nitrate, 
                Mean_Monthly_Nitrate,                 
                Mean_Monthly_Summer_Nitrate , Mean_Monthly_Upwelling_Nitrate, Min_Monthly_Anomaly_Nitrate,          
                Min_Monthly_Anomaly_Summer_Nitrate, Min_Monthly_Anomaly_Upwelling_Nitrate, Min_Monthly_Nitrate) %>%
  glimpse()

names(predictors)

length(names(predictors)) # 19

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #  460
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C

# compute the p-value of correlations --
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


## Nitrate Only 2 ----

# remove: some vars 

names(predictors)

predictors <- dat2 %>%
  dplyr::select(#Days_1N, Days_2N, Days_3N,                        
                Days_4N , #Degree_Days_1N , Degree_Days_2N ,Degree_Days_3N,                       
                #Degree_Days_4N,  
                #Max_Monthly_Anomaly_Nitrate, Max_Monthly_Anomaly_Summer_Nitrate,
                #Max_Monthly_Anomaly_Upwelling_Nitrate, 
                Max_Monthly_Nitrate, 
                #Mean_Monthly_Nitrate,                 
                Mean_Monthly_Summer_Nitrate , 
                #Mean_Monthly_Upwelling_Nitrate, Min_Monthly_Anomaly_Nitrate,          
                #Min_Monthly_Anomaly_Summer_Nitrate, 
                Min_Monthly_Anomaly_Upwelling_Nitrate, Min_Monthly_Nitrate) %>%
  glimpse()

length(names(predictors)) # 43

length(which(is.na(predictors))) # 18
predictors <- na.omit(predictors)
nrow(predictors) #  460
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##



##

## Temperature - distribution of variables  ----

names(df.sc)
names(dat2)

# nitrate --
var.preds <- c("Days_22C", "Max_Monthly_Anomaly_Summer_Temp","Max_Monthly_Anomaly_Temp" ,            
               "Max_Monthly_Anomaly_Upwelling_Temp","Max_Monthly_Temp","Mean_Monthly_Summer_Temp" ,            
               "Mean_Monthly_Temp" , "Mean_Monthly_Upwelling_Temp", "MHW_Days"  ,                           
               "MHW_Intensity", "MHW_Summer_Days", "MHW_Summer_Intensity",                 
               "MHW_Upwelling_Days" ,"MHW_Upwelling_Intensity", "Min_Monthly_Anomaly_Summer_Temp"  ,    
               "Min_Monthly_Anomaly_Temp" ,"Min_Monthly_Anomaly_Upwelling_Temp", "Min_Monthly_Temp" ) # 19 preds

predictors <- dat2 %>%
  dplyr::select(var.preds) %>% glimpse()

pred.var <- dat2 %>%
  dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = Days_22C:Min_Monthly_Temp, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()


pvar <- pred.var %>%
  #mutate(text = fct_reorder(Variable, Value)) %>%
  ggplot( aes(x=Value, color=Variable, fill=Variable)) +
  #geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
  geom_histogram(position="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~Variable, scales = 'free') +
  #theme_ipsum() +
  # theme(
  #   legend.position="none",
  #   panel.spacing = unit(0.1, "lines"),
  #   strip.text.x = element_text(size = 8)
  # ) +
  xlab("") +
  ylab("Frequency") + 
  theme(legend.position="bottom") 

pvar

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in var.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}


##

## Temperature distributions with transformations ----

glimpse(predictors)

pred.trans <- predictors %>%
  mutate(log_Days_22C = log(Days_22C + 1),
         log_Max_Monthly_Anomaly_Summer_Temp = log(Max_Monthly_Anomaly_Summer_Temp + 1),
         log_Max_Monthly_Anomaly_Temp = log(Max_Monthly_Anomaly_Temp + 1),
         log_MHW_Days = log(MHW_Days + 1),
         log_MHW_Intensity = log(MHW_Intensity + 1),
         log_MHW_Summer_Days = log(MHW_Summer_Days + 1),
         log_MHW_Summer_Intensity = log(MHW_Summer_Intensity),
         log_MHW_Upwelling_Days = log(MHW_Upwelling_Days + 1),
         log_MHW_Upwelling_Intensity = log(MHW_Upwelling_Intensity + 1)) %>%
  dplyr::select(-c(Days_22C, Max_Monthly_Anomaly_Summer_Temp,
                   Max_Monthly_Anomaly_Temp, MHW_Days,
                   MHW_Intensity, MHW_Summer_Days,
                   MHW_Summer_Intensity, MHW_Upwelling_Days,
                   MHW_Upwelling_Intensity)) %>%
  glimpse()
         

pred.var <- pred.trans %>%
  #dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = Max_Monthly_Anomaly_Upwelling_Temp:log_MHW_Upwelling_Intensity, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()


pvar <- pred.var %>%
  #mutate(text = fct_reorder(Variable, Value)) %>%
  ggplot( aes(x=Value, color=Variable, fill=Variable)) +
  #geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
  geom_histogram(position="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~Variable, scales = 'free') +
  #theme_ipsum() +
  # theme(
  #   legend.position="none",
  #   panel.spacing = unit(0.1, "lines"),
  #   strip.text.x = element_text(size = 8)
  # ) +
  xlab("") +
  ylab("Frequency") + 
  theme(legend.position="bottom") 

pvar



## Temperature Only 1 ----



names(dat2)

# temperature --
var.preds <- c("Days_22C", "Max_Monthly_Anomaly_Summer_Temp","Max_Monthly_Anomaly_Temp" ,            
               "Max_Monthly_Anomaly_Upwelling_Temp","Max_Monthly_Temp","Mean_Monthly_Summer_Temp" ,            
               "Mean_Monthly_Temp" , "Mean_Monthly_Upwelling_Temp", "MHW_Days"  ,                           
               "MHW_Intensity", "MHW_Summer_Days", "MHW_Summer_Intensity",                 
               "MHW_Upwelling_Days" ,"MHW_Upwelling_Intensity", "Min_Monthly_Anomaly_Summer_Temp"  ,    
               "Min_Monthly_Anomaly_Temp" ,"Min_Monthly_Anomaly_Upwelling_Temp", "Min_Monthly_Temp" ) # 19 preds

predictors <- dat2 %>%
  dplyr::select(var.preds) %>% glimpse()

names(predictors)

length(names(predictors)) # 18

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #   4088
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C

# compute the p-value of correlations --
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


## Temperature Only 2 ----

# remove: Days and Min Anomaly Temp

# temperature --
var.preds <- c("Days_22C", 
               "Max_Monthly_Anomaly_Summer_Temp",
               #"Max_Monthly_Anomaly_Temp" ,            
               "Max_Monthly_Anomaly_Upwelling_Temp",
               #"Max_Monthly_Temp",
               #"Mean_Monthly_Summer_Temp" ,            
               "Mean_Monthly_Temp" , 
               #"Mean_Monthly_Upwelling_Temp", 
               "MHW_Days"  ,                           
               #"MHW_Intensity", "MHW_Summer_Days", 
               #"MHW_Summer_Intensity",                 
               "MHW_Upwelling_Days" ,
               #"MHW_Upwelling_Intensity", 
               #"Min_Monthly_Anomaly_Summer_Temp"  ,    
               "Min_Monthly_Anomaly_Temp" ,
               "Min_Monthly_Anomaly_Upwelling_Temp" 
               #"Min_Monthly_Temp" 
               ) # 19 preds

predictors <- dat2 %>%
  dplyr::select(var.preds) %>% glimpse()

length(names(predictors)) # 8

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #  4088
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


## Temperature Only 3 ----

# remove: Days and Min Anomaly Temp

var.preds <- c("Days_22C", 
               "Max_Monthly_Anomaly_Summer_Temp",
               #"Max_Monthly_Anomaly_Temp" ,            
               #"Max_Monthly_Anomaly_Upwelling_Temp",
               #"Max_Monthly_Temp",
               #"Mean_Monthly_Summer_Temp" ,            
               "Mean_Monthly_Temp" , 
               #"Mean_Monthly_Upwelling_Temp", 
               #"MHW_Days"  ,                           
               #"MHW_Intensity", "MHW_Summer_Days", 
               #"MHW_Summer_Intensity",                 
               "MHW_Upwelling_Days" 
               #"MHW_Upwelling_Intensity", 
               #"Min_Monthly_Anomaly_Summer_Temp"  ,    
               #"Min_Monthly_Anomaly_Temp" ,
               #"Min_Monthly_Anomaly_Upwelling_Temp" 
               #"Min_Monthly_Temp" 
) # 19 preds

predictors <- dat2 %>%
  dplyr::select(var.preds) %>% glimpse()

length(names(predictors)) # 4

length(which(is.na(predictors))) # 18
predictors <- na.omit(predictors)
nrow(predictors) #  460
names(predictors2)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


## Bio - distribution of variables  ----

names(df.sc)
names(dat2)

# bio --
var.preds <- c("den_MACPYRAD" , "den_MESFRAAD" ,                        
               "den_PANINT", "den_STRPURAD" , "den_PYCHEL",                           
               "den_HALRUF", "seaotter_dens_sm") # 19 preds

predictors <- dat2 %>%
  dplyr::select(var.preds) %>% glimpse()

pred.var <- dat2 %>%
  dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = den_MACPYRAD:seaotter_dens_sm, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()


pvar <- pred.var %>%
  #mutate(text = fct_reorder(Variable, Value)) %>%
  ggplot( aes(x=Value, color=Variable, fill=Variable)) +
  #geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
  geom_histogram(position="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~Variable, scales = 'free') +
  #theme_ipsum() +
  # theme(
  #   legend.position="none",
  #   panel.spacing = unit(0.1, "lines"),
  #   strip.text.x = element_text(size = 8)
  # ) +
  xlab("") +
  ylab("Frequency") + 
  theme(legend.position="bottom") 

pvar

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in var.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}


##

## Bio distributions with transformations ----

glimpse(predictors)

pred.trans <- predictors %>%
  mutate(log_den_MACPYRAD = log(den_MACPYRAD + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_PANINT = log(den_PANINT + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_den_PYCHEL = log(den_PYCHEL + 1),
         log_den_HALRUF = log(den_HALRUF + 1)) %>%
  dplyr::select(-c(den_MACPYRAD,
                   den_MESFRAAD,
                   den_PANINT,
                   den_STRPURAD,
                   den_PYCHEL,
                   den_HALRUF)) %>%
  glimpse()


pred.var <- pred.trans %>%
  #dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = seaotter_dens_sm:log_den_HALRUF, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()


pvar <- pred.var %>%
  #mutate(text = fct_reorder(Variable, Value)) %>%
  ggplot( aes(x=Value, color=Variable, fill=Variable)) +
  #geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
  geom_histogram(position="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~Variable, scales = 'free') +
  #theme_ipsum() +
  # theme(
  #   legend.position="none",
  #   panel.spacing = unit(0.1, "lines"),
  #   strip.text.x = element_text(size = 8)
  # ) +
  xlab("") +
  ylab("Frequency") + 
  theme(legend.position="bottom") 

pvar



## Bio Only 1 ----


# bio --
var.preds <- c("den_MACPYRAD" , "den_MESFRAAD" ,                        
               #"den_PANINT", 
               "den_STRPURAD"  
               #"den_PYCHEL",                           
               #"den_HALRUF", 
               #"seaotter_dens_sm"
               ) # 19 preds

predictors <- dat2 %>%
  dplyr::select(var.preds) %>% glimpse()

pred.trans <- predictors %>%
  mutate(log_den_MACPYRAD = log(den_MACPYRAD + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         #log_den_PANINT = log(den_PANINT + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1)
         #log_den_PYCHEL = log(den_PYCHEL + 1),
         #log_den_HALRUF = log(den_HALRUF + 1)
         ) %>%
  dplyr::select(-c(den_MACPYRAD,
                   den_MESFRAAD,
                   #den_PANINT,
                   den_STRPURAD
                   #den_PYCHEL,
                   #den_HALRUF
                   )) %>%
  glimpse()

predictors <- pred.trans %>%
  #dplyr::select(var.preds) %>% #glimpse()
  #pivot_longer(cols = seaotter_dens_sm:log_den_STRPURAD, names_to = 'Variable', values_to = 'Value') %>%
  glimpse()

names(predictors)

length(names(predictors)) # 18

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #   4088
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C

# compute the p-value of correlations --
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:3])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##

##


## All Env Vars 1 ----

names(dat2)

predictors <- dat2 %>%
  dplyr::select(# Substrate vars
    sd_prob_of_rock, min_prob_of_rock, mean_prob_of_rock, 
    max_depth, min_depth, depth_mean, 
    mean_slope, min_vrm,
    # Nitrate vars 
    Max_Monthly_Nitrate, Mean_Monthly_Summer_Nitrate, Days_3N,
    Min_Monthly_Anomaly_Upwelling_Nitrate,
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
  glimpse() # Rows: 4,088
    
names(predictors)

length(names(predictors)) # 24

length(which(is.na(predictors))) # 
predictors <- na.omit(predictors)
nrow(predictors) #   3971
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C

# compute the p-value of correlations --
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.7, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.7, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


###


## All Env Vars 2 - tranformed ----

names(dat2)

predictors <- dat2 %>%
  dplyr::select(# Substrate vars
    sd_prob_of_rock, min_prob_of_rock, mean_prob_of_rock, 
    max_depth, min_depth, depth_mean, 
    mean_slope, min_vrm,
    # Nitrate vars 
    Max_Monthly_Nitrate, Mean_Monthly_Summer_Nitrate, Days_3N,
    Min_Monthly_Anomaly_Upwelling_Nitrate,
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
         log_Days_21C = log(Days_21C + 1),
         log_MHW_Upwelling_Days = log(MHW_Upwelling_Days + 1),
         log_den_MACPYRAD = log(den_MACPYRAD + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(min_vrm, Days_21C, MHW_Upwelling_Days,
                   den_MACPYRAD, den_MESFRAAD, den_STRPURAD)) %>%
  glimpse() # Rows: 4,088

names(predictors)

length(names(predictors)) # 24

length(which(is.na(predictors))) # 
predictors <- na.omit(predictors)
nrow(predictors) #   3971
names(predictors)

C <- cor(predictors, method = "pearson")
head(round(C,1))
C

# compute the p-value of correlations --
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.7, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.7, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


###



