###
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
ncsites <- years %>%
  dplyr::select(-c(X, latitude, longitude)) %>%
  mutate_at(vars(site_campus_unique_ID), list(as.factor)) %>%
  # get north coast 
  dplyr::filter(region == 'nc') %>%
  # get only sites with PRE MHW data 
  dplyr::filter(preMHW > 0) %>%
  # get sites with 3 or more total survey years
  dplyr::filter(no.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 16

# gets lats and lons ----
nclatlons <- years %>%
  dplyr::select(-c(X)) %>%
  mutate_at(vars(site_campus_unique_ID), list(as.factor)) %>%
  # get north coast 
  dplyr::filter(region == 'nc') %>%
  # get only sites with PRE MHW data 
  dplyr::filter(preMHW > 0) %>%
  # get sites with 3 or more total survey years
  dplyr::filter(no.years > 2) %>%
  dplyr::select(site_campus_unique_ID, latitude, longitude) %>%
  droplevels() %>%
  glimpse() # Rows: 16


## Load transect data ----

df <- read.csv(paste(merged.data, "swath_upc_fish_terr_clim_env_depth_seaotter_temp_nitrate_transect.csv", sep ='/')) %>%
  mutate_at(vars(site_campus_unique_ID), list(as.factor)) %>%
  glimpse()


## get the sites for North Coast model ----
df.nc <- df %>%
  dplyr::select(-c(X.x, Y.x)) %>%
  right_join(ncsites, by = c('site_campus_unique_ID', 'mlpa_region', 'region')) %>%
  rename(X = X.y, Y = Y.y) %>%
  droplevels() %>%
  relocate(c(no.years, preMHW, duringMHW, postMHW), .after = longitude) %>%
  relocate(c(X, Y), .after = longitude) %>%
  glimpse() # Rows: 478

length(levels(df.nc$site_campus_unique_ID)) # 16



## Remove variables that are not necessary ----

dat <- df.nc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_MACPYRAD,
    den_PANINT,
    pct_cov_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    seaotter_dens_sm:seaotter_trend5yr,
    Days_15C, Days_16C, Days_18C:Degree_Days_16C, Degree_Days_18C:Degree_Days_23C,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index
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
    Days_17C, Degree_Days_17C,
    # the same as Max_Monthly_Anomaly_Temp
    Max_Monthly_Anomaly_Summer_Temp,
    # the same as Min_Monthly_Anomaly_Temp
    Min_Monthly_Anomaly_Upwelling_Temp,
    # MHW Days are the same as Intensity
    MHW_Intensity, MHW_Summer_Intensity, MHW_Upwelling_Intensity,
    # remove other correlated vars
    mei_mean, pdo_mean, 
    Max_Monthly_Anomaly_Temp,
    sd_depth, sd_prob_of_rock
  )) %>% # 
  glimpse()



## Nitrate - distribution of variables  ----

names(dat2)

predictors <- dat2[,c(35:77)]

# nitrate --
nit.preds <- c("Days_10N", "Max_Monthly_Anomaly_Nitrate" , "Max_Monthly_Anomaly_Summer_Nitrate",
               "Max_Monthly_Anomaly_Upwelling_Nitrate", 
               "Max_Monthly_Nitrate"  , "Mean_Monthly_Nitrate" ,
               "Mean_Monthly_Summer_Nitrate" ,          "Mean_Monthly_Upwelling_Nitrate" ,      
               "Min_Monthly_Anomaly_Nitrate" ,          "Min_Monthly_Anomaly_Summer_Nitrate" ,   
               "Min_Monthly_Anomaly_Upwelling_Nitrate", "Min_Monthly_Nitrate"  ,
               "Min_Monthly_Temp_Nitrate" )

pred.nit <- predictors %>%
  dplyr::select(nit.preds) %>% #glimpse()
  pivot_longer(cols = Days_10N:Min_Monthly_Temp_Nitrate, names_to = 'Variable', values_to = 'Value') %>%
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
 
## Nitrate Only 1 ----

names(dat2)

predictors <- dat2[,c(65:67,69:77)]
names(predictors)

length(names(predictors)) # 43

length(which(is.na(predictors))) # 18
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
         number.cex = 0.6, # text labels size
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

predictors <- predictors[,c(4,5,2,6,9,1,3)]
names(predictors)

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


## Nitrate Only 3 ----

# remove: some vars 

names(predictors)

predictors <- predictors[,c(2:6)]
names(predictors)

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


## Nitrate Only 4 ----

# remove: some vars 

names(predictors)

predictors <- predictors[,c(1,2,4)]
names(predictors)

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

## Temperature - distribution of variables  ----

names(df.nc)
## Remove variables that are not necessary ----

dat <- df.nc %>%
  dplyr::select(-c(
    campus, transect,
    region:postMHW,
    den_MACPYRAD,
    den_PANINT,
    pct_cov_MACPYRHF:pct_cov_TURF,
    count_SPUL:count_TURF,
    seaotter_dens_sm:seaotter_trend5yr,
    Days_18C:Degree_Days_16C, Degree_Days_18C:Degree_Days_23C,
    mean_depth,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index
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

names(dat2)

predictors <- dat2[,c(23:41)]

# nitrate --
var.preds <- c("Days_17C"  , "Degree_Days_17C" ,                     
               "Max_Monthly_Anomaly_Summer_Temp" , "Max_Monthly_Anomaly_Temp" ,  "Max_Monthly_Anomaly_Upwelling_Temp" ,  
               "Max_Monthly_Temp", "Mean_Monthly_Summer_Temp" , "Mean_Monthly_Temp",                    
               "Mean_Monthly_Upwelling_Temp", "MHW_Days" , "MHW_Intensity"  ,                      
               "MHW_Summer_Days",  "MHW_Summer_Intensity" , "MHW_Upwelling_Days",                   
               "MHW_Upwelling_Intensity", "Min_Monthly_Anomaly_Summer_Temp", "Min_Monthly_Anomaly_Temp",             
               "Min_Monthly_Anomaly_Upwelling_Temp","Min_Monthly_Temp") # 19 preds

pred.var <- predictors %>%
  dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = Days_17C:Min_Monthly_Temp, names_to = 'Variable', values_to = 'Value') %>%
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




## Temperature Only 1 ----

## remove correlated substrate terrain variables ----


names(dat2)


## remove Days 17C and Degree Days 17C because they are just zero for the sites and years we have
## Same for Days 16C
## So I included Days 15 C
predictors <- dat2[,c(23,27:42)]
names(predictors)

length(names(predictors)) # 43

length(which(is.na(predictors))) # 18
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
         number.cex = 0.6, # text labels size
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

names(predictors)

predictors2 <- predictors[,c(1:9,11,13,15,17)]
names(predictors2)

length(names(predictors2)) # 43

length(which(is.na(predictors2))) # 18
predictors <- na.omit(predictors2)
nrow(predictors2) #  460
names(predictors2)

C <- cor(predictors2, method = "pearson")
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

names(predictors)

predictors2 <- predictors[,c(1:5,7:13)]
names(predictors2)

length(names(predictors2)) # 43

length(which(is.na(predictors2))) # 18
predictors <- na.omit(predictors2)
nrow(predictors2) #  460
names(predictors2)

C <- cor(predictors2, method = "pearson")
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




## All Env Vars 1 ----

names(dat2)

predictors <- dat2[,c(15,16,19,20,23,27,30,31,36,40,28,32,34,38,42,29,
                      51,52,53,79,77,80,75,84,76)]
names(predictors)

length(names(predictors)) # 25

length(which(is.na(predictors))) # 18
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
         number.cex = 0.45, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


###

###     ###     ###     ###     ###

###



