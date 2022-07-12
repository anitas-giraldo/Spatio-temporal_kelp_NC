###
### Script modifed from one created by : Anita Giraldo on 21 March 2022
### Script last updated by : Anita Giraldo on 11 July 2022  when received orb vel and NPP data from Tom

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

## Remove variables that are not necessary ----
names(df.nc)


dat <- df.nc %>%
  dplyr::select(-c(Days_18C:Days_23C, 
                   Degree_Days_23C:Degree_Days_23C, 
                   Days_11N:Days_7N,
                   Degree_Days_11N:Degree_Days_7N,
    Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    Max_Monthly_Nitrate_Index, Min_Monthly_Temp_Nitrate
  )) %>%
  glimpse() # 708



###


## Nitrate - distribution of variables  ----


names(dat)

#predictors <- dat2[,c(35:77)]

# nitrate --
nit.preds <- c("Days_8N", "Days_9N", "Days_10N",
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
               "Min_Monthly_Nitrate" )

predictors <- dat %>%
  dplyr::select(nit.preds) %>% 
  glimpse() # 14

pred.nit <- dat %>%
  dplyr::select(nit.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(nit.preds), names_to = 'Variable', values_to = 'Value') %>%
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

dev.off()




###




## Check correlations ----
 
## Nitrate Only 1 ----

names(predictors)

length(names(predictors)) # 14

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

predictors <- dat %>%
  dplyr::select(Days_10N, 
                Min_Monthly_Nitrate, 
                Max_Monthly_Nitrate,
                Mean_Monthly_Nitrate,
                Mean_Monthly_Upwelling_Nitrate,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate) %>%
  glimpse() # 7

length(names(predictors)) # 7

length(which(is.na(predictors))) # 
predictors <- na.omit(predictors)
nrow(predictors) #  708
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

names(df.nc)
names(dat)

# nitrate --
var.preds <- c("Days_15C", "Days_16C" , "Days_17C",
               "Max_Monthly_Anomaly_Summer_Temp" ,      "Max_Monthly_Anomaly_Temp" ,            
               "Max_Monthly_Anomaly_Upwelling_Temp",    "Max_Monthly_Temp" ,                     "Mean_Monthly_Summer_Temp",             
               "Mean_Monthly_Temp" ,                    "Mean_Monthly_Upwelling_Temp" ,          "MHW_Days",                             
               "MHW_Intensity",                         "MHW_Summer_Days"  ,                     "MHW_Summer_Intensity" ,                
               "MHW_Upwelling_Days"  ,                  "MHW_Upwelling_Intensity"  ,             "Min_Monthly_Anomaly_Summer_Temp" ,     
               "Min_Monthly_Anomaly_Temp",              "Min_Monthly_Anomaly_Upwelling_Temp" ,   "Min_Monthly_Temp" )

predictors <- dat %>%
  dplyr::select(var.preds) %>% glimpse() # 20

pred.var <- dat %>%
  dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(var.preds), names_to = 'Variable', values_to = 'Value') %>%
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
  mutate(log_Days_15C = log(Days_15C + 1),
         log_MHW_Intensity = log(MHW_Intensity + 1),
         sqrt_MHW_Summer_Days = sqrt(MHW_Summer_Days),
         sqrt_MHW_Summer_Intensity = sqrt(MHW_Summer_Intensity),
         sqrt_MHW_Upwelling_Days = sqrt(MHW_Upwelling_Days),
         sqrt_MHW_Upwelling_Intensity = sqrt(MHW_Upwelling_Intensity),
         sqrt_MHW_Days = sqrt(MHW_Days),
         log_Max_Monthly_Anomaly_Summer_Temp = log(Max_Monthly_Anomaly_Summer_Temp + 1),
         log_Max_Monthly_Anomaly_Temp = log(Max_Monthly_Anomaly_Temp + 1)) %>%
  dplyr::select(-c(Days_15C,
                   MHW_Intensity, 
                   MHW_Summer_Days,
                   MHW_Summer_Intensity, 
                   MHW_Days,
                   Max_Monthly_Anomaly_Summer_Temp,
                   Max_Monthly_Anomaly_Temp,
                   MHW_Upwelling_Intensity,
                   MHW_Upwelling_Days
                   )) %>%
  glimpse()
         

pred.var <- pred.trans %>%
  #dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(names(predictors)), names_to = 'Variable', values_to = 'Value') %>%
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



names(predictors)

# temperature --


names(predictors)

length(names(predictors)) # 20

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #   708
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
var.preds <- c("Days_16C" ,
               "Mean_Monthly_Temp" ,
               "Mean_Monthly_Summer_Temp",
               "MHW_Upwelling_Days"  , 
               "Min_Monthly_Anomaly_Temp",
               "Max_Monthly_Anomaly_Upwelling_Temp",
               "Min_Monthly_Temp", 
               "Mean_Monthly_Upwelling_Temp" )

predictors <- dat %>%
  dplyr::select(var.preds) %>% glimpse()

length(names(predictors)) # 8

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #  708
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


## Bio - distribution of variables  ----

names(df.nc)
names(dat)

# bio --
var.preds <- c("den_NERLUE" , "den_MESFRAAD" , "den_STRPURAD" , "den_PYCHEL",                           
               "den_HALRUF") # 19 preds

predictors <- dat %>%
  dplyr::select(var.preds) %>% glimpse()

pred.var <- dat %>%
  dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(var.preds), names_to = 'Variable', values_to = 'Value') %>%
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
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_den_PYCHEL = log(den_PYCHEL + 1),
         log_den_HALRUF = log(den_HALRUF + 1)) %>%
  dplyr::select(-c(den_NERLUE,
                   den_MESFRAAD,
                   den_STRPURAD,
                   den_PYCHEL,
                   den_HALRUF)) %>%
  glimpse()


pred.var <- pred.trans %>%
  #dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(var.preds), names_to = 'Variable', values_to = 'Value') %>%
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
names(predictors)


length(names(predictors)) # 5

length(which(is.na(predictors))) # 0
predictors <- na.omit(predictors)
nrow(predictors) #   686
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
         number.cex = 0.7, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##

##


## Other Env Vars 1 ----

names(dat)

var.preds <- c("npp.mean",
               "orb.95"  , "orb.max" ,                             
               "wh.95" ,   "wh.max",
               "pdo_mean" , "npgo_mean" ,"mei_mean" )

predictors <- dat %>%
  dplyr::select(var.preds) %>%
  glimpse() # Rows: 8
    
names(predictors)

length(names(predictors)) # 8

## Env - distribution of variables  ----

names(df.nc)
names(dat)


pred.var <- dat %>%
  dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(var.preds), names_to = 'Variable', values_to = 'Value') %>%
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

dev.off()


# transformations --
names(predictors)

predictors <- predictors %>%
  mutate(Log_npp.mean = log(npp.mean +1)) %>%
  dplyr::select(-npp.mean)%>%
  glimpse()


# plot 
pred.var <- predictors %>%
  #dplyr::select(var.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(var.preds), names_to = 'Variable', values_to = 'Value') %>%
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


## correlations ---

length(which(is.na(predictors))) # 
predictors <- na.omit(predictors)
nrow(predictors) #   528
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
         cl.cex =  0.7, # legend text size
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


###

###


## All Vars 1 - tranformed ----

names(dat)

predictors <- dat %>%
  dplyr::select(# Bio vars
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
    npp.mean,                             
    wh.95 ,   wh.max,
    npgo_mean , mei_mean) %>%
  # Bio transformations
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_den_PYCHEL = log(den_PYCHEL + 1),
         log_den_HALRUF = log(den_HALRUF + 1)) %>%
  dplyr::select(-c(den_NERLUE,
                   den_MESFRAAD,
                   den_STRPURAD,
                   den_PYCHEL,
                   den_HALRUF)) %>%
  # Temperature transformations
  mutate(log_Days_16C = log(Days_16C + 1)) %>%
  dplyr::select(-c(Days_16C)) %>%
  # Other Env transformations
  mutate(log_npp.mean = log(npp.mean + 1)) %>%
  dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708

names(predictors)

length(names(predictors)) # 25

length(which(is.na(predictors))) # 325
predictors <- na.omit(predictors)
nrow(predictors) #   507
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


## WAVES - distribution of variables  ----


names(dat)

#predictors <- dat2[,c(35:77)]

# nitrate --
w.preds <- c("wh_mean" , "wh_max" ,  "wh_95prc" ,  "wh_99prc" ,                            
             "mean_waveyear" ,  "max_waveyear" ,   "wh_95prc_wy",                          
             "wh_99prc_wy")

predictors <- dat %>%
  dplyr::select(w.preds) %>% 
  glimpse() # 8

pred.w <- dat %>%
  dplyr::select(w.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(w.preds), names_to = 'Variable', values_to = 'Value') %>%
  glimpse() # Rows: 5,664


pw <- pred.w %>%
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

pw

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in w.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

dev.off()




###




## Check correlations ----

## WAVES Only 1 ----

names(predictors)

length(names(predictors)) # 8

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
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##


## WAVES Only 2 ----

# remove: some vars 

names(predictors)

predictors <- dat %>%
  dplyr::select(Days_10N, 
                Min_Monthly_Nitrate, 
                Max_Monthly_Nitrate,
                Mean_Monthly_Nitrate,
                Mean_Monthly_Upwelling_Nitrate,
                Max_Monthly_Anomaly_Nitrate,
                Mean_Monthly_Summer_Nitrate) %>%
  glimpse() # 7

length(names(predictors)) # 7

length(which(is.na(predictors))) # 
predictors <- na.omit(predictors)
nrow(predictors) #  708
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





###

###

## ORB VELOCITY - distribution of variables  ----


names(dat)

#predictors <- dat2[,c(35:77)]

# nitrate --
w.preds <- c("UBR_Max" , "UBR_Mean",                             
             "UBRYear_Max" , "UBRYear_Mean" , "Wave_Max" ,                            
             "Wave_Mean" , "WaveYear_Max" , "WaveYear_Mean")

predictors <- dat %>%
  dplyr::select(w.preds) %>% 
  glimpse() # 8

pred.w <- dat %>%
  dplyr::select(w.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(w.preds), names_to = 'Variable', values_to = 'Value') %>%
  glimpse() # Rows: 5,664


pw <- pred.w %>%
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

pw

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in w.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

dev.off()

## 

## Transform variables  ----

pred.w2 <- predictors %>%
  dplyr::select(UBR_Max, UBR_Mean, UBRYear_Max, UBRYear_Mean) %>%
  mutate(sqrt.UBR_Max = sqrt(UBR_Max),
         sqrt.UBR_Mean = sqrt(UBR_Mean),
         sqrt.UBRYear_Mean = sqrt(UBRYear_Mean),
         log.UBRYear_Max = sqrt(UBRYear_Max)) %>%
  dplyr::select(-c(UBR_Max, UBR_Mean, UBRYear_Max, UBRYear_Mean)) %>%
  glimpse()

pred.w3 <- pred.w2 %>%
  #dplyr::select(w.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(pred.w2), names_to = 'Variable', values_to = 'Value') %>%
  glimpse() # Rows: 5,664


pw <- pred.w3 %>%
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

pw


###


## Check correlations ----

## ORB VELOCITY Only 1 ----

names(predictors)

length(names(predictors)) # 8

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
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##



###

###

## NPP - distribution of variables  ----


names(dat)

#predictors <- dat2[,c(35:77)]

# nitrate --
w.preds <- c("Max_Monthly_NPP" ,  "Max_Monthly_NPP_Summer" ,"Max_Monthly_NPP_Upwelling",            
             "Mean_Monthly_NPP", "Mean_Monthly_NPP_Summer", "Mean_Monthly_NPP_Upwelling",           
             "Min_Monthly_NPP", "Min_Monthly_NPP_Summer", "Min_Monthly_NPP_Upwelling")

predictors <- dat %>%
  dplyr::select(w.preds) %>% 
  glimpse() # 8

pred.w <- dat %>%
  dplyr::select(w.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(w.preds), names_to = 'Variable', values_to = 'Value') %>%
  glimpse() # Rows: 5,664


pw <- pred.w %>%
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

pw

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in w.preds) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

dev.off()

## 

## Transform variables  ----

pred.w2 <- predictors %>%
  mutate(log.Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling),
         log.Min_Monthly_NPP = log(Min_Monthly_NPP),
         log.Min_Monthly_NPP_Summer = log(Min_Monthly_NPP_Summer),
         log.Min_Monthly_NPP_Upwelling = log(Min_Monthly_NPP_Upwelling)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling, Min_Monthly_NPP, Min_Monthly_NPP_Summer, Min_Monthly_NPP_Upwelling)) %>%
  glimpse()

pred.w3 <- pred.w2 %>%
  #dplyr::select(w.preds) %>% #glimpse()
  pivot_longer(cols = 1:length(pred.w2), names_to = 'Variable', values_to = 'Value') %>%
  glimpse() # Rows: 5,664


pw <- pred.w3 %>%
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

pw


###


## Check correlations ----

## NPP Only 1 ----

names(predictors)

length(names(predictors)) # 8

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
         number.cex = 0.6, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


##

### 

### Correlations of all selected variables ----

dat1 <- dat %>%
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
    UBR_Mean,
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
  mutate(log_UBR_Mean = log(UBR_Mean + 1)) %>%
  dplyr::select(-c(UBR_Mean)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
  glimpse() # Rows: 708


## correlations --

names(dat1)

predictors <- dat1[,7:ncol(dat1)]
names(predictors)

length(names(predictors)) # 35

length(which(is.na(predictors))) 
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
         tl.cex = 0.5, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.5, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
