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
d.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.ls.outputs"
o.dir <- here("outputs_nc_ls")
p.dir <- paste(o.dir, 'plots', sep ='/')



## Load data ----

df <- read.csv(paste(d.dir, "NC_subsamp_Landsat_99-21_merged_vars_noNAs.csv", sep = '/')) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  mutate(depth = depth*(-1)) %>% # make depth positive
  glimpse()



## Remove variables that are not necessary ----
names(df)

df.nc <- df %>%
  dplyr::select(-c(# temp variables
    Days_15C, Days_19C, Days_20C, Days_21C, Days_22C, Days_23C,
                   Degree_Days_15C, Degree_Days_19C, Degree_Days_20C, Degree_Days_21C, Degree_Days_22C, Degree_Days_23C,
                   Max_Monthly_Temp_Index, Min_Monthly_Temp_Index,
    # nitrate variables
                   Days_4N, Days_5N, Days_6N, Days_7N, Days_8N, Days_9N,
                   Days_10N, Days_11N, Days_12N, Days_13N, Days_14N, Days_15N,
                   Degree_Days_4N, Degree_Days_5N, Degree_Days_6N,
                   Degree_Days_7N, Degree_Days_8N, Degree_Days_9N,
                   Degree_Days_10N, Degree_Days_11N, Degree_Days_12N,
                   Degree_Days_13N, Degree_Days_14N, Degree_Days_15N,
                   Max_Monthly_Nitrate_Index, Min_Monthly_Temp_Nitrate,
    # others
    ID, latlon.ci, latitude, longitude)) %>%
  glimpse() # Rows: 3,151 - Columns: 52


###

###


# Substrate - distribution of variables  ----

names(df.nc)

predictors <- df.nc[,c(6:9)]

# substrate --
sub.preds <- c("depth", "slope", "vrm", "prob_of_rock")

pred.sub <- predictors %>%
  dplyr::select(sub.preds) %>% #glimpse()
  #pivot_longer(cols = depth_mean:prop_map_vrm, names_to = 'Variable', values_to = 'Value') %>%
  pivot_longer(cols = 1:4, names_to = 'Variable', values_to = 'Value') %>%
  mutate_at(vars(Variable), list(as.factor)) %>%
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
for (i in 1:length(predictors)) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(sub.preds[i]))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

dev.off()


## Correlations Substrate ----

names(predictors)

length(names(predictors)) # 4


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

#corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)
#corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)


# transformed ----

head(predictors)

predictors.t <- predictors %>%
  mutate(depth_sqrt = sqrt(depth),
         slope_log = log(slope + 1),
         vrm_sqrt = sqrt(vrm),
         prob_of_rock_log = log(prob_of_rock + 1)) %>%
  dplyr::select(-c(depth, slope, vrm, prob_of_rock)) %>%
  glimpse()

## correlations --
C <- cor(predictors.t, method = "pearson")
head(round(C,1))
C

p.mat <- cor.mtest(predictors.t)

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

##




###


## Nitrate - distribution of variables  ----


names(df.nc)

predictors <- df.nc[,c(36:52)]

# nitrate --
nit.preds <- c("Days_1N" ,"Days_2N", "Days_3N" , "Degree_Days_1N" ,"Degree_Days_2N", "Degree_Days_3N", 
               "Max_Monthly_Anomaly_Nitrate" ,"Max_Monthly_Anomaly_Summer_Nitrate", 
               "Max_Monthly_Anomaly_Upwelling_Nitrate", "Max_Monthly_Nitrate", "Mean_Monthly_Nitrate", 
               "Mean_Monthly_Summer_Nitrate" ,"Mean_Monthly_Upwelling_Nitrate",       
               "Min_Monthly_Anomaly_Nitrate" , "Min_Monthly_Anomaly_Summer_Nitrate" ,
               "Min_Monthly_Anomaly_Upwelling_Nitrate", "Min_Monthly_Nitrate" )


preds <- df.nc %>%
  dplyr::select(nit.preds) %>% 
  glimpse()

pred.nit <- preds %>%
  dplyr::select(nit.preds) %>% #glimpse()
  pivot_longer(cols = 1:17, names_to = 'Variable', values_to = 'Value') %>%
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


## Correlations Nitrate Only 1 ----

names(predictors)

length(names(predictors)) # 17

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


##


## Nitrate Only 2 ----

# remove: some vars 

names(predictors)

to.remove <- c('Days_1N', 'Days_2N', #'Days_3N',
               'Degree_Days_1N',
               'Min_Monthly_Anomaly_Summer_Nitrate',
               'Min_Monthly_Nitrate')

predictors2 <- predictors %>%
  dplyr::select(-c(to.remove)) %>%
  glimpse()

length(names(predictors2)) # 11

C <- cor(predictors2, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors2)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.8, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.8, # legend text size
         number.cex = 0.55, # text labels size
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


# temp --
preds <- c("Days_16C" , "Days_17C" ,  "Days_18C",                             
           "Degree_Days_16C", "Degree_Days_17C", "Degree_Days_18C",                      
           "Max_Monthly_Anomaly_Summer_Temp", "Max_Monthly_Anomaly_Temp", "Max_Monthly_Anomaly_Upwelling_Temp",   
           "Max_Monthly_Temp","Mean_Monthly_Summer_Temp", "Mean_Monthly_Temp",                    
           "Mean_Monthly_Upwelling_Temp", "MHW_Days", "MHW_Intensity",                        
           "MHW_Summer_Days", "MHW_Summer_Intensity" , "MHW_Upwelling_Days",                   
           "MHW_Upwelling_Intensity", "Min_Monthly_Anomaly_Summer_Temp", "Min_Monthly_Anomaly_Temp",             
           "Min_Monthly_Anomaly_Upwelling_Temp", "Min_Monthly_Temp") # 19 preds

predictors <- df.nc %>%
  dplyr::select(preds) %>% 
  glimpse()

pred.var <- predictors %>%
  #dplyr::select(preds) %>% #glimpse()
  pivot_longer(cols = 1:23, names_to = 'Variable', values_to = 'Value') %>%
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
for (i in 1:length(predictors)) {
  x <- predictors[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(preds[i]))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

dev.off()


##

## Temperature distributions with transformations ----

glimpse(predictors)

length(names(predictors)) # 23

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


##


## Temperature Only 2 ----

# remove: Days and Min Anomaly Temp
names(predictors)

# temperature --
to.remove <- c('Days_18C', 'Degree_Days_18C', 'Max_Monthly_Anomaly_Summer_Temp',
               'Max_Monthly_Temp', 'Max_Monthly_Anomaly_Temp', 'MHW_Upwelling_Days',
               'MHW_Upwelling_Intensity',
               'MHW_Summer_Days', 'MHW_Summer_Intensity',
               'Min_Monthly_Anomaly_Upwelling_Temp') # 

predictors2 <- predictors %>%
  dplyr::select(-c(to.remove)) %>% 
  glimpse()

length(names(predictors2)) # 12

C <- cor(predictors2, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors2)
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

to.remove <- c('Degree_Days_17C',
               'Days_17C',
               'Degree_Days_16C', 'Max_Monthly_Anomaly_Upwelling_Temp',
               'Min_Monthly_Temp') # 

predictors3 <- predictors2 %>%
  dplyr::select(-c(to.remove)) %>% glimpse()

length(names(predictors3)) # 4

C <- cor(predictors3, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors3)
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


## Correlations All variables ----

names(df.nc)

pred.names <- c('depth', 'slope', 'vrm', 'prob_of_rock',
                'Max_Monthly_Anomaly_Summer_Nitrate',
                'Mean_Monthly_Summer_Nitrate',
                'Min_Monthly_Anomaly_Nitrate',
                'Min_Monthly_Anomaly_Upwelling_Nitrate',
                'Max_Monthly_Anomaly_Nitrate',
                'Days_3N',
                'Max_Monthly_Anomaly_Upwelling_Nitrate', 
                'Max_Monthly_Nitrate',
                'Mean_Monthly_Upwelling_Nitrate',
                'Days_16C', 'Mean_Monthly_Upwelling_Temp',
                'MHW_Intensity',
                'Mean_Monthly_Summer_Temp',
                'Mean_Monthly_Temp',
                'Min_Monthly_Anomaly_Summer_Temp',
                'Min_Monthly_Anomaly_Temp',
                'pdo_mean', 'npgo_mean', 'mei_mean')

predictors <- df.nc %>%
  dplyr::select(c(pred.names)) %>% glimpse()

length(names(predictors)) # 23

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
         tl.cex = 0.6, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.6, # legend text size
         number.cex = 0.45, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

###

## Al variables 2 ----

pred.names <- c('Max_Monthly_Anomaly_Summer_Nitrate',
                'Mean_Monthly_Summer_Nitrate',
                'Min_Monthly_Anomaly_Nitrate',
                'Min_Monthly_Anomaly_Upwelling_Nitrate',
                'Max_Monthly_Anomaly_Nitrate',
                'Days_3N',
                'Max_Monthly_Anomaly_Upwelling_Nitrate', 
                'Max_Monthly_Nitrate',
                'Mean_Monthly_Upwelling_Nitrate',
                'Days_16C', 'Mean_Monthly_Upwelling_Temp',
                'MHW_Intensity',
                'Mean_Monthly_Summer_Temp',
                'Mean_Monthly_Temp',
                'Min_Monthly_Anomaly_Summer_Temp',
                'Min_Monthly_Anomaly_Temp',
                'pdo_mean', 'npgo_mean', 'mei_mean')

predictors2 <- df.nc %>%
  dplyr::select(c(pred.names)) %>% glimpse()

length(names(predictors2)) # 23

C <- cor(predictors2, method = "pearson")
head(round(C,1))
C


# matrix of the p-value of the correlation
p.mat <- cor.mtest(predictors2)
head(p.mat[, 1:5])

# customize correlogram --
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(100),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.6, tl.col="black", tl.srt= 45, # Text label size, color and rotation 
         cl.cex =  0.6, # legend text size
         number.cex = 0.5, # text labels size
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

corrplot(C, method="color", type = "lower", order="hclust", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)

corrplot(C, type = "lower", tl.cex = 0.7, tl.col = 'black', tl.srt= 45)