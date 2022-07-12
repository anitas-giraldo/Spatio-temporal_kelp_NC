
# by Anita Giraldo, 12 April 2022
# last modified, Anita Giraldo, 11 July, 2022, added orb vel and NPP from TOM


## Check correlations between predictors and dependent ----


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
library(ggpmisc)



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
o.dir <- here("outputs_nc_rcca")
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
p.dir <- paste(o.dir, 'plots', sep ='/')

# load data ----
df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep = '/')) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  rename(wh_mean = mean,
         wh_max = max) %>%
  #dplyr::select(-c(group, site_to_match_to, min_distance, unique_id.x,  unique_id.y, x, y)) %>%
  #mutate(depth = depth*(-1)) %>% # make depth positive
  glimpse()


###

## ENV PREDICTORS ----

names(df)

# define the predictor names ---
ynames <- names(df[16:23])
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

names(df)

for(i in 1:length(ynames)){
  
  df.plot <- df[c(12, i+15)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density of N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "env_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=2, ncol=2)
dev.off()

###

###

## TEMPERATURE PREDICTORS ----

names(df)

names.temp <- names(df)[c(24:26, 43,44,46:58,60)]
length(names.temp) # 19


predictors <- df %>%
  dplyr::select(den_NERLUE, names.temp) %>%
  glimpse()



# define the predictor names ---
ynames <- names(df)[c(24:26, 43,44,46:58,60)]
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###
names(df)

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- predictors[c(1, i+1)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density of N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "temp_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=3, ncol=2)
dev.off()


###

###

## NITRATE PREDICTORS ----

names(df)


names.nit <- names(df)[c(61,74,75,91:93,95:102)]
length(names.nit) # 14


predictors <- df %>%
  dplyr::select(den_NERLUE, names.nit) %>%
  glimpse()



# define the predictor names ---
ynames <- names(df)[c(61,74,75,91:93,95:102)]
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- predictors[c(1, i+1)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "nit_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=2, ncol=2)
dev.off()



###

###

## BIO PREDICTORS ----

names(df)


names.bio <- names(df)[c(8:11)]
length(names.bio) # 4


predictors <- df %>%
  dplyr::select(den_NERLUE, names.bio) %>%
  glimpse()



# define the predictor names ---
ynames <- names(df)[c(8:11)]
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- predictors[c(1, i+1)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "bio_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=2, ncol=2)
dev.off()


###


## WAVE PREDICTORS ----

names(df)


names.bio <- names(df)[c(128:135)]
length(names.bio) # 8


predictors <- df %>%
  dplyr::select(den_NERLUE, names.bio) %>%
  glimpse()



# define the predictor names ---
ynames <- names(df)[c(128:135)]
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- predictors[c(1, i+1)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "wave_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=2, ncol=2)
dev.off()




###


## ORBITAL VELOCITY PREDICTORS ----

names(df)


names.bio <- names(df)[c(136:143)]
length(names.bio) # 8


predictors <- df %>%
  dplyr::select(den_NERLUE, names.bio) %>%
  glimpse()



# define the predictor names ---
ynames <- names(df)[c(136:143)]
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- predictors[c(1, i+1)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "orbvel_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=2, ncol=2)
dev.off()



###


## NPP PREDICTORS ----

names(df)


names.bio <- names(df)[c(144:152)]
length(names.bio) # 8


predictors <- df %>%
  dplyr::select(den_NERLUE, names.bio) %>%
  glimpse()



# define the predictor names ---
ynames <- names(df)[c(144:152)]
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- predictors[c(1, i+1)] #%>% glimpse()
  name.pred <- names(df.plot[2])
  colnames(df.plot)[1] <- "kelp.area"
  colnames(df.plot)[2] <- "predictor"
  
  # ggplot(df.plot, aes(x = log_den_MACPYRAD, y = predictor)) +        
  #   geom_point() +
  #   labs(x= name.pred, y='Log M. ppyrifera density') +
  #   theme_bw()
  
  p <- ggplot(df.plot, aes(x = predictor, y = kelp.area)) +   
    geom_smooth(method = "lm", se=FALSE, color="blue", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), color = "blue", parse = TRUE) +         
    geom_point(size = 0.5) +
    labs(x= name.pred, y='Density N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "npp_corrs-test.pdf", sep ='/'))
#par(mfrow=c(2,3), mar=c(9,4,3,1))
#do.call(grid.arrange(nrow=4, ncol=2), plot_list)
marrangeGrob(plot_list, nrow=2, ncol=2)
dev.off()

