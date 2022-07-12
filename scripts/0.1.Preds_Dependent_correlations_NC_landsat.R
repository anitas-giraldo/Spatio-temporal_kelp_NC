
# by Anita Giraldo, 12 April 2022
# last modified, Anita Giraldo, 12 April, 2022


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
d.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.ls.outputs"
#raw.dir <- here("data_raw")
#plots.dir <- here("plots")
o.dir <- here("outputs_nc_ls")

# load data ----
df <- read.csv(paste(d.dir, "NC_subsamp_Landsat_99-21_merged_vars_noNAs.csv", sep = '/')) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  mutate(depth = depth*(-1)) %>% # make depth positive
  glimpse()

## plot correlation ----
p.dir <- paste(o.dir, 'plots', sep ='/')



## SUBSTRATE PREDICTORS ----

names(df)

# define the predictor names ---
ynames <- names(df[10:13])

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- df[c(9, i+9)] #%>% glimpse()
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
    labs(x= name.pred, y='Area of N. luetkeana') +
    theme_bw()
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
# https://stackoverflow.com/questions/50371265/multiple-r-ggplots-on-same-pdf-page
# https://stackoverflow.com/questions/69067773/saving-multiple-ggplots-to-a-pdf-multiple-on-one-page-in-a-for-loop
pdf(paste(p.dir, "terrain_corrs-test.pdf", sep ='/'))
par(mfrow=c(2,2), mar=c(9,4,3,1))
do.call(grid.arrange, plot_list)
dev.off()


###

###

## ENV PREDICTORS ----

names(df)

# define the predictor names ---
ynames <- names(df[14:16])
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- df[c(9, i+13)] #%>% glimpse()
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
    labs(x= name.pred, y='Area of N. luetkeana') +
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

# remove vars
df2 <- df %>%
  dplyr::select(-c(Days_15C, Days_19C, Days_20C, Days_21C, Days_22C, Days_23C,
                   Degree_Days_15C, Degree_Days_19C, Degree_Days_20C, Degree_Days_21C, Degree_Days_22C, Degree_Days_23C,
                   Max_Monthly_Temp_Index, Min_Monthly_Temp_Index)) %>%
  glimpse()

names(df2)

# define the predictor names ---
ynames <- names(df[c(17:39)])
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- df2[c(9, i+16)] #%>% glimpse()
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
    labs(x= name.pred, y='Area of N. luetkeana') +
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

# remove vars
df3 <- df %>%
  dplyr::select(-c(Days_4N, Days_5N, Days_6N, Days_7N, Days_8N, Days_9N,
                   Days_10N, Days_11N, Days_12N, Days_13N, Days_14N, Days_15N,
                   Degree_Days_4N, Degree_Days_5N, Degree_Days_6N,
                   Degree_Days_7N, Degree_Days_8N, Degree_Days_9N,
                   Degree_Days_10N, Degree_Days_11N, Degree_Days_12N,
                   Degree_Days_13N, Degree_Days_14N, Degree_Days_15N,
                   Max_Monthly_Nitrate_Index, Min_Monthly_Temp_Nitrate)) %>%
  glimpse()

names(df3)

# define the predictor names ---
ynames <- names(df3[54:70])
ynames
length(ynames)

# define formula for ploting ---
my.formula <- y ~ x


###

library(gridExtra)
plot_list = list()

for(i in 1:length(ynames)){
  
  df.plot <- df3[c(9, i+53)] #%>% glimpse()
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
    labs(x= name.pred, y='Area N. luetkeana') +
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


