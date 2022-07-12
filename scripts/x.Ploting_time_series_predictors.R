
### By Anita Giraldo - 13 Jan 2022
### Last modified - 26 April 2022

## Libraries ----
library(here)
library(dplyr)
library(stringr)
library(googlesheets4) 
library(reshape2)
library(tidyr)
library(plot.matrix)
library(viridisLite)
library(raster)
library(terra)
library(rgdal)
library(rgeos)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(stringr)


# Clear environment ----
rm(list=ls())


# 1. Get rasters of predictors for every year ----

# set directory of rasters --

re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data"

# check variables ---
nereo.mod$formula


# load susbtrate predictors ----


## Prob of rock ----

# set extent
e.nc <- extent(-124.60, -120.6956, 37.62184, 42.00996)
# bathy
sp.dir <- here('spatial')
nc.mask <- rast(paste(sp.dir, "mask_NC.tif", sep ='/'))

sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

prob_rock <- rast(paste(sub.dir, "prob_rock_nc.all_300m_wInterp.tif", sep ='/'))
prob_rock
plot(prob_rock)

## aggregate to 2000 m --
fact <- 2000/300
prob_rock2000 <- aggregate(prob_rock, fact = fact, fun = mean, na.rm = T)
plot(prob_rock2000)

# make df --
df.var <- as.data.frame(prob_rock2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  group_by(y) %>%
  summarise(mean_value = mean(prob_rock_nc.all_30m_wInterp, na.rm = T)) %>%
  mutate(dummy = 0) %>%
  glimpse()



# plot

p <- ggplot(df.var2, aes(x=dummy, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'B', direction = -1) +
  labs(y = "Latitude", title = "Probability of rock") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())

p


###

## VRM ----

# set extent

sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

vrm <- rast(paste(sub.dir, "vrm_nc_all-mapped_300m_wInterp.tif", sep ='/'))
vrm
plot(vrm)

## aggregate to 2000 m --
fact <- 2000/300
vrm2000 <- aggregate(vrm, fact = fact, fun = mean, na.rm = T)
plot(vrm2000)

# make df --
df.var <- as.data.frame(vrm2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  group_by(y) %>%
  summarise(mean_value = mean(vrm_nc_all.mapped_30m, na.rm = T)) %>%
  mutate(dummy = 0) %>%
  glimpse()



# plot

p <- ggplot(df.var2, aes(x=dummy, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'B', direction = -1) +
  labs(y = "Latitude", title = "VRM") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())

p




# load temperature predictors ----

# Days_16C ----

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

# aggregate to 2000 m --
log_Days_16C2000 <- aggregate(Days_16C4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(log_Days_16C2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("Days_16C", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 8)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'H') +
  labs(y = "Latitude", title = "Days_16C in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p

###

# MHW  ----

mhw_up <- rast(paste(re.dir, 'Temperature', "MHW_Upwelling_Days.tif", sep ='/'))
mhw_up

# resample predictors to bathy ----
mhw_up3 <- resample(mhw_up, prob_rock)

# mask predictors to bathy ----
mhw_up4 <- mask(mhw_up3, prob_rock)
plot(mhw_up4)

# aggregate to 2000 m --
mhw_up4_2000 <- aggregate(mhw_up4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(mhw_up4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("MHW_Upwelling_Days", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 18)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'H') +
  labs(y = "Latitude", title = "MHW_Upwelling_Days in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p

###


# Max Anomaly Upwelling  ----

max_an <- rast(paste(re.dir, 'Temperature', "Max_Monthly_Anomaly_Upwelling_Temp.tif", sep ='/'))
max_an

# resample predictors to bathy ----
max_an3 <- resample(max_an, prob_rock)

# mask predictors to bathy ----
max_an4 <- mask(max_an3, prob_rock)
plot(max_an4)

# aggregate to 2000 m --
max_an4_2000 <- aggregate(max_an4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(max_an4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("Max_Monthly_Anomaly_Upwelling_Temp", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 34)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'H') +
  labs(y = "Latitude", title = "Max_Monthly_Anomaly_Upwelling_Temp in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p

###






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

# aggregate to 2000 m --
max_nit4_2000 <- aggregate(max_nit4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(max_nit4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("Max_Monthly_Nitrate", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 19)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'H') +
  labs(y = "Latitude", title = "Max_Monthly_Nitrate in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p

###

## Days_10N ----

Days_10N <- rast(paste(re.dir, "Nitrate", "Days_10N.tif", sep ='/'))
Days_10N

# # crop to NC --
# mean_sum_nit2 <- crop(mean_sum_nit, extent(e.nc))
# plot(mean_sum_nit2[[1]])

# resample predictors to bathy ----
Days_10N3 <- resample(Days_10N, prob_rock)

# mask predictors to bathy ----
Days_10N4 <- mask(Days_10N3, prob_rock)
plot(Days_10N4)

# aggregate to 2000 m --
Days_10N4_2000 <- aggregate(Days_10N4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(Days_10N4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("Days_10N", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 8)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'H') +
  labs(y = "Latitude", title = "Days_10N in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p

###


###

## Min_Nit ----

min_nit <- rast(paste(re.dir, "Nitrate", "Min_Monthly_Nitrate.tif", sep ='/'))
min_nit

# # crop to NC --
# mean_sum_nit2 <- crop(mean_sum_nit, extent(e.nc))
# plot(mean_sum_nit2[[1]])

# resample predictors to bathy ----
min_nit3 <- resample(min_nit, prob_rock)

# mask predictors to bathy ----
min_nit4 <- mask(min_nit3, prob_rock)
plot(min_nit4)

# aggregate to 2000 m --
min_nit4_2000 <- aggregate(min_nit4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(min_nit4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("Min_Monthly_Nitrate", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 19)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'H') +
  labs(y = "Latitude", title = "Min_Monthly_Nitrate in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p

###



## URCHINS ----

# load raster data --
u.dir <- paste(m.dir, "outputs_nc_rcca", "gam_urchins2", "sp_predictions", sep ='/')

n.files <- dir(u.dir)
# list files in source --
n.files <- list.files(u.dir, pattern = 'NC.tif', full.names = TRUE)
n.files
length(n.files)
# list names to load onto the Environment --
names.list <- list.files(u.dir, pattern = 'NC.tif')
names.list <- str_replace(names.list, "NC.tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)

# transform 
preds.stack2 <- exp(preds.stack)


# resample predictors to bathy ----
preds.stack3 <- resample(preds.stack2, prob_rock)

# mask predictors to bathy ----
preds.stack4 <- mask(preds.stack3, prob_rock)
plot(preds.stack4)

# aggregate to 2000 m --
preds.stack4_2000 <- aggregate(preds.stack4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(preds.stack4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("log_urchin_density", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 18)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'B', direction = -1) +
  labs(y = "Latitude", title = "Density of S. purpuratus in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p


### 


## KELP DENSITY ----


# load raster data --
k.dir <- paste(m.dir, "outputs_nc_rcca", "gam_V3", "sp_predictions", sep ='/')

n.files <- dir(k.dir)
# list files in source --
n.files <- list.files(k.dir, pattern = 'NC.tif', full.names = TRUE)
n.files
length(n.files)
# list names to load onto the Environment --
names.list <- list.files(k.dir, pattern = 'NC.tif')
names.list <- str_replace(names.list, "NC.tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)


# resample predictors to bathy ----
preds.stack3 <- resample(preds.stack, prob_rock)

# mask predictors to bathy ----
preds.stack4 <- mask(preds.stack3, prob_rock)
plot(preds.stack4)

# aggregate to 2000 m --
preds.stack4_2000 <- aggregate(preds.stack4, fact = fact, fun = mean, na.rm = T)


# make df --
df.var <- as.data.frame(preds.stack4_2000, xy = T)
head(df.var)
any(is.na(df.var))
names(df.var)

new.names.year <- paste("Nereo_density", 1998:2021, sep ='_')
new.names <- c('x', 'y', new.names.year)

names(df.var) <- new.names

# aggregate by latitude and calculate mean --
df.var2 <- df.var %>%
  pivot_longer(cols = 3:26, names_to = "variable_year", values_to = "value") %>%
  mutate_at(vars(x, y, variable_year), list(as.factor)) %>%
  group_by(y, variable_year) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  mutate(year = str_sub(variable_year, -4, -1)) %>%
  mutate(variable = substr(variable_year, 1, 18)) %>%
  mutate_at(vars(year), list(as.factor)) %>%
  glimpse()

levels(df.var2$variable_year)

# plot

p <- ggplot(df.var2, aes(x=year, y = y, fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis(option = 'D', direction = -1) +
  labs(y = "Latitude", title = "Density of N. luetkeana in a year") +
  theme(axis.text.x = element_text(angle = 45, h = 1, size = 10),
        axis.text.y = element_blank())

p


