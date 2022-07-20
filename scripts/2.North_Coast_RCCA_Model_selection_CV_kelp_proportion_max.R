# 

# Script by Anita Giraldo - 2 May 2022
# last modified by Anita Giraldo - 19 July 2022 - new as suggested by Tom


## This script takes the best models obtained from FSSgam and selects the best one using CROSS VALIDATION 

# resources:
# https://drsimonj.svbtle.com/k-fold-cross-validation-with-modelr-and-broom
# https://www.p8105.com/cross_validation.html
# https://rpubs.com/dgrtwo/cv-modelr
# https://www.r-bloggers.com/2016/11/easy-cross-validation-in-r-with-modelr/



# libraries ----
library(dplyr)
library(purrr)
library(modelr)
library(broom)
library(tidyr)
library(here)
library(mgcv)
library(ggplot2)
library(viridis)
library(beepr)


###

###


# ~ A. JUST RUN THE LOOP ----


## just to run the loop to get summary stats of best models --

# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('data')
o.dir <- here('outputs_nc_rcca')
k.dir <- paste(o.dir, "gam_prop_max_V1", sep ='/') # kelp model results
#u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
#rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
#dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"




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


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "prop_max_best_models.csv", sep ='/')) 
head(best_mods)
nrow(best_mods) # 25


### 2.2. Set submodel names ----
#no_submodels <- nrow(best_mods)
no_submodels <- nrow(best_mods)

#submodels <- c('3.0', '3.1', '3.2', '3.3')
submodels <- paste(1:no_submodels)

submodel_names <- paste("bm", submodels, sep = '')


### 2.3. Set number of folds ----
no_folds <- 15


### 2.4. Set data frame to save summary stats ----
overall.summary <- data.frame()


# 3. LOOP ----

for (i in 1:no_submodels){
  
  # 3.1. Get submodel name --
  bm_name <- best_mods$modname[i]
  
  # get submodel variables --
  bm_vars <- strsplit(bm_name, split = "+", fixed = T)
  
  # 3.2. Set model formula --
  
  # set dependent variable --
  dep <- 'prop.max_NERLUE'
  
  # set predictors --
  preds <- bm_vars[[1]]
  
  # get variables with smoothers --
  
  preds2 <- character()
  
  for (j in 1:length(preds)) {
    pred.x <- preds[j]
    pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
    preds2 <- append(preds2, pred_form)
  }
  
  preds2
  
  
  # get random variables --
  random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"
  
  
  # set model formula --
  bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
  bm_form
  
  # 4. RUN CROSS VALIDATION --
  
  # 4.1. Divide data in k folds --
  
  folds <- crossv_kfold(dat2, k = no_folds)
  folds
  
  
  # 4.2. Run model on each set of training data --
  
  cv_mods <- folds %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model = purrr::map(train, ~ gam(bm_form, family = betar, data = .))) %>%
    glimpse()
  
  
  
  # 4.3. Add the null model  --
  cv_mods <- cv_mods %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model.null = purrr::map(train, ~ gam(prop.max_NERLUE ~ 1, family = betar, data = .))) %>%
    glimpse()
  
  
  
  # 4.4. Calculate % deviance explained --
  
  # get empty df
  dev.exp.df <- data.frame()
  
  # run loop 
  for (m in 1:no_folds){
    dev.null <- deviance(cv_mods$model.null[[m]])
    dev.bm <- deviance(cv_mods$model[[m]])
    dev.exp <- (dev.null-dev.bm)/dev.null
    k.fold <- paste(m)
    dev.row <- cbind(k.fold, dev.null, dev.bm, dev.exp)
    dev.exp.df <- rbind(dev.exp.df, dev.row)
  }
  
  dev.exp.df
  
  # diagnostics --
  print("step 4.4 completed")
  
  # 4.5. Calulate residuals --
  
  # get empty dfs
  resid.df <- data.frame()
  all.resid.df <- data.frame()
  
  for (l in 1:no_folds){
    
    # get predicted and observed
    pred_values <- cv_mods$model[[l]] %>% fitted()
    obs_values <- cv_mods$train[[l]] %>% select(prop.max_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$prop.max_NERLUE, all.res)
    all.resid.df <- rbind(all.resid.df, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$prop.max_NERLUE)
    
    # calculate R-squared per k fold model
    all.rsq <- var(pred_values)/var(obs_values)
    
    resid.df1 <- cbind(mean.res, all.rsq)
    resid.df <- rbind(resid.df, resid.df1)
    
    
  }
  
  resid.df
  all.resid.df
  names(all.resid.df) <- c("pred_values",    "obs_values",     "residuals")
  names(resid.df) <- c('mean_residuals', 'mean_R.Sqr')
  row.names(resid.df) <- paste(1:no_folds)
  
  # diagnostics --
  print("step 4.5 completed")
  
  
  # 4.6. Merge dev.exp and residuals --
  
  summary.df <- cbind(dev.exp.df, resid.df)
  summary.df
  
  
  # 5. PREDICT ON TEST DATA --
  
  cv_mods <- cv_mods %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
    glimpse()
  
  # extract relevant information from these predicted results --
  predicted <- cv_mods %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
    unnest(predicted)
  predicted
  
  names(predicted)
  
  resid.df.test <- data.frame()
  all.resid.df.test <- data.frame()
  
  # 5.1. Get residuals and R squared --
  
  for (n in 1:no_folds){
    
    # get predicted and observed
    pred_values <- cv_mods$predicted[[n]] %>% select(.fitted)
    obs_values <- cv_mods$predicted[[n]] %>% select(prop.max_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$prop.max_NERLUE, all.res)
    all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$prop.max_NERLUE)
    
    # calculate R-squared per k fold model
    all.rsq <- var(pred_values)/var(obs_values)
    
    resid.df1 <- cbind(mean.res, all.rsq)
    resid.df.test <- rbind(resid.df.test, resid.df1)
    
    
  }
  
  head(resid.df.test)
  names(resid.df.test) <- c("mean_residuals_test", "mean_R.sqr_test")
  row.names(resid.df.test) <- paste(1:no_folds)
  
  head(all.resid.df.test)
  names(all.resid.df.test) <- c("pred_values_test",    "obs_values_test",     "residuals_test")
  
  
  
  # 5.2. Merge residuals from predicted and deviance explained --
  
  summary.df <- cbind(summary.df, resid.df.test)
  summary.df
  
  
  
  
  # 5.3. Calculae RMSE --
  
  rmse.folds <- cv_mods %>%
    mutate(rmse = map2_dbl(model, test, rmse)) %>%
    select(.id, rmse) %>%
    glimpse()
  
  
  # 5.4. Merge RMSE with deviance explained --
  
  summary.df <- cbind(summary.df, rmse.folds)
  summary.df
  
  
  
  # 6. SUMMARY SUBMODEL STATISTICS --
  
  
  
  glimpse(summary.df)
  
  summary.df2 <- summary.df %>%
    rename(id = .id) %>%
    relocate(id) %>%
    dplyr::select(-c(k.fold, dev.null, dev.bm)) %>%
    mutate_at(vars(dev.exp), list(as.numeric)) %>%
    pivot_longer(cols = dev.exp:rmse, names_to = "summary_stat", values_to = "values") %>%
    mutate_at(vars(id, summary_stat), list(as.factor)) %>%
    group_by(summary_stat) %>%
    summarise(mean_stat = mean(values)) %>%
    mutate(submodel = submodel_names[i]) %>%
    mutate_at(vars(submodel), list(as.factor)) %>%
    glimpse()
  
  
  
  # 6.1. add to overall summary --
  
  overall.summary <- rbind(overall.summary, summary.df2)
  overall.summary
  
  
}

beep()

# 4. SAVE SUMMARY ----
head(overall.summary)
nrow(overall.summary)

overall.summary <- na.omit(overall.summary)

model_no <- "V1"

write.csv(overall.summary, paste(k.dir, "best_models_CV_summary_stats.csv", sep ='/'))



# 5. PLOT SUMMARY STATS FOR BEST MODELS ----

p <- ggplot(overall.summary %>% dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                                                summary_stat != 'mean_residuals'), 
            aes(x = submodel, y = mean_stat, color = submodel)) +
  geom_point(size = 5) +
  facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

p


# save plot --

ggsave(plot = p, filename = "best_models_CV_summary_stats.png", device = "png", path = paste(k.dir, sep ='/'))


# 6. Compute best statistics ----

sum.sum <- overall.summary %>%
  dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                  summary_stat != 'mean_residuals') %>%
  pivot_wider(id_cols = submodel, names_from = summary_stat, values_from = mean_stat) %>%
  mutate(sum_stat = (dev.exp + mean_R.Sqr)*10 - (mean_residuals_test + rmse)) %>%
  mutate(sum_max = (dev.exp + mean_R.Sqr)) %>%
  mutate(sum_min = (mean_residuals_test + rmse)) %>%
  glimpse()

ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_stat, color = sum_stat)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_max, color = sum_max)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_min, color = sum_min)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


###

###

###


# ~ C. COSTUME LOOP ----


# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('data')
o.dir <- here('outputs_nc_rcca')
k.dir <- paste(o.dir, "gam_prop_max_V1", sep ='/') # kelp model results
#u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
#rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
#dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



## 1. Prepare the data ----

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


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "prop_max_all.mod.fits.csv", sep ='/')) 
head(best_mods)


### 2.2. Get submodels ----

# get submodel 2.1 --
bm_name <- best_mods$modname[2]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)
bm_vars

# submodel 2.2 --
bm_vars[[2]] <- c(bm_vars[[1]], "log_UBR_Mean")

# submodel 2.3 --
bm_vars[[3]] <- c(bm_vars[[1]], "log_UBR_Mean", "mean_prob_of_rock")


# get submodel 10 --
bm_name <- best_mods$modname[10]

# get submodel variables --
bm_vars10 <- strsplit(bm_name, split = "+", fixed = T)
bm_vars10

# submodel 10.1 --
bm_vars[[4]] <- c(bm_vars10[[1]])

# submodel 10.2 --
bm_vars[[5]] <- c(bm_vars10[[1]], "log_UBR_Mean")

# submodel 10.3 --
bm_vars[[6]] <- c(bm_vars10[[1]], "log_UBR_Mean", "mean_prob_of_rock")




# get submodel 18 --
bm_name <- best_mods$modname[18]

# get submodel variables --
bm_vars18 <- strsplit(bm_name, split = "+", fixed = T)
bm_vars18

# submodel 18.1 --
bm_vars[[7]] <- c(bm_vars18[[1]])

# submodel 18.2 --
bm_vars[[8]] <- c(bm_vars18[[1]], "log_UBR_Mean")

# submodel 10.3 --
bm_vars[[9]] <- c(bm_vars18[[1]], "log_UBR_Mean", "mean_prob_of_rock")



# get submodel 24 --
bm_name <- best_mods$modname[24]

# get submodel variables --
bm_vars24 <- strsplit(bm_name, split = "+", fixed = T)
bm_vars24

# submodel 24.1 --
bm_vars[[10]] <- c(bm_vars24[[1]])

# submodel 24.2 --
bm_vars[[11]] <- c(bm_vars24[[1]], "log_UBR_Mean")

# submodel 24.3 --
bm_vars[[12]] <- c(bm_vars24[[1]], "log_UBR_Mean", "mean_prob_of_rock")







### 2.3. Set submodel names ----

submodels <- c('2.1', '2.2', '2.3', '10.1', '10.2', '10.3', '18.1', '18.2', '18.3', '24.1', '24.2', '24.3')

no_submodels <- length(submodels)

submodel_names <- paste("bm", submodels, sep = '')


### 2.4. Set number of folds ----
no_folds <- 15


### 2.5. Set data frame to save summary stats ----
overall.summary <- data.frame()


# 3. LOOP ----

for (i in 1:no_submodels){
  
  # 3.1. Get submodel name --
  #bm_name <- best_mods$modname[i]
  
  # get submodel variables --
  #bm_vars <- strsplit(bm_name, split = "+", fixed = T)
  
  # 3.2. Set model formula --
  
  # set dependent variable --
  dep <- 'prop.max_NERLUE'
  
  # set predictors --
  preds <- bm_vars[[i]]
  
  # get variables with smoothers --
  
  preds2 <- character()
  
  for (j in 1:length(preds)) {
    pred.x <- preds[j]
    pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
    preds2 <- append(preds2, pred_form)
  }
  
  # diagnostics --
  print("step 3.2 completed")
  preds2
  
  
  # get random variables --
  random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"
  
  
  # set model formula --
  bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
  bm_form
  
  # 4. RUN CROSS VALIDATION --
  
  # 4.1. Divide data in k folds --
  
  folds <- crossv_kfold(dat2, k = no_folds)
  folds
  
  # diagnostics --
  print("step 4.1 completed")
  
  # 4.2. Run model on each set of training data --
  
  cv_mods <- folds %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model = purrr::map(train, ~ gam(bm_form, family = betar, data = .))) %>%
    glimpse()
  
  # diagnostics --
  print("step 4.2 completed")
  
  # 4.3. Add the null model  --
  cv_mods <- cv_mods %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model.null = purrr::map(train, ~ gam(prop.max_NERLUE ~ 1, family = betar, data = .))) %>%
    glimpse()
  
  # diagnostics --
  print("step 4.3 completed")
  
  # 4.4. Calculate % deviance explained --
  
  # get empty df
  dev.exp.df <- data.frame()
  
  # run loop 
  for (m in 1:no_folds){
    dev.null <- deviance(cv_mods$model.null[[m]])
    dev.bm <- deviance(cv_mods$model[[m]])
    dev.exp <- (dev.null-dev.bm)/dev.null
    k.fold <- paste(m)
    dev.row <- cbind(k.fold, dev.null, dev.bm, dev.exp)
    dev.exp.df <- rbind(dev.exp.df, dev.row)
  }
  
  dev.exp.df
  
  # diagnostics --
  print("step 4.4 completed")
  
  # 4.5. Calulate residuals --
  
  # get empty dfs
  resid.df <- data.frame()
  all.resid.df <- data.frame()
  
  for (l in 1:no_folds){
    
    # get predicted and observed
    pred_values <- cv_mods$model[[l]] %>% fitted()
    obs_values <- cv_mods$train[[l]] %>% dplyr::select(prop.max_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$prop.max_NERLUE, all.res)
    all.resid.df <- rbind(all.resid.df, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$prop.max_NERLUE)
    
    # calculate R-squared per k fold model
    all.rsq <- var(pred_values)/var(obs_values)
    
    resid.df1 <- cbind(mean.res, all.rsq)
    resid.df <- rbind(resid.df, resid.df1)
    
    
  }
  
  resid.df
  all.resid.df
  names(all.resid.df) <- c("pred_values",    "obs_values",     "residuals")
  names(resid.df) <- c('mean_residuals', 'mean_R.Sqr')
  row.names(resid.df) <- paste(1:no_folds)
  
  # diagnostics --
  print("step 4.5 completed")
  
  
  # 4.6. Merge dev.exp and residuals --
  
  summary.df <- cbind(dev.exp.df, resid.df)
  summary.df
  
  
  # 5. PREDICT ON TEST DATA --
  
  cv_mods <- cv_mods %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
    glimpse()
  
  # extract relevant information from these predicted results --
  predicted <- cv_mods %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
    unnest(predicted)
  predicted
  
  names(predicted)
  
  resid.df.test <- data.frame()
  all.resid.df.test <- data.frame()
  
  # diagnostics --
  print("step 5 completed")
  
  # 5.1. Get residuals and R squared --
  
  for (n in 1:no_folds){
    
    # get predicted and observed
    pred_values <- cv_mods$predicted[[n]] %>% dplyr::select(.fitted)
    obs_values <- cv_mods$predicted[[n]] %>% dplyr::select(prop.max_NERLUE)
    
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$prop.max_NERLUE, all.res)
    all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$prop.max_NERLUE)
    
    # calculate R-squared per k fold model
    all.rsq <- var(pred_values)/var(obs_values)
    
    resid.df1 <- cbind(mean.res, all.rsq)
    resid.df.test <- rbind(resid.df.test, resid.df1)
    
    
  }
  
  head(resid.df.test)
  names(resid.df.test) <- c("mean_residuals_test", "mean_R.sqr_test")
  row.names(resid.df.test) <- paste(1:no_folds)
  
  head(all.resid.df.test)
  names(all.resid.df.test) <- c("pred_values_test",    "obs_values_test",     "residuals_test")
  
  
  
  # 5.2. Merge residuals from predicted and deviance explained --
  
  summary.df <- cbind(summary.df, resid.df.test)
  summary.df
  
  
  
  
  # 5.3. Calculae RMSE --
  
  rmse.folds <- cv_mods %>%
    mutate(rmse = map2_dbl(model, test, rmse)) %>%
    dplyr::select(.id, rmse) %>%
    glimpse()
  
  
  # 5.4. Merge RMSE with deviance explained --
  
  summary.df <- cbind(summary.df, rmse.folds)
  summary.df
  
  
  
  # 6. SUMMARY SUBMODEL STATISTICS --
  
  
  
  glimpse(summary.df)
  
  summary.df2 <- summary.df %>%
    rename(id = .id) %>%
    relocate(id) %>%
    dplyr::select(-c(k.fold, dev.null, dev.bm)) %>%
    mutate_at(vars(dev.exp), list(as.numeric)) %>%
    pivot_longer(cols = dev.exp:rmse, names_to = "summary_stat", values_to = "values") %>%
    mutate_at(vars(id, summary_stat), list(as.factor)) %>%
    group_by(summary_stat) %>%
    summarise(mean_stat = mean(values)) %>%
    mutate(submodel = submodel_names[i]) %>%
    mutate_at(vars(submodel), list(as.factor)) %>%
    glimpse()
  
  
  
  # 6.1. add to overall summary --
  
  overall.summary <- rbind(overall.summary, summary.df2)
  overall.summary
  
  
}

beep()

# 4. SAVE SUMMARY ----
head(overall.summary)
nrow(overall.summary)


#write.csv(overall.summary, paste(k.dir, "best_models_CV_summary_stats_added_vars.csv", sep ='/'))



# 5. PLOT SUMMARY STATS FOR BEST MODELS ----

p <- ggplot(overall.summary %>% dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                                                summary_stat != 'mean_residuals'), 
            aes(x = submodel, y = mean_stat, color = submodel)) +
  geom_point(size = 5) +
  facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
p


# save plot --

#ggsave(plot = p, filename = "best_models_summary_stats_added_vars.png", device = "png", path = paste(d.dir, "gam_V4", sep ='/'))


# 6. Compute best statistics ----

sum.sum <- overall.summary %>%
  dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                  summary_stat != 'mean_residuals') %>%
  pivot_wider(id_cols = submodel, names_from = summary_stat, values_from = mean_stat) %>%
  mutate(sum_stat = (dev.exp + mean_R.Sqr) - (mean_residuals_test + rmse)/10) %>%
  mutate(sum_max = (dev.exp + mean_R.Sqr)) %>%
  mutate(sum_min = (mean_residuals_test + rmse)) %>%
  glimpse()


ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_stat, color = sum_stat)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_max, color = sum_max)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_min, color = sum_min)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))


##

bm_name2 <- best_mods$modname[2]
bm_name2

bm_name10 <- best_mods$modname[10]
bm_name10

bm_name18 <- best_mods$modname[18]
bm_name18

bm_name24 <- best_mods$modname[24]
bm_name24

###

###

## ~D. TEST MANUALLY ----

library(visreg)
library(caret)
library(ggpmisc)

# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('data')
o.dir <- here('outputs_nc_rcca')
k.dir <- paste(o.dir, "gam_prop_max_V1", sep ='/') # kelp model results
#u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
#rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
#dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



## 1. Prepare the data ----

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


# 2. Divide data into train and test ----

inTraining <- createDataPartition(dat2$prop.max_NERLUE, p = 0.75, list = FALSE)
train.gam <- dat2[ inTraining,]
test.gam  <- dat2[-inTraining,]



## BM2.1 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[1]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam2.1 <- gam(bm_form, family = betar, data = train.gam)
summary(gam2.1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2.1)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                wh_95prc,
                mean_waveyear,
                log_UBR_Max, 
                site_name, zone, year)

fits <- predict.gam(gam2.1, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM2.2 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[2]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam2.2 <- gam(bm_form, family = betar, data = train.gam)
summary(gam2.2)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2.2)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                wh_95prc,
                mean_waveyear,
                log_UBR_Mean, 
                site_name, zone, year)

fits <- predict.gam(gam2.2, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM2.3 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[3]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam2.3 <- gam(bm_form, family = betar, data = train.gam)
summary(gam2.3)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2.3)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                wh_95prc,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam2.3, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM10.1 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[4]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam10.1 <- gam(bm_form, family = betar, data = train.gam)
summary(gam10.1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam10.1)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam10.1, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM10.2 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[5]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam10.2 <- gam(bm_form, family = betar, data = train.gam)
summary(gam10.2)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam10.2)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam10.2, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM10.3 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[6]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam10.3 <- gam(bm_form, family = betar, data = train.gam)
summary(gam10.3)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam10.3)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam10.3, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM18.1 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[7]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam18.1 <- gam(bm_form, family = betar, data = train.gam)
summary(gam18.1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam18.1)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam18.1, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM18.2 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[8]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam18.2 <- gam(bm_form, family = betar, data = train.gam)
summary(gam18.2)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam18.2)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam18.2, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM18.3 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[9]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam18.3 <- gam(bm_form, family = betar, data = train.gam)
summary(gam18.3)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam18.3)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam18.3, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM24.1 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[10]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam24.1 <- gam(bm_form, family = betar, data = train.gam)
summary(gam24.1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam24.1)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam24.1, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


##

## BM24.2 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[11]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam24.2 <- gam(bm_form, family = betar, data = train.gam)
summary(gam24.2)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam24.2)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam24.2, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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



##

## BM24.3 ----

dep <- 'prop.max_NERLUE'

# set predictors --
preds <- bm_vars[[12]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2

# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"

# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form


gam24.3 <- gam(bm_form, family = betar, data = train.gam)
summary(gam24.3)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam24.3)
dev.off()


## Check on Test ----

testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                log_UBR_Max, 
                mean_prob_of_rock,
                site_name, zone, year)

fits <- predict.gam(gam24.3, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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


###


### BEST MODEL ----


testdata <- test.gam %>%
  dplyr::select(prop.max_NERLUE, 
                log_den_STRPURAD, 
                Mean_Monthly_Summer_Temp,
                Min_Monthly_Temp,
                Max_Monthly_Nitrate,  
                Mean_Monthly_NPP,
                log_Min_Monthly_NPP,
                log_Mean_Monthly_NPP_Upwelling,
                Mean_Monthly_Upwelling_Temp, 
                Max_Monthly_NPP_Upwelling,
                wh_95prc,
                wh_max,
                mean_waveyear,
                log_UBR_Mean, 
                log_UBR_Max, 
                mean_prob_of_rock,
                log_den_PYCHEL,
                site_name, zone, year)

## BM 2.1 ----

summary(gam2.1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2.1)
dev.off()

fits <- predict.gam(gam2.1, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

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

gam2.1 <- gam(prop.max_NERLUE ~ s(log_den_STRPURAD, k = 6, bs = "cr") +
                 #s(log_den_PYCHEL, k = 4, bs = "cr") +
                 s(log_Mean_Monthly_NPP_Upwelling, k = 6, bs = "cr") + 
                 s(log_Min_Monthly_NPP, k = 6, bs = "cr") + 
                 s(Max_Monthly_Nitrate, k = 5, bs = "cr") + 
                 s(Mean_Monthly_Summer_Temp,  k = 5, bs = "cr") + 
                 s(wh_95prc, k = 3, bs = "cr") + 
                 #s(wh_max, k = 4, bs = "cr") + 
                 #s(mean_prob_of_rock, k = 4, bs = "cr") +
                 #s(log_UBR_Max, k = 6, bs = "cr") +
                 s(site_name, zone, bs = "re") + s(year, bs = "re"),
                 family = betar(), data = train.gam)

summary(gam2.1)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(gam2.1)
dev.off()

fits <- predict.gam(gam2.1, newdata=testdata, type='response', se.fit=T)

## plot pred.vs obs ----

predicts.all <-  testdata %>% data.frame(fits)%>%
  #group_by(survey_year)%>% #only change here
  #summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>%
  glimpse()

my.formula <- y ~ x
p <- ggplot(predicts.all, aes(x = fit, y = prop.max_NERLUE)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1, color = 'red'), show.legend = F) +
  labs(x='Predicted', y='Observed', title='Proportion of maximum density N. luetkeana') +
  theme_bw()
p


## bar plot ----

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
