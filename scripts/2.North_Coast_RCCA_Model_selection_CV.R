# 

# Script by Anita Giraldo - 2 May 2022
# last modified by Anita Giraldo - 24 May 2022


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


# ~ A. JUST RUN THE LOOP ####

## just to run the loop to get summary stats of best models --

# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('outputs_nc_rcca')
k.dir <- paste(d.dir, "gam_V4", sep ='/') # kelp model results
u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



# 1. PREPARE DATA ####

### 1.1. Load info on years RCCA ----
years <- read.csv(paste(rcca.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


### 1.2. Get sites with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(dd.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse() # Rows: 1,154


### 1.4. Get the sites for North Coast model ----
# sites selected in 'ncsites'
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


### 1.5. Choose variables and transform needed ----
# as per the FSSgam
names(df.nc)

dat1 <- df.nc %>%
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
    # Other vars
    #npp.mean,                             
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc
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
  # Other Env transformations
  #mutate(log_npp.mean = log(npp.mean + 1)) %>%
  #dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "all.mod.fits.csv", sep ='/')) 
head(best_mods)


### 2.2. Set submodel names ----
#no_submodels <- nrow(best_mods)
no_submodels <- 15

#submodels <- c('3.0', '3.1', '3.2', '3.3')
submodels <- paste(1:15)

submodel_names <- paste("bm", submodels, sep = '')


### 2.3. Set number of folds ----
no_folds <- 20


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
  dep <- 'log_den_NERLUE'
  
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
    mutate(model = purrr::map(train, ~ gam(bm_form, family = tw, data = .))) %>%
    glimpse()
  
  
  
  # 4.3. Add the null model  --
  cv_mods <- cv_mods %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model.null = purrr::map(train, ~ gam(log_den_NERLUE ~ 1, family = tw, data = .))) %>%
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
    pred_values <- cv_mods$model[[i]] %>% fitted()
    obs_values <- cv_mods$train[[i]] %>% select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df <- rbind(all.resid.df, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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
    obs_values <- cv_mods$predicted[[n]] %>% select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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

overall.summary <- na.omit(overall.summary)

#write.csv(overall.summary, paste(d.dir, "gam_V4", "best_models_summary_stats.csv", sep ='/'))



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

#ggsave(plot = p, filename = "best_models_summary_stats.png", device = "png", path = paste(d.dir, "gam_V4", sep ='/'))


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


##

best_mods$formula[1]
best_mods$formula[12]

best_mods$formula[11]
###

###

###



# ~ B. TO EXPLORE EACH SUBMODEL AND PLOT  ####

## gets more plots for each step of the process and is not a loop --


# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('outputs_nc_rcca')
k.dir <- paste(d.dir, "gam_V4", sep ='/') # kelp model results
u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"




# 1. PREPARE DATA ####

### 1.1. Load info on years RCCA ----
years <- read.csv(paste(rcca.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


### 1.2. Get sites with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(dd.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse() # Rows: 1,154


### 1.4. Get the sites for North Coast model ----
# sites selected in 'ncsites'
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


### 1.5. Choose variables and transform needed ----
# as per the FSSgam
names(df.nc)

dat1 <- df.nc %>%
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
    # Other vars
    #npp.mean,                             
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc
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
  # Other Env transformations
  #mutate(log_npp.mean = log(npp.mean + 1)) %>%
  #dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "best_models.csv", sep ='/')) 
head(best_mods)

# how many best models --
length(best_mods$X) # 4

### 2.2. Get variables from first model ----
bm_form <- best_mods$formula[3]
bm_name <- best_mods$modname[3]

class(bm_form)
class(bm_name)

bm_vars <- strsplit(bm_name, split = "+", fixed = T)
bm_vars
bm_vars[[1]]
length(bm_vars[[1]]) # 6

## add other variables if needed 

bm_vars[[1]] <- c(bm_vars[[1]], "mean_depth", "mean_prob_of_rock", "mean_slope")


### 2.3. Extract variables from the data ----
names(dat2)

dat3 <- dat2 %>%
  dplyr::select(site_name, year, transect, zone, log_den_NERLUE, bm_vars[[1]]) %>%
  glimpse()


### 2.4. Set model formula ----

# set dependent variable --
dep <- 'log_den_NERLUE'

# set predictors --
preds <- bm_vars[[1]]

preds2 <- character()

# get variables with smoothers --
for (i in 1:length(preds)) {
  pred.x <- preds[i]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2


# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form

gam.test <- gam(bm_form, family = tw, data = dat3)

summary(gam.test)
gam.check(gam.test)

par(mfrow=c(4,2),mar=c(2,4,3,1))
plot(gam.test)
dev.off()

gam.check(gam.test)



# 3. RUN CROSS VALIDATION ####

### 3.1. Divide data in k folds ----

# set number of folds --
no_folds <- 10

# split data into folds --
folds <- crossv_kfold(dat3, k = no_folds)
folds

# We can now run a model on the data referenced by each train object, 
# and validate the model results on each corresponding partition in test


### 3.2. Run model on each set of training data ----

cv_mods <- folds %>% 
  mutate(train = map(train, as_tibble)) %>%
  mutate(model = purrr::map(train, ~ gam(bm_form, family = tw, data = .))) %>%
  glimpse()



### 3.3. Add the null model  ----
cv_mods <- cv_mods %>% 
  mutate(train = map(train, as_tibble)) %>%
  mutate(model.null = purrr::map(train, ~ gam(log_den_NERLUE ~ 1, family = tw, data = .))) %>%
  glimpse()

# The result is a new model column containing fitted regression models based on each of the train data 
# (i.e., the whole data set excluding each partition).

# For example, the model fitted to our first set of training data is:
cv_mods$model[[1]] %>% summary()
cv_mods$model[[1]] %>% deviance()
cv_mods$model[[1]] %>% AIC()
cv_mods$model[[1]] %>% fitted()


cv_mods$model[[1]] %>% deviance()
cv_mods$model.null[[1]] %>% deviance()


### 3.4. Calculate % deviance explained ----
# https://stats.stackexchange.com/questions/497317/how-to-calculate-percent-partial-deviance-explained-by-each-predictor-variable-i
# % dev explained = dev(null model) - dev(best model) / dev(null model)

((cv_mods$model.null[[1]] %>% deviance())-(cv_mods$model[[1]] %>% deviance()))/ cv_mods$model.null[[1]] %>% deviance()

dev.exp.df <- data.frame()

for (i in 1:no_folds){
  dev.null <- deviance(cv_mods$model.null[[i]])
  dev.bm <- deviance(cv_mods$model[[i]])
  dev.exp <- (dev.null-dev.bm)/dev.null
  k.fold <- paste(i)
  dev.row <- cbind(k.fold, dev.null, dev.bm, dev.exp)
  dev.exp.df <- rbind(dev.exp.df, dev.row)
}

dev.exp.df

#### plot deviance explained ----

dev.exp.df

dev.exp.df %>% 
  mutate_at(vars(dev.null, dev.bm, dev.exp), list(as.numeric)) %>%
  ggplot(aes(x = k.fold, y=dev.exp, color = k.fold)) +
  geom_point(size = 5) +
  geom_hline(aes(yintercept = mean(dev.exp))) +
  coord_flip() +
  theme_bw()


### 3.5. Calulate residuals ----
# https://stackoverflow.com/questions/52754951/how-to-manually-calculate-the-residuals-of-linear-model-in-r


resid.df <- data.frame()
all.resid.df <- data.frame()

for (i in 1:no_folds){
  
  # get predicted and observed
  pred_values <- cv_mods$model[[i]] %>% fitted()
  obs_values <- cv_mods$train[[i]] %>% select(log_den_NERLUE)
  
  # calculate residuals
  all.res <- obs_values - pred_values # residuals
  
  # merge in data frame for ploting later
  all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
  all.resid.df <- rbind(all.resid.df, all.resid.df1)
  
  # calculate mean per k fold model
  mean.res <- mean(all.res$log_den_NERLUE)
  
  # calculate R-squared per k fold model
  all.rsq <- var(pred_values)/var(obs_values)
  
  resid.df1 <- cbind(mean.res, all.rsq)
  resid.df <- rbind(resid.df, resid.df1)
  
  
}

resid.df
all.resid.df
names(all.resid.df) <- c("pred_values",    "obs_values",     "residuals")
names(resid.df) <- c('mean_residuals', 'mean_R.Sqr')
row.names(resid.df) <- paste(1:10)

#### plot residuals ----
all.resid.df %>%
  ggplot(aes(obs_values, residuals)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  stat_smooth(method = "loess") +
  theme_minimal()


### 3.6. Merge dev.exp and residuals ----

summary.df <- cbind(dev.exp.df, resid.df)
summary.df

#### plot mean residuals for each k fold ----
summary.df %>%
  ggplot(aes(k.fold, mean_residuals, color = k.fold)) +
  #geom_hline(aes(yintercept = mean(mean_residuals))) +
  geom_point(size = 5) +
  theme_minimal()

#### plot R squared for each k fold ----
summary.df %>%
  ggplot(aes(k.fold, mean_R.Sqr, color = k.fold)) +
  geom_hline(aes(yintercept = mean(mean_R.Sqr))) +
  geom_point(size = 5) +
  theme_minimal()



### 3.7. Predict on test data ----

cv_mods <- cv_mods %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
  glimpse()

# extract relevant information from these predicted results --
predicted <- cv_mods %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
  unnest(predicted)
predicted

names(predicted)

resid.df.test <- data.frame()
all.resid.df.test <- data.frame()

for (i in 1:no_folds){
  
  # get predicted and observed
  pred_values <- cv_mods$predicted[[i]] %>% select(.fitted)
  obs_values <- cv_mods$predicted[[i]] %>% select(log_den_NERLUE)
  
  # calculate residuals
  all.res <- obs_values - pred_values # residuals
  
  # merge in data frame for ploting later
  all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
  all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
  
  # calculate mean per k fold model
  mean.res <- mean(all.res$log_den_NERLUE)
  
  # calculate R-squared per k fold model
  all.rsq <- var(pred_values)/var(obs_values)
  
  resid.df1 <- cbind(mean.res, all.rsq)
  resid.df.test <- rbind(resid.df.test, resid.df1)

  
}

head(resid.df.test)
names(resid.df.test) <- c("mean_residuals_test", "mean_R.sqr_test")
row.names(resid.df.test) <- paste(1:10)

head(all.resid.df.test)
names(all.resid.df.test) <- c("pred_values_test",    "obs_values_test",     "residuals_test")

#### plot predicted residuals --
library(ggplot2)
all.resid.df.test %>%
  ggplot(aes(obs_values_test, residuals_test)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  stat_smooth(method = "loess") +
  theme_minimal()


all.resid.df.test %>%
  ggplot(aes(obs_values_test, pred_values_test)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  stat_smooth(method = "loess") +
  theme_minimal()


### 3.8. Merge residuals from predicted and deviance explained ----

summary.df <- cbind(summary.df, resid.df.test)
summary.df

#### plot mean residuals for each k fold ----
summary.df %>%
  ggplot(aes(k.fold, mean_residuals_test, color = k.fold)) +
  #geom_hline(yintercept = 0) +
  geom_point(size = 5) +
  theme_minimal()



### 3.9. Calculae RMSE ----
# https://www.r-bloggers.com/2016/11/easy-cross-validation-in-r-with-modelr/


rmse.folds <- cv_mods %>%
  mutate(rmse = map2_dbl(model, test, rmse)) %>%
  select(.id, rmse) %>%
  glimpse()


### 3.10. Merge RMSE with deviance explained ----

summary.df <- cbind(summary.df, rmse.folds)
summary.df

#### plot mean residuals for each k fold ----
summary.df %>%
  ggplot(aes(.id, rmse, color = .id)) +
  #geom_hline(yintercept = 0) +
  geom_point(size = 5) +
  theme_minimal()




# 4. SUMMARY SUBMODEL STATISTICS ----

overall.summary <- data.frame()

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
  mutate(submodel = "bm1") %>%
  mutate_at(vars(submodel), list(as.factor)) %>%
  glimpse()
         
         

#### add to overall summary ----

overall.summary <- rbind(overall.summary, summary.df2)
overall.summary



###

###

###


# get variables of each submodel ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[4]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)
bm_vars[[1]]


###

###

###


# ~ C. COSTUME LOOP ----

# A. JUST RUN THE LOOP ####

## just to run the loop to get summary stats of best models --

# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('outputs_nc_rcca')
k.dir <- paste(d.dir, "gam_V4", sep ='/') # kelp model results
u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



# 1. PREPARE DATA ####

### 1.1. Load info on years RCCA ----
years <- read.csv(paste(rcca.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


### 1.2. Get sites with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(dd.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse() # Rows: 1,154


### 1.4. Get the sites for North Coast model ----
# sites selected in 'ncsites'
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


### 1.5. Choose variables and transform needed ----
# as per the FSSgam
names(df.nc)

dat1 <- df.nc %>%
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
    # Other vars
    #npp.mean,                             
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc
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
  # Other Env transformations
  #mutate(log_npp.mean = log(npp.mean + 1)) %>%
  #dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "all.mod.fits.csv", sep ='/')) 
head(best_mods)


### 2.2. Get submodels ----

# get submodel 5 --
bm_name <- best_mods$modname[5]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)
bm_vars

# submodel 5.1 --
bm_vars[[2]] <- c(bm_vars[[1]], "Days_10N", "log_den_STRPURAD", "log_den_STRPURAD")

# submodel 5.2 --
bm_vars[[3]] <- c(bm_vars[[1]], "Days_10N", "log_den_STRPURAD", "log_den_STRPURAD", "log_Days_16C")


# get submodel 12 --
bm_name <- best_mods$modname[12]

# get submodel variables --
bm_vars12 <- strsplit(bm_name, split = "+", fixed = T)
bm_vars12

# submodel 12 --
bm_vars[[4]] <- c(bm_vars12[[1]])

# submodel 12.1 --
bm_vars[[5]] <- c(bm_vars12[[1]], "mean_depth", "mean_prob_of_rock", "mean_slope")

# submodel 12.2 --
bm_vars[[6]] <- c(bm_vars12[[1]], "mean_depth", "mean_slope", "log_mean_vrm")

# submodel 12.3 --
bm_vars[[7]] <- c(bm_vars12[[1]], "mean_depth", "mean_slope")

# submodel 12.4 --
bm_vars[[8]] <- c(bm_vars12[[1]], "mean_depth", "log_mean_vrm")


bm_vars

### 2.3. Set submodel names ----
no_submodels <- length(bm_vars)

submodels <- c('5', '5.1', '5.2', '12', '12.1', '12.2', '12.3', '12.4')

submodel_names <- paste("bm", submodels, sep = '')


### 2.4. Set number of folds ----
no_folds <- 20


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
  dep <- 'log_den_NERLUE'
  
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
    mutate(model = purrr::map(train, ~ gam(bm_form, family = tw, data = .))) %>%
    glimpse()
  
  # diagnostics --
  print("step 4.2 completed")
  
  # 4.3. Add the null model  --
  cv_mods <- cv_mods %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model.null = purrr::map(train, ~ gam(log_den_NERLUE ~ 1, family = tw, data = .))) %>%
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
    obs_values <- cv_mods$train[[l]] %>% dplyr::select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df <- rbind(all.resid.df, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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
    obs_values <- cv_mods$predicted[[n]] %>% dplyr::select(log_den_NERLUE)

    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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

#write.csv(overall.summary, paste(d.dir, "gam_V4", "best_models_summary_stats_added_vars.csv", sep ='/'))



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

bm_vars[6]



###

###

###

# ~ D. JUST RUN THE LOOP ####

## just to run the loop to get summary stats of best models --

# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('outputs_nc_rcca')
k.dir <- paste(d.dir, "gam_V4", sep ='/') # kelp model results
u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



# 1. PREPARE DATA ####

### 1.1. Load info on years RCCA ----
years <- read.csv(paste(rcca.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


### 1.2. Get sites with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(dd.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse() # Rows: 1,154


### 1.4. Get the sites for North Coast model ----
# sites selected in 'ncsites'
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


### 1.5. Choose variables and transform needed ----
# as per the FSSgam
names(df.nc)

dat1 <- df.nc %>%
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
    # Other vars
    #npp.mean,                             
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc
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
  # Other Env transformations
  #mutate(log_npp.mean = log(npp.mean + 1)) %>%
  #dplyr::select(-c(npp.mean)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "all.mod.fits.csv", sep ='/')) 
head(best_mods)
nrow(best_mods)


### 2.2. Set submodel names ----

# take the best 15 models --
no_submodels <- 15

#submodels <- c('3.0', '3.1', '3.2', '3.3')
submodels <- paste(1:no_submodels)

submodel_names <- paste("bm", submodels, sep = '')


### 2.3. Set number of folds ----
no_folds <- 20


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
  dep <- 'log_den_NERLUE'
  
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
    mutate(model = purrr::map(train, ~ gam(bm_form, family = tw, data = .))) %>%
    glimpse()
  
  
  
  # 4.3. Add the null model  --
  cv_mods <- cv_mods %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model.null = purrr::map(train, ~ gam(log_den_NERLUE ~ 1, family = tw, data = .))) %>%
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
    pred_values <- cv_mods$model[[i]] %>% fitted()
    obs_values <- cv_mods$train[[i]] %>% select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df <- rbind(all.resid.df, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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
    obs_values <- cv_mods$predicted[[n]] %>% select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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

#write.csv(overall.summary, paste(d.dir, "gam_V4", "best_models_summary_stats.csv", sep ='/'))



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

#ggsave(plot = p, filename = "best_models_summary_stats.png", device = "png", path = paste(d.dir, "gam_V4", sep ='/'))



# 6. Check the variables in models that perform better ----
# best submodels 3, 5, 6, 9, 10, 11, 12, 13
# best models that include substrate variables 5, 11, 13

bm_name <- best_mods$modname[6]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)
bm_vars


### Plot just best models ----
p <- ggplot(overall.summary %>% dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                                                summary_stat != 'mean_residuals') %>%
              filter(submodel == 'bm5'|
                       submodel == 'bm6'|
                       submodel == 'bm9'|
                       submodel == 'bm10'|
                       submodel == 'bm11'|
                       submodel == 'bm12'|
                       submodel == 'bm13'), 
            aes(x = submodel, y = mean_stat, color = submodel)) +
  geom_point(size = 5) +
  facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
p


p <- ggplot(overall.summary %>% dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                                                summary_stat != 'mean_residuals') %>%
              filter( submodel == 'bm5'|
                       submodel == 'bm11'|
                       submodel == 'bm13'), 
            aes(x = submodel, y = mean_stat, color = submodel)) +
  geom_point(size = 5) +
  facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
p


###

###

###




# ~ C. COSTUME LOOP ----

## just to run the loop to get summary stats of best models --


### 2.1. Load best model csv ----
best_mods <- read.csv(paste(k.dir, "all.mod.fits.csv", sep ='/')) 
head(best_mods)
nrow(best_mods)

### 2.2. Get submodels ----

# get submodel 3 --
bm_name <- best_mods$modname[5]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)
bm_vars

bm_vars[[2]] <- c(bm_vars[[1]], "log_den_STRPURAD")

# submodel 3.1 --
bm_name <- best_mods$modname[11]
bm_vars2 <- strsplit(bm_name, split = "+", fixed = T)
bm_vars2

bm_vars[[3]] <- bm_vars2[[1]]

bm_vars[[4]] <- c(bm_vars2[[1]], "Mean_Monthly_Upwelling_Temp")

bm_vars[[5]] <- c(bm_vars2[[1]], "log_Days_16C")

# submodel 3.2 --
bm_name <- best_mods$modname[13]
bm_vars3 <- strsplit(bm_name, split = "+", fixed = T)
bm_vars3


bm_vars[[6]] <- bm_vars3[[1]]

# submodel 3.2 --
bm_vars[[7]] <- c(bm_vars3[[1]] , "log_den_STRPURAD")


### 2.3. Set submodel names ----
no_submodels <- length(bm_vars)

#submodels <- paste(1:nrow(best_mods))
submodels <- c(5, 5.1, 11, 11.1, 11.2, 13, 13.1)

submodel_names <- paste("bm", submodels, sep = '')


### 2.4. Set number of folds ----
no_folds <- 20


### 2.5. Set data frame to save summary stats ----
overall.summary2 <- data.frame()


# 3. LOOP ----

for (i in 1:no_submodels){
  
  # 3.1. Get submodel name --
  #bm_name <- best_mods$modname[i]
  
  # get submodel variables --
  #bm_vars <- strsplit(bm_name, split = "+", fixed = T)
  
  # 3.2. Set model formula --
  
  # set dependent variable --
  dep <- 'log_den_NERLUE'
  
  # set predictors --
  preds <- bm_vars[[i]]
  
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
    mutate(model = purrr::map(train, ~ gam(bm_form, family = tw, data = .))) %>%
    glimpse()
  
  
  
  # 4.3. Add the null model  --
  cv_mods <- cv_mods %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model.null = purrr::map(train, ~ gam(log_den_NERLUE ~ 1, family = tw, data = .))) %>%
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
    obs_values <- cv_mods$train[[l]] %>% select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df <- rbind(all.resid.df, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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
    obs_values <- cv_mods$predicted[[n]] %>% select(log_den_NERLUE)
    
    # calculate residuals
    all.res <- obs_values - pred_values # residuals
    
    # merge in data frame for ploting later
    all.resid.df1 <- cbind(pred_values, obs_values = obs_values$log_den_NERLUE, all.res)
    all.resid.df.test <- rbind(all.resid.df.test, all.resid.df1)
    
    # calculate mean per k fold model
    mean.res <- mean(all.res$log_den_NERLUE)
    
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
  
  overall.summary2 <- rbind(overall.summary2, summary.df2)
  overall.summary2
  
  
}

beep()

# 4. SAVE SUMMARY ----
head(overall.summary2)

#write.csv(overall.summary, paste(d.dir, "gam_V4", "best_models_summary_stats.csv", sep ='/'))



# 5. PLOT SUMMARY STATS FOR BEST MODELS ----

p <- ggplot(overall.summary2 %>% dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                                                summary_stat != 'mean_residuals'), 
            aes(x = submodel, y = mean_stat, color = submodel)) +
  geom_point(size = 5) +
  facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
p


# save plot --

#ggsave(plot = p, filename = "best_models_summary_stats.png", device = "png", path = paste(d.dir, "gam_V4", sep ='/'))