# 

# Script by Anita Giraldo - 2 May 2022
# last modified by Anita Giraldo - 12 April 2023


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


# ~ C. COSTUME LOOP ----


## just to run the loop to get summary stats of best models --

# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('data')
o.dir <- here('outputs_nc_rcca')
k.dir <- paste(o.dir, "gam_V7", sep ='/') # kelp model results
cv.dir <- paste(o.dir, "new_cvs", sep ='/')
#u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
#rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
#dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



# 1. PREPARE DATA ####

### 1.1. Load info on years RCCA ----
years <- read.csv(paste(d.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


### 1.2. Get sites with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(total.years > 2) %>%
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  mutate(zone_new = case_when(
    transect == '1' ~ 'OUTER',
    transect == '2' ~ 'OUTER', 
    transect == '3' ~ 'OUTER', 
    transect == '4' ~ 'INNER',
    transect == '5' ~ 'INNER',
    transect == '6' ~ 'INNER')) %>%
  dplyr::select(-zone) %>%
  rename(zone = zone_new) %>%
  mutate_at(vars(zone), list(as.factor)) %>%
  relocate(zone, .after = transect) %>%
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
    #wh.95 ,   wh.max,
    #npgo_mean , mei_mean,
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
  dplyr::select(-c(UBR_Mean, UBR_Max)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)


# 2. LOAD BEST MODELS ####

### 2.1. Load best model csv ----
best_mods <- read.csv(paste(cv.dir, "best_models_new_cvs.csv", sep ='/')) 
head(best_mods)

names(best_mods)
names(best_mods) <- c("bm_id"  , "modname")
names(best_mods)



### 2.2. Set submodel names ----
no_submodels <- nrow(best_mods)
no_submodels


#submodels <- c('3.0', '3.1', '3.2', '3.3')
submodels <- paste(1:no_submodels)

submodel_names <- paste("bm", submodels, sep = '')


### 2.3. Set number of folds ----
no_folds <- 15


### 2.4. Set data frame to save summary stats ----
overall.summary <- data.frame()


### Standard error function ----
std.error <- function(x) sd(x)/sqrt(length(x))



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
    summarise(mean_stat = mean(values),
              se_stat = std.error(values)) %>%
    mutate(submodel = submodel_names[i]) %>%
    mutate_at(vars(submodel), list(as.factor)) %>%
    glimpse()
  
  
  
  # 6.1. add to overall summary --
  
  overall.summary <- rbind(overall.summary, summary.df2)
  overall.summary
  
  
}

nobeep()

# 4. SAVE SUMMARY ----
head(overall.summary)

overall.summary <- na.omit(overall.summary)

write.csv(overall.summary, paste(cv.dir, "best_models_CV_summary_stats.csv", sep ='/'))



# 5. PLOT SUMMARY STATS FOR BEST MODELS ----

p <- ggplot(overall.summary %>% dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                                                summary_stat != 'mean_residuals'), 
            aes(x = submodel, y = mean_stat, color = submodel)) +
  geom_errorbar(aes(ymin=mean_stat-se_stat, ymax=mean_stat+se_stat)) +
  geom_point(size = 3) +
  facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, h = 1))
p


# save plot --

ggsave(plot = p, filename = "best_models_CV_summary_stats.png", device = "png", path = cv.dir)


# 6. Compute best statistics ----

sum.sum <- overall.summary %>%
  dplyr::filter(summary_stat != 'mean_R.sqr_test' &
                  summary_stat != 'mean_residuals') %>%
  pivot_wider(id_cols = submodel, names_from = summary_stat, values_from = mean_stat) %>%
  mutate(sum_stat = (dev.exp + mean_R.Sqr)*10 - (mean_residuals_test + rmse)) %>%
  mutate(sum_max = (dev.exp + mean_R.Sqr)) %>%
  mutate(sum_min = (mean_residuals_test + rmse)) %>%
  glimpse()

p1 <- ggplot(sum.sum %>%
         dplyr::select(sum_stat, sum_max, sum_min, submodel),
       aes(x = submodel, y = sum_stat, color = sum_stat)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

p1


ggsave(plot = p1, filename = "best_models_CV_summary_stats_sum.png", device = "png", path = cv.dir)


##

p2 <- ggplot(sum.sum %>%
               dplyr::select(sum_stat, sum_max, sum_min, submodel),
             aes(x = submodel, y = sum_max, color = sum_stat)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

p2


ggsave(plot = p2, filename = "best_models_CV_summary_stats_max.png", device = "png", path = cv.dir)

##

p3 <- ggplot(sum.sum %>%
               dplyr::select(sum_stat, sum_max, sum_min, submodel),
             aes(x = submodel, y = sum_min, color = sum_stat)) +
  geom_point(size = 5) +
  theme_bw() +
  #facet_wrap(~ summary_stat, scales = "free") +
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 45, h = 1))

p3


ggsave(plot = p3, filename = "best_models_CV_summary_stats_min.png", device = "png", path = cv.dir)


##

best_mods$formula[1]
best_mods$formula[2]

###

###

###

##

##


# BM 1 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[1]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm1 <- gam(bm_form, family = tw(), data = dat2, method = "GCV.Cp")
bm1 <- gam(bm_form, family = tw(), data = dat2, method = "REML")



# check model ----

bm1$aic
summary(bm1)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()


##

##


# BM 2 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[2]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm2 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm2$aic
summary(bm2)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()


##

##


# BM 3 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[3]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm3 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm3$aic
summary(bm3)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()


##

##


# BM 4 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[4]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm4 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm4$aic
summary(bm4)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()

##

##


# BM 5 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[5]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm5 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm5$aic
summary(bm5)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()


##

##


# BM 6 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[6]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm6 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm6$aic
summary(bm6)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()


##

##


# BM 7 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[7]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm7 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm7$aic
summary(bm7)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()

##

##


# BM 8 ----
# 3.1. Get submodel name --
bm_name <- best_mods$modname[8]

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
random_vars <- "s(site_campus_unique_ID, zone, bs = 're') + s(survey_year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
bm8 <- gam(bm_form, family = tw(), data = dat2, method = "REML")


# check model ----

bm8$aic
summary(bm8)
gam.check(bm)

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()



