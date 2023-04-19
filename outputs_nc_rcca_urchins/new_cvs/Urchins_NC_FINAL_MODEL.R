# 

# Script by Anita Giraldo - 2 May 2022
# last modified by Anita Giraldo - 12 Abril 2023


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
o.dir <- here('outputs_nc_rcca_urchins')
#k.dir <- paste(o.dir, "gam_urchins4", sep ='/')
cv.dir <- paste(o.dir, "new_cvs", sep ='/')
#k.dir <- paste(o.dir, "gam_urchins5", sep ='/') # kelp model results
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
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
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
    wh_max, wh_mean, mean_waveyear, wh_95prc,
    # Orb vel
    UBR_Mean, UBR_Max,
    # NPP
    Mean_Monthly_NPP, Max_Monthly_NPP_Upwelling, Mean_Monthly_NPP_Upwelling, Min_Monthly_NPP) %>%
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
  dplyr::select(-c(UBR_Mean,
                   UBR_Max)) %>%
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
names(dat2)


# 2. LOAD BEST MODELS ####


### 2.1. Load best model csv ----
best_mods <- read.csv(paste(cv.dir, "best_model.csv", sep ='/')) 
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



# 3.1. Get submodel name --
bm_name <- best_mods$modname[1]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)

# 3.2. Set model formula --

# set dependent variable --
dep <- 'log_den_STRPURAD'

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



# Fit model bm 3 using all data ----
#bm1 <- gam(bm_form, family = tw(), data = dat2, method = "GCV.Cp")
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

# Plot response curves ----


library(visreg)
library(tidymv)
library(gratia)
library(patchwork)

length(preds)

theme_set(theme_bw())

gratia::draw(bm, residuals = TRUE,
             select = c("s(wh_mean)", "s(mean_slope)"),
             smooth_col = "black",
             ci_col = "steelblue3",
             ci_alpha = 0.2,
             resid_col = 'gray70',
             resid_size = 0.5)

# evaluate the smooths
sm <- smooth_estimates(bm1) %>%
  add_confint()
sm

# add partial residuals to data
eg1 <- dat2 %>%
  add_partial_residuals(bm1)

names(eg1)

theme_set(theme_bw())


# p1 ----
vars <- bm_vars[[1]]
var <- vars[1]
svar <- paste('s(', var, ')', sep ='')
var


p1 <- sm %>%
  filter(smooth == "s(Max_Monthly_Anomaly_Upwelling_Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Max_Monthly_Anomaly_Upwelling_Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Mean_Monthly_Upwelling_Temp),
  #             alpha = 0.2) +
  geom_point(aes(x = Max_Monthly_Anomaly_Upwelling_Temp, y = `s(Max_Monthly_Anomaly_Upwelling_Temp)`),
             data = eg1, cex = 1.9, colour = "steelblue3", shape = 16
             , alpha = 0.3
  ) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Max_Monthly_Anomaly_Upwelling_Temp),
              alpha = 0.5) +
  geom_line(aes(x = Max_Monthly_Anomaly_Upwelling_Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", x = "Max temperature anomaly - Upwelling"
       #, title = "s(log_UBR_Max)"
  ) +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
p1


# p2 ----
vars <- bm_vars[[1]]
var <- vars[2]
svar <- paste('s(', var, ')', sep ='')
var


p2 <- sm %>%
  filter(smooth == "s(Days_10N)") %>%
  ggplot() +
  geom_rug(aes(x = Days_10N),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Max_Monthly_Nitrate),
  #             alpha = 0.2) +
  geom_point(aes(x = Days_10N, y = `s(Days_10N)`),
             data = eg1, cex = 1.9, colour = "steelblue3", shape = 16
             , alpha = 0.3
  ) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Days_10N),
              alpha = 0.5) +
  geom_line(aes(x = Days_10N, y = est), lwd = 1.2) +
  labs(y = "Partial effect", x = "Days_10N"
       #, title = "s(log_UBR_Max)"
  ) +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(0),
        axis.text = element_text(size = 10))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
p2


# p3 ----
vars <- bm_vars[[1]]
var <- vars[3]
svar <- paste('s(', var, ')', sep ='')
var


p3 <- sm %>%
  filter(smooth == "s(log_UBR_Max)") %>%
  ggplot() +
  geom_rug(aes(x = log_UBR_Max),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = depth_mean),
  #             alpha = 0.2) +
  geom_point(aes(x = log_UBR_Max, y = `s(log_UBR_Max)`),
             data = eg1, cex = 1.9, colour = "steelblue3", shape = 16
             , alpha = 0.3
  ) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = log_UBR_Max),
              alpha = 0.5) +
  geom_line(aes(x = log_UBR_Max, y = est), lwd = 1.2) +
  labs(y = "Partial effect", x ="Max orbital velocity (log)"
       #, title = "s(log_UBR_Max)"
  ) +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
p3



# p4 ----
vars <- bm_vars[[1]]
var <- vars[4]
svar <- paste('s(', var, ')', sep ='')
var


p4 <- sm %>%
  filter(smooth == "s(mean_depth)") %>%
  ggplot() +
  geom_rug(aes(x = mean_depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = wh_mean),
  #             alpha = 0.2) +
  geom_point(aes(x = mean_depth, y = `s(mean_depth)`),
             data = eg1, cex = 1.9, colour = "steelblue3", shape = 16
             , alpha = 0.3
  ) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = mean_depth),
              alpha = 0.5) +
  geom_line(aes(x = mean_depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", x = "Mean depth",
       #, title = "s(log_UBR_Max)"
  ) +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
p4




# p5 ----
vars <- bm_vars[[1]]
var <- vars[5]
svar <- paste('s(', var, ')', sep ='')
var


p5 <- sm %>%
  filter(smooth == "s(log_Mean_Monthly_NPP_Upwelling)") %>%
  ggplot() +
  geom_rug(aes(x = log_Mean_Monthly_NPP_Upwelling),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = log_UBR_Max),
  #             alpha = 0.2) +
  geom_point(aes(x = log_Mean_Monthly_NPP_Upwelling, y = `s(log_Mean_Monthly_NPP_Upwelling)`),
             data = eg1, cex = 1.9, colour = "steelblue3", shape = 16
             , alpha = 0.3
  ) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = log_Mean_Monthly_NPP_Upwelling),
              alpha = 0.5) +
  geom_line(aes(x = log_Mean_Monthly_NPP_Upwelling, y = est), lwd = 1.2) +
  labs(y = "Partial effect", x = "Mean NPP - Upwelling (log)", 
       #, title = "s(log_UBR_Max)"
  ) +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
p5


# p6 ----
vars <- bm_vars[[1]]
var <- vars[6]
svar <- paste('s(', var, ')', sep ='')
var


p6 <- sm %>%
  filter(smooth == "s(log_Min_Monthly_NPP)") %>%
  ggplot() +
  geom_rug(aes(x = log_Min_Monthly_NPP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Mean_Monthly_NPP),
  #             alpha = 0.2) +
  geom_point(aes(x = log_Min_Monthly_NPP, y = `s(log_Min_Monthly_NPP)`),
             data = eg1, cex = 1.9, colour = "steelblue3", shape = 16
             , alpha = 0.3
  ) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = log_Min_Monthly_NPP),
              alpha = 0.5) +
  geom_line(aes(x = log_Min_Monthly_NPP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", x = "Min NPP (log)"
       #, title = "s(log_UBR_Max)"
  ) +
  theme(axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
p6





# plot all ----
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)

library(gridExtra)
library(ggpubr)

# ar <- grid.arrange(p1,p2,p3,p4, ncol = 2)
# 
# ggsave(
#   filename = "response_curves1.tiff", 
#   path = cv.dir, 
#   plot = marrangeGrob(p1,p2,p3,p4, nrow=2, ncol=2), 
#   width = 15, height = 9
# )

one.page <- ggarrange(p1,p2,p3,p4,p5, p6, 
                      labels = c('a', 'b', 'c', 'd', 'e', 'f'),
                      nrow=2, ncol=3,
                      #width = 1200, height = 1200,
                      res = 600) 
one.page[[1]]

ggexport(one.page[[1]], filename=paste(cv.dir, "response_curves1.tiff", sep ='/'),
         width = 2600, height = 1600,
         res = 300)

ggexport(one.page[[1]], filename=paste(cv.dir, "response_curves1.png", sep ='/'),
         width = 2600, height = 1600,
         res = 300)
