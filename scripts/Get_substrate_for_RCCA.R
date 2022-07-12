###

## To get susbtrate variables for RCCA data ----

## aggragate with threshold ----
mean.thresh <- function(x) {ifelse(length(which(is.na(x) == T)) >= 70, NA, mean(x, na.rm = T))} 

sd.thresh <- function(x) {ifelse(length(which(is.na(x) == T)) >= 70, NA, sd(x, na.rm = T))}
