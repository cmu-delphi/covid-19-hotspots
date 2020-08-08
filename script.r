library(ggplot2)
library(covidcast)
library(tidyverse)
library(tidyr) ## for pivot_wider() ?
library(devtools)
library(glmnet)
#load_all("/home/shyun/repos/covidcast/R-packages/covidcast") ## Sorry for the
## absolute path,
## you need to load
## the package
## directly
load_all("~/Desktop/CMU/Projects/Delphi-Covid-19/delphi_repos/covidcast/R-packages/covidcast") ## Sorry for the
source("helpers.r")

lags = 21
n_ahead = 21
response = "confirmed_7dav_incidence_prop"
data_sources = c("indicator-combination", 
                 "fb-survey", 
                 "fb-survey", 
                 "fb-survey", 
                 "fb-survey")
signals = c("confirmed_7dav_incidence_prop", 
            "smoothed_cli", 
            "smoothed_nohh_cmnty_cli", 
            "smoothed_wcli", 
            "smoothed_hh_cmnty_cli")
start_day = as.Date("2020-05-01")
end_day = as.Date("2020-08-04")
validation_days = seq(end_day-30, end_day, by = "days")
signals = data.frame(data_sources = data_sources, signals = signals)
#suppressMessages({
mat = covidcast_signals(signals,
                        start_day = start_day, end_day = end_day)
#})
mat <- mat %>% select(geo_value, time_value, signal, data_source, value) 
mat <- add_NAval_missing_dates(mat) ## SH: This is going away? NLO: I think it's safer to keep it here
t0 <- Sys.time()
df_model <- ready_to_model(mat, lags, n_ahead, response)
Sys.time()-t0
df_traintest <- df_model %>% filter(!(time_value %in% validation_days))
df_validation <- df_model %>% filter(time_value %in% validation_days)

## model with lagged responses only
##########
predictions_onlylaggedresponse <- fit_predict_models(df_traintest %>% select(geo_value, time_value, resp, contains(response)), 
                                                     lags = lags, n_ahead = n_ahead)
plot_adapted_roc(predictions_onlylaggedresponse)


## model with lagged responses + acebook signals
##########
predictions_laggedandfacebook <- fit_predict_models(df_traintest, lags = lags, n_ahead = n_ahead)
plot_adapted_roc(predictions_laggedandfacebook)
