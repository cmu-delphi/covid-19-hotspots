library(tidyverse)
library(tidyr) ## for pivot_wider() ?
library(lubridate)
library(devtools)
library(glmnet)
library(ggplot2)
library(covidcast)
load_all("C:/Users/dzhao/Documents/Carnegie Mellon 2018-2023/RESEARCH - Delphi/covidcast/R-packages/covidcast")
setwd("C:/Users/dzhao/Documents/Carnegie Mellon 2018-2023/RESEARCH - Delphi/covid-19-hotspots")
source("helpers.r")
source("helpers_20200916.r")



#########################################################################################
###  Load data from API
#########################################################################################
geo_type = "county"
start_day = as.Date("2020-05-01")
end_day = as.Date("2020-09-15")

### Note: FB survey also has raw and weighted CLI / ILI, and community minus household
### Note: Doctor visits also has day-of-week effect removed
### Note: Safegraph also has full time work (>6 hours), part time work (3-6 hours), median away from home time
data_sources = c("indicator-combination", 
                 "fb-survey",
                 "doctor-visits",
                 "ght",
                 "quidel",
                 "safegraph"
                 )
signals = c("confirmed_7dav_incidence_prop", 
            "smoothed_hh_cmnty_cli",
            "smoothed_cli",
            "smoothed_search",
            "covid_ag_smoothed_pct_positive",
            "completely_home_prop"
            )
signals = data.frame(data_sources = data_sources, signals = signals)
suppressMessages({
  mat = covidcast_signals(signals,
                          start_day = start_day, end_day = end_day, geo_type = geo_type)
})
mat <- mat %>% select(geo_value, time_value, signal, data_source, value)



#########################################################################################
###  restrict for now to counties with at least 20 cases on some date  (for now, May 14)
#########################################################################################
cutoff = 20
suppressMessages({
  cases = covidcast_signals(data.frame(data_sources = c("indicator-combination"), signals = c("confirmed_incidence_num")),
                            start_day = as.Date("2020-05-14"), end_day = as.Date("2020-05-14"), geo_type = geo_type)
})
cases <- cases %>% select(geo_value, time_value, value) %>% filter(value >= cutoff)
keep_geos <- pull(cases, geo_value)
mat <- mat %>% filter(geo_value %in% keep_geos)



#########################################################################################
###  get all features
#########################################################################################
###  all features summed over a week
mat <- mat %>% group_by(geo_value, week=epiweek(time_value), signal, data_source) %>% summarize(value = sum(value))

### looking 1-4 weeks back, including slopes
lags = 4
slope_lags = 4
slope = TRUE
features <- mat %>% plyr::ddply(c("signal", "data_source", "geo_value"),
                                lagged_weekly_features_onegeo,
                                lags=lags, slope_lags=slope_lags, slope=slope) %>% na.omit()



#########################################################################################
###  get cases ahead
#########################################################################################
response = "confirmed_7dav_incidence_prop"
n_ahead = 4

### looking 1-4 weeks ahead
responses <- mat %>% filter(signal == response) %>%
  plyr::ddply(c("signal", "data_source", "geo_value"), cases_ahead_weekly_onegeo,
              n_ahead = n_ahead) # %>% na.omit()

## transform the dataframe in a wide format, with one row per geo_value and date
names_to_pivot <- colnames(features %>% select(-geo_value, -week, -signal, -data_source))
features <- pivot_wider(features, id_cols = c("geo_value", "week"), names_from = c("signal", "data_source"),
                        values_from = all_of(names_to_pivot)) %>% ungroup

## join features and response
df_model <- inner_join(features, responses %>% select(-signal, -data_source),
                       by = c("geo_value", "week")) # %>% na.omit()
df_model %>% names() %>% print()

#save(df_model, file=file.path("df_model_20.Rdata"))
#save(df_model200, file=file.path("df_model_200.Rdata"))
load("df_model_20.Rdata")



#########################################################################################
###  define response variable
#########################################################################################

######## current response function used by Justin and Natalia
response_name = "feature_lag0_confirmed_7dav_incidence_prop_indicator-combination"
fn_response = response_diff_avg_1week_min20
fn_response_name = "response_diff_avg_1week_min20"
n_ahead = 4
threshold = 0.25

## add response variable to matrix
df_model1 <- fn_response(df_model, n_ahead, threshold, response_name, fn_response_name)


