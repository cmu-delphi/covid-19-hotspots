
library(ggplot2)
library(covidcast)
library(tidyverse)
library(tidyr) ## for pivot_wider() ?
library(devtools)
library(glmnet)
source("helpers.r")


## Location of your covidcast R package.
## load_all("/home/shyun/repos/covidcast/R-packages/covidcast")
load_all("~/Desktop/CMU/Projects/Delphi-Covid-19/delphi_repos/covidcast/R-packages/covidcast")
lags = 28

geo_type = "state"
slope = TRUE
response = "confirmed_7dav_incidence_prop"
fn_response = response_diff_avg_1week
fn_response_name = "response_diff_avg_1week"

# n_ahead = 28
# threshold = .40

## re-run the following every time geo_type or response changes!
if(FALSE){
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
  end_day = as.Date("2020-08-10")
  validation_days = seq(end_day-30, end_day, by = "days")
  signals = data.frame(data_sources = data_sources, signals = signals)
  suppressMessages({
    mat = covidcast_signals(signals,
                            start_day = start_day, end_day = end_day, geo_type = geo_type)
  })
  mat <- mat %>% select(geo_value, time_value, signal, data_source, value) 
}


for(n_ahead in c(28,21,14)){
  for(threshold in c(.25,.40)){
    to_file <- paste("\n\n\nTime: ", Sys.time(), "\nSpecifications: ", geo_type, ", lags = ", lags, " n_ahead = ", n_ahead, ", slope = ", slope, ", \nresponse = ", response, ", response function = ", fn_response_name, " , threshold = ", threshold, sep = "")
    cat(to_file)
    write(to_file, file = "counts.txt", append = TRUE)
    
    
    t0 <- Sys.time()
    df_model <- ready_to_model(mat, lags, n_ahead, response, slope, fn_response, threshold)
    Sys.time()-t0
    ## add census features (currently only population)
    # df_model <- add_geoinfo(df_model, geo_type)
    ## divide data into train, test, and validation sets
    df_traintest <- df_model %>% filter(!(time_value %in% validation_days))
    df_validation <- df_model %>% filter(time_value %in% validation_days)
    #splitted <- sample_split_date(df_traintest, pct_test=0.3)
    splitted <- list(df_train = df_traintest, df_test = df_validation)
    
    to_file <- paste("\n\tTraining set: ",nrow(splitted$df_train), " observations. 1's:", sum(splitted$df_train$resp), ", 0's:", sum(1-splitted$df_train$resp), "\n\tTest set: ",nrow(splitted$df_test), " observations. 1's:", sum(splitted$df_test$resp), ", 0's:", sum(1-splitted$df_test$resp),  sep = "")
    cat(to_file)
    write(to_file, file = "counts.txt", append = TRUE)
    
    
    ######################################
    ## model with lagged responses only ##
    ###################################### 
    predictions_onlylaggedresponse <- fit_predict_models(splitted$df_train %>% select(geo_value, time_value, resp, contains(response)), 
                                                         splitted$df_test %>% select(geo_value, time_value, resp, contains(response)),
                                                         lags = lags, n_ahead = n_ahead)
    a = plot_adapted_roc(predictions_onlylaggedresponse, geo_type = geo_type)
    a
    
    ####################################################
    ## model with lagged responses + facebook signals ##
    ####################################################
    predictions_laggedandfacebook <- fit_predict_models(splitted$df_train, splitted$df_test, lags = lags, n_ahead = n_ahead)
    b = plot_adapted_roc(predictions_laggedandfacebook, add=TRUE, df_plot_existing=a, geo_type = geo_type)
    b
    # ggsave(plot = b, filename = paste("figures/", toupper(geo_type), "precrecall_lag", lags,"_nahead", n_ahead, ".png", sep = ""), width = 12, height = 8, dpi = 200)
    ggsave(plot = b, filename = paste("figures/", fn_response_name,"/", geo_type, "_resp", threshold*100, "_lag", lags,"_nahead", n_ahead, "_slope", slope, ".png", sep = ""), width = 12, height = 8, dpi = 200) 
  
    a = plot_roc(predictions_onlylaggedresponse, geo_type = geo_type, popweighted = FALSE)
    a
    b = plot_roc(predictions_laggedandfacebook, add=TRUE, df_plot_existing=a, geo_type = geo_type, popweighted = FALSE)
    b
    ggsave(plot = b, filename = paste("figures/", fn_response_name,"/", geo_type, "_resp", threshold*100, "_lag", lags,"_nahead", n_ahead, "_slope", slope, "ROC.png", sep = ""), width = 12, height = 8, dpi = 200) 
    
  }
}

