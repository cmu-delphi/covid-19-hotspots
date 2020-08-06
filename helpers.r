library(dplyr)
library(lubridate)
library(plotROC)

##' adds NA to value if there is no signal for a particular day
##' 
##' @param df data frame with columns time_value, value, data_source, signal, geo_value
##' @return inputted dataframe with added rows for the dates that are missing, with NA for value
add_NAval_missing_dates <- function(df){
  df$time_value <- as.Date(df$time_value)
  ## for each source-signal-geo_value, add missing dates with NA's in the value column
  df %>% group_by(data_source, signal, geo_value) %>% group_modify(function(df_one,...){
    
    ## if there is no datapoints for this combination of source-signal-geo_value, return empty dataset
    if(nrow(df_one)==0) return(data.frame())
    
    ## selects dates that do not have available signal value
    seq_days <- seq(min(df_one$time_value), max(df_one$time_value), "days")
    temp_dates <- seq_days[!(seq_days %in% unique(df_one$time_value))]
    
    ## if all dates are available, return input dataframe
    if(length(temp_dates) == 0) return(df_one)
  
    ## join inputed data frame with a dataframe only with time column, where the times are the ones that are missing from the inputted dataset
    return(full_join(df_one, data.frame(time_value = as.Date(temp_dates)), c("time_value")) %>% arrange(time_value))
  })
}

##' creates a dataframe with lagged featuresfor one geo_value
##' assumes that time column has points for all dates!!!
##' 
##' @param df dataframe with ONE geo_value and ONE feature (which will be lagged), columns val, time
##' @param lags number of past values to include in the data frame; for time t, dataframe will have in one row X_t until X_{t-lag}
##' @param name variable name that will be used for the lagged features
##' @return dataframe with lagged features for one geo_value
lagged_features_onegeo <- function(df, lags, name = "feature"){
  signal <- df$value
  timestamp <- df$time_value
  
  ## if you want more lags than available points, returns empty dataframe
  len <- nrow(df)
  if(len<=lags){
    return(data.frame())
  }
  
  ## make sure timestamp is in date formate
  if(!is.Date(timestamp)) timestamp <- as.Date(timestamp)
  
  ## TODO if there's NA in signal, treat it BEFORE creating lagged matrix. needs helper function. 
  ## it's reasonable to interpolate the TS as long as there are not many sequential missing obs
  
  out <- data.frame(time_value = timestamp[(lags+1):len])
  ## adding lagged feature from t-0, t-1, t-2, until t-lags
  for(i in 0:lags){
    out <- suppressMessages(bind_cols(out, signal[(lags+1-i):(len-i)]))
  }
  names(out) = c("time_value", paste(name, "_lag", 0:lags, sep = ""))
  
  return(out)
}

##' constructs response variable for one geo_value using the provided function
##' 
##' @param df dataframe with ONE geo_value and ONE feature (which will be used to construct the binary response), columns val, time
##' @param n_ahead number of days ahead that response will be computed
##' @param fn_response logic for computing response, based on a provided response vector whose all points will be used for this computation
##' @return inputted dataframe with addedmresp colum, which is a binary variable indicating if there is a hotspot n days ahead of the date variable
response_onegeo <- function(df, n_ahead, fn_response = response_diff, ...){
  signal <- df$value
  timestamp <- df$time_value
  
  ## we can only determine a hotspot n_ahead days if that day is available
  len <- nrow(df)
  stopifnot(n_ahead <= len)
  
  ## make sure timestamp is in date format
  if(!is.Date(timestamp)) timestamp <- as.Date(timestamp)
  
  ## TODO if there's NA in signal, treat it BEFORE creating response. needs helper function. 
  ## it's reasonable to interpolate the TS as long as there are not many sequential missing obs
  
  out <- data.frame(time_value = timestamp[1:(len-n_ahead+1)])
  out$resp <- NA
  ## for the points that have n_ahead points available, compute hotspot based on the values available 
  ## between time and time+n_ahead using the logic provided through fn_response
  for(i in 1:(len-n_ahead+1)){
    out$resp[i] <- fn_response(signal[i:(i+n_ahead)], ...)
  }
  
  return(out)
}

##' considers increase if there is a 25% increase
##' 
##' @param x vector of values that will be used to determine hotspot
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff <- function(x, threshold = .25){
  len = length(x)
  ifelse((x[len]-x[1])/x[1]>1.25, 1, 0)
}

##' considers increase if there is any 25% increase in the provided interval
##' 
##' @param x vector of values that will be used to determine hotspot
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_maxdiff <- function(x, threshold = .25){
  ifelse((max(x)-min(x))(which.max(x)>=which.min(x))/min(x)>1.25, 1, 0)
}

## TODO
#get_countyinfo <- function()

##' outputs data ready for modeling, with all lagged features and binary response
##'    uses output from API call
##' 
##' @param mat resulted from the API call with all signals we want to use for model construction
##' @param lags number of past values to include in the data frame; for time t, dataframe will have in one row X_t until X_{t-lag}
##' @param n_ahead number of days ahead that response will be computed
##' @param response name of the response variable in mat
ready_to_model <- function(mat, lags, n_ahead, response = "confirmed_7dav_incidence_num"){
  ## construct lagged features for all available signals, including lagged responses
  ## removes all NAs 
  # TODO deal with NAs
  features <- mat %>% plyr::ddply(c("signal", "data_source", "geo_value"), lagged_features_onegeo, lags = lags) %>% na.omit()
  ## construct hotspot indicator in the resp variable
  responses <- mat %>% filter(signal == response) %>% plyr::ddply(c("signal", "data_source", "geo_value"), response_onegeo, n_ahead = n_ahead) %>% na.omit()
  ## transform the dataframe in a wide format, with one row per geo_value and date
  names_to_pivot <- colnames(features %>% select(-geo_value, -time_value, -signal, -data_source))
  features <- pivot_wider(features, id_cols = c("geo_value", "time_value"), names_from = c("signal", "data_source"), 
                          values_from = all_of(names_to_pivot)) %>% ungroup
  ## join features and response
  mat_to_model <- inner_join(features, responses %>% select(-signal, -data_source), by = c("geo_value", "time_value"))
  ## TODO add census features
  return(mat_to_model)
}

##' performs sample splitting based on date: test set will always be the most recent data
##' because it depends on date, it is NOT RANDOM!!!
##' it only splits into 2 parts
##' 
##' @param df_tomodel dataset ready to model, with all lagged covariates and binary response
##' @param pct_test percentage of the points that will be in the test set
sample_split_date <- function(df_tomodel, pct_test=0.3){
  df_tomodel <- df_tomodel %>% arrange(desc(time_value)) %>% na.omit() # TODO: treat NA's properly. Maybe Dmitry's smoother?
  start_test_date <- df_tomodel[round(pct_test*nrow(df_tomodel)),"time_value"]
  df_test <- df_tomodel %>% filter(time_value >= start_test_date$time_value[1])
  df_train <- df_tomodel %>% filter(time_value < start_test_date$time_value[1])
  return(list(df_test = df_test, df_train = df_train))
}

##' fit models and produce test set predictions
##' currently: lasso, ridge
##' 
##' @param mat output from api call
##' @param lags number of past values to include in the data frame; for time t, dataframe will have in one row X_t until X_{t-lag}
##' @param n_ahead number of days ahead that response will be computed
fit_predict_models <- function(mat, lags, n_ahead, response = "confirmed_7dav_incidence_num"){
  cat("Creating lagged variables and binary response...\n")
  df_model <- ready_to_model(mat, lags, n_ahead, response)
  cat("Done! \nNow fitting models:\n")
  splitted <- sample_split_date(df_model, pct_test=0.3)
  df_train <- splitted$df_train
  df_test <- splitted$df_test
  
  predictions <- df_test %>% select(geo_value, time_value, resp)
  
  cat("\tFitting LASSO...")
  fit_lasso <- cv.glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)), y = df_train$resp, family = "binomial", alpha = 1)
  fit_lasso <- glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)), 
                      y = df_train$resp, family = "binomial", lambda = fit_lasso$lambda.1se, alpha = 1)
  predictions[[paste("lasso_lags", lags, "_nahead", n_ahead, sep = "")]] = 
    predict(fit_lasso, newx = as.matrix(df_test %>% select(-geo_value, -time_value, -resp)), type = "response")[,1]
  
  cat(" Done!\n\tFitting Ridge...")
  ### more models here!!! SVM, xgboost... add predictions as a col to predictions
  ## eg ridge:
  fit_ridge <- cv.glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)), y = df_train$resp, family = "binomial", alpha = 0)
  fit_ridge <- glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)), 
                      y = df_train$resp, family = "binomial", lambda = fit_ridge$lambda.1se, alpha = 0)
  predictions[[paste("ridge_lags", lags, "_nahead", n_ahead, sep = "")]] = predict(fit_ridge, newx = as.matrix(df_test %>% select(-geo_value, -time_value, -resp)), type = "response")[,1]
  cat(" Done!\n")
  return(predictions)
}



# signal <- c(1,3,2,5,5,2,6,3,2,3,4,7,6,8,8,5)
# start_date <- as.Date("2020-05-10")
# timestamp <- seq(start_date, start_date+length(signal)-1, 1)
# lags = 2; n_ahead = 3
# mat_test <- data.frame(geo = 1,time=timestamp, val=signal, signal="testsignal")
# ready_to_model(mat_test, lags, n_ahead, "testsignal")




