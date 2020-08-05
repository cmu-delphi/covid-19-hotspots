library(dplyr)
library(lubridate)
library(plotROC)

##' creates a dataframe with lagged featuresfor one geo_value
##' 
lagged_features_onegeo <- function(df, lags, name = "feature"){
  signal <- df$val
  timestamp <- df$time
  len <- nrow(df)
  if(len<=lags){
    return(data.frame())
  }
  
  ## make sure timestamp is in date formate
  if(!is.Date(timestamp)) timestamp <- as.Date(timestamp)
  
  ## TODO if there's NA in signal, treat it BEFORE creating lagged matrix. needs helper function. 
  ## it's reasonable to interpolate the TS as long as there are not many sequential missing obs
  
  out <- data.frame(time = timestamp[(lags+1):len])
  for(i in 0:lags){
    out <- suppressMessages(bind_cols(out, signal[(lags+1-i):(len-i)]))
  }
  names(out) = c("time", paste(name, "_lag", 0:lags, sep = ""))
  
  return(out)
}

##' constructs response variable for one geo_value using the provided function
##' 
response_onegeo <- function(df, n_ahead, fn_response = response_diff, ...){
  signal <- df$val
  timestamp <- df$time
  len <- nrow(df)
  stopifnot(n_ahead <= len)
  
  ## make sure timestamp is in date format
  if(!is.Date(timestamp)) timestamp <- as.Date(timestamp)
  
  ## TODO if there's NA in signal, treat it BEFORE creating response. needs helper function. 
  ## it's reasonable to interpolate the TS as long as there are not many sequential missing obs
  
  out <- data.frame(time = timestamp[1:(len-n_ahead+1)])
  out$resp <- NA
  for(i in 1:(len-n_ahead+1)){
    out$resp[i] <- fn_response(signal[i:(i+n_ahead)], ...)
  }
  
  return(out)
}

##' considers increase if there is a 25% increase
##' 
response_diff <- function(x, threshold = .25){
  len = length(x)
  ifelse((x[len]-x[1])/x[1]>1.25, 1, 0)
}

##' considers increase if there is any 25% increase in the provided interval
##' 
response_maxdiff <- function(x, threshold = .25){
  ifelse((max(x)-min(x))(which.max(x)>=which.min(x))/min(x)>1.25, 1, 0)
}

## TODO
#get_countyinfo <- function()

##' outputs data ready for modeling, with all lagged features and binary response
##'    uses output from API call
##' 
ready_to_model <- function(mat, lags, n_ahead, response = "confirmed_7dav_incidence_num"){
  features <- mat %>% plyr::ddply(c("signal", "geo"), lagged_features_onegeo, lags = lags) %>% na.omit()
  response <- mat %>% filter(signal == response) %>% plyr::ddply(c("signal", "geo"), response_onegeo, n_ahead = n_ahead) %>% na.omit()
  names_to_pivot <- colnames(features %>% select(-geo, -time, -signal))
  features <- pivot_wider(features, id_cols = c("geo", "time"), names_from = c("signal"), values_from = all_of(names_to_pivot))
  mat_to_model <- inner_join(features, response %>% select(-signal), by = c("geo", "time"))
  ## TODO add census features
  return(mat_to_model)
}

##' performs sample splitting based on date: test set will always be the most recent data
##' because it depends on date, it is NOT RANDOM!!!
##' 
sample_split_date <- function(df_tomodel, pct_test=0.3){
  df_tomodel <- df_tomodel %>% arrange(desc(time)) %>% na.omit() # TODO: treat NA's properly. Maybe Dmitry's smoother?
  start_test_date <- df_tomodel[round(pct_test*nrow(df_tomodel)),"time"]
  df_test <- df_tomodel %>% filter(time >= start_test_date$time[1])
  df_train <- df_tomodel %>% filter(time < start_test_date$time[1])
  return(list(df_test = df_test, df_train = df_train))
}

##' fit models and produce test set predictions
##' currently: lasso, ridge
##' 
fit_predict_models <- function(mat, lags, n_ahead){
  cat("Creating lagged variables and binary response...\n")
  df_model <- ready_to_model(mat, lags, n_ahead)
  cat("Done! \nNow fitting models:\n")
  splitted <- sample_split_date(df_model, pct_test=0.3)
  df_train <- splitted$df_train
  df_test <- splitted$df_test
  
  predictions <- df_test %>% select(geo, time, resp)
  
  cat("\tFitting LASSO...")
  fit_lasso <- cv.glmnet(x = as.matrix(df_train %>% select(-geo, -time, -resp)), y = df_train$resp, family = "binomial", alpha = 1)
  fit_lasso <- glmnet(x = as.matrix(df_train %>% select(-geo, -time, -resp)), 
                      y = df_train$resp, family = "binomial", lambda = fit_lasso$lambda.1se, alpha = 1)
  predictions[[paste("lasso_lags", lags, "_nahead", n_ahead, sep = "")]] = predict(fit_lasso, newx = as.matrix(df_test %>% select(-geo, -time, -resp)), type = "response")[,1]
  
  cat(" Done!\n\tFitting Ridge...")
  ### more models here!!! SVM, xgboost... add predictions as a col to predictions
  ## eg ridge:
  fit_ridge <- cv.glmnet(x = as.matrix(df_train %>% select(-geo, -time, -resp)), y = df_train$resp, family = "binomial", alpha = 0)
  fit_ridge <- glmnet(x = as.matrix(df_train %>% select(-geo, -time, -resp)), 
                      y = df_train$resp, family = "binomial", lambda = fit_ridge$lambda.1se, alpha = 0)
  predictions[[paste("ridge_lags", lags, "_nahead", n_ahead, sep = "")]] = predict(fit_ridge, newx = as.matrix(df_test %>% select(-geo, -time, -resp)), type = "response")[,1]
  cat(" Done!\n")
  return(predictions)
}



# signal <- c(1,3,2,5,5,2,6,3,2,3,4,7,6,8,8,5)
# start_date <- as.Date("2020-05-10")
# timestamp <- seq(start_date, start_date+length(signal)-1, 1)
# lags = 2; n_ahead = 3
# mat_test <- data.frame(geo = 1,time=timestamp, val=signal, signal="testsignal")
# ready_to_model(mat_test, lags, n_ahead, "testsignal")




