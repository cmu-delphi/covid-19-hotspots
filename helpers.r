library(dplyr)
library(lubridate)

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
    out$resp[i] <- fn_response(signal[1:(i+n_ahead)], i, ...)
  }
  
  return(out)
}

##' considers increase if there is a 25% increase
##' very simple function, only looks at today's value and target value
##' we should write a better function
##' 
##' @param x vector of values that will be used to determine hotspot, from 1 until i+n_ahead
##' @param i position of the vector x that is "today"; everything from i+1:forward is not known as features
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff <- function(x, i, threshold = .25){
  len = length(x)
  ifelse((x[len]-x[i])/x[i]>1.25, 1, 0)
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
  mat_to_model <- inner_join(features, responses %>% select(-signal, -data_source), by = c("geo_value", "time_value")) %>% na.omit()
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

#' Make |foldid| argument for covariate matrix |x| and |nfold|-fold
#' cross-validation; makes nfold consecutive time blocks
#' @param x covariate matrix
#' @param nfold nfold
make_foldid <- function(x, nfold){
  ## inds = round(seq(from = 0, to = nrow(df_train), length=nfold))
  ## indlist = lapply(1:(length(inds)-1), function(ii) (inds[ii]+1):inds[ii+1])
  ## vec = rep(NA,nrow(x))
  ## for(ii in 1:length(indlist)){vec[indlist[[ii]]] = ii}
  
  vec = rep(1:nfold, each=ceiling(nrow(x)/nfold))
  vec = vec[1:nrow(x)]
  
  return(vec)
}

##' fit models and produce test set predictions
##' currently: lasso, ridge
##' 
##' @param df_model 
##' @param lags number of past values to include in the data frame; for time t, dataframe will have in one row X_t until X_{t-lag}
##' @param n_ahead number of days ahead that response will be computed
fit_predict_models <- function(df_model, lags, n_ahead, response = "confirmed_7dav_incidence_num"){
  cat("Fitting models:\n")
  splitted <- sample_split_date(df_model, pct_test=0.3)
  df_train <- splitted$df_train
  df_test <- splitted$df_test
  cat(paste("Training set: ",nrow(df_train), " observations. \nTest set: ",nrow(df_train), " observations. \n",  sep = ""))

  predictions <- df_test %>% select(geo_value, time_value, resp)
  
  cat("\tFitting LASSO...")
  preds <- fit_logistic_regression(df_train, df_test, nfold = 10, alpha = 0)
  predictions[[paste("lasso_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  cat(" Done!\n")

  cat("\tFitting Ridge...")
  ### more models here!!! SVM, xgboost... add predictions as a col to predictions
  ## eg ridge:
  preds <- fit_logistic_regression(df_train, df_test, nfold = 10, alpha = 1)
  predictions[[paste("ridge_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  cat(" Done!\n")
  
  cat("\tFitting SVM...")
  preds <- fit_svm(df_train, df_test)
  predictions[[paste("svm_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  cat(" Done!\n")
  
  return(predictions)
}



# signal <- c(1,3,2,5,5,2,6,3,2,3,4,7,6,8,8,5)
# start_date <- as.Date("2020-05-10")
# timestamp <- seq(start_date, start_date+length(signal)-1, 1)
# lags = 2; n_ahead = 3
# mat_test <- data.frame(geo = 1,time=timestamp, val=signal, signal="testsignal")
# ready_to_model(mat_test, lags, n_ahead, "testsignal")


##' Performs CV-ed logistic lasso prediction, given training & test matrices.
##' Outputs a vector of values the same as.
##'
##' @param df_train Training matrix. Must contain columns "geo_value",
##'   "time_value", "resp", and some other columns that will be used as
##'   covariates.
##' @param df_test Test matrix. Same format as df_train.
##' @param nfold 10.
##' @param alpha 0 for lasso, or 1 for ridge regression. Used by \code{glmnet()}.
##'
##' @return Numeric vector the same length as \code{nrow(df_test)}.
fit_logistic_regression <- function(df_train, df_test, nfold = 10, alpha = 0){
  
  ## Input checks (should be common for all fit_OOOO() functions
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_train)))
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_test)))
  
  ## Input check
  stopifnot(alpha %in% c(0,1)) ## Only allow ridge or lasso for now.
  
  ## Make contiguous time blocks for CV
  foldid <- make_foldid(df_train, nfold)
  
  ## Main part of the lasso fitting and predicting
  fit_lasso <- cv.glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)),
                         y = df_train$resp,
                         family = "binomial",
                         alpha = 1,
                         foldid = foldid,
                         nfold=nfold)
  fit_lasso <- glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)),
                      y = df_train$resp, family = "binomial", lambda = fit_lasso$lambda.1se, alpha = 1)
  preds = predict(fit_lasso, newx = as.matrix(df_test %>% select(-geo_value, -time_value, -resp)), type = "response")[,1]
  
  ## Out checks (should be common for all fit_OOOO() functions)
  stopifnot(length(preds) == nrow(df_test))
  
  preds
}

##' Performs SVM prediction, given training & test matrices.
##' Outputs a vector of values the same as.
##'
##' @param df_train Training matrix. Must contain columns "geo_value",
##'   "time_value", "resp", and some other columns that will be used as
##'   covariates.
##' @param df_test Test matrix. Same format as df_train.
##' @param ... Additional functions to \code{svm()} of the \code{e1071} R
##'   package.
##'
##' @return Numeric vector the same length as \code{nrow(df_test)}.
fit_svm <- function(df_train, df_test, ...){
  
  ## Input checks (should be common for all fit_OOOO() functions
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_train)))
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_test)))
  
  ## Fit SVM and make predictions
  train_mat <- df_train %>% select(-geo_value, -time_value)
  test_mat <- df_test %>% select(-geo_value, -time_value)
  model <- e1071::svm(resp ~ ., data = train_mat, probability = TRUE, ...)
  preds <- predict(model, test_mat)
  
  ## Next: replace with faster SVM: https://cran.r-project.org/src/contrib/Archive/RSofia/
  ## ## Also might be useful: multicore, faster radial SVM with CV?
  ## library(caret)
  ## library(doMC)
  ## registerDoMC()
  ## model <-  train(Species ~ ., data = iris, method="svmRadial",
  ##     trControl=trainControl(method='cv', number=10)) ## This CV is not great..
  ## confusionMatrix(model)
  
  ## Out checks (should be common for all fit_OOOO() functions)
  stopifnot(length(preds) == nrow(df_test))
  
  preds
}


##' gets population for specific geo_type
##'
##' @param geo_type county, msa, state
##'
##' @return dataframe with columns geo_value and population
get_population <- function(geo_type){
  if(geo_type == "county"){
    county_pop = county_census %>% 
      transmute (geo_value = 1000*as.numeric(STATE) + as.numeric(COUNTY),
                 population = POPESTIMATE2019)
    county_pop$geo_value <- sprintf("%05d", county_pop$geo_value)
    return(county_pop)
  }
  if(geo_type == "state"){
    state_pop <- state_census %>% 
      mutate(geo_value = as.numeric(STATE)) %>% 
      filter(STATE != 0) %>% 
      group_by(geo_value) %>% 
      summarise(population = sum(POPESTIMATE2019))
    state_crosswalk <- state.fips %>% 
      select(abb, fips) %>% distinct() %>% mutate(abb = tolower(abb))
    state_pop <- state_pop %>% 
      inner_join(state_crosswalk, by = c("geo_value" = "fips")) %>% 
      select(geo_value = abb, population)
    return(state_pop)
  }
  if(geo_type == "msa"){
    return(msa_census %>% 
             transmute(geo_value = as.character(CBSA),
                       population = POPESTIMATE2019))
  }
}

##' computes population weighted precision and population weighted proportion of predicted 1's for different cutoffs
##' considers ONE MODEL only
##'
##' @param df_one dataframe for one model with at least columns: value, pred, resp, population
##'
##' @return dataframe with colulmns cutoff, wpred1, wprecision
adapted_roc_onemodel <- function(df_one,...){
  df_one <- df_one %>% arrange(value)
  ## changing cutoffs
  metrics <- sapply(seq(0, 1, 0.01), function(i){
    df_temp <- df_one
    ## binary predictions for the cutoff
    df_temp$pred <- ifelse(df_temp$value<=i, 0, 1)
    
    ## computing weighted precision. if numerator and denominator are 0, assign precision = 1
    precision = (df_temp %>% filter(pred == 1, resp == 1) %>% summarise(n()) %>% unlist)/(df_temp %>% filter(pred == 1) %>% summarise(n()) %>% unlist)
    if((df_temp %>% filter(pred == 1) %>% nrow) == 0) precision = 1
    wprecision = precision*(df_temp %>% filter(pred == 1, resp == 1) %>% summarise(sum(population)) %>% unlist)/(df_temp %>% filter(pred == 1) %>% summarise(sum(population)) %>% unlist)
    if((df_temp %>% filter(pred == 1) %>% summarise(sum(population)) %>% unlist) == 0) wprecision = 1
    
    ## computing weighted % of predicted 1's
    pred1 = (df_temp %>% filter(pred == 1) %>% summarise(n()) %>% unlist)/nrow(df_temp)
    wpred1 = pred1 * (df_temp %>% filter(pred == 1) %>% summarise(sum(population)) %>% unlist)/sum(df_temp$population)
    return(c(i, wpred1, wprecision))
  })
  ## transforming matrix into dataframe and naming it appropriately
  metrics <- as.data.frame(t(metrics))
  names(metrics) <- c("cutoff","wpred1", "wprecision")
  return(metrics)
}

##' computes metrics for all models and produces roc curves (our adapted roc with different metrics)
##'
##' @param predictions dataframe with cols: geo_value, time_value, resp, 
##'                    and one column per model with predicted values whose col name is the models name
##' @param geo_type county, msa, or state. will be used to get population data
##'
##' @return ggplot of model comparison curve
plot_adapted_roc <- function(predictions, geo_type = "county"){
  df_temp <- inner_join(predictions, get_population(geo_type), by = "geo_value")
  df_temp <- reshape2::melt(df_temp, id.vars = c("geo_value", "time_value", "resp", "population"))
  df_plot <- df_temp %>% rename(model = variable) %>% group_by(model) %>% group_modify(adapted_roc_onemodel)
 
  ggplot(df_plot, aes(x = wpred1, y = wprecision, color = model)) +  
    geom_line() +
    theme_bw(base_size = 18) + 
    ylab("population weighted precision") +
    xlab("population weighted % predicted hotspots") + 
    theme(legend.position = "bottom")
}
