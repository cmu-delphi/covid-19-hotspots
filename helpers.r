library(dplyr)
library(lubridate)
library(xgboost)
library(pROC)


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

##' Creates a dataframe with lagged information for one geo_value.
##' If slopes = FALSE, creates only lagged features (up to \code{lags_val} which equals 5).
##' If slopes = TRUE, creates slopes of the time tendencies over t and the past 3, 6, 9... days + adds feature value at time t
##' assumes that time column has points for all dates!!!
##'
##' @param df dataframe with ONE geo_value and ONE feature (which will be
##'   lagged), columns val, time
##' @param lags Maximum number of past values from which to calculate the SLOPE
##'   feature, to include in the data frame.
##' @param name variable name that will be used for the lagged features
##' @param slopes if TRUE, returns a dataframe with slopes based on the past
##'   feature values and if FALSE, returs raw lagged features
##' @return dataframe with lagged features for one geo_value OR slopes
lagged_features_onegeo <- function(df, lags, name = "feature",slopes = FALSE){
  ## Basic checks
  stopifnot(length(lags) == 1 & lags >= 1)
  
  ## The fixed number of time lags for the data values themselves.
  df <- df %>% arrange(time_value)
  signal <- df$value
  timestamp <- df$time_value
  
  ## if you want more lags than available points, returns empty dataframe
  len <- nrow(df)
  lags_val = 5
  if(len <= max(lags, lags_val)){
    return(data.frame())
  }
  
  ## make sure timestamp is in date formate
  if(!is.Date(timestamp)) timestamp <- as.Date(timestamp)
  
  ## TODO if there's NA in signal, treat it BEFORE creating lagged matrix. needs helper function.
  ## I think it's reasonable to interpolate the TS as long as there are not many sequential missing obs
  ## low priority
  
  out <- data.frame(time_value = timestamp)
  
  ################################################################
  ## adding lagged feature from t-0, t-1, t-2, until t-lags  #####
  ################################################################
  lags_val = 5
  for(i in 0:lags_val){
    inds = 1:(len-i)
    out <- suppressMessages(bind_cols(out, c(rep(NA,i), signal[inds])))
  }
  names(out) = c("time_value", paste(name, "_lag", 0:lags_val, sep = ""))
  
  #############################################
  ## adding slopes feature every 3 points #####
  #############################################
  if(slopes){
    how_far_back = seq(min(lags, 3), lags, by = 3)
    for(nback in how_far_back){
      onevar = make_slope_var(signal, nback)
      stopifnot(length(onevar) == nrow(out))
      out[[paste(name, "_slope", nback, sep = "")]] <- onevar
    }
  }
  return(out)
}

##' Calculate slope variable, looking back \code{nback} days.
##'
##' @param signal Data.
##' @param nback How far back to go when calculating slope. For example,
##'   \code{nback=3} means calculate the slope between indices (t, t-1, .. t-3),
##'   for all possible \code{t} i.e. \code{t} in \code{(4:len)}. All other
##'   entries are NA.
##'
##' @param return A vector containing slopes calculated from today through
##'   \code{nback} days.
make_slope_var <- function(signal, nback){
  n = length(signal)
  aux <- rep(NA, n)
  endpts = (nback+1):(n)
  for(endpt in endpts){
    aux[endpt] <- get_slope(endpt + (-nback):0, signal)
  }
  return(aux)
}

##' Get the slope of \code{signal[inds]}.
##'
##' @param inds Indices.
##' @param signal Data.
##'
##' @return The slope of a OLS linear regression on \code{signal[inds]}.
##'
get_slope <- function(inds, signal){
  stopifnot(all(inds %in% 1:length(signal)))
  signal_vec <- signal[inds] %>% na.exclude() %>% as.numeric()
  if(any(is.na(signal_vec))){
    slp = NA
  } else {
    x <- (1:length(signal_vec))
    slp = .lm.fit(cbind(1,x), signal_vec)$coef[2] ## MUCH faster version than lm()
  }
  return(slp)
}


##' constructs response variable for one geo_value using the provided function
##'
##' @param df dataframe with ONE geo_value and ONE feature (which will be used
##'   to construct the binary response), columns val, time
##' @param n_ahead number of days ahead that response will be computed
##' @param fn_response logic for computing response, based on a provided
##'   response vector whose all points will be used for this computation
##' @param onset if TRUE, then hotspot is defined as the onset of
##'   increases. Otherwise, a hotspot is defined as an increasing trend,
##'   regardless of the past.
##' @return inputted dataframe with addedmresp colum, which is a binary variable
##'   indicating if there is a hotspot n days ahead of the date variable
response_onegeo <- function(df, n_ahead, fn_response = response_diff_avg, threshold, onset = FALSE,...){
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
    out$resp[i] <- fn_response(signal[1:(i+n_ahead)], i, threshold, ...)
  }
  
  if(onset){
    ## Go back (up to 60 days), and if x% of 14 previous days were a hot spot,
    ## don't deem as hotspot. Use default of x=70%
    out$resp_new = out$resp
    track_past = 14
    for(ii in 1:(len - n_ahead + 1)){
      start = pmax(ii - track_past + 1, 1)
      hot_past = see_past_hotness(out$resp[start:ii])
      out$resp_new[ii] <- (out$resp[ii] & !hot_past)
    }
    out$resp = out$resp_new
    ## If you're coding and you need a break, see this: https://xkcd.com/2346/
  }
  
  return(out)
}

##' From a vector of 0's and 1's, see if at least \code{perc} percent are 1's.
##'
##' @param vec vector of 0's and 1's
##' @param perc percentage of 1's required
##'
##' @return 0 or 1
see_past_hotness <- function(vec, perc = 0.7){
  if(any(is.na(vec))){ vec = vec[which(!is.na(vec))]  }
  stopifnot(all(vec %in% c(0, 1)))
  val = (sum(vec) > length(vec) * perc)
  stopifnot((val %in% c(0, 1)) & length(val)==1 )
  val
}


##' considers increase if there is a 25% increase
##' looks at last weeks average and compare to average of today+1 until today+n_ahead
##'
##' @param x vector of values that will be used to determine hotspot, from 1 until i+n_ahead
##' @param i position of the vector x that is "today"; everything from i+1:forward is not known as features
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff_avg <- function(x, i, threshold){
  len = length(x)
  up = mean(x[(i+1):len])
  low = mean(x[max(1,i-6):i], na.rm=TRUE)
  ifelse(((up-low)/low)>(1+threshold), 1, 0)
}

##' considers increase if there is a 25% increase
##' looks at (last weeks average) and the (avg of the last week between today and n_ahead)
##'
##' @param x vector of values that will be used to determine hotspot, from 1 until i+n_ahead
##' @param i position of the vector x that is "today"; everything from i+1:forward is not known as features
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff_avg_1week <- function(x, i, threshold){
  len = length(x)
  up = mean(x[max(i+1,len-6):len], na.rm=TRUE)
  low = mean(x[max(1,i-6):i], na.rm=TRUE)
  #cat(paste(round(low, 3), round(up,3), round((up-low)/low,3), ifelse(((up-low)/low)>(1+threshold), 1, 0), "\n", sep = " "))
  ifelse(((up-low)/low)>(threshold), 1, 0)
  ## ifelse(((up-low)/low)>(1+threshold), 1, 0)
}

##' considers increase if there is a 25% increase and the response at t+n_ahead has to have a minimum value of 30
##' looks at (last weeks average) and the (avg of the last week between today and n_ahead)
##'
##' @param x vector of values that will be used to determine hotspot, from 1 until i+n_ahead
##' @param i position of the vector x that is "today"; everything from i+1:forward is not known as features
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff_avg_1week_min30 <- function(x, i, threshold){
  len = length(x)
  up = mean(x[max(i+1,len-6):len], na.rm=TRUE)
  low = mean(x[max(1,i-6):i], na.rm=TRUE)
  #cat(paste(round(low, 3), round(up,3), round((up-low)/low,3), ifelse(((up-low)/low)>(1+threshold), 1, 0), "\n", sep = " "))
  ifelse((((up-low)/low)>(threshold))&&(up>=30), 1, 0)
  ## ifelse((((up-low)/low)>(1+threshold))&&(up>=30), 1, 0)
}


##' considers increase if there is a 25% increase and the response at t+n_ahead has to have a minimum value of 30
##' looks at (last weeks average) and the (avg of the last week between today and n_ahead)
##'
##' @param x vector of values that will be used to determine hotspot, from 1 until i+n_ahead
##' @param i position of the vector x that is "today"; everything from i+1:forward is not known as features
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff_avg_1week_min20 <- function(x, i, threshold){
  len = length(x)
  up = mean(x[max(i+1,len-6):len], na.rm=TRUE)
  low = mean(x[max(1,i-6):i], na.rm=TRUE)
  #cat(paste(round(low, 3), round(up,3), round((up-low)/low,3), ifelse(((up-low)/low)>(1+threshold), 1, 0), "\n", sep = " "))
  ifelse((((up-low)/low)>(threshold))&&(up>=20), 1, 0)
  ## ifelse((((up-low)/low)>(1+threshold))&&(up>=30), 1, 0)
}



##' considers increase if there is a 25% increase
##' very simple function, only looks at today's value and target value, no averages
##'
##' @param x vector of values that will be used to determine hotspot, from 1 until i+n_ahead
##' @param i position of the vector x that is "today"; everything from i+1:forward is not known as features
##' @param threshold threshold on increase val to determine hotspot
##' @return 1 if hotspot, 0 if not
response_diff <- function(x, i, threshold){
  len = length(x)
  ifelse((x[len]-x[i])/x[i]>(1+threshold), 1, 0)
}

##' considers increase if there is a 25% increase
##' very simple function, only looks at today's value and target value, no averages
##'
##' TODO add more geo information to the features other than the geo_value's population
##'
##' @param df_all any dataframe that we want to add geographical-level information
##'               it should have the geo_value column!!!
##' @param geo_type one of state, msa, county
##' @return inputted dataframe with geographical location information columns
add_geoinfo <- function(df_all, geo_type){
  if(geo_type == "county"){
    county_pop = county_census %>%
      transmute (geo_value = 1000*as.numeric(STATE) + as.numeric(COUNTY),
                 population = POPESTIMATE2019)
    county_pop$geo_value <- sprintf("%05d", county_pop$geo_value)
    return(inner_join(df_all, county_pop, by = "geo_value"))
  }
  if(geo_type == "state"){
    state_pop <- state_census %>%
      mutate(geo_value = as.numeric(STATE)) %>%
      filter(STATE != 0) %>%
      group_by(geo_value) %>%
      summarise(population = sum(POPESTIMATE2019), .groups = 'drop')
    state_crosswalk <- maps::state.fips %>%
      select(abb, fips) %>% distinct() %>% mutate(abb = tolower(abb))
    state_pop <- state_pop %>%
      inner_join(state_crosswalk, by = c("geo_value" = "fips")) %>%
      select(geo_value = abb, population)
    return(inner_join(df_all, state_pop, by = "geo_value"))
  }
  if(geo_type == "msa"){
    return(inner_join(df_all, msa_census %>%
                        transmute(geo_value = as.character(CBSA),
                                  population = POPESTIMATE2019), by = "geo_value"))
  }
}


##' Outputs data ready for modeling, with all lagged features and binary
##'    response uses output from API call.
##'
##' @param mat resulted from the API call with all signals we want to use for
##'   model construction.
##' @param lags number of past values to include in the data frame; for time t,
##'   dataframe will have in one row X_t until X_{t-lag}.
##' @param n_ahead number of days ahead that response will be computed.
##' @param response name of the response variable in mat.
##' @param fn_response logic for computing response, based on a provided
##'   response vector whose all points will be used for this computation.
##' @param threshold threshold on increase val to determine hotspot.
##' @param slopes if TRUE, produces a dataframe with slopes based on the past
##'   feature values and if FALSE, produces raw lagged features.
##' @param onset if TRUE, then hotspot is defined as the onset of
##'   increases. Otherwise, a hotspot is defined as an increasing trend,
##'   regardless of the past.

##' @return dataset ready to be fed to fitting functions
ready_to_model <- function(mat, lags, n_ahead,
                           response = "confirmed_7dav_incidence_num",
                           slope = FALSE, fn_response = response_diff_avg_1week,
                           threshold = .25,
                           onset = FALSE){
  
  ## Construct lagged features for all available signals, including lagged responses
  # TODO deal with potential NAs?
  features <- mat %>% plyr::ddply(c("signal", "data_source", "geo_value"),
                                  lagged_features_onegeo, lags = lags, slope = slope) %>% na.omit()
  
  ## Construct hotspot indicator in the resp variable
  responses <- mat %>% filter(signal == response) %>%
    plyr::ddply(c("signal", "data_source", "geo_value"), response_onegeo,
                n_ahead = n_ahead, fn_response = fn_response, threshold = threshold,
                onset = onset) %>% na.omit()
  
  ## transform the dataframe in a wide format, with one row per geo_value and date
  names_to_pivot <- colnames(features %>% select(-geo_value, -time_value, -signal, -data_source))
  features <- pivot_wider(features, id_cols = c("geo_value", "time_value"), names_from = c("signal", "data_source"),
                          values_from = all_of(names_to_pivot)) %>% ungroup
  
  ## join features and response
  mat_to_model <- inner_join(features, responses %>% select(-signal, -data_source),
                             by = c("geo_value", "time_value")) %>% na.omit()
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


##' Wrapper around \code{sample_split_geo} to continue to sample until the 0/1
##' distribution of the data train and test is very close
##'
##' @param df_model Output from \code{ready_to_model()}.
##' @param pct_test Percent of data to use for test data.
##'
##' @return The same format of out put as \code{sample_split_geo()}.
stratified_sample_split_geo <- function(df_model, pct_test = 0.3, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  nsim = 1000
  for(ii in 1:nsim){
    splitted = sample_split_geo(df_model, pct_test = pct_test)
    ratio_train = splitted$df_train %>% select(resp) %>% table() %>% (function(a){a["1"]/a["0"]})
    ratio_test = splitted$df_test %>% select(resp) %>% table() %>% (function(a){a["1"]/a["0"]})
    if(abs(ratio_train - ratio_test) < 0.02) break
  }
  return(splitted)
}

##' Sample splitting by geo valuse.
##'
##' @param df_model Output from \code{ready_to_model()}.
##' @param pct_test Percent of data to use for test data.
##'
##' @return List containing \code{df_train} and \code{df_test}.
sample_split_geo <- function(df_model, pct_test = 0.3, seed=NULL){
  
  ## df_tomodel <- df_tomodel %>% arrange(desc(time_value)) %>% na.omit()
  ## start_test_date <- df_tomodel[round(pct_test*nrow(df_tomodel)),"geo_value"]
  ## start_test_date <- df_tomodel[round(pct_test*nrow(df_tomodel)),"geo_value"]
  geos = df_model %>% select(geo_value) %>% unlist() %>% unique()
  if(!is.null(seed)) set.seed(seed)
  test_ind =  sample(length(geos), length(geos) * pct_test)
  test_geos = geos[test_ind]
  train_geos = geos[-test_ind]
  df_test <- df_model %>% filter(geo_value %in% test_geos)
  df_train <- df_model %>% filter(geo_value %in% train_geos)
  
  stopifnot(intersect(
    df_test %>% select(geo_value) %>% unique() %>% unlist(),
    df_train %>% select(geo_value) %>% unique() %>% unlist()) %>% length() == 0)
  
  ## Todo: check if train and test have equal number of hot spots. Doesn't seem
  ## to be a big problem since we are naively splitting geos, but still..
  
  ## ## Helper to calculate rise from June to July
  ## howmuchrise <- function(val, time){
  ##  mean(val[time >= "2020-7-26"]) - mean(val[time <= "2020-6-06"])
  ## }
  
  ## ## Sort all the rises
  ## rises = mat %>%
  ##   group_by(geo_value) %>%
  ##   dplyr::summarise(rise = howmuchrise(value, time_value)) %>%
  ##   arrange(desc(rise))
  
  ## ## ## See the numerical summary
  ## ## rises %>% dplyr::select(rise) %>% unlist() %>% summary()
  
  ## ## See the numerical summary
  ## topn = 40
  ## geos = rises %>%  slice_head(n = topn) %>% dplyr::select(geo_value) %>% unlist()
  ## set.seed(0)
  ## train_i = sample(topn, topn/5)
  ## train_hotspot_geos = geos[train_i]
  ## test_hotspot_geos = geos[-train_i]
  
  ## ## mat %>% select(contains("geo")) %>% unlist() %>% unique() %>% length()
  ## df_no_hotspot <- df_model %>% filter(!(geo_value %in% geos))
  ## df_yes_hotspot <- df_model %>% filter(geo_value %in% geos)
  
  # Make training set and test set separately
  return(list(df_test = df_test, df_train = df_train))
}


##' Wrapper for \code{make_foldid_geo} to sample a set of CV fold IDs that have
##' a somewhat even distribution of 0s and 1s among the folds.
##'
make_stratified_foldid_geo <- function(x, nfold, seed=NULL){
  
  if(!is.null(seed)) set.seed(seed)
  nsim = 10000
  for(ii in 1:nsim){
    ## x = splitted$df_train
    foldid = make_foldid_geo(x, nfold)
    ratios = sapply(1:nfold, function(ifold){
      x[which(foldid==ifold),] %>%
        select(resp) %>%
        table() %>%
        (function(a){a["1"]/a["0"]})
    })
    ## ratios %>% round(2) %>% print()
    if(any(is.na(ratios))) next
    if(all(abs(ratios - (1/nfold)) < (1/nfold) * 0.5) ) break
  }
  return(foldid)
}

#' Make |foldid| argument for covariate matrix |x| and |nfold|-fold
#' cross-validation; makes nfold geo partitions.
#'
#' @param x covariate matrix.
#' @param nfold nfold.
#' @param seed seed number to be used for splitting geos..
#'
#' @return A numeric vector containing elements of \code{(1:nfold)} specifying
#'   row numbers of X (or entry numbers of y) to be used for each CV fold.
make_foldid_geo <- function(x, nfold, seed=NULL){
  
  ## Validation_geos.
  geos = x %>% select(geo_value) %>% unlist()
  unique_geos = geos %>% unique() %>% sort()
  if(!is.null(seed)) set.seed(seed)
  geo_blocks = split(sample(unique_geos),
                     sort(1:length(unique_geos) %% nfold))
  cv_inds = lapply(1:nfold, function(ifold){
    which(geos %in% geo_blocks[[ifold]])
  })
  
  final_inds = rep(NA, length(geos))
  for(ifold in 1:nfold){
    inds = cv_inds[[ifold]]
    final_inds[inds] = ifold
  }
  
  ## Quick test before returning
  stopifnot(all((x[which(final_inds==1),] %>%
                   select(geo_value) %>%
                   unique() %>% unlist()) %in% geo_blocks[[1]]))
  return(final_inds)
}


#' Make |foldid| argument for covariate matrix |x| and |nfold|-fold
#' cross-validation; makes nfold consecutive time blocks.
#'
#' @param x covariate matrix
#' @param nfold nfold
#'
#' @return A numeric vector containing elements of \code{(1:nfold)} specifying
#'   row numbers of X (or entry numbers of y) to be used for each CV fold.
make_foldid <- function(x, nfold){
  
  ## The fold ids need to be made temporal
  times = x %>% select(time_value) %>% unlist() %>% as_date()
  unique_sorted_times = times %>% unique() %>% sort()
  
  endpoints = round(seq(from = 0, to = length(unique_sorted_times), length = nfold+1))
  inds = Map(function(a,b){(a+1):b},
             endpoints[-length(endpoints)],
             endpoints[-1])
  time_blocks = lapply(inds, function(ind){ unique_sorted_times[ind]})
  
  cv_inds = lapply(1:nfold, function(ifold){
    which(times %in% time_blocks[[ifold]])
  })
  
  final_inds = rep(NA, length(times))
  for(ifold in 1:nfold){
    inds = cv_inds[[ifold]]
    final_inds[inds] = ifold
  }
  return(final_inds)
}


##' fit models and produce test set predictions
##' currently: lasso, ridge
##'
##' @param df_train
##' @param df_test
##' @param lags number of past values to include in the data frame; for time t, dataframe will have in one row X_t until X_{t-lag}
##' @param n_ahead number of days ahead that response will be computed
##' @param response string just for rendering plots later on
##' @return
fit_predict_models <- function(df_train, df_test, lags, n_ahead, response = "confirmed_7dav_incidence_num",
                               foldid = NULL,
                               geo_cv_split_seed = NULL,
                               stratify_cv_split = TRUE,
                               verbose = FALSE){
  if(verbose) cat("Fitting models:\n")
  
  predictions <- df_test %>% select(geo_value, time_value, resp)
  
  if(verbose) cat("\tFitting LASSO...")
  suppressMessages({
    preds <- fit_logistic_regression(df_train, df_test, nfold = 5, alpha = 1,
                                     foldid = foldid,
                                     geo_cv_split_seed = geo_cv_split_seed,
                                     stratify_cv_split = stratify_cv_split)
  })
  predictions[[paste("lasso_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  if(verbose) cat(" Done!\n")
  
  
  ## Note: only using lasso for now.
  
  if(verbose) cat("\tFitting Ridge...")
  suppressMessages({
    preds <- fit_logistic_regression(df_train, df_test, nfold = 5, alpha = 0,
                                     foldid = foldid,
                                     geo_cv_split_seed = geo_cv_split_seed,
                                     stratify_cv_split = stratify_cv_split)
  })
  predictions[[paste("ridge_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  if(verbose) cat(" Done!\n")
  
  ## cat("\tFitting SVM...")
  ## preds <- fit_svm(df_train, df_test)
  ## predictions[[paste("svm_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  ## cat(" Done!\n")
  
  if(verbose) cat("\tFitting xgboost...")
  suppressMessages({
    preds <- fit_xgb(df_train, df_test)
  })
  predictions[[paste("xgb_lags", lags, "_nahead", n_ahead, sep = "")]] = preds
  if(verbose) cat(" Done!\n")
  
  ### can add more models here!!! add \hat{y} as a col to |predictions|
  
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
##' @param nfold 5 (previously 10).
##' @param alpha 1 for lasso, or 0 for ridge regression. Used by
##'   \code{glmnet()}.
##' @param geo_cv_split_seed Random seed for geo split for cross validation by
##'   \code{glmnet::glmnet()}.
##'
##' @return Numeric vector the same length as \code{nrow(df_test)}.
fit_logistic_regression <- function(df_train, df_test, nfold = 5, alpha = 1,
                                    foldid = NULL,
                                    geo_cv_split_seed = NULL,
                                    stratify_cv_split = TRUE){
  
  ## Input checks (should be common for all fit_OOOO() functions
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_train)))
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_test)))
  
  ## Input check
  stopifnot(alpha %in% c(0,1)) ## Only allow ridge or lasso for now.
  
  ## (Not used for now) Make contiguous time blocks for CV
  ## foldid <- make_foldid(df_train, nfold)
  
  ## Split by geo
  if(!is.null(foldid)){
    if(stratify_cv_split){
      foldid <- make_stratified_foldid_geo(df_train, nfold, seed = geo_cv_split_seed)
    } else {
      foldid <- make_foldid_geo(df_train, nfold, seed = geo_cv_split_seed)
    }
  } else {
    stopifnot(all(foldid %in% 1:nfold))
  }
  stopifnot(length(foldid) == nrow(df_train))

  ## Main part of the lasso fitting and predicting (todo: use quantile lasso's forward-looking CV functionality)
  fit_lasso <- cv.glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)),
                         y = df_train$resp,
                         family = "binomial",
                         alpha = alpha,
                         foldid = foldid,
                         nfold = nfold)
  preds = predict(fit_lasso, s = "lambda.min",
                  newx = as.matrix(df_test %>% select(-geo_value, -time_value, -resp)), type = "response")[,1]
  
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

##' Performs xgboost prediction, given training & test matrices.
##' Outputs a vector of values the same as y.
##' NOTE: this does not perform CV!!! it just uses mainly default hyperparameter values. CV would take time to run and we have higher priority things to worry for now
##'
##' @param df_train Training matrix. Must contain columns "geo_value",
##'   "time_value", "resp", and some other columns that will be used as
##'   covariates.
##' @param df_test Test matrix. Same format as df_train.
##' @param ... Additional functions to \code{xgb.train()} of the \code{xgboost} R
##'   package.
##'
##' @return Numeric vector the same length as \code{nrow(df_test)}.
fit_xgb <- function(df_train, df_test){
  
  ## Input checks (should be common for all fit_OOOO() functions
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_train)))
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_test)))
  
  
  ## transforms dataframe to XGBoost data format
  dtrain <- xgb.DMatrix(as.matrix(df_train %>% dplyr::select(-resp, -geo_value, -time_value)),
                        label = as.vector(df_train$resp))
  dtest <- xgb.DMatrix(as.matrix(df_test %>% dplyr::select(-resp, -geo_value, -time_value)),
                       label = as.vector(df_test$resp))
  ## Fit xgboost and make predictions
  mod <- xgb.train(booster = "gbtree",
                   data = dtrain,
                   nthread = 5,
                   eta = 0.3,
                   gamma = 1,
                   max_depth = 6,
                   nrounds = 500,
                   ## nfold = 5,
                   objective = "binary:logistic",
                   colsample_bytree = 0.7,
                   subsample = 0.7)
  
  preds <- predict(mod, dtest)
  
  ## Out checks (should be common for all fit_OOOO() functions)
  stopifnot(length(preds) == nrow(df_test))
  
  preds
}




##' Performs random forest prediction, given training & test matrices.  Outputs
##' a vector of values the same as.
##'
##' @param df_train Training matrix. Must contain columns "geo_value",
##'   "time_value", "resp", and some other columns that will be used as
##'   covariates.
##' @param df_test Test matrix. Same format as df_train.
##' @param ... Additional functions to \code{randomForest()} of the
##'   \code{randomForest} R package.
##'
##' @return Numeric vector the same length as \code{nrow(df_test)}.
fit_random_forest <- function(df_train, df_test, ...){
  
  ## Input checks (should be common for all fit_OOOO() functions
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_train)))
  stopifnot(all(c("time_value", "geo_value", "resp") %in% colnames(df_test)))
  
  
  train_mat <- df_train %>% select(-geo_value, -time_value)
  test_mat <- df_test %>% select(-geo_value, -time_value)
  
  ## X and y matrix for randomForest() function.
  X = train_mat %>% select(-resp)
  y = train_mat %>% select(resp) %>% unlist()
  
  ## Tip: start small, and scale up slowly.
  ## rf = randomForest(X, y,sampsize=1000, ntree=5)
  ## print(rf)
  rf = randomForest(X, y, sampsize = 5000, ntree = 500)
  print(rf)
  
  
  ## (not written yet)
  test_X = test_mat %>% select(-resp)
  preds <- predict(model, newdata = test_X, type = "prob")
  
  ## This seems like a good quick R random forest guide:
  ## https://stackoverflow.com/questions/46124424/how-can-i-draw-a-roc-curve-for-a-randomforest-model-with-three-classes-in-r
  
  ## This is about how to speed it up:
  ## https://stackoverflow.com/questions/34706654/get-randomforest-regression-faster-in-r
  
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
    state_crosswalk <- maps::state.fips %>%
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

##' computes population weighted precision, population weighted recall, and population weighted proportion of predicted 1's for different cutoffs
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
    
    wprecision = (df_temp %>% filter(pred == 1, resp == 1) %>% summarise(sum(population)) %>% unlist)/(df_temp %>% filter(pred == 1) %>% summarise(sum(population)) %>% unlist)
    if((df_temp %>% filter(pred == 1) %>% summarise(sum(population)) %>% unlist) == 0) wprecision = 1
    
    wrecall = (df_temp %>% filter(pred == 1, resp == 1) %>% summarise(sum(population)) %>% unlist)/(df_temp %>% filter(resp == 1) %>% summarise(sum(population)) %>% unlist)
    if((df_temp %>% filter(resp == 1) %>% summarise(sum(population)) %>% unlist) == 0) wrecall = 1
    
    ## computing weighted % of predicted 1's
    wpred1 =(df_temp %>% filter(pred == 1) %>% summarise(sum(population)) %>% unlist)/sum(df_temp$population)
    return(c(i, wpred1, wprecision, wrecall))
  })
  ## transforming matrix into dataframe and naming it appropriately
  metrics <- as.data.frame(t(metrics))
  names(metrics) <- c("cutoff","wpred1", "wprecision", "wrecall")
  return(metrics)
}

##' computes metrics for all models and produces roc curves (our adapted roc with different metrics)
##'
##' @param predictions dataframe with cols: geo_value, time_value, resp,
##'                    and one column per model with predicted values whose col name is the models name
##' @param geo_type county, msa, or state. will be used to get population data
##' @param add if TRUE, adds current curves to an existing plot
##' @param df_plot_existing plot that will have mroe curves added to it, and these new curves will be labeled as FbFeatures
##'
##' @return ggplot of model comparison curve
plot_adapted_roc <- function(predictions, geo_type = "county", add = FALSE, df_plot_existing = NULL){
  df_temp <- inner_join(predictions, get_population(geo_type), by = "geo_value")
  df_temp <- reshape2::melt(df_temp, id.vars = c("geo_value", "time_value", "resp", "population"))
  df_plot <- df_temp %>% rename(model = variable) %>% group_by(model) %>% group_modify(adapted_roc_onemodel)
  
  precision_thresh <- df_temp %>% select(resp, population, geo_value, time_value) %>%
    distinct() %>% filter(resp == 1) %>%
    summarise(sum(population)/sum(df_temp$population), .groups = 'drop') %>%
    unlist()
  
  if(!add){
    ggplot(df_plot, aes(x = wpred1, color = model)) +
      geom_vline(xintercept = precision_thresh, size = 1.25, col = "gray30", alpha = .8) +
      geom_line(aes(y = wprecision, linetype ="wprecision", size = "LaggedResponse", alpha = "LaggedResponse")) +
      geom_line(aes(y = wrecall, linetype = "wrecall", size = "LaggedResponse", alpha = "LaggedResponse")) +
      scale_linetype_manual(name = "",
                            values = c( "wprecision" = 1, "wrecall" = 2),
                            labels = c("Precision", "Recall")) +
      scale_size_manual(name = "",
                        values = c( "LaggedResponse" = 0.5, "FbFeatures" = 1.5)) +
      scale_alpha_manual(name = "",
                         values = c( "LaggedResponse" = 1, "FbFeatures" = 0.4)) +
      ylim(0,1) +
      xlim(0,1) +
      theme_bw(base_size = 18) +
      guides(color=FALSE, size=guide_legend(nrow=2,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("population weighted precision") +
      xlab("population weighted % predicted hotspots") +
      scale_y_continuous(sec.axis = sec_axis(~., name = "population weighted recall (dashed)")) +
      theme(legend.position = "bottom") +
      scale_x_log10() +
      facet_wrap(~model)
  } else {
    df_plot_existing +
      geom_line(data=df_plot,  aes(x=wpred1, y = wprecision, linetype ="wprecision", size = "FbFeatures", alpha ="FbFeatures")) +
      geom_line(data=df_plot, aes(y=wrecall, linetype = "wrecall", size = "FbFeatures", alpha = "FbFeatures"))
  }
}


##' computes sensitivity and specificity (population weighted or not) for different cutoffs
##' considers ONE MODEL only
##'
##' @param df_one dataframe for one model with at least columns: value, pred, resp, population
##' @param popweighted indicates if metrics for ROC curve should be population weighted
##'
##' @return dataframe with colulmns cutoff, specificity, sensitivity
roc_onemodel <- function(df_one, popweighted = FALSE){
  df_one <- df_one %>% arrange(value)
  ## changing cutoffs
  metrics <- sapply(seq(0, 1, 0.005), function(i){
    df_temp <- df_one
    ## binary predictions for the cutoff
    df_temp$pred <- ifelse(df_temp$value<=i, 0, 1)
    
    ## specificity
    # pop weighted
    if(popweighted){
      specificity = (df_temp %>% filter(pred == 0, resp == 0) %>% summarise(sum(population)) %>% unlist)/(df_temp %>% filter(resp == 0) %>% summarise(sum(population)) %>% unlist)
      if((df_temp %>% filter(resp == 0) %>% summarise(sum(population)) %>% unlist) == 0) specificity = 1
    } else{ # NOT pop weighted
      specificity = (df_temp %>% filter(pred == 0, resp == 0) %>% nrow())/(df_temp %>% filter(resp == 0) %>% nrow())
      if((df_temp %>% filter(resp == 0) %>% nrow()) == 0) specificity = 1
    }
    
    ## sensitivity = recall
    # pop weighted
    if(popweighted){
      sensitivity = (df_temp %>% filter(pred == 1, resp == 1) %>% summarise(sum(population)) %>% unlist)/(df_temp %>% filter(resp == 1) %>% summarise(sum(population)) %>% unlist)
      if((df_temp %>% filter(resp == 1) %>% summarise(sum(population)) %>% unlist) == 0) sensitivity = 1
    } else{ # NOT pop weighted
      sensitivity = (df_temp %>% filter(pred == 1, resp == 1) %>% nrow())/(df_temp %>% filter(resp == 1) %>% nrow())
      if((df_temp %>% filter(resp == 1) %>% nrow()) == 0) sensitivity = 1
    }
    
    return(c(i, specificity, sensitivity))
  })
  ## transforming matrix into dataframe and naming it appropriately
  metrics <- as.data.frame(t(metrics))
  names(metrics) <- c("cutoff", "specificity", "sensitivity")
  return(metrics)
}





##' computes metrics for all models and produces roc curves (our adapted roc with different metrics)
##'
##' @param predictions dataframe with cols: geo_value, time_value, resp,
##'                    and one column per model with predicted values whose col name is the models name
##' @param geo_type county, msa, or state. will be used to get population data
##' @param add if TRUE, adds current curves to an existing plot
##' @param df_plot_existing plot that will have mroe curves added to it, and these new curves will be labeled as FbFeatures
##' @param popweighted indicates if metrics for ROC curve should be population weighted
##'
##' @return ggplot of model comparison curve
plot_roc <- function(predictions, geo_type = "county", add = FALSE, df_plot_existing = NULL, popweighted = FALSE, only_return_auc=FALSE){
  
  df_temp <- inner_join(predictions, get_population(geo_type), by = "geo_value")
  df_temp <- reshape2::melt(df_temp, id.vars = c("geo_value", "time_value", "resp", "population"))
  df_auc <- df_temp %>% rename(model = variable) %>% group_by(model) %>% group_modify(function(df, ...){data.frame(auc = round(pROC::auc(response = df$resp, predictor = df$value)[1], 3))})
  if(only_return_auc) return(df_auc)
  df_plot <- df_temp %>% rename(model = variable) %>% group_by(model) %>% group_modify(~roc_onemodel(.x, popweighted = popweighted))
  
  if(!add){
    ggplot(df_plot, aes(x = 1-specificity,y=sensitivity, color = model, size = "LaggedResponse", alpha = "LaggedResponse")) +
      geom_line() +
      geom_abline(slope = 1, intercept = 0, size = 1.25, col = "gray30", alpha = .8)  +
      geom_text(data = df_auc, mapping = aes(x = rep(.5,nrow(df_auc)), y = seq(.2, .2-0.05*(nrow(df_auc)-1), -0.05), color = model, label = auc), size = 5, show.legend = FALSE) +
      annotate(geom="text", x=.5, y=.25, label="AUC") +
      scale_size_manual(name = "",
                        values = c( "LaggedResponse" = 0.5, "FbFeatures" = 1.5)) +
      scale_alpha_manual(name = "",
                         values = c( "LaggedResponse" = 1, "FbFeatures" = 0.4)) +
      ylim(0,1) +
      xlim(0,1) +
      theme_bw(base_size = 18) +
      guides(color=guide_legend(nrow=2,byrow=TRUE), size=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("sensitivity") +
      xlab("1-specificity") +
      theme(legend.position = "bottom")
  } else {
    df_plot_existing +
      geom_line(data=df_plot,  aes(x = 1-specificity,y=sensitivity, size = "FbFeatures", alpha ="FbFeatures")) +
      geom_text(data = df_auc, mapping = aes(x = rep(.65,nrow(df_auc)), y = seq(.2, .2-0.05*(nrow(df_auc)-1), -0.05), color = model, label = auc), size = 5, show.legend = FALSE) +
      annotate(geom="text", x=.65, y=.25, label="FbAUC")
  }
}


make_preds <- function(destin = "figures", splitted, lags, n_ahead, geo_type,
                       response,
                       fn_response_name, threshold, slope, onset,
                       geo_cv_split_seed = NULL,
                       include_fb = TRUE){
  
  ######################################
  ## Model with lagged responses only ##
  ######################################
  if(!include_fb){
    predictions_onlylaggedresponse <- fit_predict_models(splitted$df_train %>% select(geo_value, time_value, resp, contains(response)),
                                                         splitted$df_test %>% select(geo_value, time_value, resp, contains(response)),
                                                         lags = lags, n_ahead = n_ahead,
                                                         geo_cv_split_seed = geo_cv_split_seed)
    preds = predictions_onlylaggedresponse
  }
  ####################################################
  ## Model with lagged responses + facebook signals ##
  ####################################################
  if(include_fb){
    predictions_laggedandfacebook <- fit_predict_models(splitted$df_train, splitted$df_test,
                                                        lags = lags, n_ahead = n_ahead,
                                                        geo_cv_split_seed = geo_cv_split_seed)
    preds = predictions_laggedandfacebook
  }
  return(preds)
}

##' Fit various models to train data, and make plots for test data.
##'
##' @param destin Where to save plots.
##' @param splitted a list containing train and test data.
##'
##' @return
make_plots <- function(destin = "figures", splitted, lags, n_ahead, geo_type,
                       response,
                       fn_response_name, threshold, slope, onset,
                       geo_cv_split_seed = NULL){
  
  ######################################
  ## Model with lagged responses only ##
  ######################################
  predictions_onlylaggedresponse <- fit_predict_models(splitted$df_train %>% select(geo_value, time_value, resp, contains(response)),
                                                       splitted$df_test %>% select(geo_value, time_value, resp, contains(response)),
                                                       lags = lags, n_ahead = n_ahead,
                                                       geo_cv_split_seed = geo_cv_split_seed)
  a = plot_adapted_roc(predictions_onlylaggedresponse, geo_type = geo_type)
  
  ####################################################
  ## Model with lagged responses + facebook signals ##
  ####################################################
  predictions_laggedandfacebook <- fit_predict_models(splitted$df_train, splitted$df_test,
                                                      lags = lags, n_ahead = n_ahead,
                                                      geo_cv_split_seed = geo_cv_split_seed)
  
  b = plot_adapted_roc(predictions_laggedandfacebook, add=TRUE, df_plot_existing=a, geo_type = geo_type)
  plot_adj_roc = b
  plotname_root = paste0(geo_type, "__resp_", threshold*100, "__lag_", lags,"__nahead_",
                         n_ahead, "__slope_", slope,
                         "__onset_", onset)
  plotname_adj_roc = paste0(plotname_root, ".png")
  plotname_roc = paste0(plotname_root, "_ROC", ".png")
  ggsave(plot = b,
         filename = file.path(destin, fn_response_name, plotname_adj_roc),
         width = 12, height = 8, dpi = 200)
  
  ## Also plot regular ROC urves
  a = plot_roc(predictions_onlylaggedresponse, geo_type = geo_type, popweighted = FALSE)
  b = plot_roc(predictions_laggedandfacebook, add=TRUE, df_plot_existing=a, geo_type = geo_type, popweighted = FALSE)
  ggsave(plot = b,
         filename = file.path(destin, fn_response_name,  plotname_roc),
         width = 12, height = 8, dpi = 200)
  plot_roc = b
  return(list(adj_roc = plot_adj_roc,
              roc = plot_roc))
}


##' Fit various models to train data, and make plots for test data.
##'
##' @param destin Where to save plots.
##' @param splitted a list containing train and test data.
##'
##' @return
calc_auc <- function(destin = "figures", splitted, lags, n_ahead, geo_type,
                     fn_response_name, threshold, slope,  onset,
                     foldid = NULL,
                     geo_cv_split_seed = NULL){
  
  ######################################
  ## Model with lagged responses only ##
  ######################################
  predictions_onlylaggedresponse <- fit_predict_models(splitted$df_train %>% select(geo_value, time_value, resp, contains(response)),
                                                       splitted$df_test %>% select(geo_value, time_value, resp, contains(response)),
                                                       lags = lags, n_ahead = n_ahead,
                                                       foldid = foldid,
                                                       geo_cv_split_seed = geo_cv_split_seed,
                                                       verbose = TRUE)
  df_auc_no_fb = plot_roc(predictions_onlylaggedresponse, geo_type = geo_type, popweighted = FALSE, only_return_auc = TRUE)
  
  ####################################################
  ## Model with lagged responses + facebook signals ##
  ####################################################
  predictions_laggedandfacebook <-
    fit_predict_models(splitted$df_train, splitted$df_test, lags = lags, n_ahead = n_ahead,
                       foldid = foldid, geo_cv_split_seed = geo_cv_split_seed)
  df_auc_yes_fb = plot_roc(predictions_laggedandfacebook,
                           geo_type = geo_type,
                           popweighted = FALSE,
                           only_return_auc = TRUE)
  
  ## column_names = paste0(c("no_fb", "yes_fb"), paste0("_n_ahead", n_ahead))
  column_names = c("no_fb", "yes_fb")
  auc_df = df_auc_no_fb %>% bind_cols(df_auc_yes_fb, .id = "model")
  auc_df = df_auc_no_fb %>% full_join(df_auc_yes_fb, by="model")
  colnames(auc_df)[2:3] =  column_names
  return(auc_df)
}


##' Main function to get a data frame of the form y|X.
get_data <- function(geo_type = "state", lags = 28, n_ahead = 28, threshold = 0.25,
                     response = "confirmed_7dav_incidence_prop",
                     fn_response = response_diff_avg_1week_min20,
                     fn_response_name = "response_diff_avg_1week_min20",
                     slope = TRUE,
                     onset = FALSE,
                     start_day = as.Date("2020-05-01"),
                     end_day = as.Date("2020-08-30")
){

  ########################
  ## Read in data once ###
  ########################
  data_sources = c("indicator-combination",
                   "fb-survey")
  signals = c("confirmed_7dav_incidence_prop",
              "smoothed_hh_cmnty_cli")
  response = "confirmed_7dav_incidence_prop"
  signals = data.frame(data_sources = data_sources, signals = signals)
  suppressMessages({
    mat = covidcast_signals(signals,
                            start_day = start_day, end_day = end_day, geo_type = geo_type)
  })
  mat <- mat %>% select(geo_value, time_value, signal, data_source, value)
  

  ##############################################################################
  ## only includes counties that have at most 5% of missing fb survey signals ##
  ##############################################################################
  mat$time_value <- as.Date(mat$time_value)
  geo_pct_missing <- mat %>% filter(data_source == 'fb-survey') %>%
    plyr::ddply("geo_value", function(df){
      data.frame(prop_missing = 1 - length(df$time_value)/as.numeric(end_day - start_day + 1))
  })
  geo_include <- geo_pct_missing %>% filter(prop_missing <= 0.05) %>% select(geo_value) %>% unlist()
  mat <- mat %>% filter(geo_value %in% geo_include)

  ##########################
  ## Form the y|X matrix ###
  ##########################
  df_model <- ready_to_model(mat, lags, n_ahead, response, slope, fn_response, threshold, onset)
  return(list(mat = mat,
              df_model = df_model,
              geo_type = geo_type,
              lags = lags,
              n_ahead = n_ahead,
              threshold = threshold,
              response = response,
              fn_response = fn_response,
              fn_response_name = fn_response_name,
              slope = slope,
              onset = onset,
              start_day = start_day,
              end_day = end_day))
}




##' Renumber any collection of (1,..,6) e.g. (4,4,2,2,3,3,3,1,6,6,6) to
##' (1,1,2,2,3,3,3,4,5,5,5).
##'
##' @param foldid CV Fold IDs.
##' @param nfold Original number of folds intended
##'
##' @return Renumbered integer vector of the same length.
renumber_foldid <- function(foldid, nfold){

  ## Basic check
  ## foldid = c(1,1,1,4,4,4,2,2,2,3,3,3,6,6,6)
  stopifnot(all(foldid %in% 1:(nfold+1)))

  ## Renumber foldid
  renumbered_foldid = foldid
  all_ids = unique(foldid)
  for(ii in 1:length(all_ids)){
    old_id = all_ids[ii]
    renumbered_foldid[which(foldid == old_id)] = ii
    ## print(paste0("Changing ", old_id))
    ## print(renumbered_foldid)
  }

  ## Basic checks before returning
  stopifnot(length(renumbered_foldid) == length(foldid))
  stopifnot(all(renumbered_foldid %in% 1:nfold))
  return(renumbered_foldid)
}



##' Creates (nfold) splits, and one at a time, designates one fold as test data,
##' and the rest as training data.
##'
##' @param df_model Data matrix.
##'
##' @param nfold How many splits to make
##'
##' @return A (nfold) lengthed list containing the splitted training/test data
##'   sets from \code{df_model}.
##'
outer_split <- function(df_model, nfold = 4){
  foldid = make_stratified_foldid_geo(df_model, nfold)
  splitted_list <- lapply(1:nfold, function(ifold){
    splitted = list(df_test = df_model[which(foldid == ifold), ],
                    df_train = df_model[which(foldid != ifold), ])
    return(splitted)
  })
  return(splitted_list)
}
