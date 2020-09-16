

#########################################################################################
# make lagged features
#########################################################################################
lagged_weekly_features_onegeo <- function(df, lags, slope_lags, name = "feature",slopes = FALSE){
  
  ## The fixed number of time lags for the data values themselves.
  df <- df %>% arrange(week)
  signal <- df$value
  timestamp <- df$week
  
  ## if you want more lags than available points, returns empty dataframe
  len <- nrow(df)
  if(len <= max(lags, slope_lags)){
    return(data.frame())
  }
  
  out <- data.frame(time_value = timestamp[(lags+1):len])
  
  ################################################################
  ## adding lagged feature from t-0, t-1, t-2, until t-lags  #####
  ################################################################
  for(i in 0:lags){
    inds = (lags+1-i):(len-i)
    out <- suppressMessages(bind_cols(out, signal[inds]))
  }
  names(out) = c("week", paste(name, "_lag", 0:lags, sep = ""))
  
  #############################################
  ## adding slopes feature every 1 points #####
  #############################################
  if(slopes){
    npoints = slope_lags+1
    nfeats = floor(npoints/1) ## 1 is a magic number. will construct a new feature (new slope) every 1 points
    limits_lm = round(seq(slope_lags+1, 1, length.out = nfeats+1))
    ## out[[paste(name, "_lag0", sep = "")]] = signal[(slope_lags+1):(len)] ## also adds lag0 to the slopes feature matrix
    for(j in 1:nfeats){
      aux <- rep(NA, nrow(out))
      row_pos <- 1
      for(i in (slope_lags+1):(len)){
        signal_vec <- signal[i:(limits_lm[j+1]+row_pos-1)] %>% na.exclude() %>% as.numeric()
        if(any(is.na(signal_vec))){
          aux[row_pos] <- NA
        } else{
          x <- (1:length(signal_vec))
          aux[row_pos] <- .lm.fit(cbind(1,x), signal_vec)$coef[2] ## MUCH faster version than lm()
        }
        row_pos <- row_pos + 1
      }
      out[[paste(name, "_slope", j, sep = "")]] <- aux
    }
  }
  
  return(out)
}



#########################################################################################
# make response columns
#########################################################################################
cases_ahead_weekly_onegeo <- function(df, n_ahead){
  signal <- df$value
  timestamp <- df$week
  
  ## we can only determine a hotspot n_ahead days if that day is available
  len <- nrow(df)
  stopifnot(n_ahead <= len)
  
  out <- data.frame(week = timestamp[1:len])
  for (n in 1:n_ahead){
    out[, paste0("cases_ahead", n)] <- NA
    ## weekly response n weeks ahead
    for(i in 1:(len-n)){
      out[i, paste0("cases_ahead", n)] <- signal[i+n]
    }
  }
  
  return(out)
}



#########################################################################################
# response functions
#########################################################################################
response_diff_avg_1week_min20 <- function(df, n_ahead, threshold, response_name, fn_response_name){
  future = df[, paste0("cases_ahead", n_ahead)]
  current = df[, response_name]
  df[fn_response_name] = ((future-current)/current > threshold) * (future >= 20)
  return(df)
}


