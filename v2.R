## A mirror of the script in v1.Rmd
outputdir = "./figures"

##' Parse command line arguments and assigns the values of them. |args| is meant
##' to just be additional command line arguments. Something like "Rscript somefile
##' @param args Argument.
##' @examples \dontrun{parse_args(args = commandArgs(trailingOnly = TRUE))}
parse_args <- function(args){
  args = sapply(args, strsplit, "=")
  for(arg in args){
    assign(arg[1], as.numeric(arg[2]), inherits=TRUE)
  }
}
parse_args(args = commandArgs(trailingOnly = TRUE))##, verbose=TRUE)


## Load libraries
library(ggplot2)
library(tidyverse)
library(tidyr) ## for pivot_wider() ?
library(devtools)
library(glmnet)
library(parallel)
library(lubridate)
source("/scratch/sangwonh/repos/covid-19-hotspots/helpers.r")
## source("~/repos/covid-19-hotspots/helpers.r")
load_all("/scratch/sangwonh/repos/covidcast/R-packages/covidcast")
## load_all("~/repos/covidcast/R-packages/covidcast")

## See how many cores this machine has.
mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                           ",")[[1]]
mycommand = "sinfo -o '%40N %c' --Node | grep "
## mycommand = "sinfo -o '%40N %c' --Node --long | grep "
mystring = system(paste0(mycommand, mynodename), intern=TRUE)
mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])/3


#######################
## "read-data" block ##
#######################

if(FALSE){
## Setup
lags = 28
geo_type = "state"   ## geo_type = "county"
response = "confirmed_7dav_incidence_prop"

## Read in data once
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
end_day = as.Date("2020-08-25")
validation_days = seq(end_day-30, end_day, by = "days") ## Old! to be retired soon

signals = data.frame(data_sources = data_sources, signals = signals)
suppressMessages({
  mat = covidcast_signals(signals,
                          start_day = start_day, end_day = end_day, geo_type = geo_type)
})
mat <- mat %>% select(geo_value, time_value, signal, data_source, value)

## New: instead of validation_days, use validation_geos
geos = mat %>% select(geo_value) %>% unlist() %>% unique()
set.seed(1000)
pct_validation = 0.3
validation_ind = sample(length(geos), length(geos) * pct_test)
validation_geos = geos[validation_ind]
validation_geos = c()
}

## save(list=ls(), file=file.path(outputdir, "data-county.Rdata")
## save(list=ls(), file=file.path(outputdir, "data-state.Rdata")
if(geo==1)load(file=file.path(outputdir, "data-county.Rdata"))
if(geo==2)load(file=file.path(outputdir, "data-state.Rdata"))
validation_geos = c() ## No validation set for now.
source("/scratch/sangwonh/repos/covid-19-hotspots/helpers.r")
## source("~/repos/covid-19-hotspots/helpers.r")

#####################
## "analyze" block ##
#####################
n_ahead = 28
threshold = 0.25
split_type = "geo"
onset = FALSE
slope = TRUE

## nn = 20
## n_ahead_list = round(seq(from=10, to=30, length=nn))
n_ahead_list = 10:30##round(seq(from=10, to=30, length=nn))
nn = length(n_ahead_list)
## nnlist = 1:(nn/2)
## nnlist = (nn/2+1):nn
nnlist = 1:nn
preload_df_model = FALSE
auc_list <- mclapply(nnlist, function(ii){

  n_ahead = n_ahead_list[ii]


  ## ## Load configuration
  ## configs[irow,] %>% list2env(globalenv())
  ## cat(## fn_response_name,
  ##     n_ahead, threshold, split_type, slope, onset, fill=TRUE) ## check that it's working
  fn_response <- match.fun("response_diff_avg_1week_min20")
  fn_response_name <- "response_diff_avg_1week_min20"

  ## ## Temporary
  ## slope = FALSE
  ## ## End of temporary

  ## Make a y|X model matrix ready to model
  if(preload_df_model){
    df_model <- ready_to_model(mat, lags, n_ahead, response, slope, fn_response, threshold, onset)
    save(df_model, file = file.path(outputdir, "v2-files", paste0(geo_type, "_n_ahead_", n_ahead, ".pdf")))
    return()
  }

  ## Instead of the above block, simply load the data.
  load_df_model <- function(n_ahead){
    load(file = file.path(outputdir, "v2-files", paste0(geo_type, "_n_ahead_", n_ahead, ".pdf")))
  }
  load_df_model(n_ahead)

  ## Split by time or geo
  if(split_type == "time"){
    df_traintest <- df_model %>% filter(!(time_value %in% validation_days))
    df_validation <- df_model %>% filter(time_value %in% validation_days)
    splitted <- sample_split_date(df_traintest, pct_test = 0.3)
  } else {
    df_traintest <- df_model %>% filter(!(geo_value %in% validation_geos))
    df_validation <- df_model %>% filter(geo_value %in% validation_geos)
    splitted <- sample_split_geo(df_traintest, pct_test = 0.3)

    ## Temporary check: show how many 1's exist
    df_model %>% select(resp) %>% table()
    df_traintest %>% select(resp) %>% table()
    splitted$df_train %>% select(resp) %>% table()
    splitted$df_test %>% select(resp) %>% table()
    df_validation %>% select(resp) %>% table()
    nfold = 5
    foldid <- make_foldid(splitted$df_train, nfold)
    for(ifold in 1:nfold){
      splitted$df_train[which(foldid==ifold),] %>% select(resp) %>% table() %>% print()
    }
  }
  ## source("~/repos/covid-19-hotspots/helpers.r")
  auc_df = calc_auc(destin = outputdir, splitted, lags, n_ahead, geo_type,
                    fn_response_name, threshold, slope, split_type, onset)
  save(auc_df, file = file.path(outputdir, paste0(geo_type, "_auc_n_ahead", n_ahead, ".Rdata")))
  return(auc_df)
}, mc.cores = mc.cores)##pmax(length(nnlist), mc.cores))

names(auc_list) = n_ahead_list
save(auc_list, file = file.path(outputdir, paste0(geo_type, "_auc_all_n_ahead", ".Rdata")))



## Plotting
if(FALSE){

  ## State
  geo_type = "state"
  load(file = file.path(outputdir, paste0(geo_type, "_auc_all_n_ahead", ".Rdata")))

  ## County
  geo_type = "county"
  n_ahead_list = 10:30##round(seq(from=10, to=30, length=nn))
  auc_list = list()
  for(ii in 1:length(n_ahead_list)){
    n_ahead = n_ahead_list[ii]
    load(file = file.path(outputdir, paste0(geo_type, "_auc_n_ahead", n_ahead, ".Rdata")))
    auc_list[[ii]] = auc_df
  }
  par(mfrow=c(1,5))

  ## Plot the difference
  n_ahead_list = 10:30##round(seq(from=10, to=30, length=nn))
  ## n_ahead_list = 10:19##round(seq(from=10, to=30, length=nn))
  auc_mat = do.call(rbind, lapply(auc_list, function(a){ unlist(a[,3] - a[,2])}))
  auc_mat %>% matplot(x=n_ahead_list, lwd=c(2,1,2,2), lty=c(1,2,1,1), type='l', xlab="n_ahead")
    abline(h=seq(from=0.1,to=1, by=0.1), col='grey80', lty=2)
  abline(h=0)
  legend("topleft", col=1:4, lty=1, legend=c("lasso", "ridge", "svm", "xgb"))

  ## Plot the two lines superimposed
  yes_fb = do.call(rbind, lapply(auc_list, function(a){ a[,3] %>% unlist()}))

  ## %>%
  ##   matplot(x=n_ahead_list, lwd=c(2,1,2,2)*3, lty=c(1,2,1,1), type='l', xlab="n_ahead", ylim=c(0.5,1))
  ## abline(h=0)
  ## title(main="yes fb")

  no_fb = do.call(rbind, lapply(auc_list, function(a){ a[,2] %>% unlist() }))
  plot_names = c("lasso", "ridge", "svm", "xgb")
  for(ii in 1:4){
    matplot(cbind(yes_fb[,ii], no_fb[,ii]), col=ii, lwd=3, lty=c(1,2), ylim=c(0.5,1), type='l')
    legend("topright", col=ii, lwd=2, lty=c(1,2), legend=c("With FB", "Without FB"))
    title(main=plot_names[ii])
    abline(h=seq(from=0.1,to=1, by=0.1), col='grey80', lty=2)
  }
  ## %>%
  ##   matlines(x=n_ahead_list, lwd=c(2,1,2,2), lty=c(1,2,1,1), type='l', xlab="n_ahead")

}
