# A mirror of the script in v1.Rmd
outputdir = "./figures"
mc.cores = 1

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


if(FALSE){
  for(geo in 1:2){
    print("geo")
    print(geo)

    ## Setup
    lags = 28
    ## n_ahead = 21 ## 28
    threshold = 0.25
    if(geo==1) geo_type = "county"## "state" ## or "county"
    if(geo==2) geo_type = "state"## "state" ## or "county"
    response = "confirmed_7dav_incidence_prop"
    fn_response = response_diff_avg_1week_min20
    fn_response_name = "response_diff_avg_1week_min20"
    slope = TRUE
    onset = FALSE
    split_type = "geo"

    ## Read in data once
    data_sources = c("indicator-combination",
                     "fb-survey")
    signals = c("confirmed_7dav_incidence_prop",
                "smoothed_hh_cmnty_cli")
    start_day = as.Date("2020-05-01")
    end_day = as.Date("2020-08-30")
    signals = data.frame(data_sources = data_sources, signals = signals)
    suppressMessages({
      mat = covidcast_signals(signals,
                              start_day = start_day, end_day = end_day, geo_type = geo_type)
    })
    mat <- mat %>% select(geo_value, time_value, signal, data_source, value)

    ## ## Form the y|X matrix
    ## df_model <- ready_to_model(mat, lags, n_ahead, response, slope, fn_response, threshold, onset)

    ## ## save(list=ls(), file=file.path(outputdir, "data-county.Rdata")
    if(geo==1) filename = "data-county.Rdata"
    if(geo==2) filename = "data-state.Rdata"
    save(lags, threshold, geo_type, response, fn_response, fn_response_name,
         slope, onset, split_type, mat, data_sources, signals, start_day, end_day, signals, mat,
         file=file.path(outputdir, filename))
    print("Saved to")
    print(file.path(outputdir, filename))
  }
}

if(geo==1)load(file=file.path(outputdir, "data-county.Rdata"))
if(geo==2)load(file=file.path(outputdir, "data-state.Rdata"))

n_ahead_list = 10:30
if(blocknum==1) n_ahead_list = 10:20
if(blocknum==2) n_ahead_list = 21:30
nn = length(n_ahead_list)
if(is.null(mc.cores)) mc.cores = 1
auc_list <- mclapply(1:nn, function(ii){

  n_ahead = n_ahead_list[ii]
  printprogress(n_ahead, 30, "auc_list", fill=TRUE)
  fn_response <- match.fun("response_diff_avg_1week_min20")
  fn_response_name <- "response_diff_avg_1week_min20"

  ## ## Make a y|X model matrix ready to model, and save it.
  df_model <- ready_to_model(mat, lags, n_ahead, response, slope, fn_response, threshold, onset)

  ## ## Get AUC for 5 different splits
  ## nsim = 5
  ## list_of_auc_df = lapply(1:nsim, function(isim){
  ##   splitted <- stratified_sample_split_geo(df_model, pct_test = 0.3)
  ##   auc_df = calc_auc(destin = outputdir, splitted, lags, n_ahead, geo_type,
  ##                     fn_response_name, threshold, slope, split_type, onset)
  ##   auc_df
  ## })

  ## NEW: Get AUC for 5 different splits, but cycling over five test data
  ## splits, instead of five random splits.
  set.seed(12345)
  obj <- outer_split(df_model, nfold)
  list_of_auc_df = lapply(1:(nfold+1), function(ii){
    splitted <- obj$splitted_list[[ii]
    auc_df = calc_auc(destin = outputdir, splitted, lags, n_ahead, geo_type,
                      fn_response_name, threshold, slope, split_type, onset)
    auc_df
  })


  ## Take an average
  avg_auc_df = Reduce("+", list_of_auc_df) / length(list_of_auc_df)
  avg_auc_df[["model"]] = list_of_auc_df[[1]][["model"]]
  auc_df = avg_auc_df
  save(list_of_auc_df, auc_df, file = file.path(outputdir, paste0(geo_type, "_auc_n_ahead",
                                                                  n_ahead, ".Rdata")))
  cat(fill=TRUE)
  return(auc_df)
}, mc.cores = mc.cores)
names(auc_list) = n_ahead_list
save(auc_list, n_ahead_list, file = file.path(outputdir, paste0(geo_type, "_auc_all_n_ahead", ".Rdata")))


## Plotting
if(FALSE){

  ## ## State
  ## geo_type = "state"
  ## load(file = file.path(outputdir, paste0(geo_type, "_auc_all_n_ahead", ".Rdata")))

  ## County
  for(geo_type in c("state", "county")){
  n_ahead_list = 10:30##round(seq(from=10, to=30, length=nn))
  auc_list = list()
  for(ii in 1:length(n_ahead_list)){
    n_ahead = n_ahead_list[ii]
    load(file = file.path(outputdir, paste0(geo_type, "_auc_n_ahead", n_ahead, ".Rdata")), verbose = TRUE)
    auc_list[[ii]] = auc_df
  }
  par(mfrow=c(1,4))
  par(oma = c(0,0,2,0))

  ## Plot the difference
  n_ahead_list = 10:30##round(seq(from=10, to=30, length=nn))
  ## n_ahead_list = 10:19##round(seq(from=10, to=30, length=nn))
  auc_mat = do.call(rbind, lapply(auc_list, function(a){ unlist(a[,3] - a[,2])}))
  auc_mat %>% matplot(x=n_ahead_list, lwd=c(2,1,2,2), lty=c(1,2,1,1), type='l', xlab="n_ahead")
    abline(h=seq(from=0.1,to=1, by=0.1), col='grey80', lty=2)
  abline(h=0, col='grey80', lwd=2)
  legend("topleft", col=1:4, lty=1, legend=c("lasso", "ridge", "svm", "xgb"))
  title(main="AUC(yes FB) - AUC(no FB)")

  ## Plot the two lines superimposed
  yes_fb = do.call(rbind, lapply(auc_list, function(a){ a[,3] %>% unlist()}))

  no_fb = do.call(rbind, lapply(auc_list, function(a){ a[,2] %>% unlist() }))
  plot_names = paste0("AUC, ", c("lasso", "ridge",  "xgb"))
  for(ii in 1:length(plot_names)){
    matplot(y=cbind(yes_fb[,ii], no_fb[,ii]), x=n_ahead_list, col=ii, lwd=3, lty=c(1,2), ylim=c(0.5,1), type='l',
            ylab = "AUC", xlab="n_ahead")
    abline(h=seq(from=0.1,to=1, by=0.1), col='grey80', lty=2)
    legend("topright", col=ii, lwd=2, lty=c(1,2), legend=c("With FB", "Without FB"), bg="white")
    title(main=plot_names[ii])
  }
  mtext(outer=TRUE, text=bquote(bold(.(toupper(geo_type)))), side=3)
  }

}
