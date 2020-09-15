setwd("C:/Users/dzhao/Documents/Carnegie Mellon 2018-2023/RESEARCH - Delphi/covid-19-hotspots")
load("testdata_mat.RData")
load("testdata_df_model.RData")
source("helpers.r")
library(testthat)

lags <- 28
geo_type <- "state"  # "county"
response <- "confirmed_7dav_incidence_prop"
fn_response <- response_diff_avg_1week_min20
fn_response_name <- "response_diff_avg_1week_min20"
start_day <- as.Date("2020-05-01")
end_day <- as.Date("2020-09-04")


test_that("check raw df", {
  expect_equal(names(mat), c("geo_value", "time_value", "signal", "data_source", "value"))
  expect_equal(n_distinct(mat$geo_value), 52)
  expect_equal(min(mat$time_value), start_day)
  expect_equal(max(mat$time_value), end_day)
  expect_equal(any(is.na(mat)), FALSE)
  expect_true(max(mat$value) < 1e6)
  expect_true(min(mat$value) > -1e6)
})


test_that("check ready to model df", {
  expect_true(as.logical(prod(c("geo_value", "time_value", "resp") %in% names(df_model))))
  expect_equal(n_distinct(df_model$geo_value), 52)
  expect_equal(min(df_model$time_value), start_day + 5)
  expect_equal(any(is.na(df_model)), FALSE)
  expect_true(max(df_model$resp) == 1)
  expect_true(min(df_model$resp) == 0)
})


test_that("check lasso results", {
  
})


test_that("check response visualization", {
  
})


test_that("check ROC plots", {
  
})


test_that("check AUC plots", {
  
  # find one example of a model, compute AUC and check it by eye
  # can use pROC:auc
  
  # find one model, test its ROC separately (with code I wrote, see that answers match)
  
  # can further modularize it
  # 
  
})

