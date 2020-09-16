setwd("C:/Users/dzhao/Documents/Carnegie Mellon 2018-2023/RESEARCH - Delphi/covid-19-hotspots")
load("testdata_mat.RData")
load("testdata_df_model.RData")
load("testdata_splitted.RData")
load("testdata_fit_lasso.Rdata")
load("testdata_metrics.Rdata")
load("testdata_real_labels.Rdata")
source("helpers.r")
library(pROC)
library(testthat)

lags = 28
lags_val = 5
n_ahead = 21 ## 28
threshold = 0.25
geo_type = "state" ## or "county"
response = "confirmed_7dav_incidence_prop"
fn_response = response_diff_avg_1week_min20
fn_response_name = "response_diff_avg_1week_min20"
slope = TRUE
onset = FALSE
split_type = "geo"
data_sources = c("indicator-combination", 
                 "fb-survey")
signals = c("confirmed_7dav_incidence_prop", 
            "smoothed_hh_cmnty_cli")
start_day = as.Date("2020-05-01")
end_day = as.Date("2020-08-30")


test_that("check raw df", {
  expect_equal(names(mat), c("geo_value", "time_value", "signal", "data_source", "value"))
  expect_equal(n_distinct(mat$geo_value), 52)
  expect_equal(min(mat$time_value), start_day)
  expect_equal(max(mat$time_value), end_day)
  expect_equal(any(is.na(mat)), FALSE)
  expect_true(max(mat$value) < 1e6)
  expect_true(min(mat$value) > -1e6)
})


df_train = splitted$df_train
df_test = splitted$df_test
test_that("check ready to model df", {
  expect_true(as.logical(prod(c("geo_value", "time_value", "resp") %in% names(df_model))))
  expect_equal(n_distinct(df_model$geo_value), 52)
  expect_equal(min(df_model$time_value), start_day + lags_val)
  expect_true(max(df_model$time_value) <= end_day)
  expect_equal(any(is.na(df_model)), FALSE)
  expect_equal(max(df_model$resp), 1)
  expect_equal(min(df_model$resp), 0)
  expect_equal(dim(df_train)[1] + dim(df_test)[1], dim(df_model)[1])
  expect_equal(dim(df_train)[2], dim(df_test)[2], dim(df_model)[2])
})


library(glmnet)
preds = predict(fit_lasso, s = "lambda.min",
                newx = as.matrix(df_test %>% select(-geo_value, -time_value, -resp)),
                type = "response")[,1]
test_that("check lasso results", {
  expect_true(max(preds) <= 1)
  expect_true(min(preds) >= 0)
  expect_equal(length(preds), dim(df_test)[1])
  expect_equal(any(is.na(preds)), FALSE)
  expect_true(all(diff(fit_lasso$lambda) <= 0))
})


# NOTE: AUC is the probability that a true 1 will have a higher predicted probability than a true 0
combinations <- expand.grid(positiveProbs=preds[real_labels == 1L], 
                            negativeProbs=preds[real_labels == 0L])
test_that("check ROC and AUC", {
  expect_true(all(diff(metrics$fpr) <= 0))
  expect_true(all(diff(metrics$fnr) <= 0))
  expect_true(max(metrics$fpr) <= 1)
  expect_true(max(metrics$fnr) <= 1)
  expect_true(min(metrics$fpr) >= 0)
  expect_true(min(metrics$fnr) >= 0)
  expect_equal(mean(combinations$positiveProbs > combinations$negativeProbs),
               as.numeric(auc(real_labels, preds)))
})

