scramble_some_training_columns <- function(nn, df_train, df_test, foldid){

  ## ## One diagnostic: take some columns, permute the values (killing the signal in that column)
  nr = df_train %>% nrow()
  nc = df_train %>% ncol()
  if(nn > nc) stop(paste("Can't scramble more than", nc, "columns!"))
  icols = sample(3:(nc-1), nn, replace=FALSE)
  for(icol in icols){
    df_train[,icol] = df_train[sample(1:nr, replace=FALSE), icol]
  }

  ## Main part of the lasso fitting and predicting
  fit_lasso <- cv.glmnet(x = as.matrix(df_train %>% select(-geo_value, -time_value, -resp)),
                         y = df_train$resp,
                         family = "binomial",
                         alpha = 1,
                         foldid = foldid,
                         nfold = nfold)

  ## ## See the fitted model
  ## coef(fit_lasso, s = "lambda.min") %>% print()

  ## Obtain the predictions
  preds = predict(fit_lasso, s = "lambda.min",
                  newx = as.matrix(df_test %>% select(-geo_value, -time_value, -resp)),
                  type = "response")[,1]
  suppressMessages({
    auc = pROC::auc(response = df_test$resp, predictor = preds)[1]
  })
  return(auc)
}
