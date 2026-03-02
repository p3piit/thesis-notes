if (sub_id) {
  cluster_sub <- count(data, id) %>%
    filter(n > threshold_subset)
  data_model <- data_model %>%
    filter(id %in% cluster_sub$id)
}

data_split <- mlml::split_gmert_data(data_model, 
                                    split_by_cluster = split_by_cluster)

fit_gmert <- mlml::fit_gmert_small(data_split$train,
  id = args_gmert$id,
  target = args_gmert$target,
  random_effect = args_gmert$random_effect,
  max_iter_inn = args_gmert$max_iter_inn,
  max_iter_out = args_gmert$max_iter_out,
  tol = args_gmert$tol,
  cp = args_gmert$cp,
  minsplit = args_gmert$minsplit,
  minbucket = args_gmert$minbucket,
  maxdepth = args_gmert$maxdepth,
  xval = args_gmert$xval,
  sanity_check = args_gmert$sanity_check)

fit_gmerf <- mlml::fit_gmerf_small(
  data_split$train,
  id = args_gmerf$id,
  target = args_gmerf$target,
  random_effect = args_gmerf$random_effect,
  max_iter_inn = args_gmerf$max_iter_inn,
  max_iter_out = args_gmerf$max_iter_out,
  tol = args_gmerf$tol,
  ntrees = args_gmerf$ntrees,
  mtry = args_gmerf$mtry,
  min_node_size = args_gmerf$min_node_size,
  max.depth = args_gmerf$max.depth,
  seed = args_gmerf$seed,
  num.threads = args_gmerf$num.threads,
  sanity_check = args_gmerf$sanity_check
)

data_train_noid <- data_split$train %>%
  select(-id)

form <- as.formula(
  paste(target, "~ .")
)

fit_cart <- rpart::rpart(form,
                         data = data_train_noid,
                         method = "class",
                         control = ctrl)

fit_rf <- ranger::ranger(
  formula = form,
  data = data_train_noid,
  num.trees = args_rn$ntrees,
  mtry = args_rn$mtry,
  min.node.size = args_rn$min_node_size,
  max.depth = args_rn$max.depth,
  seed = args_rn$seed,
  num.threads = args_rn$num.threads,
  probability = FALSE
)

pred_gmert <- mlml::predict_gmert(fit_gmert,
                                  data_split$test,
                                  random_effect = random_effect)

pred_gmerf <- mlml::predict_gmerf(fit_gmerf,
                                  data_split$test,
                                  random_effect = random_effect)

pred_cart <- as.numeric(as.vector(predict(fit_cart,
                    data_split$test,
                    type = "class")))

pred_rf <- ifelse(predict(fit_rf,
                   data_split$test)$predictions > .5,
                    1, 0)
  
  

## Confusion matrices

actual_obs <- unlist(data_split$test[target])

cm_gmert <- table(predicted = pred_gmert, actual = actual_obs)
cm_gmerf <- table(predicted = pred_gmerf, actual = actual_obs)
cm_cart  <- table(predicted = pred_cart, actual = actual_obs)
cm_rf    <- table(predicted = pred_rf, actual = actual_obs)

## Accuracy

acc_gmert <- mean(pred_gmert == actual_obs)
acc_gmerf <- mean(pred_gmerf == actual_obs)
acc_cart  <- mean(pred_cart == actual_obs)
acc_rf    <- mean(pred_rf == actual_obs)

## F1

f1_gmert_maj <- mlml::f1_fun(cm_gmert)
f1_gmerf_maj <- mlml::f1_fun(cm_gmerf)
f1_cart_maj  <- mlml::f1_fun(cm_cart)
f1_rf_maj    <- mlml::f1_fun(cm_rf)

f1_gmert_min <- mlml::f1_fun(cm_gmert, majority = FALSE)
f1_gmerf_min <- mlml::f1_fun(cm_gmerf, majority = FALSE)
f1_cart_min  <- mlml::f1_fun(cm_cart, majority = FALSE)
f1_rf_min    <- mlml::f1_fun(cm_rf, majority = FALSE)

## Bias

  #Bias
  bias_gmert <- sum(pred_gmert)/sum(actual_obs) - 1
  bias_gmerf <- sum(pred_gmerf)/sum(actual_obs) - 1
  bias_cart  <- sum(pred_cart)/sum(actual_obs) - 1
  bias_rf    <- sum(pred_rf)/sum(actual_obs) - 1
  

summary <- data.frame(
  model   = c("GMERT", "GMERF", "CART", "RF"),
  acc = c(acc_gmert, acc_gmerf, acc_cart, acc_rf),
  f1_mag  = c(f1_gmert_maj, f1_gmerf_maj, f1_cart_maj, f1_rf_maj),
  f1_min  = c(f1_gmert_min, f1_gmerf_min, f1_cart_min, f1_rf_min),
  bias = c( bias_gmert, bias_gmerf, bias_cart, bias_rf )
  )




print(summary)


