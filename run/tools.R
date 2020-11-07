library(dplyr)
library(ggplot2)
library(caret)
library(pROC)


############################## FUNCTION TO TRANSFORM TARGET
## INPUT: 
# y - target values
## RETURNS:
# y - transformed target values
## Returns transformed target values
##############################
density_transform <- function(y, mean_dens=NA){
  if (is.na(mean_dens)){
    mean_dens <- mean(y)
    }
  return(list(mean_dens, c(log2(mean_dens + y))))
}



############################## FUNCTION TO GET STRATIFIED 10-FOLD CV
## INPUT: 
# data - original dataframe 
# target_col - column name of target
# start_pos_col- column name of start position of window
# chr_col - column name of chromosome
# feature_cols - vector of features names 
# train_ratio - ratio of data to place at train dataset
# seed - random seed for reproducibility
## RETURNS: all_data - list of: 
# x_train - train features dataset 
# y_train - train target (factor)
# x_test - test features dataset
# y_test - test target (factor)
## Returns data split stratified by target, chromosome and position in chromosome (divided on 4 parts)
##############################
get_train_test_split <- function(data, target_col, start_pos_col, chr_col, feature_cols, n_folds, seed){
  
  chrs <- c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22"
  )
  data[chr_col] <- lapply(data[chr_col], as.character)
  
  # get position bins
  all_pos_map <- data.frame()
  
  for (chr in chrs){
    starts <- data[data[chr_col] == chr, start_pos_col]
    starts <- starts[order(starts)]
    bin_name <- cut(starts, breaks = 4, labels = c("first_bin", "second_bin", "third_bin", "last_bin"),
                    include.lowest = TRUE, right = FALSE)
    pos_map <- data.frame(start = starts, pos_bin = bin_name, chr = chr)
    names(pos_map) <- c(start_pos_col, "pos_bin", chr_col)
    all_pos_map <- rbind(all_pos_map, pos_map)
  }
  all_pos_map[chr_col] <- lapply(all_pos_map[chr_col], as.character)
  data <- data %>%
    inner_join(all_pos_map, by = c(chr_col, start_pos_col))
  
  # create strata column by chromosome and quantile of position in chromosome
  data['strata'] <- paste0(data[, chr_col], "_", data[, 'pos_bin'])
  rownames(data) <- NULL

  # split by strata column
  set.seed(seed)
  ind <- createFolds(y=data$strata, k=n_folds, list = TRUE)
  all_data_resamples <- list()
  for (fold in names(ind)){
    data_train <- data[-ind[[fold]], ]
    data_test  <- data[ind[[fold]], ]
    
    x_train <- data_train[feature_cols]
    x_test <- data_test[feature_cols]
    
    # transform y
    ans <- density_transform(data_train[, target_col], mean_dens=NA)
    mean_dens <- ans[[1]]
    y_train <- ans[[2]]
    ans <- density_transform(data_test[, target_col], mean_dens=mean_dens)
    y_test <- ans[[2]]
    
    all_data <- list()
    all_data[[1]] <- x_train
    all_data[[2]] <- y_train
    all_data[[3]] <- x_test
    all_data[[4]] <- y_test    
    names(all_data) <- c("x_train", "y_train", "x_test", "y_test")
    all_data_resamples[[fold]] <- all_data
  }
  return(all_data_resamples)
}



############################## FUNCTION TO GET QUALITY METRICS FOR TRAIN AND TEST
## INPUT: 
# train_pred - dataframe with columns "X1"(predicted probability of class "1") and "target" 
# with levels "X0" and "X1" (class "1")
# test_pred - dataframe with columns "X1"(predicted probability of class "1") and "target" 
# with levels "X0" and "X1" (class "1")
# model - caret model for which varImp function is applicable
## RETURNS: results - list of: 
# roc_auc - ROC AUC on train and test 
# importance - feature importance from model
# recall - recall data for different probability quantiles 
##############################
get_model_quality <- function(train_pred, test_pred){

  res_train <- postResample(pred=train_pred$pred, obs=train_pred$target)
  res_test <- postResample(pred=test_pred$pred, obs=test_pred$target)
  out_df <- data.frame(
    res_train[1],
    res_train[2],
    res_test[1],
    res_test[2]
  )
  names(out_df) <- c("rmse_train", "r_squared_train", "rmse_test", "r_squared_test")
  rownames(out_df) <- NULL
  return(out_df)
}


######### publication plots theme
theme_light_custom <- function (scale_fill=TRUE) {
  base_size <- 9
  base_family <- "sans"
  base_line_size <- base_size / 22
  base_rect_size <- base_size / 22
  half_line <- base_size / 2
  th <- theme_light(base_size, base_family, base_line_size, base_rect_size) %+replace%
    theme(
      text = element_text(family = base_family, face = "plain",
                          colour = "black", size = base_size,
                          lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
                          margin = margin(), debug = FALSE, inherit.blank = TRUE
      ),
      axis.title.x = element_text(margin = margin(t = 2.75),
                                  vjust = 1, size = 10, inherit.blank = T),
      axis.title.y = element_text(angle = 90, margin = margin(r = 2.75),
                                  vjust = 1, size = 10),
      legend.position = "bottom",
      legend.title = element_text(hjust = 0, size = 10),
    )
  if (scale_fill) {
    return(list(th, scale_fill_jco()))
  }
  else return(th)
  
}


############################## FUNCTION TO GET HIGHER-LEVEL OF AGGREGATION FEATURES
## INPUT: 
# data - original dataframe
# win_len_upper - level of aggregation for which to take new (higher-level) features
# path_to_upper_data - full path to datafrom higher level
# features_cols - vector of all columns names which to take from higher-level data

## RETURNS: data - dataframe: initial dataframe + new features
##############################
get_higher_level_features <- function(data, win_len_upper, path_to_upper_data, features_cols){
  
  data_upper <- read.csv(path_to_upper_data, stringsAsFactors = F)
  data_upper$to_window <- ceiling(data_upper$to / win_len_upper)
  
  data_upper <- data_upper %>%
    select(c("chr", "to_window", features_cols))
  
  nm_feats <- names(data_upper)[3:length(names(data_upper))]
  nm_feats <- paste0("upper_", nm_feats)
  names(data_upper)[3:length(names(data_upper))] <- nm_feats
  
  data$to_window <- ceiling(data$to / win_len_upper)
  new_features <- data %>%
    left_join(data_upper, by=c("chr", "to_window")) %>%
    select(-to_window)
  
  for (col in nm_feats){
    
    new_features <- new_features %>% 
      mutate(!!as.name(col) := if_else(is.na(!!as.name(col)), 0, !!as.name(col)))    
  }
  
  return(new_features)
}



############################## FUNCTION TO GET STRATIFIED TRAIN/TEST SPLIT
## INPUT: 
# data - original dataframe 
# target_col - column name of target
# start_pos_col- column name of start position of window
# chr_col - column name of chromosome
# feature_cols - vector of features names 
# train_ratio - ratio of data to place at train dataset
# seed - random seed for reproducibility
## RETURNS: all_data - list of: 
# x_train - train features dataset 
# y_train - train target (factor)
# x_test - test features dataset
# y_test - test target (factor)
## Returns data split stratified by target, chromosome and position in chromosome (divided on 4 parts)
##############################
get_train_test_split_classification <- function(data, target_col, start_pos_col, chr_col, feature_cols,
                                 train_ratio, seed, dens_col=NULL){
  
  chrs <- c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
  data[chr_col] <- lapply(data[chr_col], as.character)
  
  # get position bins
  all_pos_map <- data.frame()
  
  for (chr in chrs){
    starts <- data[data[chr_col] == chr, start_pos_col]
    starts <- starts[order(starts)]
    bin_name <- cut(starts, breaks = 4, labels = c("first_bin", "second_bin", "third_bin", "last_bin"),
                    include.lowest = TRUE, right = FALSE)
    pos_map <- data.frame(start = starts, pos_bin = bin_name, chr = chr)
    names(pos_map) <- c(start_pos_col, "pos_bin", chr_col)
    all_pos_map <- rbind(all_pos_map, pos_map)
  }
  all_pos_map[chr_col] <- lapply(all_pos_map[chr_col], as.character)
  data <- data %>%
    inner_join(all_pos_map, by = c(chr_col, start_pos_col))
  
  # create strata column by chromosome and quantile of position in chromosome
  data['strata'] <- paste0(data[chr_col], "_", data['pos_bin'])
  
  # split data on examples with positive and negative target
  rownames(data) <- NULL
  pos_label_data <- data[data[target_col] == 1, ]
  neg_label_data <- data[data[target_col] == 0, ]
  
  # split in train_ratio proportion stratified by strata column
  set.seed(seed)
  tr_neg <- createDataPartition(neg_label_data$strata, p = train_ratio, list = FALSE)
  data_train <- neg_label_data[tr_neg,]
  data_test  <- neg_label_data[-tr_neg,]
  
  set.seed(seed)
  tr_pos <- createDataPartition(pos_label_data$strata, p = train_ratio, list = FALSE)
  data_train <- rbind(data_train, pos_label_data[tr_pos, ])
  data_test <- rbind(data_test, pos_label_data[-tr_pos, ])
  
  x_train <- data_train[feature_cols]
  x_test <- data_test[feature_cols]
  
  y_train <- as.factor(data_train[, target_col])
  levels(y_train) <- c("X0", "X1")
  y_test <- as.factor(data_test[, target_col])
  levels(y_test) <- c("X0", "X1")
  
  all_data <- list()
  all_data[[1]] <- x_train
  all_data[[2]] <- y_train
  all_data[[3]] <- x_test
  all_data[[4]] <- y_test
  
  if (is.null(dens_col)) {
    names(all_data) <- c("x_train", "y_train", "x_test", "y_test")
  } else {
    train_dens <- data_train[dens_col]
    test_dens <- data_test[dens_col]
    all_data[[5]] <- train_dens
    all_data[[6]] <- test_dens    
    names(all_data) <- c("x_train", "y_train", "x_test", "y_test", "density_train", "density_test")
  }
  
  return(all_data)
}



############################## FUNCTION TO GET QUALITY METRICS FOR TRAIN AND TEST
## INPUT: 
# train_pred - dataframe with columns "X1"(predicted probability of class "1") and "target" 
# with levels "X0" and "X1" (class "1")
# test_pred - dataframe with columns "X1"(predicted probability of class "1") and "target" 
# with levels "X0" and "X1" (class "1")
# model - caret model for which varImp function is applicable
## RETURNS: results - list of: 
# roc_auc - ROC AUC on train and test 
# importance - feature importance from model
# recall - recall data for different probability quantiles 
##############################
get_model_quality_classification <- function(train_pred, test_pred, model, quant=seq(0.01, 0.1, 0.01)){
  
  # ROC AUC
  tr_auc <- as.numeric(auc(cases = train_pred$X1[train_pred$target == "X1"], controls = train_pred$X1[train_pred$target == "X0"],
                           direction = "<"))
  te_auc <- as.numeric(auc(cases = test_pred$X1[test_pred$target == "X1"], controls = test_pred$X1[test_pred$target == "X0"],
                           direction = "<"))
  out_df <- data.frame(tr_auc = tr_auc, te_auc = te_auc)
  
  # importance
  imp <- varImp(model, scale = FALSE)$importance
  names(imp) <- "importance"
  imp$feature <- rownames(imp)  
  rownames(imp) <- NULL
  
  # recall
  test_pred <- test_pred[order(test_pred$X1, decreasing=TRUE), ]
  test_pred$target <- unlist(lapply(test_pred$target, as.character))
  test_pred$target <- gsub(pattern = "X", replacement = "", x = test_pred$target)
  test_pred$target <- as.numeric(test_pred$target)
  test_pos <- sum(test_pred$target)
  
  # pr auc
  test_pred$target_num <- as.numeric(gsub(x = test_pred$target, pattern = "X", replacement = ""))
  roc <- pr.curve(scores.class0 = test_pred$X1, weights.class0 = test_pred$target_num)
  out_df$pr_auc <- roc$auc.integral
  
  # efron r2
  mu <- mean(test_pred$target_num)
  out_df$efron_r2 <- 1 - sum((test_pred$target_num - test_pred$X1)^2)/sum((test_pred$target_num - mu)^2)
  
  df_metr <- data.frame(quantile = quant, n_selected = round(nrow(test_pred) * quant))
  df_metr1 <- apply(
    df_metr, 
    MARGIN = 1, 
    function(x) {
      data.frame(
        quantile = x['quantile'],
        n_selected = x['n_selected'],
        recall = sum(test_pred$target[1:x['n_selected']]) / test_pos,
        tp = sum(test_pred$target[1:x['n_selected']]),
        precision = sum(test_pred$target[1:x['n_selected']]) / x['n_selected']
      )
    }
  )
  
  ratio_data <- data.frame()
  for (df in df_metr1){
    ratio_data <- rbind(ratio_data, df)
  }
  rownames(ratio_data) <- NULL
  ratio_data$n_test_pos <- test_pos
  ratio_data$lift_recall <- ratio_data$recall / ratio_data$quantile
  ratio_data$f1 <- ratio_data$recall * ratio_data$precision * 2 / (ratio_data$recall + ratio_data$precision)
  ratio_data$f1[is.na(ratio_data$f1)] <- 0
  
  results <- list()
  results[['roc_auc']] <- out_df
  results[['importance']] <- imp
  results[['recall']] <- ratio_data
  
  return(results)
}

get_sample_size <- function(n_pos, n_neg){
  
  bal <- n_pos / n_neg
  if (bal < 0.2){
    # if balance < 20% then make it to be equal 20% by taking all positives 
    # and sampling negatives
    bal_vect <- c(X0 = floor(n_pos * 5), X1 = n_pos)
  } else if (bal > 1){
    # if there are more (target) positive examples than negative examples
    # then make balance equal to perfect 1 (same number of positives and negatives)
    bal_vect <- c(X0 = n_neg, X1 = n_neg)
  } else {
    # if balance is between 0.2 and 1 then consider it as a rather natural balance and 
    # take it as is
    bal_vect <- c(X0 = n_neg, X1 = n_pos)
  }
  return(bal_vect)
}



rf_fit <- function(x_train, y_train, trCtrl, n_pos, n_neg, mtryGrid){
  
  model <- train(
    x = x_train, 
    y = y_train,
    preProcess = c("zv", "YeoJohnson", "center", "scale"),
    method = "rf",
    metric="ROC",   
    maximize = TRUE,
    trControl = trCtrl,
    ntree = 500,
    nodesize = 30,
    maxnodes = 3,
    sampsize = get_sample_size(n_pos = n_pos, n_neg = n_neg),
    tuneGrid = mtryGrid
  )
  
  return(model)
}

