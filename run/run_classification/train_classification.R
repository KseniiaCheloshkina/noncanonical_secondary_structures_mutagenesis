script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)
library(openxlsx)
library(wesanderson)

library(progress)
library(doParallel)

library(caret)
library(randomForest)
library(pROC)
library(PRROC)

source("tools.R")

n_cores <- 3
registerDoParallel(n_cores)

set.seed(7)

win_len <- 500000
win_len_upper <- win_len * 2

output_path <- "../data/output/classification_restricted/" 

# load data
data <- read.csv("../data/hsp_dataset.csv")
data$X <- NULL

hsp_cols <- grep(x = names(data), pattern = "hsp", value = TRUE)
conserved_features <- c(
  "A_Phased_Repeat", "Direct_Repeat", 
  "Inverted_Repeat", "Mirror_Repeat", "Short_Tandem_Repeat", 
  "Z_DNA_Motif", "G_quadruplex", 
  "stemloops_16_50", "stemloops_6_15", 
  "X3UTR", "X5UTR", "codingExons", "downstream", "introns", "promoters", "WholeGenes", 
  "tad_boundaries_liver", "tad_boundaries_ovary", "tad_boundaries_pancreatic",
  
  "cancer_liver_ATF3.human", "cancer_liver_CTCF.human", "cancer_liver_EGR1.human",
  "cancer_liver_FOXA1.human", "cancer_liver_FOXA2.human", "cancer_liver_GABPA.human",
  "cancer_liver_HNF4A.human", "cancer_liver_HNF4G.human", "cancer_liver_JUND.human",
  "cancer_liver_MAX.human", "cancer_liver_NR2F2.human", "cancer_liver_REST.human",
  "cancer_liver_RXRA.human", "cancer_liver_SP1.human", "cancer_liver_YY1.human",
  "cancer_liver_ZBTB33.human"
)
conserved_features <- intersect(conserved_features, names(data))
tissue_spec_feats <- setdiff(names(data), c("chr","from","to", hsp_cols, conserved_features))

ss_cols <- c ("chr", "from", "to")
features_cols <- c(conserved_features, tissue_spec_feats)

all_data <- data

all_tissue_spec_feats <- vector()
for (feat in features_cols){
  for (col in tissue_spec_feats){
    all_tissue_spec_feats <- c(all_tissue_spec_feats, grep(x = feat, pattern = col, value = TRUE))  
  }
}
all_conserved_feats <- setdiff(features_cols, all_tissue_spec_feats)

# train/test split
train_ratio <- 0.7
n_repeats <- 30

# for storing results
df_roc_auc_all <- data.frame()
df_imp_all <- data.frame()
df_recall_all <- data.frame()

# select quantiles for quality assessment
recall_quantiles <- c(0.001, 0.005, 0.002,  0.003, 0.015, 0.025, seq(0.01, 0.05, 0.01), seq(0.1, 0.9, 0.05))

# set progress bar
n_iterations <- length(hsp_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in hsp_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats,
                                                   pattern = cancer_type, value = TRUE))
  
  # models parameters
  trCtrl <- trainControl(
    method="none",
    verboseIter = TRUE,
    classProbs=TRUE,
    summaryFunction = twoClassSummary,
    seeds = seq(1, n_repeats),
    allowParallel = TRUE
  )
  # set "seeds" parameters to make results reproducible
  if (length(all_features_cols) < 5){
    mtryGrid <- expand.grid(mtry = length(all_features_cols))
  } else {
    mtryGrid <- expand.grid(mtry = 5)
  }
  
  repeats_res <- foreach(i=seq(1, n_repeats)) %dopar% {
    
    # train/test split
    splitted_dataset <- get_train_test_split_classification(data=all_data, target_col=target_column, 
                                             start_pos_col="from",
                                             chr_col="chr", feature_cols=all_features_cols,
                                             train_ratio=train_ratio, 
                                             seed=i)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  
    
    n_pos <- length(y_train[y_train == "X1"])
    n_neg <- length(y_train[y_train == "X0"])
    
    # fit models
    model <- rf_fit(x_train = x_train[all_features_cols], y_train = y_train, trCtrl = trCtrl,
                    n_pos = n_pos, n_neg = n_neg, mtryGrid = mtryGrid)
    
    # prediction
    train_pred <- predict(model, newdata = x_train[all_features_cols], type = "prob")
    test_pred <- predict(model, newdata = x_test[all_features_cols], type = "prob")
    train_pred$target <- y_train
    test_pred$target <- y_test
    
    # model quality
    model_qual <- get_model_quality_classification(train_pred, test_pred, model, recall_quantiles)
    return(model_qual)
  }
  
  
  # save results
  for (split_iter in 1:length(repeats_res)){
    
    res_iter <- repeats_res[[split_iter]]
    
    # extract specific datasets
    df_roc_auc <- res_iter[['roc_auc']]
    df_imp <- res_iter[['importance']]
    df_recall <- res_iter[['recall']]
    
    # add general info
    df_roc_auc <- df_roc_auc %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        win_len = format(win_len, scientific = F)
      )
    
    df_imp <- df_imp %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        win_len = format(win_len, scientific = F)
      )
    
    df_recall <- df_recall %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        win_len = format(win_len, scientific = F)
      )
    
    df_roc_auc_all <- rbind(df_roc_auc_all, df_roc_auc)
    df_imp_all <- rbind(df_imp_all, df_imp)
    df_recall_all <- rbind(df_recall_all, df_recall)
  }
  
  write.csv(df_roc_auc_all, file = paste0(output_path, "result_roc_auc.csv"))
  write.csv(df_imp_all, file = paste0(output_path, "result_imp.csv"))
  write.csv(df_recall_all, file = paste0(output_path, "result_recall.csv"))
  pb$tick()
}
