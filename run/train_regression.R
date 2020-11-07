script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)

library(progress)
library(doParallel)

library(caret)
library(randomForest)
library(pROC)

source("run/tools.R")

n_cores <- 3
registerDoParallel(n_cores)
seed <- 7
set.seed(seed)
output_path <- "data/output/" 
n_folds <- 10


win_len <- 500000
win_len_upper <- win_len * 2

# load data
data <- read.csv("data/dataset.csv")
data$X <- NULL

dens_cols <- grep(x = names(data), pattern = "dens", value = TRUE)
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
tissue_spec_feats <- setdiff(names(data), c("chr","from","to", dens_cols, conserved_features))

ss_cols <- c ("chr", "from", "to")
features_cols <- c(conserved_features, tissue_spec_feats)

# data <- get_higher_level_features(data=data, features_cols = features_cols, 
#                                       win_len_upper = win_len_upper,
#                                       path_to_upper_data = paste0(data_path, "dataset_",
#                                                                   format(win_len_upper, scientific = FALSE), ".csv"))
# 
# high_features <- setdiff(names(data), c(features_cols, ss_cols, dens_cols))
# features_cols <- c(features_cols, high_features)

all_tissue_spec_feats <- vector()
for (feat in features_cols){
  for (col in tissue_spec_feats){
    all_tissue_spec_feats <- c(all_tissue_spec_feats, grep(x = feat, pattern = col, value = TRUE))  
  }
}
all_conserved_feats <- setdiff(features_cols, all_tissue_spec_feats)

# for storing results
df_r2_all <- data.frame()


# set progress bar
n_iterations <- length(dens_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in dens_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][2]
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats,
                                                   pattern = cancer_type, value = TRUE))
  all_resamples <- get_train_test_split(data, target_col=target_column, 
                                        start_pos_col="from", chr_col="chr", 
                                        feature_cols=all_features_cols, 
                                        n_folds=n_folds, seed=seed)
  # models parameters
  trCtrl <- trainControl(
    method="none",
    verboseIter = TRUE,
    seeds = seq(1, n_folds),
    allowParallel = TRUE
  )
  # set "seeds" parameters to make results reproducible
  if (length(all_features_cols) < 5){
    mtryGrid <- expand.grid(mtry = length(all_features_cols))
  } else {
    mtryGrid <- expand.grid(mtry = 5)
  }
  
  repeats_res <- foreach(i=seq(1, n_folds), .combine=rbind) %dopar% {
    
    # train/test split
    splitted_dataset <- all_resamples[[i]]
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  
    
    # fit models
    model <- train(
      x = x_train, 
      y = y_train,
      preProcess = c("center", "scale"),
      method = "rf",
      trControl = trCtrl,
      tuneGrid = mtryGrid,
      tuneLength = 1
    )
    
    # prediction
    train_pred <- data.frame(pred=predict(model, newdata = x_train[all_features_cols]))
    test_pred <- data.frame(pred=predict(model, newdata = x_test[all_features_cols]))
    train_pred$target <- y_train
    test_pred$target <- y_test
    
    # model quality
    model_qual <- get_model_quality(train_pred, test_pred)
    return(model_qual)
  }
  
  # add general info
  repeats_res <- repeats_res %>%
      mutate(
        iter = seq(1, n_folds),
        cancer_type = cancer_type
      )
  df_r2_all <- rbind(df_r2_all, repeats_res)
  write.csv(df_r2_all, file = paste0(output_path, "regression_restricted.csv"))
  pb$tick()
}
