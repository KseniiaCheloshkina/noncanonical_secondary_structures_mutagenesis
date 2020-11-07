script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)

output_folder <- "data/datasets/"


##### Read target

df_bkpt <- read.csv("data/target/final/all_data_5e+05.csv", stringsAsFactors = F)
df_bkpt <- df_bkpt[df_bkpt$chr != "X", ]
identity_cols <- c("chr", "from", "to")
density_cols <- grep(x = names(df_bkpt), pattern = "density", value = T)

# log transform with 1
for (col in density_cols){
  par(mfrow=c(1,2))
  hist(df_bkpt[, c(col)], main="original")
  hist(log2(1 + df_bkpt[, c(col)]), main="log-transformed")
}

# log transform with mean 
mean_dens <- mean(df_bkpt$density_BRCA)
for (col in density_cols){
  par(mfrow=c(1,2))
  hist(df_bkpt[, c(col)], main="original")
  hist(log2(mean_dens + df_bkpt[, c(col)]), main="log-transformed")
}

df_bkpt <- df_bkpt[, c(identity_cols, density_cols)]


##### Join features
features_path <- "data/features/"
all_features_paths <- list.files(features_path)
all_features_paths <- grep(x=all_features_paths, pattern = "500000", value=T)


# collect conserved features - non-B DNA motifs
tissue_spec_path <- c(
  "DNase_seq_500000.csv",
  "H3K4me1-human_500000.csv", "H3K4me3-human_500000.csv", "H3K9me3-human_500000.csv",
  "H3K27ac-human_500000.csv", "H3K27me3-human_500000.csv", "H3K36me3-human_500000.csv"
  )
all_full_conserved_features_path <- setdiff(all_features_paths, tissue_spec_path)

all_data <- df_bkpt

for (feat_path in all_full_conserved_features_path){
  
  feature_data <- read.csv(paste0(features_path, feat_path), stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_data <- all_data  %>%
    left_join(feature_data, by = c("chr", "from", "to"))
  
  feat_col <- names(feature_data)[ncol(feature_data)]
  all_data[is.na(all_data[feat_col]), feat_col] <- 0
}


# collect tisssue-specific features

# Missing from our study:
# "ESAD" - Esophagus
# "GACA" - Stomach 
# "RECA" - kidney

map_abb_to_cancer_type <- data.frame(
  abb = c("BRCA", "LIRI", "MALY", "OV", "PACA", "PBCA"),
  cancer_type = c("breast", "liver", "blood", "ovary", "pancreatic", "brain")
)

for (feat_path in tissue_spec_path){
  
  feature_data <- read.csv(paste0(features_path, feat_path), stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_cancer_types <- unique(feature_data$cancer_type)
  
  for (canc_type in all_cancer_types){
    
    feature_data_part <- feature_data[feature_data$cancer_type == canc_type, ]
    feature_data_part$cancer_type <- NULL
    feat_name <- names(feature_data_part)[ncol(feature_data_part)]
    new_canc_type <- as.character(map_abb_to_cancer_type[map_abb_to_cancer_type$cancer_type == canc_type,
                                                         "abb"])
    if (length(new_canc_type) != 0){
      new_feat_name <- paste0("cancer_", new_canc_type, "_", feat_name)
      feature_data_part <- feature_data_part %>%
        rename(!!new_feat_name := feat_name)
      
      all_data <- all_data  %>%
        left_join(feature_data_part, by=c("chr", "from", "to"))
      all_data[is.na(all_data[new_feat_name]), new_feat_name] <- 0        
    }
  }
}


# log transform
feature_names <- setdiff(names(all_data), c(identity_cols, density_cols))
cancer_spec_features <- grep(x = feature_names, pattern = "cancer", value = T)
conserved_feat <- setdiff(feature_names, cancer_spec_features)

# check log transform with 1
for (col in conserved_feat){
  print(col)
  par(mfrow=c(1,2))
  hist(all_data[, c(col)], main="original")
  hist(log2(1 + all_data[, c(col)]), main="log-transformed")
}

for (col in cancer_spec_features){
  print(col)
  par(mfrow=c(1,2))
  hist(all_data[, c(col)], main="original")
  hist(log2(1 + all_data[, c(col)]), main="log-transformed")
}

# check log transform with mean 
for (col in cancer_spec_features){
  print(col)
  mean_dens <- mean(all_data[, c(col)])
  par(mfrow=c(1,2))
  hist(all_data[, c(col)], main="original")
  hist(log2(mean_dens + all_data[, c(col)]), main="log-transformed")
}

# log transform with 1
for (col in feature_names){
  all_data[, c(col)] <- log2(1 + all_data[, c(col)])
}

write.csv(all_data, file = "data/dataset.csv", row.names = FALSE)
