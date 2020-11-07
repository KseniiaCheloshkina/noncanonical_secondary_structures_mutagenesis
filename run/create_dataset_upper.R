script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)

features_path <- "data/features/"

all_features_paths <- list.files(features_path)
all_features_paths <- grep(x=all_features_paths, pattern = "1000000", value=T)

# collect conserved features - non-B DNA motifs
tissue_spec_path <- c(
  "DNase_seq_1000000.csv",
  "H3K4me1-human_1000000.csv", "H3K4me3-human_1000000.csv", "H3K9me3-human_1000000.csv",
  "H3K27ac-human_1000000.csv", "H3K27me3-human_1000000.csv", "H3K36me3-human_1000000.csv"
)
all_full_conserved_features_path <- setdiff(all_features_paths, tissue_spec_path)

for (i in 1:length(all_full_conserved_features_path)){
  
  feature_data <- read.csv(paste0(features_path, all_full_conserved_features_path[i]),
                           stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  if (i == 1){
    all_data <- feature_data
  } else {
    all_data <- all_data  %>%
      left_join(feature_data, by = c("chr", "from", "to"))    
  }
  feat_col <- names(feature_data)[ncol(feature_data)]
  all_data[is.na(all_data[feat_col]), feat_col] <- 0
}

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

# collect other conserved features - folders
conserved_folders <- c("genome_regions", "tad")
all_full_conserved_features_path <- vector()
for (folder in conserved_folders){
  full_path <- paste0(features_path, folder, "/")
  fl <- list.files(full_path)
  all_full_conserved_features_path <- append(
    all_full_conserved_features_path, 
    paste0(
      full_path, 
      grep(x = fl, pattern = paste0("_", format(1000000, scientific = FALSE), ".csv"), 
           value = TRUE)
    )
  )
}

for (feat_path in all_full_conserved_features_path){
  feature_data <- read.csv(feat_path, stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_data$chr <- as.character(all_data$chr)
  all_data <- all_data  %>%
    left_join(feature_data, by = c("chr", "from", "to"))
  feat_col <- names(feature_data)[ncol(feature_data)]
  all_data[is.na(all_data[feat_col]), feat_col] <- 0
}

# collect other conserved features
tissue_spec_features_path <- c("DNA_methylation", "TF")
all_full_spec_features_path <- vector()

for (folder in tissue_spec_features_path){
  full_path <- paste0(features_path, folder, "/")
  fl <- list.files(full_path)
  all_full_spec_features_path <- append(
    all_full_spec_features_path, 
    paste0(
      full_path, 
      grep(x = fl, pattern = paste0("_", format(1000000, scientific = FALSE), ".csv"), 
           value = TRUE)
    )
  )
}
map_abb_to_cancer_type <- data.frame(
  abb = c("BRCA", "LIRI", "MALY", "OV", "PACA", "PBCA"),
  cancer_type = c("breast", "liver", "blood", "ovary", "pancreatic", "brain")
)
for (feat_path in all_full_spec_features_path){
  
  feature_data <- read.csv(feat_path, stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_cancer_types <- unique(feature_data$cancer_type)
  all_cancer_types <- all_cancer_types[all_cancer_types != "pancreas"]
  for (canc_type in all_cancer_types){
    feature_data_part <- feature_data[feature_data$cancer_type == canc_type, ]
    feature_data_part$cancer_type <- NULL
    feat_name <- names(feature_data_part)[ncol(feature_data_part)]
    new_canc_type <- as.character(map_abb_to_cancer_type[map_abb_to_cancer_type$cancer_type == canc_type,
                                                         "abb"])
    new_feat_name <- paste0("cancer_", new_canc_type, "_", feat_name)
    if (length(new_canc_type) != 0){
      feature_data_part <- feature_data_part %>%
        rename(!!new_feat_name := feat_name)
      
      all_data <- all_data  %>%
        left_join(feature_data_part, by=c("chr", "from", "to"))
      all_data[is.na(all_data[new_feat_name]), new_feat_name] <- 0  
    }
    
  }
}


# log transform with 1
identity_cols <- c("chr", "from", "to")
feature_names <- setdiff(names(all_data), identity_cols)
for (col in feature_names){
  all_data[, c(col)] <- log2(1 + all_data[, c(col)])
}

write.csv(all_data, file = "data/upper_dataset.csv", row.names = FALSE)

