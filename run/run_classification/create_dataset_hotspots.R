library(dplyr)

bkpt_data <- read.csv("data/target/final/all_data_5e+05.csv")


#### Create hotspots
q <- 0.99
density_cols <- grep(x = names(bkpt_data), pattern = "density", value = TRUE)

for (col in density_cols){
  q_value <- quantile(x = bkpt_data[, col], probs = q)
  q_vars <- names(q_value)
  for (q_var_cur in q_vars){
    new_col_name <- paste0("hsp_", format(q_var_cur, scientific = FALSE), "_", strsplit(x = col, split = "_")[[1]][2])
    bkpt_data <- bkpt_data %>%
      mutate(!!new_col_name := ifelse(!!as.name(col) > q_value[q_var_cur], 1, 0))    
  }
}

hsp_cols <- grep(x = names(bkpt_data), pattern = "hsp", value = TRUE)
select_cols <- c("chr", "from", "to", hsp_cols)
bkpt_data <- bkpt_data[select_cols]



##### Join features (small set)
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

all_data <- bkpt_data

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

write.csv(all_data, file = "data/hsp_dataset.csv", row.names = FALSE)




##### Join features (full set)
win_len <- "500000"
features_path <- "data/features/"
all_features_paths <- list.files(features_path)
all_features_paths <- grep(x=all_features_paths, pattern = win_len, value=T)

all_data <- bkpt_data

# collect conserved features - non-B DNA motifs
tissue_spec_path <- c(
  "DNase_seq_", "H3K4me1-human_", "H3K4me3-human_", "H3K9me3-human_", "H3K27ac-human_", 
  "H3K27me3-human_", "H3K36me3-human_")
tissue_spec_path <- paste0(tissue_spec_path,  win_len, ".csv")

all_full_conserved_features_path <- setdiff(all_features_paths, tissue_spec_path)

for (i in 1:length(all_full_conserved_features_path)){
  
  feature_data <- read.csv(paste0(features_path, all_full_conserved_features_path[i]),
                           stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_data <- all_data  %>%
    left_join(feature_data, by = c("chr", "from", "to"))    
  
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
      grep(x = fl, pattern = paste0("_", win_len, ".csv"), 
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
      grep(x = fl, pattern = paste0("_", win_len, ".csv"), 
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

write.csv(all_data, file = "data/full_hsp_dataset.csv", row.names = FALSE)
