script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd("..")

library(dplyr)
library(ggplot2)
library(reshape2)

source("../cbp_data/get_intersections.R")



## UNION ALL BREAKPOINTS IN ONE FILE

data_path <- "data/target/"

fl <- list.files(data_path)
all_fl <- data.frame()
for (f in fl){
  nm <- gsub(".bed", "", gsub("All_rearrangements_", "", f))
  nm_common <- strsplit(x = nm, split = "_")[[1]][1]
  all_fl <- rbind.data.frame(
    all_fl,
    data.frame(nm=nm, nm_common=nm_common)
  )
}
all_fl

df_all <- data.frame()
for (i in 1:nrow(all_fl)){
  nm <- paste0("All_rearrangements_", all_fl[i, 'nm'], ".bed")
  df <- read.table(paste0(data_path, nm))
  df$cancer_type <- all_fl[i, 'nm_common']
  df_all <- rbind.data.frame(df_all, df)
}
names(df_all) <- c("chr",  "start", "end", "cancer_type")
df_all$chr <- as.character(df_all$chr)
# remove Y chromosome
df_all <- df_all[df_all$chr != "chrY", ]
nrow(df_all)
# 211 803
nrow(unique(df_all))
# 209 966
df_all <- unique(df_all)

df_all$chr <- gsub("chr", "", df_all$chr)




## INTERSECT WITH GOOD BINS

# set window size
win_len <- 500000

# path_to_data_folder <- "../cbp_data/data/"
output_path <- "data/target_preprocessed/"

# read bins
bins_path <- "data/other/"
bins <- read.csv(paste0(bins_path, "bins_", format(win_len, scientific = FALSE), ".csv"), 
                 stringsAsFactors = FALSE, row.names = 1) 
all_bkpt_counts <- bins
all_bkpt_counts <- all_bkpt_counts %>%
  select(chr, start, end) %>%
  setNames(c("chr", "from", "to"))  

# get intersection with good bins
# for each cancer type
for (canc in unique(df_all$cancer_type)){
  data_cancer <- df_all[df_all$cancer_type == canc, ]
  data_cancer$cancer_type <- NULL
  data_intersected <- get_intersection_intervals(bins, data_cancer)
  data_intersected <- data_intersected %>%
    group_by(chr, start, end) %>%
    summarize(bkpt_in_window = n())
  
  data_intersected <- data_intersected %>%
    select(chr, start, end, bkpt_in_window) %>%
    setNames(c("chr", "from", "to", paste0("bkpt_in_window_", canc)))
  
  all_bkpt_counts <- all_bkpt_counts %>%
    left_join(data_intersected, by = c('chr', "from", 'to'))
}

all_bkpt_counts[is.na(all_bkpt_counts)] <- 0
unique(all_bkpt_counts$chr)
all_bkpt_counts <- all_bkpt_counts[all_bkpt_counts$chr != "Y", ]



## CALCULATE DENSITY

quantity_names <- setdiff(names(all_bkpt_counts), c("chr", "from", "to"))

total_counts <- all_bkpt_counts %>%
  group_by(chr) %>%
  summarise_at(quantity_names, sum) %>%
  ungroup() %>%
  rename_at(quantity_names, function(x) paste0("total_", x))

all_bkpt_counts <- all_bkpt_counts %>%
  inner_join(total_counts, by=c("chr"))


for (col in quantity_names){
  new_col_name <- paste0("density_", tail(strsplit(x = col, split = "_")[[1]], 1))
  all_bkpt_counts <- all_bkpt_counts %>%
    mutate(!!new_col_name := !!as.name(col) / !!as.name(paste0("total_", col)))
}


densities <- all_bkpt_counts %>%
  select_at(vars(contains("density")))

densities_id <- all_bkpt_counts %>%
  select(c("chr", "from", "to", names(densities)))

data_path <- "data/target/final/"

write.csv(all_bkpt_counts, file = paste0(data_path, "all_data_", win_len, ".csv"), 
          row.names = FALSE)

