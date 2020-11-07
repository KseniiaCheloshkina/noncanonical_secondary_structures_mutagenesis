script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd("..")

library(dplyr)
library(reshape2)

#### breast

df <- read.table("data/target/breast.rearrangements.n560.highSpec.noInvArtefact.brassII.01.04", sep="\t")
df_start <- df[, c("V1", "V3", "V4")]
df_start %>% nrow()
# 77 695
df_start %>% unique() %>% nrow()
# 77 585
names(df_start) <- c("chr", "start", "end")

df_end <- df[, c("V5", "V7", "V8")]
df_end  %>% unique() %>% nrow()
# 77 603
names(df_end) <- c("chr", "start", "end")

full_df <- rbind.data.frame(
  df_start,
  df_end
  )
full_df %>% nrow()
# 155390
full_df %>% unique() %>% nrow()
# 155 143

full_df$chr <- paste0("chr", full_df$chr)
write.table(x = full_df, file = "data/target/bed/All_rearrangements_BRCA.bed", col.names = F, row.names = F, quote = F)
