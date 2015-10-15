require(dplyr)

samples <- read.table("PACE-samples.tsv", sep="\t", stringsAsFactors = F, header=T) %>% tbl_df
allreads.lengths <- read.table("map.pb/allreads.lengths.txt", header=F, sep="\t", col.names=c("Sample.Abbrev", "Read.Length")) %>% tbl_df
allreads <- left_join(samples, allreads.lengths)
allreads$Condition <- factor(allreads$Condition, c("Initial", "Low_Mut", "Med_Mut", "Low_WT", "High_WT"))

# Output files for figure construction
allreads %>% group_by(Read.Length) %>% summarize(Read.Count = n()) %>%
  write.table("output/PacBio.ReadLength.ReadCount.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
allreads %>% group_by(Time) %>% summarize(Read.Count = n()) %>% 
  write.table("output/PacBio.Timepoint.ReadCount.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
allreads %>% group_by(Time, Read.Length) %>% summarize(Read.Count = n()) %>%
  write.table("output/PacBio.Timepoint.ReadLength.ReadCount.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
