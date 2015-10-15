require(dplyr)
require(ggplot2)
require(scales)
require(seqinr)

vcf <- read.table("freebayes.20150717122736.data.tsv", sep="\t", header = T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(POS = POS - 2832)

# Table for SNP quality figure
vcf %>% select(POS, QUAL) %>%
  write.table("output/Ill.Position.SNPQual.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

# What does the QUAL distribution look like?
vcf %>% ggplot(aes(POS, QUAL)) + geom_point() + 
  xlab("CDS position") + ylab("QUAL") + scale_y_log10() + annotation_logticks(sides = "l") + 
  scale_x_continuous(breaks=c(0,500,1000,1500,2000,2500)) + theme_bw()

vcf.hq <- vcf %>% filter(QUAL > 1000)
alleles <- vcf.hq$ALT %>% strsplit(",")
numalleles <- alleles %>% lapply(length) %>% unlist
refs <- c()
pos <- c()
quals <- c()
for (i in 1:length(alleles)) {
  refs <- c(refs, rep(vcf.hq$REF[i], numalleles[i])) %>% as.character
  pos <- c(pos, rep(vcf.hq$POS[i], numalleles[i])) %>% as.integer
  quals <- c(quals, rep(vcf.hq$QUAL[i], numalleles[i])) %>% as.numeric
}
alleles <- alleles %>% unlist %>% as.character

timecourse <- data.frame(Sample=c(), Position=c(), Ref=c(), Allele=c(), QUAL=c(), DP=c(), RO=c(), AO=c()) %>% tbl_df
for (sample in names(vcf.hq)[10:length(vcf.hq)]) {
  smp <- vcf.hq[[sample]] %>% strsplit(":") %>% as.data.frame %>% t %>% as.data.frame(stringsAsFactors = F) %>% tbl_df
  names(smp) <- c("GT","DP","RO","QR","AO","QA","GL")
  dps <- c()
  ros <- c()
  aos <- c()
  for (i in 1:length(numalleles)) {
    dps <- c(dps, rep(smp$DP[i], numalleles[i])) %>% as.integer
    ros <- c(ros, rep(smp$RO[i], numalleles[i])) %>% as.integer
    aos <- c(aos, unlist(strsplit(smp$AO[i], ","))) %>% as.integer
  }
  timecourse <- rbind(timecourse, data.frame(Sample=substr(sample,8,10), Position=pos, Ref=refs, Allele=alleles, 
                                             QUAL=quals, DP=dps, RO=ros, AO=aos)) %>% tbl_df
}

# Which ones are complex?
timecourse %>% mutate(freq = AO/DP) %>% group_by(Position, Ref, Allele) %>% dplyr::summarize(max=max(freq)) %>% View

# HERE - Manually correct down to single-allele SNPs, save as Ill.SNP.Ref.Mut.tsv

# Figure out AA changes
SNPs <- read.table("Ill.SNP.Ref.Mut.tsv", sep="\t", header=T, stringsAsFactors = F) %>% tbl_df
FASTA <- read.table("SP055-rpoZ-cMyc-Cry1Ac1-d123.fasta", header = T, stringsAsFactors = F) %>% tbl_df
seq.wt <- FASTA$X.SP055.rpoZ.cMyc.Cry1Ac1.d123 %>% paste(collapse="")
cds <- seq.wt %>% substr(2833, 4971) %>% tolower
cds.wt <- strsplit(cds, "")[[1]]
cds.mut <- strsplit(cds, "")[[1]]
for (i in (1:nrow(SNPs))) {
  ntpos <- SNPs$Position[i]
  cds.mut[ntpos] <- tolower(SNPs$Allele[i])[1]
}
pep.wt <- translate(cds.wt)
pep.mut <- translate(cds.mut)

for (i in (1:nrow(SNPs))) {
  aapos <- ceiling(SNPs$Position[i]/3)
  SNPs$Mutation[i] <- paste0(pep.wt[aapos], aapos, pep.mut[aapos])
  SNPs$Synonymous[i] <- ifelse(pep.wt[aapos] == pep.mut[aapos], "S", "N")
}

write.table(SNPs, "Ill.SNP.Ref.Mut.tsv", sep="\t", quote = F, row.names = F, col.names = T)
