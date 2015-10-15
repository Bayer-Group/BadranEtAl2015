require(dplyr)
require(reshape2)
require(igraph)
require(stringdist)
require(scales)
require(ggplot2)

# Read in SNP data, correct AA pos for Cry1Ac indexing, figure out initial oligotype
SNPs <- read.table("Ill.SNP.Ref.Mut.tsv", sep="\t", header=T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(mutpos = as.integer(gsub("[A-Z]", "", Mutation)) - 103, mut = gsub("[0-9]", "", Mutation)) %>%
  mutate(Mutation = paste0(substr(mut, 1, 1), mutpos, substr(mut, 2, 2))) %>% select(-mutpos, -mut)
initial <- SNPs$Ref %>% paste(collapse="")

# Read in observed oligotypes
matrix.percent <- read.table("map.pb/all_reads-sc25-s1-a1.0-A0-M0/MATRIX-PERCENT.txt", header = T, stringsAsFactors = F) %>%
  melt(variable.name = "Oligotype", value.name = "Prevalence") %>% tbl_df %>% mutate(Prevalence = Prevalence/100) %>%
  mutate(Timepoint = as.integer(gsub("t", "", samples)))
oligotypes <- unique(matrix.percent$Oligotype) %>% as.character
# Correct oligotypes resulting from misalignments
rights <- oligotypes[grep("\\.", oligotypes, invert = T)]
for (wrong in oligotypes[grep("\\.", oligotypes)]) {
  right <- rights[which(rank(stringdist(gsub("X", "", wrong), rights, method = "hamming"), ties.method = "min") == 1)]
  if (length(right) == 1) {
    matrix.percent$Prevalence[which(matrix.percent$Oligotype == right)] <-
      matrix.percent$Prevalence[which(matrix.percent$Oligotype == right)] + 
      matrix.percent$Prevalence[which(matrix.percent$Oligotype == wrong)]
  }
  matrix.percent <- matrix.percent %>% filter(Oligotype != wrong)
}
oligotypes <- matrix.percent %>% filter(Prevalence > 0.01) %>% group_by(Oligotype) %>% summarize(first = min(Timepoint)) %>%
  arrange(first) %>% select(Oligotype) %>% unlist %>% as.vector

# Determine which oligotypes have coding changes and label the mutations
mutants <- data.frame(oligotype = oligotypes, coding = "WT", stringsAsFactors = F) %>% tbl_df
SNPs$row <- (1:nrow(SNPs))
SNPs.nonsyn <- SNPs %>% filter(Synonymous == "N")
for (i in (1:nrow(SNPs.nonsyn))) {
  snp <- SNPs.nonsyn$row[i]
  mutants$coding[which(substr(mutants$oligotype, snp, snp) == SNPs.nonsyn$Allele[i])] <-
    paste(mutants$coding[which(substr(mutants$oligotype, snp, snp) == SNPs.nonsyn$Allele[i])], SNPs.nonsyn$Mutation[i], sep=",")
}
mutants$coding <- gsub("WT,", "", mutants$coding)

# Place the oligotypes in order in a pedigree igraph
oligotypes <- mutants$oligotype[order(stringdist(initial, mutants$oligotype, method = "hamming"))]
pedigree <- graph.empty() + vertex(initial, coding = "WT")

# Check all oligotypes not currently placed in the pedigree
for (oligotype in setdiff(oligotypes, V(pedigree)$name)) {
  
  # Figure out how far away the oligotype is from all vertices in the pedigree
  hammingdists <- stringdist(oligotype, V(pedigree)$name, method = "hamming")
  
  # The closest vertex(es) are the parents (multiple parents indicates recombination)
  #   PROBLEM - what if a strain with 5 mutations recombines with a strain with 1?
  #   This method doesn't check that we haven't seen the new mutation before, it'll just call the recombinant a derivative of the 5-mutation strain
  parents <- V(pedigree)$name[which(rank(hammingdists) == min(rank(hammingdists)))]
  
  # Add the oligotype as a node to the pedigree
  pedigree <- pedigree + vertex(oligotype, coding = mutants$coding[which(mutants$oligotype==oligotype)])
  
  # Add edges linking the oligotype to its parent(s)
  for (parent in parents) {
    pedigree <- pedigree + edge(parent, oligotype, relationship = "Parent", distance=min(hammingdists))
  }
}

# Plot pedigree as tree
plot(pedigree, vertex.color = hcl(seq(15, 375, length = length(oligotypes)+1), l = 65, c = 100), vertex.label.color = "black", 
     vertex.label = V(pedigree)$name, edge.label = E(pedigree)$distance, edge.label.color = "black", edge.arrow.size = 0.25,
     edge.width = 3, layout = layout.reingold.tilford(pedigree, root = 1), vertex.size = 8, vertex.label.cex = 0.75, edge.label.cex = 0.75)

# Edit pedigree by hand, reload nodes and edges
nodes <- read.table("PACEpb-nodes.csv", sep = ",", header = T, stringsAsFactors = F) %>% tbl_df %>% rename(Oligotype = name)
# nodes <- read.table("/import/transfer/user/kturn1/node.csv", sep = ",", header = T, stringsAsFactors = F) %>% tbl_df %>% rename(Oligotype = name)
edges <- read.table("PACEpb-edges.csv", sep = ",", header = T, stringsAsFactors = F) %>% tbl_df
newpedigree <- graph.empty(directed = T) + vertices(nodes$Index) +
  edges(strsplit(paste(edges$Index1, edges$Index2, sep=".", collapse="."), split = "\\.")[[1]])
newpedigree <- newpedigree %>% igraph::set.vertex.attribute("coding", value = nodes$coding) %>%
  igraph::set.vertex.attribute("oligotype", value = nodes$Oligotype) %>%
  igraph::set.edge.attribute("mutation", value = edges$Mutation) %>%
  igraph::set.edge.attribute("type", value = edges$Type) %>%
  igraph::set.edge.attribute("curvature", value = edges$Curvature)

# Plot corrected pedigree
l <- read.table("PACEpb-layout.tsv", sep = "\t", header = F) %>% as.matrix
svg("pedigree.svg", width = 6, height = 10)
plot(newpedigree, vertex.color = hcl(seq(15, 375, length = length(oligotypes)+1), l = 65, c = 100), vertex.label.color = "black", 
     vertex.label = V(newpedigree)$name, edge.label = E(newpedigree)$mutation, edge.label.color = "black", edge.arrow.size = 0.5,
     edge.color = ifelse(E(newpedigree)$type == "Recombination", "#8da0cb", ifelse(E(newpedigree)$type == "Coding", "#66c2a5", "#fc8d62")),
     edge.width = 1.5, vertex.size = 8, vertex.label.cex = 0.75, edge.label.cex = 0.75, asp = F,
     layout = l, edge.curved = E(newpedigree)$curvature)
dev.off()

# Order oligotypes by position in pedigree, determine ymin and ymax for Muller plot
toplot <- matrix.percent %>% slice(match(V(newpedigree)$oligotype, Oligotype)) %>% arrange(Timepoint) %>% group_by(Timepoint) %>%
  mutate(ymax = cumsum(Prevalence), ymin = cumsum(Prevalence) - Prevalence) %>% inner_join(nodes)

# Muller plot
lab <- read.table("PACEpb-labels.tsv", sep = "\t", header = T) %>% tbl_df
p <- toplot %>% ggplot(aes(Timepoint, ymin = ymin, ymax = ymax, fill = factor(Index))) + geom_ribbon() + theme_bw() +
  scale_y_continuous("Prevalence", labels = percent) + xlab("Time (h)") + scale_fill_discrete(guide = F) +
  geom_vline(xintercept = 156) + geom_vline(xintercept = 300) + geom_vline(xintercept = 408) +
  annotate("text", lab$x, lab$y, label = lab$label, size = 4)
ggsave("mullerplot.svg", p, width = 8, height = 6)
