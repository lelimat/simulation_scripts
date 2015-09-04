# simulate_CNVs_RSVSim.R
# C: Sep 2, 2015
# M: Sep 2, 2015
# A: Leandro Lima <lelimaufc@gmail.com> / <leandrol@usc.edu>


## RSVSim ##

# Install packages
# source("http://bioconductor.org/biocLite.R")
# biocLite("RSVSim")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(RSVSim)
library(BSgenome.Hsapiens.UCSC.hg19)

for (chrom in paste('chr', c(10, 15, 20), sep='')) {

cat("\n\n* * *Simulating CNVs for ", chrom, "\n\n")
hg19_chr = DNAStringSet(Hsapiens[[chrom]])
names(hg19_chr) = chrom

k = 1000
dupN    <- delN    <- 100
dupMean <- delMean <- 10*k
dupSd   <- delSd   <- 5*k
dup_lens <- abs(as.integer(rnorm(dupN, dupMean, dupSd)))
del_lens <- abs(as.integer(rnorm(delN, delMean, delSd)))
dup_lens[1] <- del_lens[1] <- 100*k
dup_lens[2] <- del_lens[2] <- 300*k
dup_lens[3] <- del_lens[3] <- 500*k
dup_lens[4] <- del_lens[4] <- 800*k
dup_lens[5] <- del_lens[5] <- 1000*k # (1mb)

cat("Haplotype 1\n")
chr_with_CNV1 = simulateSV(genome=hg19_chr, dups=dupN, sizeDups=dup_lens, maxDups=1, dels=delN, sizeDels=del_lens)

# Creating 0 copies deletions
zero_copies_range <- 91:100
del_starts <- metadata(chr_with_CNV1)$deletions$Start[zero_copies_range]
del_ends   <- metadata(chr_with_CNV1)$deletions$End[zero_copies_range]
del_sizes  <- metadata(chr_with_CNV1)$deletions$Size[zero_copies_range]
knownDels  <- GRanges(IRanges(del_starts,del_ends), seqnames=chrom)
# Creating 4 copies duplications
four_copies_range <- 81:90
dup_starts <- metadata(chr_with_CNV1)$tandemDuplications$Start[four_copies_range]
dup_ends   <- metadata(chr_with_CNV1)$tandemDuplications$End[four_copies_range]
dup_sizes  <- metadata(chr_with_CNV1)$tandemDuplications$Size[four_copies_range]
knownDups  <- GRanges(IRanges(dup_starts,dup_ends), seqnames=chrom)

cat("Haplotype 2\n")
# names(knownDeletion) = "myDeletion"
chr_with_CNV2 = simulateSV(genome=hg19_chr, random=FALSE,
                           regionsDels=knownDels, sizeDels=del_sizes, dels=length(zero_copies_range),
                           regionsDups=knownDups, sizeDups=dup_sizes, dups=length(four_copies_range))

writeXStringSet(chr_with_CNV1, paste(chrom, "_with_CNV1.fa", sep=""), format="fasta")
writeXStringSet(chr_with_CNV2, paste(chrom, "_with_CNV2.fa", sep=""), format="fasta")
write.table(metadata(chr_with_CNV1)$deletions[,-1], paste('deletions_',chrom,'.bed', sep=""), quote=FALSE, row.names=FALSE)
write.table(metadata(chr_with_CNV1)$tandemDuplications[,-1], paste('duplications_',chrom,'.bed', sep=""), quote=FALSE, row.names=FALSE)

}


# Boxplots
duplications_chr10 <- read.table('duplications_chr10.bed', header=T)
duplications_chr15 <- read.table('duplications_chr15.bed', header=T)
duplications_chr20 <- read.table('duplications_chr20.bed', header=T)

deletions_chr10 <- read.table('deletions_chr10.bed', header=T)
deletions_chr15 <- read.table('deletions_chr15.bed', header=T)
deletions_chr20 <- read.table('deletions_chr20.bed', header=T)

all_dels <- data.frame(rbind(deletions_chr10, deletions_chr15, deletions_chr20)$Size, 'del')
all_dups <- data.frame(rbind(duplications_chr10, duplications_chr15, duplications_chr20)$Size, 'dup')
colnames(all_dels) <- colnames(all_dups) <- c('size', 'type')

all_CNVs <- rbind(all_dels, all_dups)

boxplot(all_CNVs$size ~ all_CNVs$type)
boxplot(all_CNVs$size[all_CNVs$size < 100000] ~ all_CNVs$type[all_CNVs$size < 100000])


