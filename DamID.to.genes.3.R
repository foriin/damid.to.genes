rm(list=ls())
library("plyr")

assign_sig <- function(line){
  # This function will count averaged signal for domain,
  # based on tiles in which it located
  #####################################################
  
  # if gene correspond just to one gatc, get ratio of its length to gatc length
  # and multiply damid ratio by it
  if (line$start.idx == line$end.idx){
    sig <- tiles.chrom$logdamid[line$end.idx]
    return(sig)
  }
  span <- line$start.idx:line$end.idx #indices of tiles in which current domain is located
  ave.on <- length(span)
  sig <- 0
  for (j in span){
    if (is.na(tiles.chrom$logdamid[j])) {ave.on <- ave.on - 1; next}
    sig <- sig + 2^tiles.chrom$logdamid[j] # don't average in log2 scale!
  }
  return(log2(sig/ave.on))
  
}

setwd('~/IMG/DamID/piwi_pind_prof_ryaz/05.08.16 piwi huiwi/')

# Get the table with genes coords
genes <- read.table('~/Загрузки/Drosophila_melanogaster.BDGP5.77.gtf', skip=7, header = F, sep = '\t')
chroms <- c("chr2L", "chr2R", "chr3R", "chr3L", "chrX")
genes <- genes[genes[,1] %in% chroms,]
genes <- genes[, c(1,3,4,5,7,9)]
genes <- droplevels(genes) # drop unused levels of factors
genes <- genes[order(genes[,1]),] # order by chrom, bro
genes[,6] <- as.character(genes[,6]) # no need to use gene names in factor
names(genes) <- c("chr", "type", "start", "end", "strand", "misc")
(str(genes))

# Retrieve the domain data
domains <- read.table('Piwi 2xHMM domains after remapping to MA coordinates and removal of TEs(DamID Piwi HMM2 MAr)_with average Piwi scores.bed', header = F)
names(domains) <- c("chr","start","end", "score")


gene.scores <- data.frame()
for (i in chroms) {
  genes.chrom <- genes[genes$chr == i,]
  domains.chrom <- merge(domains[domains$chr == i,],
                         data.frame("chr" = i, "start" = c(0, domains[domains$chr == i,]$end),
                                    "end" = NA, "score" = NA), all.y = T, all.x = T)
  genes.chrom$start.idx <- findInterval(genes.chrom$start, domains.chrom$start) # this determines in which gatcs the gene is located
  genes.chrom$end.idx <- findInterval(genes.chrom$end, domains.chrom$start)
  genes.chrom$score <- 0
  for (j in 1:length(genes.chrom$score)){
    genes.chrom$score[j] <- max(domains.chrom$score[genes.chrom$start.idx[j]:genes.chrom$end.idx[j]], na.rm = T)
  }
  gene.scores <- rbind(gene.scores, genes.chrom)
}

gene.sc <- gene.scores[gene.scores$type == "gene",]
misc.info <- strsplit(gene.sc$misc, split = '; *')

misc.pre <- lapply(misc.info, function(x){
  a <- data.frame(strsplit(x, ' '), stringsAsFactors = F)
  colnames(a) <- a[1,]
  a <- a[-1,]
})

add.info <- do.call("rbind", misc.pre)
gene.sc <- cbind(gene.sc[,-(6:8)], add.info)

gene.sc.srt <- gene.sc[order(-gene.sc$score),]

transcripts <- genes[genes$type == "transcript",]
tss <- data.frame("chr" = transcripts$chr, "tss" = ifelse(transcripts$strand == "+", transcripts$start, transcripts$end), transcripts[, c(5,6)])

tss.scores <- data.frame()
for (i in chroms) {
  tss.chrom <- tss[tss$chr == i,]
  domains.chrom <- merge(domains[domains$chr == i,],
                         data.frame("chr" = i, "start" = c(0, domains[domains$chr == i,]$end),
                                    "end" = NA, "score" = NA), all.y = T, all.x = T)
  idx <- findInterval(tss.chrom$tss, domains.chrom$start) # this determines in which domain the tss is located
  tss.chrom$score <- domains.chrom$score[idx]
  tss.scores <- rbind(tss.scores, tss.chrom)
}

tss.scores <- tss.scores[!is.na(tss.scores$score),]

misc.info <- strsplit(tss.scores$misc, split = '; *')

misc.pre <- lapply(misc.info, function(x){
  a <- data.frame(strsplit(x, ' '), stringsAsFactors = F)
  colnames(a) <- a[1,]
  a <- a[-1,]
})

add.info <- do.call("rbind.fill", misc.pre)
tss.scores.srt <- cbind(tss.scores[, c(1,2,3,5)], add.info[,c(1,3,5,7,8,10)])
tss.scores.srt <- tss.scores.srt[order(-tss.scores.srt$score),]


write.table(gene.sc.srt, "piwi.damid.gene.enrichm.list.csv", row.names = F, sep=';')
write.table(tss.scores.srt, "piwi.damid.tss.domain.score.csv", row.names = F, sep=';')
  
  
  