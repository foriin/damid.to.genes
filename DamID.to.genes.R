# Artem Ilyin, IMG RAS, LARG 2016
# The purpose of this script is to get sorted list of genes
# which have the highest Dam-PIWI/Dam ratio

rm(list=ls())
library("plyr") #this package provides very useful function rbind.fill

get_sig <- function(line){
  # This function will count averaged signal for gene,
  # based on gatcs in which it located
  #####################################################
  
  # if gene correspond just to one gatc, get ratio of its length to gatc length
  # and multiply damid ratio by it
  if (line$start.idx == line$end.idx){
    coef <- (line$end - line$start)/(gatcs.chrom$end[line$end.idx] - gatcs.chrom$start[line$end.idx])
    sig <- coef*(2^gatcs.chrom$logdamid[line$end.idx])
    return(log2(sig))
  }
  span <- line$start.idx:line$end.idx #indices of gatcs in which current gene is located
  sig <- 0
  for (j in span){
    if (is.na(gatcs.chrom$logdamid[j])) next
    left.coef <- (gatcs.chrom$end[j] - line$start)/(gatcs.chrom$end[j] - gatcs.chrom$start[j])
    right.coef <- (line$end - gatcs.chrom$start)/(gatcs.chrom$end[j] - gatcs.chrom$start[j])
    if (left.coef > 1) left.coef <- 1 # left.coef is needed when gatc contains start of gene
    if (right.coef > 1) right.coef <- 1 # right.coef is needed when gatc contains the end of the gene
    sig <- sig + left.coef * right.coef * 2^gatcs.chrom$logdamid[j] # don't average in log2 scale!
  }
  return(log2(sig))
  
}

# get list of genes/transcripts and leave only those on significant chromosomes
genes <- read.table('~/Загрузки/Drosophila_melanogaster.BDGP5.77.gtf', skip=7, header = F, sep = '\t')
chroms <- c("2L", "2R", "3R", "3L", "X")
genes <- genes[genes[,1] %in% chroms,]
genes <- genes[, c(1,3,4,5,7,9)]
genes <- droplevels(genes) # drop unused levels of factors
genes <- genes[order(genes[,1]),] # order by chrom, bro
genes[,6] <- as.character(genes[,6]) # no need to use gene names in factor
names(genes) <- c("chr", "type", "start", "end", "strand", "misc")
(str(genes))

# get the DamID data
# Actually you could use data from the beginning of R-counts scripts when you have nothing but read
# counts, but I used bedgraph profile because there weren't many modifications to original data except logarythmisation
DamID <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/20.06.16_piwi_redone/bedgraph/PIWI.bedgraph', header=F, skip=1)
names(DamID) <- c("chr", "start", "end", "logdamid")

# get the list of gatc which was used for damid experiments
# Worth mentioning that the second coord of gatc must be the first coord of the next gatc!
gatcs <- read.table('~/IMG/DamID/GATC-coords/DmelGATCfragments-r5_AI120515.gff', header=F)
gatcs[, 5] <- gatcs[, 5] - 3
gatcs <- gatcs[, c(1, 4, 5)]
names(gatcs) <- c("chr", "start", "end")
gatcs <- merge(x = gatcs, y = DamID, by=c("chr", "start", "end"), all.x=T) # very important part, assign damid values
# to gatcs that have it in experiment, but leave all

#############################################
# the core of the script
########################################
dupe <- data.frame()
for (i in chroms){
  chrom <- paste0("chr", i)
  genes.chrom <- genes[genes$chr == i,] # leave genes and gatcs from current chromosome
  gatcs.chrom <- gatcs[gatcs$chr == chrom,]
  genes.chrom$start.idx <- findInterval(genes.chrom$start, gatcs.chrom$start) # this determines in which gatcs the gene is located
  genes.chrom$end.idx <- findInterval(genes.chrom$end, gatcs.chrom$start)
  genes.chrom$DamID.sig <- 0
  for (j in 1:length(genes.chrom$DamID.sig)){
    genes.chrom$DamID.sig[j] <- get_sig(genes.chrom[j,])
  }
  dupe <- rbind(dupe, genes.chrom)
}


# Leave only genes with normal damid ratios (not NAs or -Inf), sort by DamID ratio and grab first 1k
gene.list <- dupe[is.finite(dupe$DamID.sig),-c(7,8)]
gene.list.sorted <- gene.list[order(-gene.list$DamID.sig),]
top1000.all <- gene.list.sorted[1:1000,]
top1000.genes <- gene.list.sorted[gene.list.sorted$type == "gene",][1:1000,]

# Get information about gene name and etc from misc column and parse it to sep columns
# Do it for dataframe with genes and dataframe with everything (cds, transcripts, etc...)
misc.info <- strsplit(top1000.genes$misc, split = '; *')
misc.info.all <- strsplit(top1000.all$misc, split = '; *')

misc.pre <- lapply(misc.info, function(x){
  a <- data.frame(strsplit(x, ' '), stringsAsFactors = F)
  colnames(a) <- a[1,]
  a <- a[-1,]
})

misc.pre.all <- lapply(misc.info.all, function(x){
  a <- data.frame(strsplit(x, ' '), stringsAsFactors = F)
  colnames(a) <- a[1,]
  a <- a[-1,]
})

add.info <- do.call("rbind.fill", misc.pre)
add.info.all <- do.call("rbind.fill", misc.pre)

# Add this information to dataframes and write them to files
top1000.genes <- cbind(top1000.genes[, -6], add.info)
top1000.all <- cbind(top1000.all[, -6], add.info.all)

setwd('~/IMG/DamID/piwi_pind_prof_ryaz/top_genes/')
write.table(top1000.genes, 'top1k.genes.damid.piwi.csv', sep = ';', row.names = F)
write.table(top1000.all, 'top1k.all.seq.damid.piwi.csv', sep = ';', row.names = F)




