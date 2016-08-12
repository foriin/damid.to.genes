# Artem Ilyin, IMG RAS, LARG 2016
# The purpose of this script is to get sorted list of genes
# which have the highest Dam-PIWI/Dam ratio and are located in enriched domens,
# identified by 2 state HMM from SnapGHC package.
# Note that the last analysis was made using not GATC coords, but microarray
# tiles coords for better HMM function and comparison with DamID experiments
# on microarrays

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

# Get the table with genes coords
genes <- read.table('~/Загрузки/Drosophila_melanogaster.BDGP5.77.gtf', skip=7, header = F, sep = '\t')
chroms <- c("2L", "2R", "3R", "3L", "X")
genes <- genes[genes[,1] %in% chroms,]
genes <- genes[, c(1,3,4,5,7,9)]
genes <- droplevels(genes) # drop unused levels of factors
genes <- genes[order(genes[,1]),] # order by chrom, bro
genes[,6] <- as.character(genes[,6]) # no need to use gene names in factor
names(genes) <- c("chr", "type", "start", "end", "strand", "misc")
(str(genes))

# retrieve HMM domain data
domains <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/19.07.16_HMM2_piwi_mapped_to_MA/HMM2x_piwi_recounted_MA_tiles.bed', skip=1)
names(domains) <- c("chr", "start", "end")
# get the DamID data
# Actually you could use data from the beginning of R-counts scripts when you have nothing but read
# counts, but I used bedgraph profile because there weren't many modifications to original data except logarythmisation
DamID <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/19.07.16_HMM2_piwi_mapped_to_MA/piwi_dam.log2ratio.MA.wig', header=F, skip=1)
names(DamID) <- c("chr", "start", "end", "logdamid")

# Retrieve positions of tiles
tiles <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/MA.tile.coords.gff', sep=';', header = F)
tiles <- tiles[, c(1, 4, 5)]
names(tiles) <- c("chr", "start", "end")
all.tiles <- merge(x = tiles, y = DamID, by=c("chr", "start", "end"), all.x=T) # very important part, assign damid values

#############################################
# the core of the script
########################################
domains.score <- data.frame()
for (i in chroms){
  chrom <- paste0("chr", i)
  domains.chrom <- domains[domains$chr == chrom,] # leave genes and gatcs from current chromosome
  tiles.chrom <- all.tiles[all.tiles$chr == chrom,]
  domains.chrom$start.idx <- findInterval(domains.chrom$start, tiles.chrom$start)
  domains.chrom$end.idx <- findInterval(domains.chrom$end, tiles.chrom$start)
  
  domains.chrom$DamID.sig <- 0
  for (j in 1:length(domains.chrom$DamID.sig)){
    domains.chrom$DamID.sig[j] <- assign_sig(domains.chrom[j,])
  }
  domains.chrom <- merge(domains.chrom,
                         data.frame("chr" = chrom, "start" = c(0, domains.chrom$end), "end" = NA, "DamID.sig" = 0),
                         all.y = T, all.x = T)
  domains.score <- rbind(domains.score, domains.chrom)
}


dupe <- data.frame()
for (i in chroms){
  chrom <- paste0("chr", i)
  genes.chrom <- genes[genes$chr == i,] # leave genes and gatcs from current chromosome
  doms.chrom <- domains.score[domains.score$chr == chrom,]
  genes.chrom$start.idx <- findInterval(genes.chrom$start, doms.chrom$start) # this determines in which gatcs the gene is located
  genes.chrom$end.idx <- findInterval(genes.chrom$end, doms.chrom$start)
  genes.chrom$score <- 0
  for (j in 1:length(genes.chrom$score)){
    genes.chrom$score[j] <- max(doms.chrom$DamID.sig[genes.chrom$start.idx[j]:genes.chrom$end.idx[j]])
  }
  dupe <- rbind(dupe, genes.chrom)
}

# Leave only genes with normal damid ratios (not NAs or -Inf), sort by DamID ratio and grab first 1k
gene.list <- dupe[dupe$score > 0,-c(7,8)]
gene.list.sorted <- gene.list[order(-gene.list$score),]
top1000.all <- gene.list.sorted[1:1000,]
top1000.genes <- gene.list.sorted[gene.list.sorted$type == "gene",][1:1000,]
top.genes <- gene.list.sorted[gene.list.sorted$type == "gene",]

misc.info <- strsplit(top1000.genes$misc, split = '; *')

misc.info <- strsplit(top.genes$misc, split = '; *')

misc.pre <- lapply(misc.info, function(x){
  a <- data.frame(strsplit(x, ' '), stringsAsFactors = F)
  colnames(a) <- a[1,]
  a <- a[-1,]
})

add.info <- do.call("rbind.fill", misc.pre)
top1000.genes <- cbind(top1000.genes[, -6], add.info)
top.genes <- cbind(top.genes[, -6], add.info)

setwd('~/IMG/DamID/piwi_pind_prof_ryaz/top_genes/')
write.table(top1000.genes, 'top1k.genes.damid.piwi.domains.csv', sep = ';', row.names = F)
write.table(top.genes, 'top.genes.damid.piwi.domains.csv', sep = ';', row.names = F)
