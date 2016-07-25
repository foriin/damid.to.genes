# Artem Ilyin, IMG RAS, LARG 2016
# Get the list of domains and based on DamID data
# count their score

rm(list=ls())
#############################################
# Functions
#############################################
assign_sig <- function(line, tiles.chrom){
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

count.score <- function(domains, damid.data){
  domains.score <- data.frame()
  chroms <- c("2L", "2R", "3R", "3L", "X")
  for (i in chroms){
    chrom <- paste0("chr", i)
    domains.chrom <- domains[domains$chr == chrom,] # leave genes and gatcs from current chromosome
    tiles.chrom <- damid.data[damid.data$chr == chrom,]
    domains.chrom$start.idx <- findInterval(domains.chrom$start, tiles.chrom$start)
    domains.chrom$end.idx <- findInterval(domains.chrom$end, tiles.chrom$start)
    
    domains.chrom$DamID.sig <- 0
    for (j in 1:length(domains.chrom$DamID.sig)){
      domains.chrom$DamID.sig[j] <- assign_sig(domains.chrom[j,], tiles.chrom)
    }
    
    domains.score <- rbind(domains.score, domains.chrom)
  }
  return(domains.score)
}

#########################################
# Data preparation
#########################################

# Retrieve positions of tiles
tiles <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/MA.tile.coords.gff', sep=';', header = F)
tiles <- tiles[, c(1, 4, 5)]
names(tiles) <- c("chr", "start", "end")

# retrieve HMM domain data
domains.piwi <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/DamID.PIWI.HMM2.norep.bed', skip=1)
names(domains.piwi) <- c("chr", "start", "end")

domains.nup98 <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/NUP98.HMM2.norep.bed', skip = 2)
names(domains.nup98) <- c("chr", "start", "end")

domains.nup98.nuc <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/NUP98.nuc.HMM2.norep.bed', skip = 2)
names(domains.nup98.nuc) <- c("chr", "start", "end")

domains.nup98.npc <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/NUP98.npc.HMM2.norep.bed', skip = 2)
names(domains.nup98.npc) <- c("chr", "start", "end")

# get the piwi DamID data
# Actually you could use data from the beginning of R-counts scripts when you have nothing but read
# counts, but I used bedgraph profile because there weren't many modifications to original data except logarythmisation
DamID.piwi <- read.table('~/IMG/DamID/piwi_pind_prof_ryaz/19.07.16_HMM2_piwi_mapped_to_MA/piwi_dam.log2ratio.MA.wig', header=F, skip=1)
names(DamID.piwi) <- c("chr", "start", "end", "logdamid")
DamID.piwi <- merge(x = tiles, y = DamID.piwi, by=c("chr", "start", "end"), all.x=T) # very important part, assign damid values

# retrieve NUP98 DamID data
DamID.nup <- read.table('~/IMG/DamID/NUP98/redo HMM 27.04/Nup98_DamID_bychr.txt')
DamID.nup <- cbind(tiles, DamID.nup[, 3])
names(DamID.nup) <- c("chr", "start", "end", "logdamid")

DamID.nup.nuc <- read.table('~/IMG/DamID/NUP98/redo HMM 27.04/Nup98_DamID_nucleo_bychr.txt')
DamID.nup.nuc <- cbind(tiles, DamID.nup.nuc[, 3])
names(DamID.nup.nuc) <- c("chr", "start", "end", "logdamid")

DamID.nup.npc <- read.table('~/IMG/DamID/NUP98/redo HMM 27.04/Nup98_DamID_NPC_tet_bychr.txt')
DamID.nup.npc <- cbind(tiles, DamID.nup.npc[, 3])
names(DamID.nup.npc) <- c("chr", "start", "end", "logdamid")

########################################################
# Count score for each experiment
########################################################
domains.piwi.sc <- count.score(domains.piwi, DamID.piwi)
domains.nup98.sc <- count.score(domains.nup98, DamID.nup)
domains.nup98.nuc.sc <- count.score(domains.nup98.nuc, DamID.nup.nuc)
domains.nup98.npc.sc <- count.score(domains.nup98.npc, DamID.nup.npc)

#######################################################
# Write data to table
#######################################################

setwd('~/IMG/DamID/piwi_pind_prof_ryaz/piwi_to_MA/domains.scores/')
write.table(domains.piwi.sc[, -c(4,5)], "piwi.repeatmasker.domain.scores.csv", row.names = F, sep=';')
write.table(domains.nup98.sc[, -c(4,5)], "NUP98.repeatmasker.domain.scores.csv", row.names = F, sep=';')
write.table(domains.nup98.nuc.sc[, -c(4,5)], "NUP98.nucleoplasmic.repeatmasker.domain.scores.csv", row.names = F, sep=';')
write.table(domains.nup98.npc.sc[, -c(4,5)], "NUP98.NPC.tethered.repeatmasker.domain.scores.csv", row.names = F, sep=';')

