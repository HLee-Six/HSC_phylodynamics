# filter mpileup of caveman calls before genotyping with shearwater.

library(GenomicRanges)

useful <- read.csv('mpileup_useful_cols_198_samps.txt', sep="\t", header=T, stringsAsFactors = F)

# Only keep clonal samples
clones <- read.csv("Clonal_colonies.txt", header=F, stringsAsFactors = F)[,1]
clonals <- useful[, c(1:4, which(sapply(strsplit(colnames(useful), "_"), "[[", 1) %in% clones))]

mtr <- grep("MTR", colnames(clonals))
dep <- grep("DEP", colnames(clonals))
clonals$counts <- apply(clonals[,mtr], 1, function(row) length(which(row>0)))
hist(clonals$counts, col="grey", 100)

# I know from the HSC tree that there is an early bifurcation with at least about 20 samples on either side. 
# so no true somatic muts should be present in more than 120 or so.
# probably much less if actually does bifurcate.
cl <- clonals[clonals$counts>0 & clonals$counts<120,]

# get rid of mutations that all fall very near each other.
tt <- cl
cc <- tt[with(tt, order(Chrom, Pos)),]

uu <- data.frame()
for (i in unique(cc$Chrom)) {
  thischr <- cc[cc$Chrom==i,]
  thischr$toprev <- NA
  
  # get the distance to the previous mutation
  for (n in 2:nrow(thischr)) {
    dd <- thischr[n, "Pos"] - thischr[n-1, "Pos"]
    thischr$toprev[n] <- dd
  }
  thischr$toprev[1] <- 11 # no need to lose the first line.
  
  # put in the distance to the next one too
  thischr$tonext <- c(thischr$toprev[2:nrow(thischr)],11)
  uu <- as.data.frame(rbind(uu, thischr))
}

vv <- uu[as.numeric(uu$toprev>10) & as.numeric(uu$tonext) > 10,]

## filter out subs near indels

inds <- read.csv("all_pass_indels_198_samps", sep="\t", header=F, stringsAsFactors = F)
colnames(inds) <- c("Chrom", "Pos", "Ref", "Alt")
inds$type <- NA
inds$type[(nchar(inds$Ref))>(nchar(inds$Alt))] <- "del"
inds$type[(nchar(inds$Ref))<(nchar(inds$Alt))] <- "ins"
inds$Start <- NA
inds$End <- NA
inds$Start[inds$type=="ins"] <- inds$Pos[inds$type=="ins"]-10
inds$End[inds$type=="ins"] <- inds$Pos[inds$type=="ins"]+10
inds$Start[inds$type=="del"] <- inds$Pos[inds$type=="del"]-10
inds$End[inds$type=="del"] <- inds$Pos[inds$type=="del"] + nchar(inds$Ref[inds$type=="del"]) + 10

grinds <- GRanges(inds$Chrom, IRanges(start = inds$Start, end=inds$End))
grsubs <- GRanges(vv$Chrom, IRanges(start=as.numeric(vv$Pos), width=1))
overlaps <- findOverlaps(grsubs, grinds)
nearindels <- unique(queryHits(overlaps))

vv2 <- vv[-c(nearindels),]
vv2$toprev <- NULL
vv2$tonext <- NULL

### Now filter on depth. NB treath autosomes and sex chromosomes differently.
dep <- grep("DEP", colnames(vv2))
mtr <- grep("MTR", colnames(vv2))

aut <- vv2[!vv2$Chrom %in% c("X", "Y"),]
sex <- vv2[vv2$Chrom %in% c("X", "Y"),]

claut <- aut[,dep]
claut[claut<6] <- NA 
aut[,dep] <- claut

clsex <- sex[,dep]
clsex[clsex<3] <- NA 
sex[,dep] <- clsex

aa <- as.data.frame(rbind(aut, sex))
aa$counts <- NULL

deps <- grep("DEP", colnames(aa))
mtrs <- grep("MTR", colnames(aa))

avafs <- (aa[,mtrs])/(aa[,deps])
avafs <- as.data.frame(cbind(aa[,1:4], avafs))

# get rid of the mutations with lots of NAs.
nabymut <- apply(avafs[,5:ncol(avafs)], 1, function(row) length(which(is.na(row))))

# get rid of NAs in more than 5
xvafs <- avafs[nabymut<=5,]

# get rid of ones with no mutant counts
xvafs$counts <- apply(xvafs[,5:ncol(xvafs)], 1, function(row) length(which(row>0)))
yvafs <- xvafs[xvafs$counts %in% c(1:120),] 

# apply mapping filters
# Filter 1: must be at least 1 TP per position.
# Assume at least 1 TP per position. If there are 2 counts, meanvaf should be about 0.5/2 = 0.25. 
# If n counts and one of them is a TP then the meanvaf should be 0.5/n. 
# As there is some fluctuation in vaf of a TP due to sampling variation, allow anything above a line at 0.3/counts.

# Filter 2: get rid of sites with a lot of FPs at low vaf. Can have a cutoff on the max number of samples where have mtrs at a v. low vaf.
# If e.g. wanted to exclude sites with FP mtrs in more than 10% of the samples, the cutoff would have to be sthg like:
# TP = the number of samples that have real mutations.
# FP = the number of samples that have mtrs but are not real - so mtrs at low vaf.
# counts = TP + FP
# vaf from TPs = 0.5*TP
# vaf from FPs = 0.1*FP

# My cutoff FP depends on the FP *rate* - I want to get rid of sites that tend to get lots of FPs.
# If say want an FP rate of max 10%, i.e.
# FP rate = 0.1 * susceptible samples
# susceptible samples = total samples - TP
# so: FP rate= 0.1 * (143 - TP)

# Therefore: meanvaf=(tpvaf*TP + fpvaf*(fprate*(samples-TP)))/(TP+fprate*(samps-TP))
# where:
# TP is the number of true positive mutations
# tpvaf is the vaf of true positive mutations (e.g. 0.5)
# fpvaf is the vaf of false positive mutations (e.g. 0.1)
# samples is the total number of samples (143)

# counts  = (TP+fprate*(samps-TP))
#         = TP*(1-fprate) + fprate*samps
#     TP  = (counts - fprate*samps)/(1-fprate)

# therefore:
# meanvaf = (tpvaf*((counts-fprate*samps)/(1-fprate)) + fpvaf*(fprate*(samps-((counts-fprate*samps)/(1-fprate)))))/counts

yvafs$meanvaf <- apply(yvafs[,grep("BM", colnames(yvafs))], 1, function(row) mean(row[row>0], na.rm = T))

counts <- seq(from=1, to=120,by=1)
minvaf1 <- 0.3/counts
tpvaf <- 0.4
fpvaf <- 0.1
fprate <- 0.1
samps <- 140
minvaf2 <- (tpvaf*((counts-fprate*samps)/(1-fprate)) + fpvaf*(fprate*(samps-((counts-fprate*samps)/(1-fprate)))))/counts
mindat <- as.data.frame(cbind(counts, minvaf1, minvaf2))

plot(meanvaf~counts, data=yvafs, main="cutoffs: \n 0.3/counts \n 0.4*((counts-0.1*143)/(1-0.1)) + 0.1*(0.1*(143-((counts-0.1*143)/(1-0.1))))/counts")
points(minvaf1~counts, data=mindat, type="l", col="red")
points(minvaf2~counts, data=mindat, type="l", col="red")

#### apply these filters.
yvafs$minvaf1 <- apply(yvafs, 1, function(row, na.rm=T) (0.3/as.numeric(row["counts"])))
yvafs$minvaf2 <- apply(yvafs, 1, function(row) ((tpvaf*((as.numeric(row["counts"])-fprate*samps)/(1-fprate)) + fpvaf*(fprate*(samps-((as.numeric(row["counts"])-fprate*samps)/(1-fprate)))))/as.numeric(row["counts"])))

zvafs <- yvafs[yvafs$meanvaf>yvafs$minvaf1 & yvafs$meanvaf>yvafs$minvaf2,]
zvafs <- zvafs[with(zvafs, order(-counts)),]
zvafs$point1 <- apply(zvafs[,grep("BM", colnames(zvafs))], 1, function(row) length(which(row>0.1)))
zvafs$point2 <- apply(zvafs[,grep("BM", colnames(zvafs))], 1, function(row) length(which(row>0.2)))
zvafs <- zvafs[zvafs$point2 > 0,]

# make mutation list for shearwater
Mutlist <- zvafs[,c(1:4)]
write.table(Mutlist, "mutation_list_for_shearwater.txt", sep="\t", col.names = T, row.names=F, quote=F)