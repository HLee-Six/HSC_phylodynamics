#### HSC Approximate Bayesian Computation for stem cell numbers and generation times.

print(Sys.time())

# generate random number.
myrand <- (runif(1, min=0, max=1000)*runif(1, min=100, max=200))+ runif(1, min=0, max=1000)

print(paste0("Setting random seed: ", myrand))
# set this as the seed, and print out the seed. That way should be able to replicate any results, I hope.
set.seed(myrand)

##############

# create HSC distribution
loghscs <- runif(1, min=log(as.numeric(commandArgs(TRUE)[1])), max=log(as.numeric(commandArgs(TRUE)[2])))
# and create a genday distribution from 50 to 150.
loggendays <- runif(1, min=log(as.numeric(commandArgs(TRUE)[3])), max=log(as.numeric(commandArgs(TRUE)[4])))

tothscs <- round(exp(loghscs))
gendays <- exp(loggendays)

print(paste0("Adult stable HSC pool size is ", tothscs))
print(paste0("Adult generation time between symmetrical divisions for a single cell is ", gendays, " days"))

# make a directory for this combination.
jobdir <- as.character((commandArgs(TRUE)[5]))   
subdir <- paste0(tothscs, "_HSCs_", gendays, "_gendays_", myrand, "_randseed")
## dir.create(file.path(jobdir, subdir), showWarnings = FALSE)

# go to specific directory for this combination.
setwd(file.path(jobdir, subdir))
print(paste0("Directory: ", getwd()))

# create an id that I can append to each file that I write out so that I know which simulation combination it came from.
# also include a random string so that I can identify each repeat of a simulation
id <- paste0(tothscs, "_HSCs_", gendays, "_gentime_", myrand, "_")
print(paste0("ID for this simulation is ", id))

# pick age at which stabel population size is reach (transition point)
# from my initial simulations, this didn't seem to matter much, so go with 5 years old.
age <- 5
transpt <- age * 365 # transition point: the age at which a stable population size is reached, measured in days
tot_muts <- 1000 # we observe 1000 substitutions in our tree.
# the number of mutations that occur on each lineage in childhood is 1/10th of the total by age 60 (1000 mutations by age 60)
childmuts <- tot_muts/10 

# now calculate the number of generations that we need to run the FW model for.
genfw <- gendays*24*60*60 # this is FW generation time in seconds.
bm_age <- 60 # this is the time from conception to the BM sample, measured in years: 60 years.

# the number of generations is the age at which the bm was taken minus the age at the transition point (when FW drift starts), all divided by the generation time.
# translate everything into seconds.
bm_age_secs <- bm_age * 365 * 24 * 60 * 60
transpt_secs <- transpt * 24 * 60 * 60
fw_time_secs <- bm_age_secs - transpt_secs

fw_gens <- round(fw_time_secs/genfw)
print(paste0("The number of Fisher Wright generations in this simulation is ", fw_gens))

# calculate the number of mutations acquired in every FW generation
fw_muts_per_gen <- (tot_muts - childmuts)/fw_gens
print(paste0("The number of mutations acquired per Fisher Wright generation is ", fw_muts_per_gen))

# read in uids.
# uids are stored in a file with 1000 uids per row.
rowsofuids <- ceiling(tothscs/1000)+1

if (rowsofuids < 15019) {
  myuids <- read.csv("uids.txt", header=F, sep=" ", stringsAsFactors = F, nrows = rowsofuids)
} else {
  myuids <- read.csv("uids.txt", header=F, sep=" ", stringsAsFactors = F)
}

myuids <- as.vector(unlist(myuids))
myuids <- myuids[myuids != ""]
myuids <- myuids[!is.na(myuids)] 
print(paste0(length(myuids), " unique ids made")) 

# run the FW model for the specified number of generations.

ngen <- paste0(0, "~", myuids[1:tothscs])

# every 10 generations write it out and only keep the last generation.
for (f in 1:fw_gens) {
  print(paste0("FW generation ", f))
  print("about to sample")
  tgen <- sample(ngen, replace=T)
  print("have finished sampling")
  
  newgen <- myuids[1:length(tgen)]
  
  # add a prefix related to this bin to all the ids in this bin.
  prefnewgen <- paste0(f, "~", newgen)
  ngen <- paste0(tgen, "_", prefnewgen)
  
  # write it out every 10 generations.
  if (f %% 10 == 0) {
    print(paste0("Writing out ", f, "th iteration"))
    write(ngen, paste0(id, "generation_", f, ".txt"), ncolumns=1)
    
    # just take the last letter of the chain as the id signifying that chain.
    ngen <- prefnewgen
  }
}
file <- paste0(id, "final_gen.txt")
write(ngen, file=file, ncolumns=1)

### now stitch all the fw generations back together. 
# only need to do this if there were more than 9 generations.
if (fw_gens>9) {
  print("About to start stitching generations back together")
  
  tselect <- read.csv(list.files(pattern='final_gen.txt'), sep="_", header=F, stringsAsFactors=F)
  tsel <- as.data.frame(tselect[with(tselect, order(V1)),])
  write.table(tsel, paste0(id, "final_gen_reordered.txt"), sep="_", quote=F, col.names=F, row.names=F)
  concat <- as.data.frame(tsel)
  colnames(concat)[1] <- "V1"
  colnames(tsel)[1] <- "V1"
  if (ncol(tsel)==1) {
    tsel$V1 <- as.character(tsel$V1)
    concat$V1 <- as.character(concat$V1)
  }
  
  # work backwards from the final generation.
  for (gen in sort((1:fw_gens)[(1:fw_gens %% 10) ==0], decreasing=T)) {
    print(gen)
    tfile <- paste0(id, "generation_", gen, ".txt")
    wgen <- read.csv(tfile, sep="_", header=F, stringsAsFactors=F)
    
    print("about to start merging batches")
    merged <- merge(y=wgen, x=tsel, by.y=ncol(wgen), by.x=1, sort=F, all.x=T, all.y=F)
    print("have finished merging batches")
    
    ##########################
    
    # write out the chosen ones
    newdat <- as.data.frame(merged[,c(1, (ncol(merged) - 9):ncol(merged))]) # NB have changed this - now include the linker column.
    
    write.table(newdat, paste0(id, "generation_", gen, "_selected_from_next_gen.txt"),
                sep="_", col.names = F, row.names=F, quote=F)
    
    # now join concat and newdat together. Order them both on the same column, so that they should match up.
    colnames(newdat) <- c("x", paste0("N", 1:(ncol(newdat)-1)))
    concat <- concat[with(concat, order(V1, decreasing=F)),]
    newdat <- newdat[with(newdat, order(x, decreasing=F)),]
    concat <- as.data.frame(cbind(newdat, concat))
    print(paste0(nrow(concat[concat$x!=concat$V1,]), " discrepancies"))
    concat[,1] <- NULL
    colnames(concat) <- paste0("V", 1:ncol(concat))
    
    # join them together to form a big string.
    write.table(concat, paste0(id, "all_fw_gens_concatenated.txt"), sep="_", col.names = F, row.names = F, quote=F) # let this over-write itself as go along.
    
    # the new tselect now becomes the first column of concat, so can pick the next one.
    tsel <- concat[,1]
  }
  # read back in the final file
  ngen <- read.csv(paste0(id, "all_fw_gens_concatenated.txt"), sep="\t", header=F, stringsAsFactors = F)
  ngen <- as.vector(ngen[,1])
  
  # delete the intermediate files
  print("About to remove intermediate files")
  toremove <- list.files(pattern="_generation_")
  toremove <- c(toremove, list.files(pattern="final_gen"))
  file.remove(toremove)
} 


##### PERIPHERAL BLOOD PART
print("Moving onto peripheral blood part.")

# function to add mutations:
myconstants <- apply(combn(c(letters, LETTERS), 2), 2, function(col) paste(as.vector(col), collapse=""))

add.mutns.fn <- function(cell, mutn_rate) {
  newmuts <- myconstants[1:rpois(1, lambda=mutn_rate)]
  withmuts <- paste(paste0(cell, newmuts), collapse="-")
  return(withmuts)
}

tid <- id  
hgen <- ngen
samplenum <- 155
wgs <- hgen[sample(length(hgen),samplenum, replace = F)] # sample the stem cells to be whole genome sequenced.

write(wgs, paste0(tid, "sampled_HSCs.txt"), ncolumns=1)
bygen <- read.csv(paste0(tid, "sampled_HSCs.txt"), sep="_", header=F, stringsAsFactors = F)

### add mutations in. 
adulthood <- bygen

print("Adding mutations to adult phase")
admarkers <- (unique(unlist(adulthood)))
print(paste0(length(admarkers), " unique generation markers from adulthood"))

batchsize <- 500
batchstart <- seq(1,length(admarkers), by=batchsize)

adkey <- data.frame()
for (b in batchstart) {
  print(paste0("Up to adding mutations to adult generation ", b))
  batchend <- b + (batchsize-1)
  if (batchend > length(admarkers)) {
    batchend <- length(admarkers)
  }
  tadmarkers <- admarkers[b:batchend]
  tadmuts <- sapply(tadmarkers, function(marker) add.mutns.fn(marker, fw_muts_per_gen))
  tadout <- cbind(tadmarkers, tadmuts)
  adkey <- as.data.frame(rbind(adkey, tadout))
}
colnames(adkey) <- c("gen_marker", "mutations")
adkey$gen_marker <- as.character(adkey$gen_marker)
adkey$mutations <- as.character(adkey$mutations)

# now replace every cell in adulthood with its partner in adkey. 
ahood <- unlist(adulthood)
batchstart <- seq(1,length(ahood), by=batchsize)

adultmuts <- c()
for (b in batchstart) {
  print(paste0("Up to replacing adult gen marker ",b, " with mutations"))
  batchend <- b + (batchsize-1)
  if (batchend > length(ahood)) {
    batchend <- length(ahood)
  }
  tahood <- ahood[b:batchend]
  tadultmuts <- sapply(tahood, function(marker) adkey[adkey$gen_marker==marker, "mutations"])
  adultmuts <- c(adultmuts, tadultmuts)
}
ad <- matrix(adultmuts, nrow=nrow(adulthood), ncol=ncol(adulthood), byrow=F)

wimuts <- as.data.frame(ad)

bymuts <- apply(wimuts, 1, function(row) paste(as.character(row), collapse="-"))

asmat1 <- sapply(sapply(bymuts, function(hap) unlist(strsplit(hap, "-"))[1:1500]), "[[", 1)
# get rid of the mutation identifier
asmat2 <- sapply(asmat1, function(mut) substr(mut,1,(nchar(mut)-2))) 
bymut <- matrix(asmat2, nrow=length(bymuts), ncol=1500, byrow=T)

#####################
# Now we have our simulated tree of whole genomes.
# So do the peripheral blood draws

mymuts <- list()
allmuts <- c()
sharmuts <- c()
for (b in 1:length(bymuts)) {
  print(paste0("Up to sampled WGS clone ", b, " of bymuts when getting shared and private mutations"))
  tmuts <- unlist(strsplit(bymuts[b], "-"))
  mymuts[[b]] <- tmuts
  sharmuts <- c(sharmuts, tmuts[tmuts %in% allmuts])
  allmuts <- c(allmuts, unique(tmuts))
}

# read in input file which contains all the necessary info.
input <- read.csv("essential_input.txt", sep=" ", header=F, stringsAsFactors = F)
input <- as.vector(apply(input,1, unlist))[1:7257]
input <- as.vector(unlist(input))

realbaitshared <- as.numeric(input[1])
privnums <- as.numeric(input[2:141])

# take all the mutations that are shared (i.e. present in more than one row)
# sample a proportion of these to put in the baitset. The shared mutations are the ones higher up the tree, plus the ones also in the polyclonal samples.
# I'm only considering some of the polyclonal samples - the ones that share >10 mutations.
# only count the bait mutations past mutation 100, as these are the only ones that I am considering

if (realbaitshared>length(unique(sharmuts))) {
  simbaitshared <- unique(sharmuts)
} else {
  simbaitshared <- sample(unique(sharmuts),realbaitshared, replace=F)
}

# take some of the private mutations.
# to do this, sample 140 of the 155 lineages to result in clonal colonies. The 15 branches leading to polyclonal colonies do not have any private mutations.
# from these 140, sample the same number of mutations as we have on each branch.
clonhaps <- bymuts[sample(155,140, replace=F)]

privmuts <- c()
leftovers <- c()
for (i in 1:length(clonhaps)) {
  thisclone <- clonhaps[i]
  tc <- unlist(strsplit(thisclone, "-"))
  tpriv <- tc[!tc %in% sharmuts]
  if (length(tpriv) >= privnums[i]) {
    tbait <- sample(tpriv, privnums[i], replace=F)
  } else {
    if (length(tpriv)==0) {
      print("No private mutations on this branch")
      next
    }
    tbait <- unique(sample(tpriv, privnums[i], replace=T))
    print("Not enough private mutations on this branch to have a full baitset")
  }
  privmuts <- c(privmuts, tbait)
  leftovers <- c(leftovers, tpriv[!tpriv %in% tbait])
}

# if the resulting baitset is a bit small, because some branches had too few mutations, make up for it by randomly sampling from the other branches. Not perfect,  as longer branches will get more mutations, but will do.
extramuts <- sample(leftovers, (sum(privnums) - length(privmuts)), replace=F)
privmuts <- c(privmuts, extramuts)

print(paste0(length(unique(sharmuts)), " shared mutations"))
print(paste0(length(privmuts), " private mutations"))

simbait <- c(simbaitshared, privmuts)

write(simbait, paste0(tid, "_baitset.txt"), ncolumns=1000)

##########################################
print("Now extracting ltt summary statistics.")

# summary statistics:
# For each plot, break into bins. The mean for each bin can be a summary statistic.
# 1) lineages through time  (ltt) plot. i.e. mean number of lineages in bins through time. 

simltt <- as.data.frame(cbind(c(1:700), apply(bymut[,1:700], 2, function(col) length(unique(col)))))
colnames(simltt) <- c("moltime", "lineages")

# get the mean number of lineages in 100 mutation bins from 0 to 700.
lttmeans <- c()
for (i in c(1:7)*100) {
  tltt <- simltt[simltt$moltime<=i & simltt$moltime>(i-100),]
  lttmeans <- c(lttmeans, mean(tltt$lineages, na.rm=T))
}

##############################################
print("")
print("Now moving onto resampling grans part - v2")
print("")

# 2) resample from fisher wright generations 6 times and look at overlap in mtr counts.
# try assessing capture-recapture of different granulocyte samples.
# we split the June grans into 6 pots, and lysed them separately ("biological replicates"). 
# from the DNA quants, there are about 90000 grans in each of these bio reps.

# make an id
tid <- paste0(tothscs, "_HSCs_", round(gendays), "_gendays_", round(myrand))

# read in the whole stem cell pool.
infileopt1 <- list.files(pattern="_all_fw_gens_concatenated.txt") # this is if had to do stitching
infileopt2 <- list.files(pattern="_final_gen.txt") # this is if did no stitching
infile <- c(infileopt1, infileopt2)[1]

# read in the whole HSC population. These end in the suffix "_all_fw_gens_concatenated.txt"
print("About to read in ngen")
ngen <- read.csv(infile, sep="_", header=F, stringsAsFactors = F)
print("Have finished reading in ngen")

# read in the baitset.
baits <- as.vector(unlist(read.csv(list.files(pattern="_baitset.txt"), sep=" ", header=F, stringsAsFactors = F)))
baits <- baits[baits != ""]

# sample 90 000 stem cells that make peripheral blood WITH REPLACEMENT
contribhscs <- 90000

# work out vafs of mutations in baitset
# first sample 3952 of the mutations in the baitset.
# this was the number of very clean positions in the observed data that we kept for this analysis.
if (length(baits)>3952) {
  baits <- sample(baits, 3952, replace = F)
}

# remove the mutation identifier to leave the generation marker
simba <- substr(baits, 1, (nchar(baits)-2))  # have changed since lengthened the mutation identifier.
baitmarks <- (unique(simba)) # these are all the unique generation markers for which need to work out vafs.
baitmarks <- baitmarks[baitmarks != ""]

# make a key for adding to as loop through.
baitsmutgen <- as.data.frame(cbind(baits, simba))
colnames(baitsmutgen)[2] <- "baitmarks"

# now have to sample with each biological replicate.
for (biorep in c(1:6)) {
  print(paste0("Sampling biological replicate: ", biorep))
  
  tocontrib <- sample(nrow(ngen), contribhscs, replace=T)
  write(tocontrib, paste0(tid, "biorep", biorep, "_contributing_hsc_indices_cap_recap.txt"), ncolumns=100)
  hgen <- ngen[tocontrib,]
  
  ### work out the count of each baitset mutation in this sample.
  tcounts <- c()
  for (tbait in baitmarks) {
    gen <- as.numeric(sapply(strsplit(tbait, "~"), "[[", 1))
    tcount <- length(which(hgen[,(gen+1)]==tbait))
    tcounts <- c(tcounts, tcount)
  }
  # tcounts is then the number of cells, out of contribhscs, that have each mutation in that granulocyte sample.
  
  # so now, for each of two technical replicates, need to resample from these according to the depth.
  # NB, I probably want to resample without replacement here.
  sampcounts <- as.data.frame(cbind(baitmarks, tcounts))
  colnames(sampcounts)[2] <- paste0("cells", biorep)
  
  # NB, need to convert from genmarkers into mutations
  
  baitsmutgen <- merge(baitsmutgen, sampcounts, all.x=T)
}

# baitsmutgen now contains, for each of the biological replicates, how many cells with each mutation ended up in that bio rep.
# now need to work out how many mutant reads get in each bio rep.

######
grans <- baitsmutgen
grans[,c(3:8)] <- matrix((sapply(unlist(grans[,c(3:8)]), as.numeric)), ncol=6)

# read in a file of sequencing depths here. 
# In the observed data, some samples have more uneven coverage than other.
# so read in a file of the observed depths.
grandeps <- read.csv("Deps_per_subsample_for_sims.txt",
                     sep="\t", header=T, stringsAsFactors=F)

# now have to resample different alleles.

# can simulate a bit of sequencing error.
# for every mutation, pick an error rate from the distribution observed in the control samples.
# use this to assign a random number of errors to each sample for that mutation.
# have one with no errors, one with an error rate drawn from the controls, and one with double that error rate
# NB these error rates are generous, as I set any site with 0 mtrs in the control to 1/10000

errors <- as.numeric(unlist(read.csv("error_rate_from_controls_zero_set_to_one_in_tenthousand.txt",
                                     stringsAsFactors = F, header=F)))

bmg <- as.data.frame(cbind(grans, matrix(nrow=nrow(grans), ncol=18)))
colnames(bmg)[9:26] <- c(paste0("sample", c(1:6), "error0"), paste0("sample", c(1:6), "error1"), paste0("sample", c(1:6), "error2"))

# there is also a tiny chance that the same mutation might have occurred independently in a different cell.
# cord blood would underestimate this as fewer genuine somatic mutations.
# 1000 mutations in 3e+09 base pairs. So for every mutation there is this source of error too.
# possible that grans acquire lots of muts as they divide, but I would be surprised if they could get more than another 1000 muts.
homoplasmyrate <- 2000/(3e+09) # to be generous

# NB because cells are diploid have to double the number of non-mutant alleles!
for (mut in 1:nrow(bmg)) {
  muterror <- sample(errors, 1) # every mutation gets its error rate.
  for (biorep in c(1:6)) {
    tcounts <- as.numeric(as.character(bmg[mut,paste0("cells", biorep)]))
    talleles <- c(rep("x", tcounts), rep("y", (contribhscs*2 - tcounts)))
    tdep <- sample(grandeps[,biorep], 1)
    tseq <- sample(talleles, tdep, replace=F)
    tmtrs <- length(which(tseq=='x'))
    bmg[mut,paste0("sample", biorep, "error0")] <- tmtrs
    
    # add in error.
    bmg[mut,paste0("sample", biorep, "error1")] <- tmtrs + rbinom(n=1, size=tdep, p=muterror) + rbinom(n=1, size=tdep, p=homoplasmyrate)
    bmg[mut,paste0("sample", biorep, "error2")] <- tmtrs + rbinom(n=1, size=tdep, p=(muterror*2)) + rbinom(n=1, size=tdep, p=homoplasmyrate)
  }
}
write.table(bmg, "alleles_with_mtrs_per_biorep_cap_recap_2018.02.09.txt", sep="\t", col.names=T, row.names=F)

#####
# now extract these new summary statistics.
# for VAF cutoffs  of MTR 1 to MTR 5, how many samples have more than the specified number of mutant reads?

cutoffs <- c(1:6)
sumstats <- c()
for (errorrate in c(0,1,2)) {
  for (cutoff in cutoffs) {
    passcut <- apply(bmg[,grep(paste0("error", errorrate), colnames(bmg))], 1, function(row) length(which(row>=cutoff)))
    ttab <- as.numeric(table(factor(passcut, levels=c(0:6))))
    sumstats <- c(sumstats, ttab)
  }
}

caprecapstats <- c(tothscs, gendays, sumstats, lttmeans)

# write out
write(caprecapstats, paste0(tid, "summary_statistics_cap_recap_v2_with_ltt.txt"), ncolumns = length(caprecapstats))


# zip files.
print("about to zip files")
tozip <- list.files()
zip(paste0(tid, "output_files_drift_and_caprecap.zip"), files=tozip)
allfiles <- list.files()
toremove <- allfiles[grep("zip", allfiles, invert=T)]
file.remove(toremove)

print(Sys.time())
