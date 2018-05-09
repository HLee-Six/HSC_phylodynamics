# HSC project
# fit MCMCglmm to analyse targeted data and plot

library(ape)
library(phytools)
library(MCMCglmm)
library(reshape2)
library(plyr)
library(RColorBrewer)

# read in model input data.

flm <- read.csv("model_input.txt", sep="\t", header=T, stringsAsFactors = F)

# isolate the mtrs.
a1 <- flm[, c("key", "node", grep("MTR", colnames(flm), value=T))]

# isolate the deps
a2 <- flm[, c("key", "node", grep("DEP", colnames(flm), value=T))]

# make it such that every mutation in every sample becomes a line.
recap.df <- melt(a1, id.vars = c("key", "node"), value.name = "MTR")
temp <- melt(a2, id.vars = c("key", "node"), value.name = "DEPTH")
recap.df$DEPTH <- temp$DEPTH
recap.df$node.fact <- factor(recap.df$node)

# make a variable saying what mutation it is and what sample it's in
recap.df$key.variable <- paste(recap.df$key, recap.df$variable, sep="_") 

# make a variable saying wheter it's a test or a control sample.
recap.df$Test <- 1
recap.df$Test[recap.df$variable=="controls_MTR"] <- 0

recap.df$Test.indicator <- recap.df$key.variable
recap.df$Test.indicator[recap.df$Test == 0] <- 0    # make this variable a 0 for all the control samples

recap.df$Node.indicator <- paste(recap.df$node, recap.df$variable, sep="_") # name the node and whether it's a test or a control sample
recap.df$Node.indicator[recap.df$Test == 0] <- 0 # make it a 0 for all the control samples
recap.df$log.DEPTH <- log(recap.df$DEPTH)

#### Now fit the model.
# set the prior 
prior <- list(R = list(V = 0.1, nu=10000), 
              B = list(mu = matrix(c(0,1),2), V=matrix(c(1e6,0,0,1e-6), nrow=2)),
              G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2), 
                       G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2),
                       G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)))

recap.mcmc.full <- MCMCglmm(MTR ~ log.DEPTH,
                            random= ~ key + Test.indicator + Node.indicator,
                            data=recap.df,
                            nitt=300000, burnin=10000, thin=1000,
                            pr=TRUE, family = "poisson", prior = prior)

pdf("MCMC all samples 300000.pdf")
plot(recap.mcmc.full)
dev.off()

celltypes <- c("jungrans", "jangrans", "novgrans", "bcells", "tcells")

# get the median for every cell type.
# the difference to previously is that I don't want to set everything that spans 0 to 0.
# Just save not only the posterior median, but also the upper and lower 90% intervals.
confidence_level <- 0.9
for(cell_type in celltypes[2:length(celltypes)]) {
  print(cell_type)
  rdf <- recap.df[grepl(cell_type, recap.df$variable),]
  tpost.median <- sapply(1:length(unique(rdf$key)), function(i,x,y) {
    median(
      x$Sol[,paste0("Node.indicator.", y$node[i], "_", cell_type, "_MTR")] +
        x$Sol[,paste0("Test.indicator.", y$key[i], "_", cell_type, "_MTR")]
    )}, 
    x=recap.mcmc.full, y=rdf)
  
  tpost.LCI <- sapply(1:length(unique(rdf$key)), function(i,x,y) 
  {quantile(x$Sol[,paste0("Node.indicator.", y$node[i], "_", cell_type, "_MTR")] +
              x$Sol[,paste0("Test.indicator.", y$key[i], "_", cell_type, "_MTR")], ((1-confidence_level)/2))}, 
  x=recap.mcmc.full, y=rdf)
  
  tpost.UCI <- sapply(1:length(unique(rdf$key)), function(i,x,y) 
  {quantile(x$Sol[,paste0("Node.indicator.", y$node[i], "_", cell_type, "_MTR")] +
              x$Sol[,paste0("Test.indicator.", y$key[i], "_", cell_type, "_MTR")], ((1+confidence_level)/2))}, 
  x=recap.mcmc.full, y=rdf)
  
  # tpost.median[tpost.90 < 0] <- 0
  
  # turn this into the vaf
  tmedvaf <- exp(tpost.median + median(recap.mcmc.full$Sol[,"(Intercept)"]))
  tucivaf <- exp(tpost.UCI + median(recap.mcmc.full$Sol[,"(Intercept)"]))
  tlcivaf <- exp(tpost.LCI + median(recap.mcmc.full$Sol[,"(Intercept)"]))
  
  flm[,paste0(cell_type, "_med.vaf")] <- tmedvaf
  flm[,paste0(cell_type, "_lci.vaf")] <- tlcivaf
  flm[,paste0(cell_type, "_uci.vaf")] <- tucivaf
}
flm$truezero <- exp(median(recap.mcmc.full$Sol[,"(Intercept)"]))

write.table(flm, "posterior_vafs_combined_model_with_node_effects_300000_with_CI.txt", sep="\t", col.names=T, row.names = T, quote=F)

# halve all x chromosome muts for plotting.
flm[(substring(flm$key,1,1) == "X"), grep("vaf", colnames(flm), value=T)] <- flm[(substring(flm$key,1,1) == "X"), grep("vaf", colnames(flm), value=T)]/2

# annotate with the timings of when branches start and stop
timings <- read.csv("branch_timings.txt", sep="\t", header=T, stringsAsFactors = F)
fla <- merge(flm, timings)



# plot the vafs down a given branch of the tree.
# function for plotting
plot.contrib.fn2 <- function(data, tree, bloodfraction, main, adjfactor, ylim=NULL, colourkey=FALSE, tiplegened=FALSE) {
  par(xpd=FALSE)
  mycols <- rep("grey", length(tree$edge.length))
  plot(tree, show.tip.label = F, edge.color = mycols, no.margin = F, direction = "downwards",
       main=main, font=2, tip.color = FALSE, cex=1.5, y.lim=ylim)
  
  # add scale bar
  axisPhylo(side=2, root.time = 0, backward=F, las=2)
  axis(side=2, labels="Molecular time (mutations)", pos=c(-20,600), at=600, tick=F)
  
  maxlen <- adjfactor # I'm not sure why I need this number
  for (i in 1:nrow(data)) {
    if (data[i, bloodfraction]=="none") {
      width <- 0.15
    } else {
      width <- 0.4
    }
    segments(y0=(maxlen - as.numeric(data[i,"muttime"])), 
             y1=(maxlen - as.numeric(data[i,"muttime"])),
             x0=as.numeric(data[i,"ycoord"]-width), x1=as.numeric(data[i,"ycoord"]+width),
             col=getcols.fn(data[i, bloodfraction]), lwd=2)
  }
}
getcols.fn <- function(combination) {
  if (combination=="none") {
    return("grey")
  }
  if (combination=="gran:B:T") {
    return("black")
  }
  if (combination=="gran:B") {
    return("purple")
  }
  if (combination=="gran:T") {
    return("tan1")
  }
  if (combination=="B:T") {
    return("cyan2")
  }
  if (combination=="gran") {
    return("firebrick3")
  }
  if (combination=="granlow") {
    return("pink") # this may still be too dark
  }
  if (combination=="B") {
    return("deepskyblue3")
  }
  if (combination=="T") {
    return("seagreen3")
  }
}

atf <- fla[,c("key", "node", "branchstart", "branchend", "ycoord", "jungrans_med.vaf", "bcells_med.vaf", "tcells_med.vaf")]
atff <- fla[,c("key", "node", "branchstart", "branchend", "ycoord", "jungrans_lci.vaf", "bcells_lci.vaf", "tcells_lci.vaf")]
atfg <- fla[,c("key", "node", "branchstart", "branchend", "ycoord", "jungrans_MTR", "bcells_MTR", "tcells_MTR")]
atf[,grep("vaf", colnames(atf), value=T)][atff[,grep("vaf", colnames(atff), value=T)] < flm$truezero[1]] <- 0
atf[,grep("vaf", colnames(atf), value=T)][atfg[,grep("MTR", colnames(atfg), value=T)]==0] <- 0

atf$granhigh <- 0
atf$granhigh[atf$jungrans_med.vaf>=(1/2000)] <- 1
atf$granlow <- 0
atf$granlow[atf$jungrans_med.vaf<(1/2000) & atf$jungrans_med.vaf>0] <- 1
atf$Bpres <- 0
atf$Bpres[atf$bcells_med.vaf>0] <- 1  
atf$Tpres <- 0
atf$Tpres[atf$tcells_med.vaf>0] <- 1  

atf$combn <- apply(atf[,grep("pres|^gran", colnames(atf), value=T)],1,function(row)
  paste(gsub("pres", "", grep("pres|^gran", colnames(atf), value=T)[row==1]), collapse=":"))
table(atf$combn)
atf$combn[atf$combn==""] <- "none"

# doesn't matter if gran low or high if has been called in another lineage
atf$combn[grep('granlow:', atf$combn)] <- gsub("granlow", "gran", atf$combn[grep('granlow:', atf$combn)])
atf$combn[grep('granhigh:', atf$combn)] <- gsub("granhigh", "gran", atf$combn[grep('granhigh:', atf$combn)])
table(atf$combn)
atf$combn[atf$combn=="granhigh"] <- "gran"

# get the new time. Reorder for this.
ordf <- data.frame()
for (node in unique(atf$node)) {
  print(node)
  tnode <- atf[atf$node==node,]
  tnode$cvaf <- rowSums(tnode[,grep("med.vaf", colnames(tnode), value=T)], na.rm=T)
  tnode <- tnode[with(tnode, order(cvaf, decreasing=T)),]
  if (nrow(tnode)==1) {
    tspacing <- (tnode$branchend[1] - tnode$branchstart[1])
  } else {
    tspacing <- (tnode$branchend[1] - tnode$branchstart[1])/(nrow(tnode) - 1)
  }
  tnode$ntimes <- seq(from=(tnode$branchstart[1] + (tspacing/2)), to=(tnode$branchend[1] - (tspacing/2)),length.out = nrow(tnode))  
  ordf <- as.data.frame(rbind(ordf, tnode))
}
ath <- merge(atf, ordf[,c("key", "ntimes")])
write.table(ath, "targeted_annotated_by_lineage_lci_cutoff_high_low_gran_vaf_only_jungrans_no_mtrs_set_to_zero.txt", sep="\t", col.names=T, row.names=F, quote=F)

ati <- ath[,c("key", "node", "branchstart", "branchend", "ntimes", "ycoord", "combn")]
colnames(ati)[colnames(ati)=="ntimes"] <- "muttime"

# read in the tree.
tree <- read.tree("hspc_tree.tree")

pdf("Figure_6a.pdf", width=10, height=7)
plot.contrib.fn2(tree=tree, data=ati, bloodfraction="combn", main="If no mtrs set to zero", adjfactor=max(nodeHeights(tree)))
par(xpd=T)
legend(x=20, y=50, legend = c("G,B,T", "G,B", "G,T", "B,T"),
       fill=sapply(c("gran:B:T", "gran:B", "gran:T", "B:T", "gran", "B", "T", "none"), getcols.fn), bty="n")
legend(x=40, y=50, legend = c("G", "G low VAF", "B", "T", "none"),
       fill=sapply(c("gran", "granlow", "B", "T", "none"), getcols.fn), bty="n")
dev.off()



# plot the vafs down the tree for figure 4.



############################################################
############################################################

# replot figure 4 with greater sensitivity.
th <- ath
min(th$jungrans_med.vaf[th$jungrans_med.vaf>0])
hist(log(th$jungrans_med.vaf[th$jungrans_med.vaf>0]), col="grey", 50)

# will need to log-transform, normalise, and correct highest vaf mutations, which are underestimated
# pick cutoff, below which we don't believe any mutations based on this histogram.
cutoff <- log(1/5000)

head(th)
aa <- th[,c("key", "node", "branchstart", "branchend", "ntimes", "ycoord", "jungrans_med.vaf")]
cell_type <- "jungrans"
aa[,paste0(cell_type, "_best.est.vaf")] <- aa[,paste0(cell_type, "_med.vaf")]
aa[,paste0(cell_type, "_best.est.vaf")][log(aa[,paste0(cell_type, "_best.est.vaf")])<cutoff] <- 0

# add on the uncorrected vafs.
ab <- merge(aa, flm[,c("key", "jungrans_VAF")])

# they tend to be a bit off for higher vaf muts, where we can trust the observed values - trying above a certain cutoff set them back to their original value
topcut <- 0.02  # beyond this none of the observed are being brought right down to zero
ab[,paste0(cell_type, "_best.est.vaf")][ab[,paste0(cell_type, "_VAF")] > topcut] <- ab[,paste0(cell_type, "_VAF")][ab[,paste0(cell_type, "_VAF")] > topcut]

# log
ab[,paste0(cell_type, "_log.best.est.vaf")] <- log(ab[,paste0(cell_type, "_best.est.vaf")])
ab[,paste0(cell_type, "_log.best.est.vaf")][ab[,paste0(cell_type, "_log.best.est.vaf")]=="-Inf"] <- cutoff

# normalise for plotting.
ab[,paste0(cell_type, "_norm.log.best.est.vaf")] <- (ab[,paste0(cell_type, "_log.best.est.vaf")] - cutoff)/(max(unlist(ab[,grep("log", colnames(ab), value=T)])) - cutoff)
write.table(ab, 'model_with_node_effects_output_corrected_and_logged_cutoff_minus8.5_topcutpt02_300000its.txt', sep="\t", col.names = T, row.names = T, quote=F)

#################
### now plot. ###
#################
# functions needed for plotting.
frbs <- ab
colnames(frbs)[colnames(frbs)=="ntimes"] <- "modeltime"

plot.mutation.fn2 <- function(data, tree, bloodfraction, main, adjfactor, ylim=NULL, xlim=NULL) {
  par(xpd=FALSE)
  mycols <- rep("grey", length(tree$edge.length))
  plot(tree, show.tip.label = F, edge.color = mycols, no.margin = F, direction = "downwards",
       main=main, font=2, cex=1.5, y.lim=ylim, x.lim=xlim)
  
  # add scale bar
  axisPhylo(side=2, root.time = 0, backward=F, las=2)
  axis(side=2, labels="Molecular time (mutations)", pos=c(-20,600), at=600, tick=F)
  
  maxlen <- adjfactor # Only need this if have tip labels
  for (i in 1:nrow(data)) {
    if (data[i, bloodfraction]==0) {
      width <- 0.15
      mylwd <- 1
    } else {
      width <- 0.4
      mylwd=3
    }
    segments(y0=(maxlen - as.numeric(data[i,"modeltime"])), 
             y1=(maxlen - as.numeric(data[i,"modeltime"])),
             x0=as.numeric(data[i,"ycoord"]-width), x1=as.numeric(data[i,"ycoord"]+width),
             col=getvafcols.fn(data[i, bloodfraction]), lwd=mylwd)
  }
}

# the colour scale I think I want is:
# <1/2500 cells (i.e. vaf < 1/500) = grey
# 1/2500-1/500 = yellow
# 1/500 - 1/100 = orange
# 1/100 - 1/20 = red
# 1/20 - 1/4 = dark red
# > 1/4 = black
# i.e. 4 shades of red - but ignore the palest as hard to see. Get all the colours and then pick from them.

# work out the transition point to the darkest colour
RColorBrewer::display.brewer.pal(9, name="YlOrRd")
heatcols <- RColorBrewer::brewer.pal(9, "YlOrRd")
plot(rep(1,4), col=heatcols[c(3,5,7,9)], pch=16, cex=7, xlim=c(0,5))

# plot the colour key:
plot(rep(1,6), col=c("grey", heatcols[c(3,5,7,9)], "black"), pch=16, cex=5, xlim=c(0,7))
dev.off()

myvafcols <- c("grey", heatcols[c(3,5,7,9)], "black")
getvafcols.fn <- function(bestvaf) {
  if (bestvaf<=(1/5000)) {
    return(myvafcols[1])
  } 
  if (bestvaf>(1/5000) & bestvaf<=(1/1000)) {
    return(myvafcols[2])
  } 
  if (bestvaf>(1/1000) & bestvaf<=(1/200)) {
    return(myvafcols[3])
  } 
  if (bestvaf>(1/200) & bestvaf<=(1/40)) {
    return(myvafcols[4])
  } 
  if (bestvaf>(1/40) & bestvaf<=(1/8)) {
    return(myvafcols[5])
  } 
  if (bestvaf>(1/8)) {
    return(myvafcols[6])
  }
}


# reorder based on vafs, so higher vaf mutations are plotted on top.
head(frbs)
saved <- frbs
frbs <- frbs[with(frbs, order(jungrans_best.est.vaf, decreasing=F)),]

pdf("figure4.pdf", width=7, height=5)
plot.mutation.fn2(data=frbs, tree=tree, bloodfraction=paste0(cell_type, "_best.est.vaf"),
                    main=cell_type, adjfactor=max(nodeHeights(tree)), ylim=NULL)
dev.off()

# plot just the embryonic muts.
plot.mutation.fn2(data=frbs, tree=tree, bloodfraction=paste0(cell_type, "_best.est.vaf"),
                  main=cell_type, adjfactor=10, ylim=c(0,10))

# plot a subpart of the tree as an inset key
plot.mutation.fn2(data=frbs, tree=tree, bloodfraction=paste0(cell_type, "_best.est.vaf"),
                  main=cell_type, adjfactor=30, ylim=c(0,30), xlim=c(5,20))
dev.off()
