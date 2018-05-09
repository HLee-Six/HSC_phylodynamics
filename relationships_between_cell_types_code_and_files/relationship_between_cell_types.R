## Script to analyse relationships between cells in tree.

library(ape)
library(phytools)

# Read in a file of every cell type.
cellkey <- read.csv("Key_for_colouring_tip_labels.txt", sep="\t", header=T, stringsAsFactors = F)

# Read tree in and make tree ultrametric.
tree <- read.tree(file="hspc_tree.tree")
plot(tree)
tree$tip.label[!tree$tip.label %in% cellkey$Sample] <-  'BMP312'

# # make tree ultrametric.
utree <- tree
for (i in 1:length(utree$tip.label)) {
  par <- utree$edge[utree$edge[,2]==i, 1]
  parheight <- nodeheight(tree, par)
  newbranlen <- max(nodeHeights(tree)) - parheight
  utree$edge.length[utree$edge[,2]==i] <- newbranlen
}

# resolve polytomies
ubtree <- multi2di(utree)
ubtree$edge.length[ubtree$edge.length<1] <- 1
utree <- ubtree
for (i in 1:length(tree$tip.label)) {
  par <- utree$edge[utree$edge[,2]==i, 1]
  parheight <- nodeheight(utree, par)
  newbranlen <- max(nodeHeights(tree)) - parheight
  utree$edge.length[utree$edge[,2]==i] <- newbranlen
}

# build a distance matrix
myvcv <- vcv(utree)
mymax <- max(nodeHeights(tree))
distmat <- matrix(sapply(myvcv, function(cell) ((mymax-cell)*2)), nrow=nrow(myvcv))
colnames(distmat) <- colnames(myvcv)
rownames(distmat) <- rownames(myvcv)
diag(distmat) <- 0

# function to perform amova.
amova.fn <- function(distmat, groupnames, cell_key) {
  groupnums <- length(groupnames)
  dw <- c()
  cellnums <- c()
  for (i in 1:groupnums) {
    tgroup <- groupnames[i]
    tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type==tgroup])
    cellnums <- c(cellnums, length(tcells))
    dw <- c(dw, sum(distmat[tcells, tcells]))
  }
  tdist <- distmat[colnames(distmat) %in% cellkey$Sample[cell_key$Cell_type %in% groupnames], rownames(distmat) %in% cell_key$Sample[cellkey$Cell_type %in% groupnames]]
  dap <- (sum(tdist) - sum(dw))/2
  
  dfAP <- length(groupnames) - 1 
  dfWP <- sum(cellnums-1)
  N <- sum(cellnums)
  SSwp <- sum(dw/(cellnums*2))
  SSap <- sum(((dw + dap)/(2*N)) - (dw/(2*cellnums)))
  
  MSwp <- SSwp/dfWP
  MSap <- SSap/dfAP
  nc <- (N - (sum(cellnums^2)/N))/dfAP
  varwp <- MSwp
  varap <- (MSap - MSwp)/nc
  obsphi <- varap/(varwp + varap)
  return(obsphi)
}

# function to randomise sample labels and repeat
randamova.fn <- function(distmat, groupnames, cell_key) {
  
  # change added 2018.01.22: when randomizing, only include the part of the distance matrix that involves the cell types being considered.
  tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames])
  #
  
  randmat <- distmat[tcells, tcells]
  colnames(randmat) <- sample(colnames(randmat))
  rownames(randmat) <- colnames(randmat)
  randphi <- amova.fn(distmat=randmat, groupnames=groupnames, cell_key=cell_key)
  return(randphi)
}

# function tying it all in together
amovapval.fn <- function(distmat, groupnames, cell_key, iterations, plottitle) {
  # calculate observed
  obsphi <- amova.fn(distmat=distmat, groupnames=groupnames, cell_key=cell_key)
  # calculate null
  randphis <- sapply(1:iterations, function(cell) randamova.fn(distmat = distmat, groupnames=groupnames, cell_key=cell_key))
  # calculate pval
  pval <- length(which(randphis>obsphi))/length(randphis) 
  
  hist(randphis, col="grey", 100, main=plottitle, xlab="Phi statistic")
  abline(v=obsphi, col="red", lwd=2)
  legend("topright", legend=paste0("Observed\n p = ", pval), lwd=2, col="red", bty="n")
}

# make separate groups for HSCs vs progenitors.
cellkey2 <- cellkey
cellkey2$Cell_type[cellkey$Cell_type %in% c("HSC", "HSC from periph bld")] <- "HSC"
cellkey2$Cell_type[cellkey$Cell_type %in% c("CMP", "MEP", "GMP")] <- "Progenitor"
table(cellkey2$Cell_type)

pdf("AMOVA_different_comparisons.pdf")
amovapval.fn(distmat=distmat, groupnames=c("HSC", "HSC from periph bld"), cell_key=cellkey, iterations = 30000, plottitle="BM vs PB HSCs")
amovapval.fn(distmat=distmat, groupnames=c("MEP", "CMP", "GMP"), cell_key=cellkey, iterations = 30000, plottitle="Progenitor types")
amovapval.fn(distmat=distmat, groupnames=unique(cellkey2$Cell_type), cell_key=cellkey2, iterations = 30000, plottitle="HSCs vs progenitors")
dev.off()



