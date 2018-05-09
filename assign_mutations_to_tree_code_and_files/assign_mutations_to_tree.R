# Assign mutations back to the tree

library(ape)
library(Rphylip)

## read in tree
tree <- read.tree("consensus_tree.tree")

## read in matrix of mutations to assign from the tree.
muts <- read.csv("Shared_mutations.txt", sep="\t", header=T, stringsAsFactors = F)

# assign mutations to the tree
# function to get all the nodes of a tree downstream of a given node.
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}
# function to get all the counts downstream of a given node
getcounts <- function(tree, node, mut) {
  descendants <- getDescendants(tree, node)
  tips <- descendants[descendants %in% c(1:length(tree$tip.label))]
  tipnames <- tree$tip.label[tips]
  counts <- mut[colnames(mut) %in% tipnames]
  return(counts)
}

# Choose FP and FN rates 
TP <- 0.99
TN <- 0.99

# make a data frame that can fill in with probabilities. Each row is a mutation, each column is a branch.
# step 1) assign the mutations that fit the tree perfectly.
# step 2) any that do not fit perfectly, assign with a particular TP and TN rate.

noleaf <- tree$edge[!(tree$edge[,2] %in% 1:length(tree$tip.label)),] # keep the edges that do not lead to leaves.
countdf <- as.data.frame(matrix(nrow=nrow(muts), ncol=(nrow(noleaf)))) # make a matrix where every row is a mutation, and every column is a node in the tree
rownames(countdf) <- muts$Pos
colnames(countdf) <- noleaf[,2] 
countdf$best <- as.numeric(0)
dim(countdf)
head(countdf)

# NB the following takes a few hours to run.
for (i in 1:nrow(muts)) {
  m <- muts[i,]
  for (j in 1:nrow(noleaf)) {  
    node <- noleaf[j,2]
    inside <- getcounts(tree=tree,node=node,mut = m) # counts in all the daughters underneath this node.
    inside <- inside[,grep("BM", colnames(inside))]
    outside <- m[,!(colnames(m) %in% colnames(inside))][,-1] # counts in all the daughters not underneath this node.
    outside <- outside[,grep("BM", colnames(outside))]
      # don't count NAs inside or outside. So NA counts as a zero in general except that when they're outside they don't count against
    if (length(inside[(inside==1 & !is.na(inside))])/length(inside)==1 & length(outside[(outside==1 & !is.na(outside))])==0) {  
      countdf[i, "best"] <- node
      break
    }
    
    TPhere = dbinom(length(inside[inside==1]), length(inside), TP)
    TNelsewhere = dbinom(length(outside[outside==0]), length(outside), TN)
    ProbHere <- TPhere * TNelsewhere
    countdf[i,j] <- ProbHere
  }
  if (countdf[i,"best"]!=0) {
    next
  }
  else {
    if (length(colnames(countdf)[which(countdf[i,1:nrow(noleaf)]==max(as.numeric(countdf[i,1:nrow(noleaf)])))]) > 1) {
      countdf[i,"best"] <- "ambig"
    } else {
      countdf[i,"best"] <- colnames(countdf)[which(countdf[i,1:nrow(noleaf)]==max(as.numeric(countdf[i,1:nrow(noleaf)])))]
    }
  }
}

write.csv(countdf, "assigning_muts_to_tree_140_HSCs_2017.03.08_two_phases_2.csv", row.names=T)


## match up the edge lengths with the downstream node for the branch.
# give all branches length zero to start off with
tree$edge.length <- rep(0, nrow(tree$edge))
# give a length of 1 to terminal branches.
tree$edge.length[tree$edge[,2] %in% c(1:length(tree$tip.label))] <- 1
for (i in 1:nrow(tree$edge)) {
  if (tree$edge[i,2] %in% names(table(countdf$best))) {
    tree$edge.length[i] <- table(countdf$best)[names(table(countdf$best)) == tree$edge[i,2]]
  }
}

### add in the leaves
shcalls <- read.csv("Shearwater_calls_FDR0.95_all_muts.txt", sep="\t", header=T, stringsAsFactors = F)
privates <- shcalls[shcalls$counts==1,]
privates[,grep("BM", colnames(privates))] <- sapply(privates[,grep("BM", colnames(privates))], function(cell) as.numeric(as.character(cell)))
privcounts <- colSums(privates[,grep("BM", colnames(privates))], na.rm = T)
# now assign these private mutations to the branches.
privcounts <- privcounts[unlist(tree$tip.label)]

# the order of tip labels is the same as the order of tips in edge, so can just do this:
tree$edge.length[tree$edge[,2] %in% c(1:length(tree$tip.label))] <- privcounts

write.tree(tree, "tree_with_leaves_shearwater_two_stages.tree")