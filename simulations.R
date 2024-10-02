##
## D.M. de Vienne
## Sept 2024
## Set of functions used for computing INDUCED TRANSFERS 
## from a species tree and a collection of gene trees
## Fonctions for conversion to Newick and other tree manipulations
## are also available.
##

## Here we peform the simulations for the article "HGTs are not SPRs"
## There are two main functions that call the base functions 
## present in the other file "InferTransfersAfterSamplingCorrectly.R": 
##
## DOTHESIMULATION(): this function takes as input a number of tips for 
## the full species tree (Ntotal), a number of tips for the sampled species tree (n)
## and a number of HGTs (Ntransf), and does the following: 
## - Simulate the HGTs (as SPRs) in the "full species tree" to get the "full gene tree"
## - Sample the "full gene tree" (remove ghosts) to get the "sampled gene tree"
## - Compute the "induced HGTs" given the full species tree, the HGTs and the sampled species tree
## - Simulate the induced HGTs as SPRs in the sampled species tree to get the "induced gene tree" 
##
## The function returns the two gene trees, the sampled one and the induced one
## 
## GO(): this small function simply calls DOTHESIMULATION() and then compute the RF
## topological distance between the two obtained trees. This function is called multiple times
## changing the different arguments to perform the simulations presented in the article.


## SOURCE FUNCTIONS
source("InferTransfersAfterSamplingCorrectly.R")


DOTHESIMULATION<-function(Ntotal=1000, N=5, Ntransf=5, birth=1, death=0, plot=TRUE) {
	## 
	if (plot)	par(mfrow=c(1,3))

	## SIMULATE TREE

	tree<-rphylo(Ntotal, birth, death, fossils=TRUE)
	#print(tree)
	tips<-sample(1:Ntotal, N)
	ghosts<-setdiff(1:Ntotal, tips)
	EdgeWithGhost<-PREPARETREE(tree, tips)

	#print("a")
	## SIMULATE TRANSFERS

	nbhgt<-0
	HGT<-matrix(ncol=3, nrow=0)
	while(nbhgt<Ntransf) {
		hgt<-TRANSF(tree)
		if (KEEPTRANSFER(tree, EdgeWithGhost, hgt)) HGT<-rbind(HGT,hgt)
		nbhgt<-nrow(HGT)
	}

	##transform transfers to give it NODE names and not BRANCh names
	HGT2<-HGT
	HGT2[,1]<- EdgeWithGhost[HGT[,1],2]
	HGT2[,2]<- EdgeWithGhost[HGT[,2],2]
	HGT2<-as.data.frame(HGT2)
	colnames(HGT2)<-c("FROM", "TO", "TIME")
	HGT2<-HGT2[order(HGT2[,3]),]
	#print("b")
	# WARNING: in the HGT matrix, the numbers correspond to BRANCHES, i.e. rows of the edge matrix tree$edge.
	# In the HGT2 matrix, the numbers in FROM and TO correspond to descendant NODES. TIME of the transfer is absolute.


	## PLOT TREE TRANSFERS AND GHOSTS

 	if (plot) 	plotAllTransf(tree, tips, EdgeWithGhost, HGT)


	## PREDICT TRANSFERS INDUCED GIVEN TREE, GHOSTS AND TRANSFERS 

	Induced<-PredictTransfers(tree, HGT2, ghosts=ghosts, plot=plot, removeUndetectable=TRUE)

	## PERFORM SPRs on the small tree (following Inuced transfers) 
	#print("c")

	#We don't care about branch length, we thus keep this simple, with ONLY edge matrix to deal with
	smalledge<-GetSmallEdgeMatrix(tree, ghosts)$edge
	NEWEDGSMALL<-SIMUTRANSF(smalledge, Induced)
	NEWTREESMALL<-Newick(NEWEDGSMALL)
	SmallTree<-read.tree(text=NEWTREESMALL)


	## PERFORM SPRs ON THE LARGE TREE AND THEN SAMPLE NON-GHOSTS
	#print("d")
	bigedge<-tree$edge
	NEWEDGBIG<-SIMUTRANSF(bigedge, HGT2)
	NEWTREEBIG<-Newick(NEWEDGBIG)
	BigTree<-read.tree(text=NEWTREEBIG)
	SampledBigTree<-drop.tip(BigTree, match(as.character(ghosts), BigTree$tip.label))

	if (plot)	cophyloplot(SmallTree, SampledBigTree, assoc=cbind(SmallTree$tip.label, SmallTree$tip.label))
	return(list(induced=SmallTree, sampled=SampledBigTree))

}

GO<-function(Ntotal=1000, N=5, Ntransf=5, birth=1, death=0, plot=TRUE) {
	TREES<-DOTHESIMULATION(Ntotal=Ntotal, N=N, Ntransf=Ntransf, birth=birth, death=death, plot=plot)
	RFDIST<-COMPUTERF(TREES)
	return(RFDIST)
}

################################################################################
###                                                                          ###
### THE SIMULATIONS ARE BELOW. THE OUTPUT IS WRITTEN IN FILE "allresults"    ###
### (uncomment to execute)                                                   ###
###                                                                          ###
################################################################################


# cpt1<-0
# for (n in c(50, 100, 500, 1000)) {
# 	for (T in seq(2,20, by=3)) {
# 		print(paste(n, T, sep=" - "))
# 		X<-NULL
# 		X<-replicate(100,GO(10000,n,T,1,0,FALSE), simplify=array)
# 		cat(c(n, T, X, "\n"), file="allresults", append=TRUE)
# 	}
# }


#################################################################################################
###                                                                                           ###
### THE PLOTING FUNCTION IS BELOW. IT READS THE FILE "allresults" AND DO THE PLOT WITH GGPLOT ###
### (uncomment to execute)                                                                    ###
###                                                                                           ###
#################################################################################################

# require(ggplot2)

## read and format results
# res<-readLines("allresults")
# RES<-sapply(res, function(x) as.numeric(strsplit(x," ")[[1]]), simplify=FALSE, USE.NAMES=FALSE)
# nbspecies<-unlist(lapply(RES, function(x) x[1]))
# nbtransf<-unlist(lapply(RES, function(x) x[2]))
# perctopodiff<-unlist(lapply(RES, function(x) sum(x[3:length(x)]>0)/100))
# df<-data.frame(nbspecies,nbtransf, perctopodiff)

## do the plot
# ggplot(df, aes(x=factor(nbspecies), y=perctopodiff, fill=factor(nbtransf)), color=factor(nbspecies)) + geom_bar(stat="identity", position=position_dodge()) + xlab("Percentage of species sampled (out of 10000)")+ylab("Proportion of simulations where\ninduced and sampled gene trees differ topologically") + ylim(0,1) + labs(fill="Number of simulated transfers\nto non-ghost recipients")+scale_fill_brewer(palette="Spectral") + scale_x_discrete(labels = c("0.5%", "1%","5%","10%"))+theme(axis.text=element_text(size=11), text=element_text(size=13))

