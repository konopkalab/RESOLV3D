##======================================================
runRESOLV3D<-function(ifn, ofp, itr, expcutoff, genm, key, backg)
	{
	# Defaults	
	#ifn <- "example_rpkm.txt"
	#ofp <- "example1"
	#itr <- 1000
	#expcutoff <- 0.5
	#genm <- "hg19"
	#key <- "geneSymbol"
	#backg <- "background_genes.txt"

	suppressWarnings(suppressMessages(library(RESOLV3D)))
	cat("===> Running package for", ifn, "with", itr, "iterations...", "\n")

	# Running data validation...
	newdt1 <- validatedata(ifn, ofp, as.numeric(expcutoff))
	#write.table(newdt1, paste(ofp, "filteredData.txt", sep="_"), row.names=T, col.names=T, quote=F, sep="\t")

	# Resolving single cells...
	newdt2 <- clustercells(newdt1, ofp, as.numeric(itr))
	save(newdt2, file="example1_main.RData")

	# Extracting cluster-specific genes
	extractgenes(newdt2)

	# Gene Ontology (GO) analysis
	#cluster1GO <- getsetgo("genes_cluster_1_postCor_genes_top50.txt", backg, ofp, genm, key)
	#save(cluster1GO, file="example1_go.RData")
	for(i in 1:max(newdt2[,1]))
		{
		cat("\tCluster ", i, "\n")
		gotemp<-getsetgo(paste("genes_cluster_", i, "_postCor_genes_top50.txt", sep=""), backg, ofp, genm, key)
		write.table(gotemp, paste("go_cluster_", i, "_genes_top50.txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		}

	# Creating violin plot and bar plot for gene of interest
	#plotgene(gene_of_interest, ofp)
	#plotgene("STMN2", ofp)
	#plotgene("VIM", ofp)
	#plotgene("DLX2", ofp)

	# Creating violin plot and bar plot for list of genes of interest
	# example for "genes_cluster_1_postCor_genes_top50.txt"
	#genelist<-scan("genes_cluster_1_postCor_genes_top50.txt", what="")
	#for(gene in genelist)
	#	{
	#	plotgene(gene, ofp)
	#	cat(gene, "\n")
	#	}
	
	# Creating animation for clustered cells in 3D space...
	animate(ofp)

	cat("===> Finished running for", ifn, "with", itr, "iterations...", "\n")
	}

##======================================================
cat("===> Checking input parameters...", "\n")
args = commandArgs(TRUE)
if(length(args) < 7)
	{
	stop("Invalid parameters. Please refer user manual or tutorial.")
	}
if(length(args) == 7)
	{
	cat("===> 7 input parameters detected...", "\n")
	#cat("===> Input file: ", args[[1]], "\n")
	#cat("===> Output prefix: ", args[[2]], "\n")
	#cat("===> Iterations: ", args[[3]], "\n")
	#cat("===> Cutoff: ", args[[4]], "\n")
	#cat("===> genome: ", args[[5]], "\n")
	#cat("===> key: ", args[[6]], "\n")
	#cat("===> background genes: ", args[[7]], "\n")
	runRESOLV3D(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]], args[[6]], args[[7]])
	}
##======================================================

