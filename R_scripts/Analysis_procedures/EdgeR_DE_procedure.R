########## My Edge R differential expression analysis procedure ##########
##########################################################################
########## This should be checked, validated, and updated as needed ######
##########################################################################

library(edgeR)

####### New table from filtered raw counts
#######

x <- genes_filtered


####### remove unnecessary rows (if any)
#######

x <- x[ -c(12:20) ]


####### group replicates
#######

group <- factor(c(1,1,1,2,2,2,2,3,3,3,3))


####### makes your deglist (This is how edgeR stores data in a simple list object)
#######

y <- DGEList(counts=x, group=group)


####### filter out low expressed genes (this is to remove genes with values less than 1 in at least 3/4 replicates)
####### CPM is used as it acounts for library size. (1 cpm = 6-7 genes expressed and 5-10 is needed for effective expression)
#######

keep <- rowSums(cpm(y)>1) >=3

y <- y[keep, , keep.lib.sizes=FALSE]


####### normalize read counts (not always needed. For samplic specific effects)
####### (please re-read info in documentation)

y <- calcNormFactors(y)


####### create design?? (no explanation but necessary)
#######

design <- model.matrix(~group)



####### Estimate dispersion (calcs the likelihood by conditioning the total count for each tag)
####### (quantile-adjusted conditional maximum likelihood (qCML))
#######

y <- estimateDisp(y,design)

et <- exactTest(y)
topTags(et)

#######  MDS-plot - distance corresponding to BCM #####
####### This is a quick way to check for similarity among replicates but could be wrong (not necessary for DE)
#######

plotMDS(y)

pch <- c(1,1,1)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=3)


########## for quasi-likelihood F-tests  ############

#################### First step #####################

####### Generalized linear Model quasi-likelihood Fitting
####### Primes the F-test
#######


fit <- glmQLFit(y, design)

#################### group 1 vs 2 ###################

####### Here we are performing the QL F-test
####### It reflects the uncertainty in estimating the dispersion for each gene
#######

qlf.I8_vs_D8 <- glmQLFTest(fit, coef = 2)


####### Quick view of DE genes (FDR < .05)
#######

summary(de <- decideTestsDGE(qlf.I8_vs_D8))


####### selecting the results and putting in an object?
#######

I8vsD8 <- topTags(qlf.I8_vs_D8, n= 16)#length(qlf.I8_vs_D8))


####### Here we are viewing the results of the test
#######

All_I8vsD8 <- I8vsD8$table

All_I8vsD8

#################### group 1 vs 3 ###################

qlf.S8_vs_D8 <- glmQLFTest(fit, coef=3)

S8vsD8 <- topTags(qlf.S8_vs_D8, n=length(qlf.S8_vs_D8))

All_S8vsD8 <- S8vsD8$table

summary(de <- decideTestsDGE(qlf.S8_vs_D8))
#################### group 2 vs 3 ###################

qlf.S8_vs_I8 <- glmQLFTest(fit, contrast=c(0,-1,1))

S8vsI8 <- topTags(qlf.S8_vs_I8, n=length(qlf.S8_vs_I8))
All_S8vsI8 <- S8vsI8$table

et <- exactTest(y, pair=c(2,3))
topTags(et, n=57)

summary(de <- decideTestsDGE(qlf.S8_vs_I8))






####### The rest is extra stuff used post analysis. (might remove later)
#######

######### for a heatmap? ############################


FDR5_S8vsD8$genes <- row.names(FDR5_S8vsD8)

FDR5_I8vsD8$genes <- row.names(FDR5_I8vsD8)




merged_DE_genes <- merge(FDR5_I8vsD8, FDR5_S8vsD8, all=TRUE, no.dups = FALSE)
#merged_DE_genes <- merged_DE_genes[-42,]
deduped.data <- unique( merged_DE_genes[ ,6] )

HM_genes <- as.data.frame(deduped.data)


write.table(merged_DE_genes, "~/projects/McKinsey_givinostat/911_merged_DEG_70.txt", sep="\t") 



#logcpm <- cpm(y, prior.count=2, log=TRUE)

#et <- exactTest(y, pair=c(1,2))
#topTags(et, n=49)
#plotMD(y, column=c(4,5,6)

       