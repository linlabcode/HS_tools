### GGplot - Boxplot ###
library(ggplot2)
library(data.table)
options(scipen=999)
######################## make vectors of genes to use & format table ########################

########## Here we just have names of genes for particular genesets -- could be enriched only or total but probably enriched only
##########

x1 <- jechlinger_mesenchymal_list
x2 <- poola_invasive_genes
x3 <- LIAN_inflam_genes

#######  when using mice genes the IDs are lowercase where as most things require captilizing all the letters 
#######

McKinsey_rnaseq_all_fpkm_means[,1] <- toupper(McKinsey_rnaseq_all_fpkm_means[,1])

####### Create a table that includes only desired samples (columns)
#######

new_means_exprs1 <- McKinsey_rnaseq_all_fpkm_means[,c(1,2,3,4)]

####### Change the name of the IDs for relevance
#######

colnames(new_means_exprs1)[colnames(new_means_exprs1)=="row.names"] <- "Genes"

####### Now create a new table only grabbing rows whos entries are selected from another table
#######

new_means_exprs2_Poola <- new_means_exprs1[is.element(new_means_exprs1$Genes, x2$V1),] 


####### No we want to make the IDs the rownames for following plots. (making a new table is optional but good for backtracking)
#######

Poola_invasive <- new_means_exprs2_Poola

Poola_invasive2 = Poola_invasive[,-1]
rownames(Poola_invasive2) <- Poola_invasive[,1]



######################################## Actual plotting ############################################

####### Create box plot of values and adjust or scale for visibility 
#######


boxplot(Poola_invasive2, ylim=c(y1=-1,y2=20),
  col=("gray50"), main="Poola Invasive Boxplots", xlab= "Sample Groups", ylab="Poola Invasive Gene Expression ")




##########################################Box plots for log2 of the values ########################################

####### calculate the log2 of values for some reason. (scaling for visibility?)

####### new table with log2 values
#######

Poola_invasive_log2 <- log2(Poola_invasive2)


####### log2(values/other_values) = fold change?
#######

dvs = log2(Poola_invasive$DOCA_8/Poola_invasive$Sham_8)
dvi = log2(Poola_invasive$DOCA_8/Poola_invasive$ITF_8)
ivs = log2(Poola_invasive$ITF_8/Poola_invasive$Sham_8)

####### The alternate log2 boxplots (not sure but leave it)
#######

boxplot(dvs,dvi,ivs,names=c('Doca vs. Sham','Doca vs. ITF', 'ITF vs. Sham'), ylim=c(y1=-2,y2=3),col=("gray50"), 
        main="Poola Invasive LOG2 Boxplots", xlab= "Sample Groups", 
        ylab="LOG2 FC of Poola Invasive Gene Expression ")


####### t.test for P_values? (just leave it for now)
#######

t.test(dvs,dvi)
t.test(dvs,ivs)
t.test(dvi,ivs) 
