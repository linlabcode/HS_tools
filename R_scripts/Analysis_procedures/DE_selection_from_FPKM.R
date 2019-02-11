#################### Using FPKM for DE expression analysis on Cufflinks output############################

####### similar to 'finding_DEG' but longer and does not create the gene list for the heatmap
#######


library("dplyr")
options=scipen(999)



#############ITF vs SHAM (8 weeks) ##########################


###### Negative values #########

####### just a new object to work with and retain the orignal 
#######

adjusted_P_I8vsS8 <- McKinsey_rnaseq_ITF_8_vs_Sham_8_exprs_matrix


####### Here we use the bonferroni and Hoffman method of adjusting P_value
####### The method may be too stringent
#######

adjusted_P_I8vsS8$P_VALUE <- p.adjust(adjusted_P_I8vsS8$P_VALUE, method = 'BH')



####### simply sort the values by Fold change
#######

adjusted_P_I8vsS8 <- arrange(adjusted_P_I8vsS8, LOG2_FOLD_CHANGE)


####### Removing ones after a cutoff 
#######

adjusted_P_I8vsS8_neg <- adjusted_P_I8vsS8[-c(150:12321),]


####### sorting the remaining values by P_values
#######

adjusted_P_I8vsS8_neg <- arrange(adjusted_P_I8vsS8_neg, P_VALUE)


####### remove genes that are over the cutoff values
#######

adjusted_P_I8vsS8_neg <- adjusted_P_I8vsS8_neg[-c(33:149),]


####### Create a new table containing resulting genes
#######

I8vsS8_neg_genes_under_.2 <- as.data.frame(adjusted_P_I8vsS8_neg$row.names) 


####### Write the table to a csv file
#######

write.csv(I8vsS8_neg_genes_under_.2, file= "~/I8vsS8_neg_under_20perc.csv")


##############################################################

###### Positive Values #########



adjusted_P_I8vsS8_pos <- arrange(adjusted_P_I8vsS8, desc(LOG2_FOLD_CHANGE)) 


adjusted_P_I8vsS8_pos <- adjusted_P_I8vsS8_pos[-c(94:12321),]


adjusted_P_I8vsS8_pos <- arrange(adjusted_P_I8vsS8_pos, P_VALUE)

adjusted_P_I8vsS8_pos <- adjusted_P_I8vsS8_pos[-c(29:93),]

I8vsS8_pos_genes_under_.2 <- as.data.frame(adjusted_P_I8vsS8_pos$row.names) 

write.csv(I8vsS8_pos_genes_under_.2, file= "~/I8vsS8_pos_under_20perc.csv")

##############################################################

############# DOCA vs ITF (8 weeks) ##########################

##### Neg values ########

adjusted_P_D8vsI8 <- McKinsey_rnaseq_DOCA_8_vs_ITF_8_exprs_matrix


adjusted_P_D8vsI8$P_VALUE <- p.adjust(adjusted_P_D8vsI8$P_VALUE, method = 'BH')


adjusted_P_D8vsI8 <- arrange(adjusted_P_D8vsI8, LOG2_FOLD_CHANGE)


adjusted_P_D8vsI8_neg <- adjusted_P_D8vsI8[-c(106:12321),]


adjusted_P_D8vsI8_neg <- arrange(adjusted_P_D8vsI8_neg, P_VALUE)

adjusted_P_D8vsI8_neg <- adjusted_P_D8vsI8_neg[-c(2:105),]

I8vsS8_neg_genes_under_.2 <- as.data.frame(adjusted_P_I8vsS8_neg$row.names) 

write.csv(I8vsS8_neg_genes_under_.2, file= "~/I8vsS8_neg_under_20perc.csv")

###### Pos values #######

adjusted_P_D8vsI8_pos <- arrange(adjusted_P_D8vsI8, desc(LOG2_FOLD_CHANGE)) 


adjusted_P_D8vsI8_pos <- adjusted_P_D8vsI8_pos[-c(77:12321),]


adjusted_P_D8vsI8_pos <- arrange(adjusted_P_D8vsI8_pos, P_VALUE)

adjusted_P_D8vsI8_pos <- adjusted_P_D8vsI8_pos[-c(),]

I8vsS8_pos_genes_under_.2 <- as.data.frame(adjusted_P_I8vsS8_pos$row.names) 

write.csv(I8vsS8_pos_genes_under_.2, file= "~/I8vsS8_pos_under_20perc.csv")

##############################################################

###########  DOCA vs SHAM (8 weeks) ##########################

##### neg values #############

adjusted_P_D8vsS8 <- McKinsey_rnaseq_DOCA_8_vs_Sham_8_exprs_matrix

adjusted_P_D8vsS8$P_VALUE <- p.adjust(adjusted_P_D8vsS8$P_VALUE, method = 'BH')

adjusted_P_D8vsS8_neg <- arrange(adjusted_P_D8vsS8, LOG2_FOLD_CHANGE)

adjusted_P_D8vsS8_neg <- adjusted_P_D8vsS8_neg[-c(154:12321),]

adjusted_P_D8vsS8_neg <- arrange(adjusted_P_D8vsS8_neg, P_VALUE)

##### pos values ############

adjusted_P_D8vsS8_pos <- arrange(adjusted_P_D8vsS8, desc(LOG2_FOLD_CHANGE))
