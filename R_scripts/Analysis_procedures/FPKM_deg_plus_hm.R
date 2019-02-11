################ For DE genes to be used as a list for a heatmap from FPKM data ###############

### This scripts feeds in data tables, adjusts the P-value using the Benjamini & Hoffman method, and filters the results for signficance  then ####


####### Here we feed in the tables and name them based on the information given. 
####### You may add tables and change there names if you like
#######

Data_D4_S4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_DOCA_4_vs_Sham_4_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_D8_D4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_DOCA_8_vs_DOCA_4_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_D8_I8 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_DOCA_8_vs_ITF_8_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_D8_S4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_DOCA_8_vs_Sham_4_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_D8_S8 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_DOCA_8_vs_Sham_8_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_I8_D4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_ITF_8_vs_DOCA_4_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_I8_S4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_ITF_8_vs_Sham_4_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_I8_S8 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_ITF_8_vs_Sham_8_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_S8_D4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_Sham_8_vs_DOCA_4_exprs_matrix.txt",header = TRUE, sep = '\t')

Data_S8_S4 = read.delim("C:/Users/Harrison/Desktop/project_folder/full_output/McKinsey_rnaseq_Sham_8_vs_Sham_4_exprs_matrix.txt",header = TRUE, sep = '\t')


####### Now we are adjusting the P-value for more accurate significance
#######

Data_D4_S4["P_Adjust"] = p.adjust(Data_D4_S4[,4], method = "BH")

Data_D8_D4["P_Adjust"] = p.adjust(Data_D8_D4[,4], method = 'BH')

Data_D8_I8["P_Adjust"] = p.adjust(Data_D8_I8[,4], method = 'BH')

Data_D8_S4["P_Adjust"] = p.adjust(Data_D8_S4[,4], method = 'BH')

Data_D8_S8["P_Adjust"] = p.adjust(Data_D8_S8[,4], method = 'BH')

Data_I8_D4["P_Adjust"] = p.adjust(Data_I8_D4[,4], method = 'BH')

Data_I8_S4["P_Adjust"] = p.adjust(Data_I8_S4[,4], method = 'BH')

Data_I8_S8["P_Adjust"] = p.adjust(Data_I8_S8[,4], method = 'BH')

Data_S8_D4["P_Adjust"] = p.adjust(Data_S8_D4[,4], method = 'BH')

Data_S8_S4["P_Adjust"] = p.adjust(Data_S8_S4[,4], method = 'BH')



####### Now we want to filter for significant genes (P_value < .05 and LOG2Fold change >|1|) 
#######

DE_D4_S4 = Data_D4_S4[Data_D4_S4$P_Adjust < .05 & abs(Data_D4_S4$LOG2_FOLD_CHANGE) >1,]

DE_D8_D4 = Data_D8_D4[Data_D8_D4$P_Adjust < .05 & abs(Data_D8_D4$LOG2_FOLD_CHANGE) >1,]

DE_D8_I8 = Data_D8_I8[Data_D8_I8$P_Adjust < .05 & abs(Data_D8_I8$LOG2_FOLD_CHANGE) >1,]

DE_D8_S4 = Data_D8_S4[Data_D8_S4$P_Adjust < .05 & abs(Data_D8_S4$LOG2_FOLD_CHANGE) >1,]

DE_D8_S8 = Data_D8_S8[Data_D8_S8$P_Adjust < .05 & abs(Data_D8_S8$LOG2_FOLD_CHANGE) >1,]

DE_I8_D4 = Data_I8_D4[Data_I8_D4$P_Adjust < .05 & abs(Data_I8_D4$LOG2_FOLD_CHANGE) >1,]

DE_I8_S4 = Data_I8_S4[Data_I8_S4$P_Adjust < .05 & abs(Data_I8_S4$LOG2_FOLD_CHANGE) >1,]

DE_I8_S8 = Data_I8_S8[Data_I8_S8$P_Adjust < .05 & abs(Data_I8_S8$LOG2_FOLD_CHANGE) >1,]

DE_S8_D4 = Data_S8_D4[Data_S8_D4$P_Adjust < .05 & abs(Data_S8_D4$LOG2_FOLD_CHANGE) >1,]

DE_S8_S4 = Data_S8_S4[Data_S8_S4$P_Adjust < .05 & abs(Data_S8_S4$LOG2_FOLD_CHANGE) >1,]



####### Here we make the row names their own column in each table if not already done.
####### Naming doesnt matter since you willl merge them together anyway 

df <- data.frame(names = row.names(DE_D4_S4), DE_D4_S4)
df1 <- data.frame(names = row.names(DE_D8_D4), DE_D8_D4)
df2 <- data.frame(names = row.names(DE_D8_I8), DE_D8_I8)
df3 <- data.frame(names = row.names(DE_D8_S4), DE_D8_S4)
df4 <- data.frame(names = row.names(DE_D8_S8), DE_D8_S8)
df5 <- data.frame(names = row.names(DE_I8_D4), DE_I8_D4)
df6 <- data.frame(names = row.names(DE_I8_S4), DE_I8_S4)
df7 <- data.frame(names = row.names(DE_I8_S8), DE_I8_S8)
df8 <- data.frame(names = row.names(DE_S8_D4), DE_S8_D4)
df9 <- data.frame(names = row.names(DE_S8_S4), DE_S8_S4)



####### Now we need to have only the gene names as heatmap input
#######

dr1 = df[-c(2:6)]
dr2 = df1[-c(2:6)]
dr3 = df2[-c(2:6)]
dr4 = df3[-c(2:6)]
dr5 = df4[-c(2:6)]
dr6 = df5[-c(2:6)]
dr7 = df6[-c(2:6)]
dr8 = df7[-c(2:6)]
dr9 = df8[-c(2:6)]
dr10 = df9[-c(2:6)]



####### Now we need to combine all DE genes into one file as input for a clustered heatmap/GSEA 
#######

m1 = merge(dr1[1], dr2[1], all=TRUE, no.dups=TRUE)
m2 = merge(m1[1], dr3[1], all=TRUE, no.dups=TRUE)
m3 = merge(m2[1], dr4[1], all=TRUE, no.dups=TRUE)
m4 = merge(m3[1], dr5[1], all=TRUE, no.dups=TRUE)
m5 = merge(m4[1], dr6[1], all=TRUE, no.dups=TRUE)
m6 = merge(m5[1], dr7[1], all=TRUE, no.dups=TRUE)
m7 = merge(m6[1], dr8[1], all=TRUE, no.dups=TRUE)
m8 = merge(m7[1], dr9[1], all=TRUE, no.dups=TRUE)
m9 = merge(m8[1], dr10[1], all=TRUE, no.dups=TRUE)

write.table(m9, 'DE_genes.txt', col.names=FALSE)


####### remember to use in terminal for noncoding RNA removed if there is any 
#######

grep -vE "(MIR|SNOR|LINC|ERCC)" m9 > m9_filtered
