#### For table merging ####

####### simple one line merge for an entire table with or without duplicates
#######

merged_D4S4 <- merge(D4S4_S4_gsea, D4S4_D4_gsea, all=TRUE, no.dups = FALSE)











#### This reads in the tables as data frames ###

control2_vs_doca1     <- read.delim('control2_vs_doca1.txt', stringsAsFactors = F, head=FALSE)

Control2_vs_control1  <- read.delim('Control2_vs_control1.txt', stringsAsFactors = F, head=FALSE)

doca1_vs_control1     <- read.delim('doca1_vs_control1.txt', stringsAsFactors = F, head=FALSE)

doca2_vs_control1     <- read.delim('doca2_vs_control1.txt', stringsAsFactors = F, head=FALSE)

ITF_vs_control1       <- read.delim('ITF_vs_control1.txt', stringsAsFactors = F, head=FALSE)

ITF_vs_doca1          <- read.delim('ITF_vs_doca1.txt', stringsAsFactors = F, head=FALSE)

### This merges the first column containing genes together into one table ###

merge1 <- merge(control2_vs_doca1[1], Control2_vs_control1[1], all=TRUE, no.dups = TRUE)

merge2 <- merge(merge1[1], doca2_vs_control1[1], all=TRUE, no.dups = TRUE)

merge3 <- merge(merge2[1], ITF_vs_control1[1], all=TRUE, no.dups = TRUE)

merge4 <- merge(merge3[1], ITF_vs_doca1[1], all=TRUE, no.dups = TRUE)

write.table(All_DE_genes, 'All_DE_genes.txt', col.names = FALSE)
