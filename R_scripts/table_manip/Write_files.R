###### one-liners for creating files ######



####### Here we convert a table into a text file while keeping it in a table format sperated by tabs
#######

write.table(merged_DE_genes, "~/projects/McKinsey_givinostat/911_merged_DEG_70.txt", sep="\t") 


####### Here we create csv that retains its table format (just another way)
#######

write.csv(I8vsS8_neg_genes_under_.2, file= "~/I8vsS8_neg_under_20perc.csv")


####### this is for making PDF files
#######

pdf('~/projects/McKinsey_givinostat/Project_folder (updated)/test_expr_bargraphs.pdf')