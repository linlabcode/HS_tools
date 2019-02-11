##################### Here we have a script that creates barplots(bargraphs?) of selected genes in different groups ################


#E2f4 bar plot

####### Loading raw table data into a new table object
#######

mean_table <- read.delim('~/projects/McKinsey_givinostat/Project_folder (updated)/output/McKinsey_rnaseq_all_fpkm_exprs_raw.txt', header=TRUE, sep = '\t')


####### If we were using normalized data 
#######

#norm_table <- read.delim('/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_all_fpkm_exprs_norm.txt', header=TRUE, sep = '\t')


####### Global Stuff (This was already here. Must just mean overall expression from the sample?)
#######

####### Here we input the genes of choice
#######

gene_name = c('Rpl34-ps1', 'Ass1','Rps27', 'H2-Q8', 'H2-Q7', 'Pou3f1', 'Psca', 'Mt2', 'Nmrk2', 'Amd1', 'H2-Q6', 'Gm1987', 'Insl3', 'Dynlt1f', 'Sap25', 'Hbb-bs', 'Top2a', 'Ptx3' )


####### here we make a pdf file of the data
#######

pdf('~/projects/McKinsey_givinostat/Project_folder (updated)/test_expr_bargraphs.pdf')


####### Not sure whats happening here
####### Par is for parameters and mfrow how something to do with splitting data
#######

par(mfrow=c(2,2))


####### Ok this function does all of the work
#######

####### For each object in the list, search for it in the 'mean' table (might need to change the name)
#######

for( gene in gene_name){ 
  gene_row = which(rownames(mean_table) == gene)
  
  
  ####### just to check
  #######
  
  print(gene_row)
  print(gene)

  
  ####### Here we create a 'container' objects for the following code
  ####### (not sure why there needs to be an if statement)
  #######
  
  if(length(gene_row) != 0){ 
     plot_DOCA = c()
     plot_SHAM =  c()
     plot_ITF =  c()
    
      
      
     ####### Here we gather the replicates into the 'containers'
     #######
     
     plot_DOCA = c(mean_table[gene_row,1], mean_table[gene_row,2], mean_table[gene_row,3])
     plot_SHAM =  c(mean_table[gene_row, 8], mean_table[gene_row,9], mean_table[gene_row, 10], mean_table[gene_row, 11])
     plot_ITF =  c(mean_table[gene_row, 4], mean_table[gene_row,5], mean_table[gene_row, 6], mean_table[gene_row, 7])
     
  
     
     ####### Now we find the mean of these seperate replicate groups
     #######
     
     mean_DOCA = mean(plot_DOCA)
     mean_SHAM = mean(plot_SHAM)
     mean_ITF = mean(plot_ITF)
     
     
     ####### Standard Deviation
     #######
  
     sd_DOCA = sd(plot_DOCA)
     sd_SHAM = sd(plot_SHAM)
     sd_ITF = sd(plot_ITF)
    
     
     ####### Standard error for error bars?
     #######
     
     se_DOCA = sd_DOCA/(sqrt(length(plot_DOCA)))
     se_SHAM = sd_SHAM/(sqrt(length(plot_SHAM)))
     se_ITF = sd_ITF/(sqrt(length(plot_ITF)))
  
  
     
     ####### No idea whats happening from here just yet but it works
     #######
     
     maxy=max(mean_DOCA,mean_SHAM,mean_ITF) + 2.5*max(sd_DOCA,sd_SHAM,sd_ITF) 
     print(maxy) 

     centers <- barplot(height = c(mean_DOCA,mean_SHAM,mean_ITF),
                     beside = true, las = 2,
                     ylim = c(0, maxy),
                     cex.names = 0.75, 
                     main = gene,
                     ylab = "Expression Signal",
                     border = "black", axes = TRUE)
  

  # 45 degree string rotation
     text(x = centers, y = par("usr")[3] - 1, srt = 45,
       adj = 1, labels = c('DOCA', 'SHAM','ITF'), xpd = TRUE)
  
     segments(centers, c(mean_DOCA,mean_SHAM,mean_ITF) - c(se_DOCA,se_SHAM, se_ITF) * 2, centers,
           c(mean_DOCA,mean_SHAM, mean_ITF) + c(se_DOCA,se_SHAM,se_ITF) * 2, lwd = 1.5)
  
     arrows(centers, c(mean_DOCA,mean_SHAM,mean_ITF) - c(se_DOCA, se_SHAM, se_ITF) * 2, centers,
         c(mean_DOCA,mean_SHAM, mean_ITF) + c(se_DOCA,se_SHAM,se_ITF) * 2, lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
  }
}

dev.off()


