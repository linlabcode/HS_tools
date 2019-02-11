############### NES plot tool (Normalized Enrichment score) ##########################

####### Here we take data from gsea output and determine
####### particular sets that are enriched in one group vs another
#######

####### Input them into new objects to retain original copy
#######

D4S4_S4_gsea <- D4S4_SHAM_gsea_report

D4S4_D4_gsea <- D4S4_DOCA_gsea_report

####### Merge into one 
#######

merged_D4S4 <- merge(D4S4_S4_gsea, D4S4_D4_gsea, all=TRUE, no.dups = FALSE)



################ select for specifc rows #################


####### Were doing this to sort and grab only data we want
####### There may be an easier way
#######

merged_inflammatory4 <- merged_D4S4[-c(25:3302),]

merged_invasive4 <- merged_D4S4[-c(54:3302),]
merged_invasive4 <- merged_invasive4[-c(1:24),]

merged_mesenchymal4 <- merged_D4S4[-c(71:3302),]
merged_mesenchymal4 <- merged_mesenchymal4[-c(1:53),]

##########################################################

############### plotting the NES plot ####################

####### Here we simply create the plot, no packages or libraries needed
####### For some reason the dashes are crossing the borders rn :/ 
#######


plot(merged_mesenchymal4$NES, merged_mesenchymal4$FDR.q.val, main="Mesenchymal_D4vsS4",
     
     xlab="NES", ylab="FDR", pch=19, col="yellow", xlim = c(-3,+3), ylim = c(-.05,+1.25))

abline(h=0.05,lty=2,col='black')

pdf(file = "~/projects/McKinsey_givinostat/ALL_D8D4.pdf")

