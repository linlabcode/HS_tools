############ this script is to simply make a table isolating only differentially expressed genes from EdgeR's analysis ###########

####### There is an easier way of doing this but this is what you have now and it works
#######

library("dplyr")
options(scipen=999) ##### this is for making the values not be in scientific notation

#### Using EdgeR #############################################

#### table formatting ######

d <- All_S8vsI8

d <- cbind(rownames(d), data.frame(d, row.names=NULL))

d2 <- d
########

g <- All_I8vsD8

g <- cbind(rownames(g), data.frame(g, row.names=NULL))

########

f <- All_S8vsD8

f <- cbind(rownames(f), data.frame(f, row.names=NULL))

#############ITF vs SHAM (8 weeks) ##########################

###### Negative values #########

FDR_I8vsS8 <- d

FDR_I8vsS8 <- FDR_I8vsS8[-c(152:12899),]

FDR_I8vsS8 <- arrange(FDR_I8vsS8, logFC)

FDR_I8vsS8_2 <- FDR_I8vsS8

FDR_I8S8_pos <- FDR_I8vsS8[-c(1:101),]

FDR_I8S8_pos <- arrange(FDR_I8S8_pos, desc(logFC))

FDR_I8S8_neg <- FDR_I8vsS8[-c(72:151),]

write.csv(FDR_I8S8_neg, file= "~/I8vsS8_neg2.csv")

write.csv(FDR_I8S8_pos, file= "~/I8vsS8_pos2.csv")


##############################################################

############# DOCA vs ITF (8 weeks) ##########################

FDR_I8vsD8 <- g

FDR_I8vsD8 <- g[-c(37:12899),]

FDR_I8vsD8 <- arrange(FDR_I8vsD8, desc(logFC))

FDR_I8vsD8_2 <- FDR_I8vsD8

FDR_I8vsD8_pos <- FDR_I8vsD8[-c(1:10),]

FDR_I8vsD8_neg <- FDR_I8vsD8_2[-c(5:17),]

write.csv(FDR_I8vsD8, file="~/I8D8_all.csv")

##############################################################

###########  DOCA vs SHAM (8 weeks) ##########################

FDR_S8vsD8 <- f

FDR_S8vsD8 <- FDR_S8vsD8[-c(42:12899),]

FDR_S8vsD8 <- arrange(FDR_S8vsD8, desc(logFC))

FDR_S8vsD8_neg <- FDR_S8vsD8[-c(5:7),]

write.csv(FDR_S8vsD8, file="~/S8D8_all.csv")


####### similarity I vs S (8 weeks)###############

I_S8_similar <- d

I_S8_similar <- I_S8_similar[-c(152:12899),]

I_S8_similar <- arrange(I_S8_similar, logFC)

I_S8_similar <- I_S8_similar[-c(1:49),]

I_S8_similar <- I_S8_similar[-c(31:80),]

write.csv(I_S8_similar, file="~/IS_similar_all.csv")
