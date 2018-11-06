library(tidyverse)
#library(knitr)
install.packages("matrixStats")
library(matrixStats)

########## Load the files ##############
setwd("B:/MauroTutino/ConcordanceRate")

# Load Sequencing data
sequencing_SNPs<-read.table("sequencing_vcfSubset_no_filterDP.txt",header = T, comment.char = '&', check.names=FALSE)

######## FILTER THE FILE ##########
sequencing_SNPs<-sequencing_SNPs[,! colnames(sequencing_SNPs) %in% c("ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]

# Make rownames with CHROM_POS
# In the sequencing ...
rownames(sequencing_SNPs) <- paste(sequencing_SNPs$`#CHROM`, sequencing_SNPs$POS, sep="_")
sequencing_SNPs<-sequencing_SNPs[,! colnames(sequencing_SNPs) %in% c("#CHROM","POS")]


# Extract the dp field for the ALT and REF allele together from the sequencing
sequencing_dp1<-sapply(sequencing_SNPs, function(x) lapply(strsplit(as.character(x), split="[:,]"), '[[', 2)) 


# Extract the dp field for the ALT and REF allele together from the sequencing
sequencing_dp2<-sapply(sequencing_SNPs, function(x) lapply(strsplit(as.character(x), split="[:,]"), '[[', 3)) 

# Extract the GT field from the sequencing 
sequencing_gt<-sapply(sequencing_SNPs, function(x) lapply(strsplit(as.character(x), ":"), '[[', 1)) 

rm(sequencing_SNPs)

Ref_dp<-list()
Alt_dp<-list()
i=0
for(columns in 1:ncol(sequencing_gt)){
  for(rows in 1:nrow(sequencing_gt)){
    if(as.character(sequencing_gt[rows,columns]) == "0/1"){
      i=i+1
      Ref_dp[i]<-sequencing_dp1[rows,columns]
      Alt_dp[i]<-sequencing_dp2[rows,columns]
    }
  }
}

mean(as.numeric(Alt_dp))
mean(as.numeric(Ref_dp))


rm(sequencing_dp1)
rm(sequencing_dp2)
rm(sequencing_gt)

test1<-do.call(rbind, Ref_dp)
test1<-as.data.frame(test1)

test1$Alt<-as.numeric(Alt_dp)

t.test(as.numeric(as.character(test1$V1)),as.numeric(test1$Alt))

###
'data:  as.numeric(as.character(test1$V1)) and as.numeric(test1$Alt)
t = 19.524, df = 2404000, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  11.33427 13.86384
sample estimates:
  mean of x mean of y 
42.84884  30.24979 
'
###