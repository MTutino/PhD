library(tidyverse)

setwd("B:/MauroTutino/VariantCalling/")
 
###### Missing genotype calls per individual #############

MissInd<-read_tsv("out.imiss",col_names = T)
MissInd<-read_tsv("out.imiss_before_mindep7_filtering.txt",col_names = T)

MissInd<- MissInd %>% arrange(F_MISS) 

MissInd$INDV<-factor(MissInd$INDV,levels = MissInd$INDV[order(-MissInd$F_MISS)])


table(MissInd$F_MISS > 0.55)
mean(MissInd$F_MISS)*100
median(MissInd$F_MISS)*100

ggplot(data=MissInd, aes(MissInd$F_MISS)) + 
  geom_histogram(breaks=seq(0.00, 1.00, by =0.05), 
                 col="red",
                 aes(fill=..count..)) +
  labs(title="Frequency of missing genotype calls",
       x="Frequency of missing genotype calls",
       y="# of samples") +
  scale_fill_gradient("# of samples") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"))
 
######### Missing genotype per SNP ####################
setwd("B:/MauroTutino/ConcordanceRate")
# Load the file
sequencing_SNPs<-read.table("sequencing_vcfSubset_no_filterDP.txt",header = T, comment.char = '&', check.names=FALSE)

# Add rownames
rownames(sequencing_SNPs) <- paste(sequencing_SNPs$`#CHROM`, sequencing_SNPs$POS, sep="_")
# Delete the columns we do not need
sequencing_SNPs<-sequencing_SNPs[,! colnames(sequencing_SNPs) %in% c("ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]
sequencing_SNPs<-sequencing_SNPs[,! colnames(sequencing_SNPs) %in% c("#CHROM","POS")]

# Extract the GT field from the sequencing 
sequencing_gt<-sapply(sequencing_SNPs, function(x) lapply(strsplit(as.character(x), ":"), '[[', 1)) 

# Extract the number of missing genotype calls
no_call<-rowSums(sequencing_gt == "./.")
no_call<-data.frame(no_call)
# Calculate franction
no_call$Frac<-(no_call$no_call/ncol(sequencing_gt))*100

# Plot the Fraction of Missing Genotype Calls Per SNP
ggplot(no_call, aes(x=Frac) ) + 
  geom_histogram(binwidth = 0.5, 
                 col="black", fill="red") + 
  labs(title="Percentage of Missing Genotype Calls Per SNP") +
  labs(x="Percentage Of Missing Genotype Calls", y="Number Of SNPs") +
  theme_bw()+ 
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  scale_x_continuous(breaks=seq(0,max(no_call$Frac)+0.5,1),
                     labels=c( "0%","1%","2%","3%",
                               "4%", "5%"))

