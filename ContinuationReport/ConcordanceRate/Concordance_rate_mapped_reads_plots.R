library(tidyverse)

# Plot before filtering sequencing vcf for mindepth 7 reads

setwd("B:/MauroTutino/QC_and_preprocessing")
mapped_reads_no_chr6<-read.table("Total_mapped_reads_no_chr6.txt", header = F)
colnames(mapped_reads_no_chr6)<-c("Sample_ID","M_reads")


setwd("B:/MauroTutino/ConcordanceRate")
Conc_rate<-read.table("Conc_rate.txt", header = T)

setwd("B:/MauroTutino/eQTL/")
IDs_link<-read.table("Sequecing_Cytokine_names.txt", header = F)
colnames(IDs_link)<-c("SeqID","Sample_ID")
mapped_reads_IDs<-merge(mapped_reads_no_chr6, IDs_link, by.x="Sample_ID", by.y="SeqID")

plot_concRate_mappedReads <- merge(x=mapped_reads_IDs, y=Conc_rate, by.x="Sample_ID.y", by.y = "sample")

#Excluding Missing genotypes
ggplot(plot_concRate_mappedReads, 
       aes(log2(plot_concRate_mappedReads$M_reads), plot_concRate_mappedReads$Non.missing, 
       color=ifelse(((plot_concRate_mappedReads$Non.missing < 85)),"A", "B") #names of levels are arbitrary, putting in alphabetical order
       )) +
  geom_point() + #drawing a scatterplot
  scale_color_manual(guide=FALSE, values=c("red", "black")) + #turn off the legend, define the colors
  labs(title="Concordance Rate and # Mapped reads correlation",
       x="Log2(# Mapped Reads) - Not Chr6",
       y="% Sequencing-Genotyping Concordance Rate") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) 

# Including missing genotypes
ggplot(plot_concRate_mappedReads, 
       aes(log2(plot_concRate_mappedReads$M_reads), plot_concRate_mappedReads$Overall, 
           color=ifelse(((plot_concRate_mappedReads$Non.missing < 85)),"A", "B") #names of levels are arbitrary, putting in alphabetical order
       )) +
  geom_point() + #drawing a scatterplot
  scale_color_manual(guide=FALSE, values=c("red", "black")) + #turn off the legend, define the colors
  labs(title="Concordance Rate and # Mapped reads correlation\n including missing genotypes ",
       x="Log2(# Mapped Reads) - Not Chr6",
       y="% Sequencing-Genotyping Concordance Rate") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"))


# Mapped reads vs missing genotypes
ggplot(plot_concRate_mappedReads, 
       aes(log2(plot_concRate_mappedReads$M_reads), plot_concRate_mappedReads$missing_genotype_First_analysis_primary_assembly_only_only_SNP_HardFiltered_AC2_biallelic_db151_ann_noHLA_PASS_MissInd_CallR_hg19, 
           color=ifelse(((plot_concRate_mappedReads$Non.missing < 85)),"A", "B") #names of levels are arbitrary, putting in alphabetical order
       )) +
  geom_point() + #drawing a scatterplot
  scale_color_manual(guide=FALSE, values=c("red", "black")) + #turn off the legend, define the colors
  labs(title="Missing Genotype calls and # Mapped reads \ncorrelation",
       x="Log2(# Mapped Reads) - Not Chr6",
       y="Sequencing Missing Genotype Calls") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"))




#################################
# # Plot after filtering sequencing vcf for mindepth 7 reads

setwd("B:/MauroTutino/QC_and_preprocessing")
mapped_reads_no_chr6<-read.table("Total_mapped_reads_no_chr6.txt", header = F)
colnames(mapped_reads_no_chr6)<-c("Sample_ID","M_reads")


setwd("B:/MauroTutino/ConcordanceRate")
Conc_rate<-read.table("Conc_rate_after_filtering_sequencing.txt", header = T)

setwd("B:/MauroTutino/eQTL/")
IDs_link<-read.table("Sequecing_Cytokine_names.txt", header = F)
colnames(IDs_link)<-c("SeqID","Sample_ID")
mapped_reads_IDs<-merge(mapped_reads_no_chr6, IDs_link, by.x="Sample_ID", by.y="SeqID")

plot_concRate_mappedReads <- merge(x=mapped_reads_IDs, y=Conc_rate, by.x="Sample_ID.y", by.y = "sample")

#Excluding Missing genotypes
ggplot(plot_concRate_mappedReads, 
       aes(log2(plot_concRate_mappedReads$M_reads), plot_concRate_mappedReads$Non.missing, 
           color=ifelse(((plot_concRate_mappedReads$Non.missing < 85)),"A", "B") #names of levels are arbitrary, putting in alphabetical order
       )) +
  geom_point() + #drawing a scatterplot
  scale_color_manual(guide=FALSE, values=c("red", "black")) + #turn off the legend, define the colors
  labs(title="Concordance Rate and # Mapped reads correlation\n excluding missing genotypes",
       x="Log2(# Mapped Reads) - Not Chr6",
       y="% Sequencing-Genotyping Concordance Rate") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) 

# Including missing genotypes
ggplot(plot_concRate_mappedReads, 
       aes(log2(plot_concRate_mappedReads$M_reads), plot_concRate_mappedReads$Overall, 
           color=ifelse(((plot_concRate_mappedReads$Non.missing < 85)),"A", "B") #names of levels are arbitrary, putting in alphabetical order
       )) +
  geom_point() + #drawing a scatterplot
  scale_color_manual(guide=FALSE, values=c("red", "black")) + #turn off the legend, define the colors
  labs(title="Concordance Rate and # Mapped reads correlation\n including missing genotypes ",
       x="Log2(# Mapped Reads) - Not Chr6",
       y="% Sequencing-Genotyping Concordance Rate") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"))

 