library(tidyverse)

########## Load the files ##############
setwd("B:/MauroTutino/ConcordanceRate")
# Load genotyping data
genotype_SNPs<-read.table("genotyping_vcfSubset.txt", header = T, comment.char = '&', check.names=FALSE)
# Load Sequencing data
sequencing_SNPs<-read.table("sequencing_vcfSubset_no_filterDP.txt",header = T, comment.char = '&', check.names=FALSE)
# Load Concordance_rate by SNP
Conc_rate<-read.table("concordance_after_fix_per_SNP.txt", header = T)

# Load Info score
info_score<-read.table("MAAS_ImputedGenotypes_12Apr2017InfoFileAllSequencedSNPs.txt", header = T)

################ FORMAT ##################
# Change rownames to Chr_pos
Conc_rate$chr_pos<-paste(Conc_rate$chr,Conc_rate$pos,sep = "_")
#Conc_rate$chr<-paste("chr",Conc_rate$chr,sep="")
rownames(Conc_rate)<-paste("chr",Conc_rate$chr_pos,sep="")

rownames(sequencing_SNPs)<-paste(sequencing_SNPs$`#CHROM`,sequencing_SNPs$POS, sep="_")


# Find duplicated SNPs in the Genotyping Info file adn remove them
duplicated_snps<-info_score$Merge[duplicated(info_score$Merge)]
info_score<-info_score[!info_score$Merge %in% duplicated_snps,]
# Rename the rownames to match the other files
rownames(info_score)<-paste("chr",info_score$Merge,sep="")

# Filter sequencing SNPs to have the same as Conc_rate
sequencing_SNPs<-sequencing_SNPs[rownames(sequencing_SNPs) %in% rownames(Conc_rate),]

############  GET AF AND DP FROM SEQUENCING DATA ##############
# Get the INFO column from sequencing
Conc_rate$info<-sequencing_SNPs$INFO
# Get the allele frequency
Conc_rate$AF<-lapply(strsplit(as.character(Conc_rate$info), "[=,;]"), '[[', 4)
# Get the overall depth
Conc_rate$DP<-lapply(strsplit(as.character(Conc_rate$info), "[=,;]"), '[[', 13)
#Get the mapping quality (MQ)
Conc_rate$MQ<-lapply(strsplit(as.character(Conc_rate$info), "[=,;]"), '[[', 19)
#Get the mapping quality (QD)
Conc_rate$QD<-lapply(strsplit(as.character(Conc_rate$info), "[=,;]"), '[[', 23)

# Add the imputation info scores
Conc_rate<-left_join(Conc_rate, info_score, by = c("chr_pos" = "Merge"))


# Format and re-order for plotting
c_rate<-as.data.frame(Conc_rate)
c_rate$AF<-as.numeric(c_rate$AF)
c_rate<-as.data.frame(lapply(c_rate, unlist))
c_rate<-c_rate[order(c_rate$AF),]
c_rate$AF<-round(as.numeric(as.character(c_rate$AF)),4)
c_rate$AF<-factor(c_rate$AF,levels = unique(c_rate$AF))

c_rate$Conc_rate<-as.numeric(as.character(c_rate$Conc_rate))
c_rate$Conc_rate[c_rate$Conc_rate == 1] <- NA



# Plot allele frequency vs concordance rate
subset(c_rate,!is.na(c_rate$Conc_rate)) %>%
ggplot(aes(x=as.numeric(as.character(AF)), y=Conc_rate, colour = info.y)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept =0.05, col="RED", lty=4, size=1) +
  geom_vline(xintercept =0.95, col="RED", lty=4, size=1) +
  geom_vline(xintercept =0.01, col="ORANGE", lty=4, size=1) +
  geom_vline(xintercept =0.99, col="ORANGE", lty=4, size=1) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  scale_x_continuous(breaks = c(0,0.01,0.05,0.1,0.2,0.3,0.7,0.8,0.9,0.95,0.99,1.0)) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title="Sequencing-Genotyping Concordance Rate By \nSNP Allele Frequency",
       x="SNP Allele Frequency",
       y="Percent Concordance") 


# Plot Allele frequency vs overall depth
c_rate<-c_rate[order(as.numeric(as.character(c_rate$DP))),]
c_rate$DP<-factor(c_rate$DP,levels = unique(c_rate$DP))

ggplot(c_rate, aes(x=as.numeric(as.character(AF)), y=log10(as.numeric(as.character(DP))/965))) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept =0.05, col="RED", lty=2, size=1) +
  geom_vline(xintercept =0.95, col="RED", lty=2, size=1) +
  geom_vline(xintercept =0.01, col="ORANGE", lty=4, size=1) +
  geom_vline(xintercept =0.99, col="ORANGE", lty=4, size=1) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title="SNP Allele Fequency By \nSequencing Depth",
       x="SNP Allele Fequency",
       y="Log10(Average Sequencing Depth Per SNP Per Sample)")


# Plot Info score vs concordance rate for imputed data(N = 3635)
pdf("Concordance_rate_at_het_sites_by_imputationScore_and_AlleleFreq.pdf", width = 14, height = 11, )

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0) %>%
ggplot(aes(x=as.numeric(as.character(info.y)), y=Conc_rate, colour= log10(as.numeric(as.character(DP))/965))) +
  geom_point() +
  #geom_smooth() +
  scale_color_gradientn(breaks = c(0,0.01,0.05,0.1,0.2,0.3,0.7,0.8,0.9,0.95,0.99,1.0), 
                        colours = c("red","yellow","green","lightblue","darkblue",
                                    "darkblue","lightblue","green","yellow","red"),
                        values=c(0,0.01,0.05,0.1,0.2,0.3,0.7,0.8,0.9,0.95,0.99,1.0),
                        name = "Allele frequency") +
  xlim(0.69, 1.0) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title=paste("Sequencing-Genotyping Concordance Rate At Het. sites",
                   "(N= ", nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0)), ")",
"\nBy Imputation Info Score and Allele Frequency"),
       x="Imputation Rsqr",
       y="Percent Concordance") +
  guides(colour = guide_colorbar(barwidth = 2, barheight = 30))

dev.off()


# Plot Info score vs concordance rate for imputed data coloured by seq. depth(N = 3635)
pdf("Concordance_rate_at_het_sites_by_imputationScore_and_DP.pdf", width = 14, height = 11, )

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0) %>%
  ggplot(aes(x=as.numeric(as.character(info.y)), y=Conc_rate, colour= log10(as.numeric(as.character(DP))/965))) +
  geom_point() +
  #geom_smooth() +
  scale_color_gradientn(breaks = c(0,0.6,0.84,1.17,1.47,2,2.5,3),
                        colours = c("red","yellow","green","lightblue","darkblue"),
                        name = "Average Read Depth",
                        labels = c(round(10^0),round(10^0.6), round(10^0.84), round(10^1.17), round(10^1.47), 10^2.0, round(10^2.5),round(10^3.0))) +
  xlim(0.69, 1.0) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title=paste("Sequencing-Genotyping Concordance Rate At Het. sites",
                   "(N= ", nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0)), ")",
                   "\nBy Imputation Info Score and Sequencing Depth"),
       x="Imputation Rsqr",
       y="Percent Concordance") +
  guides(colour = guide_colorbar(barwidth = 2, barheight = 30))

dev.off()


# Only look at genotyped SNPs (N = 585)
pdf("ConcordanceRate_at_genotyped_HetSites_by_alleleFreq_and_Seq_depth.pdf", width = 18, height = 11, )

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2) %>%
  ggplot(aes(x=as.numeric(as.character(AF)), y=Conc_rate, colour= log10(as.numeric(as.character(DP))/965))) +
  geom_point() +
  scale_color_gradientn(breaks = c(0,0.6,0.84,1.17,1.47,2,2.5,3),
                        colours = c("red","yellow","green","lightblue","darkblue"),
                        name = "Average Read Depth",
                        labels = c(round(10^0),round(10^0.6), round(10^0.84), round(10^1.17), round(10^1.47), 10^2.0, round(10^2.5),round(10^3.0))) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  scale_x_continuous(breaks = c(0,0.05,0.1,0.2,0.3,0.7,0.8,0.9,0.95,1.0)) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title=paste("Sequencing-Genotyping Concordance Rate At Genotyped Het. sites (N =",
                   nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2)),")",
"\nBy Allele Frequency And Sequencing Depth"),
       x="Allele Frequency",
       y="Percent Concordance") +
  guides(colour = guide_colorbar(barwidth = 2, barheight = 30))

dev.off()


#Calculate the mean concordance rate at genotyped SNPs, sequenced at depth > 30X
mean(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2 & as.numeric(as.character(DP))/965 > 30)$Conc_rate)
#Calculate the mean concordance rate at imputed SNPs, sequenced at depth > 30X
mean(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0 & as.numeric(as.character(DP))/965 > 30)$Conc_rate)

# Calculate the average concordance rate at different depths
# for gentoyped SNPs
mat_gen <- matrix(, nrow = length(seq(10,100,by=10)) , ncol = 3)
i=0
for(depth in seq(10,100,by=10)){
  if(depth < 100){
    i=i+1
    mat_gen[i,1]<-depth
    mat_gen[i,2]<-nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2 & as.numeric(as.character(DP))/965 > depth & as.numeric(as.character(DP))/965 < depth+10))
    mat_gen[i,3]<-mean(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2 & as.numeric(as.character(DP))/965 > depth & as.numeric(as.character(DP))/965 < depth+10)$Conc_rate)
  }
  else {
  i=i+1
  mat_gen[i,1]<-depth
  mat_gen[i,2]<-nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2 & as.numeric(as.character(DP))/965 > depth ))
  mat_gen[i,3]<-mean(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2 & as.numeric(as.character(DP))/965 > depth )$Conc_rate)
  }
}
# Plot
ggplot(as.data.frame(mat_gen), aes(x=V1, y=V3))+
  geom_line() +
  geom_point()+
  geom_text(aes(label=V2),hjust=-0.4, vjust=1.5) +
  scale_y_continuous(breaks = c(0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00)) +
  scale_x_continuous(breaks = c(seq(10,100, by = 10)),
                     labels=c( "10X","20X","30X",
                               "40X", "50X", "60X", "70X",
                               "80X", "90X", ">100X")) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title="Average per-SNP concordance for each bin depth \nat genotyped heterozygous sites",
       x="Sequencing Depth of Coverage",
       y="Per-SNP Concordance") 



# Calculate the average concordance rate at different depths
# For imputed SNPs
mat_imp <- matrix(, nrow = length(seq(10,100,by=10)) , ncol = 3)
i=0
for(depth in seq(10,100,by=10)){
  if(depth < 100){
    i=i+1
    mat_imp[i,1]<-depth
    mat_imp[i,2]<-nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0 & as.numeric(as.character(DP))/965 > depth & as.numeric(as.character(DP))/965 < depth+10))
    mat_imp[i,3]<-mean(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0 & as.numeric(as.character(DP))/965 > depth & as.numeric(as.character(DP))/965 < depth+10)$Conc_rate)
  }
  else {
    i=i+1
    mat_imp[i,1]<-depth
    mat_imp[i,2]<-nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0 & as.numeric(as.character(DP))/965 > depth ))
    mat_imp[i,3]<-mean(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0 & as.numeric(as.character(DP))/965 > depth )$Conc_rate)
  }
}
# Plot
ggplot(as.data.frame(mat_imp), aes(x=V1, y=V3))+
  geom_line() +
  geom_point()+
  geom_text(aes(label=V2),hjust=-0.4, vjust=1.5) +
  scale_y_continuous(breaks = c(0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00)) +
  scale_x_continuous(breaks = c(seq(10,100, by = 10)),
                     labels=c( "10X","20X","30X",
                               "40X", "50X", "60X", "70X",
                               "80X", "90X", ">100X")) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title="Average per-SNP concordance for each bin depth \nat imputed heterozygous sites",
       x="Sequencing Depth of Coverage",
       y="Per-SNP Concordance") 




# Plot Info score vs concordance rate but filter 
# by AF 5% and colour by read depth (N = 3220)

pdf("ConcordanceRate_at_Imputed_HetSites_by_Rsqr_and_Seq_depth_AF05.pdf", width = 18, height = 11, )

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0 & as.numeric(as.character(AF)) > 0.05 & as.numeric(as.character(AF)) < 0.95) %>%
  ggplot(aes(x=as.numeric(as.character(info.y)), y=Conc_rate, colour= log10(as.numeric(as.character(DP))/965))) +
  geom_point() +
  #geom_smooth() +
  scale_color_gradientn(breaks = c(0,0.6,0.84,1.17,1.47,2,2.5,3),
                        colours = c("red","yellow","green","lightblue","darkblue"),
                        name = "Average Read Depth",
                        labels = c(round(10^0),round(10^0.6), round(10^0.84), round(10^1.17), round(10^1.47), 10^2.0, round(10^2.5),round(10^3.0))) +
  xlim(0.69, 1.0) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title=paste("Sequencing-Genotyping Concordance Rate At Het. sites (AF > 0.05)",
                   "(N= ", nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & as.numeric(as.character(AF)) > 0.05 & as.numeric(as.character(AF)) < 0.95)), ")",
                   "\nBy Imputation Info Score and Read Depth"),
       x="Imputation Rsqr",
       y="Percent Concordance") +
  guides(colour = guide_colorbar(barwidth = 2, barheight = 30))

dev.off()


# Plot of sequencing frequency vs imputation expected frequency

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0) %>%
ggplot(aes(x=as.numeric(as.character(AF)), y = exp_freq_a1, colour= log10(as.numeric(as.character(DP))/965))) +
  geom_point() +
  scale_color_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                        name = "log10(Sequencing Depth)",
                        breaks = c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5))
subset(c_rate,!is.na(c_rate$Conc_rate) & type == 0) %>%
ggplot(aes(x=as.numeric(as.character(AF)), y = exp_freq_a1, colour= as.numeric(as.character(info.y)))) +
  geom_point() +
  scale_color_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                        name = "Info Score")



# Plot Conc rate vs allele frequency colour by MQ
# For genotyped SNPs (N = 585)
pdf("ConcordanceRate_at_genotyped_HetSites_by_alleleFreq_and_MQ.pdf", width = 18, height = 11, )

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2) %>%
  ggplot(aes(x=as.numeric(as.character(AF)), y=Conc_rate, colour= MQ)) +
  geom_point() +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  scale_x_continuous(breaks = c(0,0.05,0.1,0.2,0.3,0.7,0.8,0.9,0.95,1.0)) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title=paste("Sequencing-Genotyping Concordance Rate At Genotyped Het. \nSites (N =",
                   nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2)),")",
                   "By Allele Frequency And Mapping Quality"),
       x="Allele Frequency",
       y="Percent Concordance") 

dev.off()

# Plot Conc rate vs allele frequency colour by QD
# For genotyped SNPs (N = 585)
pdf("ConcordanceRate_at_genotyped_HetSites_by_alleleFreq_and_QD.pdf", width = 18, height = 11, )

subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2) %>%
  ggplot(aes(x=as.numeric(as.character(AF)), y=Conc_rate, colour= as.numeric(as.character(QD)))) +
  geom_point() +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  scale_x_continuous(breaks = c(0,0.05,0.1,0.2,0.3,0.7,0.8,0.9,0.95,1.0)) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title=paste("Sequencing-Genotyping Concordance Rate At Genotyped Het. \nSites (N =",
                   nrow(subset(c_rate,!is.na(c_rate$Conc_rate) & type == 2)),")",
                   "By Allele Frequency And Quality_by_Depth"),
       x="Allele Frequency",
       y="Percent Concordance") +
  scale_color_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                        name = "QD score")

dev.off()
