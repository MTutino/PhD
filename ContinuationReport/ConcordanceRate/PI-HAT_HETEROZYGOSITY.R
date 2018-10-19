library(tidyverse)

#Set working directory
setwd("B:/MauroTutino/ConcordanceRate/")

#Load IBD data
ibd<- read.table("ibd.genome", header=T)

#Make the plot
ggheatmap <- ggplot(ibd, aes(x=ibd$FID1, y=ibd$FID2, fill=ibd$PI_HAT)) + 
  geom_tile(color="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="LD Correlation\nR2") +
  scale_x_discrete(position="top") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_fixed()

print(ggheatmap)

#Modify the plot before printing
p<-ggheatmap + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.2),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 14, barheight = 3,
                               title.position = "top", title.hjust = 0.5)) 


# Plot the distribution of specific pi_hat
ibd[ibd$FID1 == "14350" | ibd$FID2 == 14350,] %>%
ggplot(aes(x=PI_HAT)) + geom_histogram(binwidth = .005)

ibd[ibd$FID1 == "11265" | ibd$FID2 == 11265,] %>%
  ggplot(aes(x=PI_HAT)) + geom_histogram(binwidth = .005)

nrow(ibd[ibd$FID1 == 11265 | ibd$FID2 == 11265 & ibd$PI_HAT > 0,])/nrow(ibd[ibd$FID1 == 11265 | ibd$FID2 == 11265,])
nrow(ibd[ibd$FID1 == 14350 | ibd$FID2 == 14350 & ibd$PI_HAT > 0,])/nrow(ibd[ibd$FID1 == 14350 | ibd$FID2 == 14350,])



# Plot the distribution of all pi_hat > 0.5
ggplot(ibd[ibd$PI_HAT > 0,],aes(x=PI_HAT)) + 
  geom_histogram(binwidth = .01) +
  labs(title="PI-HAT count",
       x="PI-HAT",
       y="Count") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) 

# Create a data frame containing only the maximum pi_hats per sample
mylist<-list()
i=1
for (smpl in unique(ibd$FID1)){
  mylist[[i]]<-ibd[ibd$FID1 == smpl & ibd$PI_HAT == max(ibd$PI_HAT[ibd$FID1 == smpl]),][1,1:ncol(ibd)]
  i=i+1
  }
df <- do.call("rbind",mylist) #combine all vectors into a matrix
df<-as.data.frame(df)
colnames(df)<- colnames(ibd)

# Get the pi-hat of only the samples that have a match in genotype 
df_05<-df[df$FID1[df$PI_HAT > 0] %in% Conc_rate$sample,]
df_05<-df_05[df_05$FID1 %in% Conc_rate$sample,]

# Plot the distribution of the maximum pi_hat
ggplot(df_05,aes(x=df_05$PI_HAT)) + geom_histogram(binwidth = .01)
# Plot the distribution of the maximum pi_hat, excluding those with no match
ggplot(df_05[df_05$PI_HAT > 0,],aes(x=PI_HAT)) + 
  geom_histogram(binwidth = .01) +
  labs(title="Max PI-HAT count",
       x="Max PI-HAT/Sample",
       y="Count") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) 

df[df$FID1 == 35798, ]
df[df$FID1 == 11265, ]
df[df$FID1 == "11265G", ]
df[df$FID1 == 14350, ]
df[df$FID1 == "14350G", ]

setwd("B:/MauroTutino/ConcordanceRate")
Conc_rate<-read.table("Conc_rate.txt", header = T)



write.table(df_05, "Sequencing_pi_hat.txt", sep = "\t", col.names = T, row.names = F)

##### PI-HAT vs Mapped Readsn ######################

# Load mapped reads counts
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

# Load maximum PI-HATs
setwd("B:/MauroTutino/ConcordanceRate")
PI_HAT<-read.table("Sequencing_pi_hat.txt", header = T)
plot_Pi_Hat_mappedReads <-merge(mapped_reads_IDs, PI_HAT, by.x="Sample_ID.y", by.y="FID1")

# PLot PI-HAT vs Mapped reads
plot_Pi_Hat_mappedReads %>%
ggplot(aes(x=log10(M_reads), y=PI_HAT)) +
  geom_point()

################ Heterozygosity #######################
#Set working directory
setwd("B:/MauroTutino/ConcordanceRate/")

#Load the heterozygosity
Het<- read.table("Sequencing_heterozygosity.txt", header=T)

#Subset to only the sequencing samples
Het<-Het[!grepl("G", Het$FID),]

# Calculate mean heterozygosity (N-O)/N : 
#(Total number of non missing genotypes-Observed homozygous calls)/Total number of non missing genotypes
Het$meanHet = (Het$N.NM. - Het$O.HOM.)/Het$N.NM.

#Plot Heterozygosity
ggplot(Het,aes(x=F, y=..count..)) + #names of levels are arbitrary, putting in alphabetical order  
               geom_histogram(binwidth = .01)

#Load missing data info
setwd("B:/MauroTutino/VariantCalling/")
MissInd<-read_tsv("out.imiss",col_names = T)
MissInd<- MissInd %>% arrange(F_MISS) 
MissInd$INDV<-factor(MissInd$INDV,levels = MissInd$INDV[order(-MissInd$F_MISS)])

# Plot Mean heterozygosity vs Missing genotypes
colors  <- densCols(MissInd$F_MISS,Het$meanHet)
plot(log10(MissInd$F_MISS+0.000001),Het$meanHet, col=colors, xlim=c(-3,0),ylim=c(0,0.5),pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
abline(h=mean(Het$meanHet)-(3*sd(Het$meanHet)),col="RED",lty=2)
abline(h=mean(Het$meanHet)+(3*sd(Het$meanHet)),col="RED",lty=2)
abline(v=-1.522879, col="RED", lty=2)

# Plot distribution of the Mean heterozygosity +/- 3 SD
ggplot(Het,aes(x=Het$meanHet, y=..count..)) + #names of levels are arbitrary, putting in alphabetical order  
  geom_histogram(binwidth = .001) +
  geom_vline(xintercept =mean(Het$meanHet)-(3*sd(Het$meanHet)), col="RED", lty=2) +
  geom_vline(xintercept=mean(Het$meanHet)+(3*sd(Het$meanHet)), col="RED", lty=2) +
labs(title="Mean heterozygosity Distribution",
     x="Mean Heterozygosity",
     y="Count") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) 



