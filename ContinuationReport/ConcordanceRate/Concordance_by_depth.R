#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")

library(tidyverse)
#library(knitr)
library(matrixStats)

########## Load the files ##############
setwd("B:/MauroTutino/ConcordanceRate")
# Load genotyping data
genotype_SNPs<-read.table("genotyping_vcfSubset.txt", header = T, comment.char = '&', check.names=FALSE)
# Load Sequencing data
sequencing_SNPs<-read.table("sequencing_vcfSubset_no_filterDP.txt",header = T, comment.char = '&', check.names=FALSE)

######## FILTER THE FILES TO BE THE SAME ##########

# Filter out the columns we do not want
genotype_SNPs<-genotype_SNPs[,! colnames(genotype_SNPs) %in% c("ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]
sequencing_SNPs<-sequencing_SNPs[,! colnames(sequencing_SNPs) %in% c("ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]

# Filter the files to have the same columns
sequencing_SNPs<-sequencing_SNPs[, colnames(sequencing_SNPs) %in% colnames(genotype_SNPs)]
genotype_SNPs<-genotype_SNPs[, colnames(genotype_SNPs) %in% colnames(sequencing_SNPs)]

# Make rownames with CHROM_POS
# In the sequencing ...
rownames(sequencing_SNPs) <- paste(sequencing_SNPs$`#CHROM`, sequencing_SNPs$POS, sep="_")
sequencing_SNPs<-sequencing_SNPs[,! colnames(sequencing_SNPs) %in% c("#CHROM","POS")]
# And in the genotyping
rownames(genotype_SNPs) <- paste(genotype_SNPs$`#CHROM`, genotype_SNPs$POS, sep="_")
genotype_SNPs<-genotype_SNPs[,! colnames(genotype_SNPs) %in% c("#CHROM","POS")]

# Filter the files to have the same rows
sequencing_SNPs<-sequencing_SNPs[rownames(sequencing_SNPs) %in% rownames(genotype_SNPs),]
genotype_SNPs<-genotype_SNPs[rownames(genotype_SNPs) %in% rownames(sequencing_SNPs),]

# Make sure that the files are in the same order
sequencing_SNPs<-sequencing_SNPs[,match(colnames(genotype_SNPs), colnames(sequencing_SNPs))]


######### GET DEPTH AND GENOTYPE FROM SEQUENCING FILE ##############

# Extract the dp field for the ALT and REF allele together from the sequencing
sequencing_dp<-sapply(sequencing_SNPs, function(x) lapply(strsplit(as.character(x), split=":"), '[[', 3)) 

# Extract the GT field from the sequencing 
sequencing_gt<-sapply(sequencing_SNPs, function(x) lapply(strsplit(as.character(x), ":"), '[[', 1)) 

######### PERFORM THE CALCULATIONS #################################

# Remove this matrix to save space in memory
rm(sequencing_SNPs)

# Calculate the concordance-by-depth per sample
# mat <- matrix(, nrow = 20, ncol = ncol(sequencing_gt))
# for(columns in 1:ncol(sequencing_gt)){
#   for(dp in 1:20){
#     matching=0
#     mismatching=0
#     for(rows in 1:nrow(sequencing_gt)){
#       if(as.character(genotype_SNPs[rows,columns]) == "0/1"){
#         if(as.character(genotype_SNPs[rows,columns]) == as.character(sequencing_gt[rows,columns]) & as.numeric(sequencing_dp[rows,columns]) == as.numeric(dp)){
#           matching=matching+1
#         } 
#         else if(as.character(genotype_SNPs[rows,columns]) != as.character(sequencing_gt[rows,columns]) & as.character(sequencing_gt[rows,columns]) != "./." & as.numeric(sequencing_dp[rows,columns]) == as.numeric(dp)){
#           mismatching=mismatching+1
#         }
#       }
#     }
#     print(matching)
#     print(mismatching)
#     perc_matching=matching/(matching+mismatching)
#     print(perc_matching)
#     mat[dp, columns] = perc_matching
#   }
# }


# Calculate the concordance-by-depth per sample
# This code does the same as the one above but it reduces the
# computational time by A LOT!
# For each column it returns a logical vector for matching and mismatching
# Then I can use sum to sum up all the TRUE values which are = 1
mat <- matrix(, nrow = 20, ncol = ncol(sequencing_gt))
for(columns in 1:ncol(sequencing_gt)){
  for(dp in 1:20){
    matching=0
    mismatching=0
    
    # Compare the column of the matrices and 
    # Return a logical vector of the same length as the column
    # For the matching
    matchmat <- as.character(genotype_SNPs[,columns]) == "0/1" & genotype_SNPs[, columns] == sequencing_gt[, columns] & sequencing_dp[, columns] == dp
    # And mismatching
    mismatchmat <- as.character(genotype_SNPs[,columns]) == "0/1" & genotype_SNPs[, columns] != sequencing_gt[, columns] & sequencing_gt[, columns] != "./." & sequencing_dp[, columns] == dp     
    # Calculate the percentage of matching
    # When we sum a logical vector TRUE = 1 and FALSE = 0
    perc_matching=sum(matchmat)/(sum(matchmat)+sum(mismatchmat))
    print(perc_matching)
    mat[dp, columns] = perc_matching
  }
  
}

# Sum up all the samples for each depth
# mat2 <- matrix(, nrow = 20, ncol = 1)
# for(dp in 1:20){
#   print(paste("depth is:", dp))
#   matching=0
#   mismatching=0
#   for(columns in 1:ncol(sequencing_gt)){
#     for(rows in 1:nrow(sequencing_gt)){
#       if(as.character(genotype_SNPs[rows,columns]) == "0/1"){
#         if(as.character(genotype_SNPs[rows,columns]) == as.character(sequencing_gt[rows,columns]) & as.numeric(sequencing_dp[rows,columns]) == as.numeric(dp)){
#           matching=matching+1
#           print(paste("row is:",rows))
#           print(paste("column is:",columns))
#         } 
#         else if(as.character(genotype_SNPs[rows,columns]) != as.character(sequencing_gt[rows,columns]) & as.character(sequencing_gt[rows,columns]) != "./." & as.numeric(sequencing_dp[rows,columns]) == as.numeric(dp)){
#           mismatching=mismatching+1
#         }
#       }
#     }
#   }
#   print(matching)
#   print(mismatching)
#   perc_matching=matching/(matching+mismatching)
#   print(perc_matching)
#   mat2[dp, 1] = perc_matching
# }


# Calculate the total number of heterozygous sites at each depth
# mat_het <- matrix(, nrow = 20, ncol = 1)
# for(dp in 1:20){
#   Het=0
#   for(columns in 1:ncol(genotype_SNPs)){
#     for(rows in 1:nrow(genotype_SNPs)){
#       if(as.character(genotype_SNPs[rows,columns]) == "0/1" & as.numeric(sequencing_dp[rows,columns]) == as.numeric(dp)){
#         Het=Het+1
#         print(paste("Depth = ", dp," Het Count:", Het))
#       }
#     }
#     
#   }
#   mat_het[dp, 1] = Het/ncol(genotype_SNPs)
# }

# Calculate the total number of heterozygous sites at each depth
# This code does the same as the one above but it is much quicker
mat_het <- matrix(, nrow = 20, ncol = 1)
for(dp in 1:20){
  
  matching<-genotype_SNPs == "0/1" & sequencing_dp == dp
  mat_het[dp, 1] = sum(matching)/ncol(genotype_SNPs)
  print(paste("Depth = ", dp," Het Count:", sum(matching)))
  }

### TRANSFORM TO DATAFRAME AND WRITE TO FILE

#write.table(mat_altONLY_per_sample,"Concordance_by_depth_ALT_ONLY_per_sample.txt", col.names = F, sep = "\t")

#Transform the matrix in dataframe for plotting
df<-as.data.frame(mat)
#Add row means
df$RMean<-rowMeans(mat, na.rm = T)
# Add column depth of sequencing
df$dp<-rownames(df)
df$dp<-factor(df$dp,levels = df$dp)
df$SD<-rowSds(mat, na.rm = T)
# Change the order of the columns
df<-df[,c(ncol(df),ncol(df)-1,ncol(df)-2,1:(ncol(df)-3))]
colnames(df)[4:ncol(df)]<-colnames(genotype_SNPs)
df$HetCount<-mat_het[,1]
df<-df[,c(ncol(df),1:(ncol(df)-1))]
# Write the data to file
write.table(df,"Concordance_by_depth_ALT_and_REF_per_sample.txt", col.names = T, row.names = F, sep = "\t")

### PLOTTING

# Concordance rate by Sequencing Depth
# Plot the data except the depth == 1
df[-c(1),] %>%
ggplot(aes(x=dp, y=RMean,group=1)) +
  geom_line(colour="#000099") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=pmax(RMean-SD,0), ymax=pmin(RMean+SD,1)), width=1,
                position=position_dodge(0.05), colour="#000099") +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels=c( "0%","10%","20%","30%",
                               "40%", "50%", "60%", "70%",
                               "80%", "90%", "100%")) +
  theme_bw() +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title="Sequencing-Genotyping Concordance Rate By \nSequencing Depth For Heterozygous Sites",
       x="Sequencing Read Depth",
       y="Percent Concordance")

# Heterozygous count by Sequencing Depth
# Plot the data except the depth == 1
df[-c(1),] %>%
  ggplot(aes(x=dp, y=HetCount,group=1)) +
  geom_line(colour="#000099") +
  geom_point(size=3)  +
  theme_bw() +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  labs(title="Average Number Of Heterozygous Sites By \nSequencing Depth",
       x="Sequencing Read Depth",
       y="Average Number Of Heterozygous sites")


