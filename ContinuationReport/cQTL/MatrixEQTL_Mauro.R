#install.packages("MatrixEQTL")
#install.packages("tidyverse")
library(MatrixEQTL)
library(tidyverse)

setwd("C:/Users/mdxasmt5/Dropbox/Systems-Immunology/MauroPhD/Lijing_cytokines/nonImputed/")

#Location of files
base.dir=getwd()


#########################################################################
# Decide which genotype file to use

# Genotype file name: RS6967330 ONLY
SNP_file_name = paste(base.dir, "/CDHR3_Dosage_MatrixEQTL.txt", sep="");
snps_location_file_name = paste(base.dir, "/SNP_location_matrixEQTL.txt", sep="");

# Genotype file name: CDHR3 +/- 1KB 
SNP_file_name = paste(base.dir, "/DosageFile_cdhr3_MatrixEQTL.txt", sep="");
snps_location_file_name = paste(base.dir, "/SNP_location_cdhr3_MatrixEQTL.txt", sep="");

# Genotype file name: IL6 +/- 1KB 
#SNP_file_name = paste(base.dir, "/DosageFile_il6_MatrixEQTL.txt", sep="");
#snps_location_file_name = paste(base.dir, "/SNP_location_il6_MatrixEQTL.txt", sep="");

# Genotype file name: All SNPs
SNP_file_name = paste(base.dir, "/DosageFile_AllSNPs_MatrixEQTL.txt", sep="");
snps_location_file_name = paste(base.dir, "/SNP_location_AllSNPs_MatrixEQTL.txt", sep="");

#########################################################################
# IF READING A NEW FILE

# Read file
#rawData<-read.csv("cyto_stim_pairs_year11_raw_LowResponseReplacedWithHalfDtm.csv", header=F)

#Rename columns and rows
#colnames(rawData)<-paste(as.matrix(rawData[1,]),as.matrix(rawData[2,]), sep="_")
#rownames(rawData)<-rawData[,1]

# Get Media
#rawDataMedia<-rawData[,rawData[2,] == "media"]
# Get Stimuli
#rawDataStimulus<-rawData[,!rawData[2,] == "media"]

#Get rid of rows an columns we do not need
#rawDataMedia<-rawDataMedia[-c(1:3),]
#rawDataStimulus<-rawDataStimulus[-c(1:3),-1]

# Transpose
#rawDataMedia<-t(rawDataMedia)
#rownames(rawDataMedia)<-gsub("_media","",rownames(rawDataMedia))
#rawDataStimulus<-t(rawDataStimulus)

# Write in file
#write.csv(rawDataMedia, "Media_year11_raw_LowResponseReplacedWithHalfDtm.csv")
#write.csv(rawDataStimulus, "Stimuli_year11_raw_LowResponseReplacedWithHalfDtm.csv")



###########################################################################
###########################################################################
# Decide which cytokine file to use

#Gene expression for media non imputed Log2
#expression_file_name = paste(base.dir, "/Media_log_nonImputed_MatrixEQTL.csv", sep="");
#gene_location_file_name = paste(base.dir, "/Cytokine_location_matrixEQTL.csv", sep="");

# Gene expression for stimulus non imputed media normalised
expression_file_name = paste(base.dir, "/AllCytokineStimuliPairs_Year11_MediaNormalised_nonImputed_MatrixEQTL.csv", sep="");
gene_location_file_name = paste(base.dir, "/Cytokine_stimulus_location_matrixEQTL.csv", sep="");

#Gene expression for media non imputed Log2 neg0 (NOT USEFULL)
#expression_file_name = paste(base.dir, "/Media_log_nonImputed_Neg0_MatrixEQTL.csv", sep="");
#gene_location_file_name = paste(base.dir, "/Cytokine_location_matrixEQTL.csv", sep="");

### Gene expression for stimulus non imputed Log2 neg0
expression_file_name = paste(base.dir, "/AllCytokineStimuliPairs_Year11_MediaNormalised_nonImputed_Neg0_MatrixEQTL.csv", sep="");
gene_location_file_name = paste(base.dir, "/Cytokine_stimulus_location_matrixEQTL.csv", sep="");

# Gene expression Raw media with 0s
#expression_file_name = paste(base.dir, "/Raw_nonImputed_Media_cytokine.csv", sep="");
#gene_location_file_name = paste(base.dir, "/Cytokine_location_matrixEQTL.csv", sep="");

# Gene expression Raw stimulus with 0s
#expression_file_name = paste(base.dir, "/Raw_nonImputed_stimulus_cytokine.csv", sep="");
#gene_location_file_name = paste(base.dir, "/Cytokine_stimulus_location_matrixEQTL.csv", sep="");

### Gene expression Raw media 0s converted to 1/2 of detection limit of Batch2
expression_file_name = paste(base.dir, "/Media_year11_raw_LowResponseReplacedWithHalfDtm.csv", sep="");
gene_location_file_name = paste(base.dir, "/Cytokine_location_matrixEQTL.csv", sep="");

### Gene expression Raw stimulus 0s converted to 1/2 of detection limit of Batch2
expression_file_name = paste(base.dir, "/Stimuli_year11_raw_LowResponseReplacedWithHalfDtm.csv", sep="");
gene_location_file_name = paste(base.dir, "/Cytokine_stimulus_location_matrixEQTL.csv", sep="");



############################################################################
# LOAD THE CYTOKINES EXPRESSION FILE
# AND FILTER OUT THE CYTOKINES WE ARE NOT INTERESTED IN
expression_cytokines<-read_csv(expression_file_name, col_names = T)
colnames(expression_cytokines)[1]<-"ID"

# EXCLUDE CYTOKINES IL-25 and IL-33
expression_cytokines<-expression_cytokines[!grepl("IL-25", expression_cytokines$ID) & !grepl("IL-33", expression_cytokines$ID),]

# IFN-beta was always 0 for Media, exclude from analysis
expression_cytokines<-expression_cytokines[!grepl("IFN-beta", expression_cytokines$ID),]

# FOR MEDIA-NORMALISED STIMULI ONLY 
#CONVERT NEGATIVE VALUES TO 0
expression_cytokines[expression_cytokines<0]<-0

##### FOR STIMULI #####################
# ALSO EXCLUDE THE NON-SIGNIFICANT CYTOKINES-STIMULUS PAIRS 
list_pairs<-read.csv("List_of_cytokines_to_exclude.csv")
expression_cytokines<-expression_cytokines[!expression_cytokines$ID %in% list_pairs$List.of.all.NO,]

# BATERIAL OR VIRAL CYTOKINE LIST
list_bact<- read.csv("List_of_bacteria_cytokines_pairs.csv")
list_vir<- read.csv("List_of_virus_cytokines_pairs.csv")

# Select the appropriate list
# BACTERIA
expression_cytokines<-expression_cytokines[expression_cytokines$ID %in% list_bact$pro.inf.and.bacteria.pairs,]

# VIRUS
# expression_cytokines<-expression_cytokines[expression_cytokines$ID %in% list_vir$IFN.and.virus.pairs,]

#########################################################################
# Normalisation methods for raw data

#Normalisation is applied on the columns
# We need to transpose the file first
expression_cytokines<-t(expression_cytokines)
colnames(expression_cytokines)<-as.matrix(expression_cytokines[1,])
expression_cytokines<-expression_cytokines[-1,]

# Convert class to numeric
class(expression_cytokines) <- "numeric"

# NORMALISATIONS
# Autoscaling
#expression_cytokines<-as.data.frame(scale(expression_cytokines,center=T,scale=T))

# Paretoscaling
#sample_classes <- NULL
#x <- expression_cytokines
#x.centered <- apply(x, 2, function(x) x - mean(x, na.rm = T))
#x.sc <- apply(x.centered, 2, function(x) x/sqrt(sd(x, na.rm = T)))
#x.sc <- cbind(sample_classes, x.sc)
#expression_cytokines<-as.data.frame(x.sc)

# Log2 transformation
expression_cytokines<-log2(expression_cytokines)


# Flip back and re-format
expression_cytokines<-t(expression_cytokines)
expression_cytokines<-as.data.frame(expression_cytokines)
expression_cytokines$ID <- rownames(expression_cytokines)
expression_cytokines<- expression_cytokines[,c(ncol(expression_cytokines),1:ncol(expression_cytokines)-1)]

#########################################################################
# MODIFY THE SNP file
#
# FORMAT SNPS AND CYTOKINES TO BE THE SAME
# Load the file
SNPs_cytokines<-read_tsv(SNP_file_name)
colnames(SNPs_cytokines)[1]<- "ID"

####  Exclude SNPs in strong Linkage Disequilibrium R2 > 0.8
LD<-read.csv("LD_Exclude.csv", header = F)
SNPs_cytokines<-SNPs_cytokines[! SNPs_cytokines$ID %in% LD$V1,]

# Make sure that the SNPs and Cytokine files have the same columns
SNPs_cytokines<-SNPs_cytokines[colnames(SNPs_cytokines) %in% colnames(expression_cytokines)]
# One file "14055" is not in the sequencing. Filter it out
expression_cytokines<-expression_cytokines[colnames(expression_cytokines) %in% colnames(SNPs_cytokines)]

#################################################################################
####### MERGE SINGLETON HOMOZIGUS WITH HETEROZIGUS ######################

# Flip the dataframe to have SNPs as columns
# This will also convert the dataframe into matrix
SNPs_cytokines<-t(SNPs_cytokines)
#Rename the columns names
colnames(SNPs_cytokines)<-as.matrix(SNPs_cytokines[1,])
colnames(SNPs_cytokines)<-gsub("-","_", colnames(SNPs_cytokines))
SNPs_cytokines<-SNPs_cytokines[-1,]

# This loop will merge the homozigus and heterozigus in the columns that 
# only contain one homozigus  
class(SNPs_cytokines)<- "numeric"
# loop through the columns
for(i in 1:ncol(SNPs_cytokines)){
  # for each column count the number of homozigus reference and alternative
  homalt=sum(as.integer(as.character(SNPs_cytokines[,i])) == 2, na.rm = T);
  homref=sum(as.integer(as.character(SNPs_cytokines[,i])) == 0, na.rm = T); 
  # if only one homozigus is present, then merge it with the heterozigus
  if(homalt==1){
    SNPs_cytokines[SNPs_cytokines[,i]==2, i]<- 1
  }
  # The Hom. Ref. can be 1 in case the ref and alt were flipped in hg38
  if(homref==1){
    SNPs_cytokines[SNPs_cytokines[,i]==0, i]<- 1
  }
}

# Flip and convert back to data.frame
SNPs_cytokines<-t(SNPs_cytokines)
SNPs_cytokines<-as.data.frame(SNPs_cytokines)


# Create the column "ID" and move it to the first position
SNPs_cytokines[,1:ncol(SNPs_cytokines)] <- lapply(SNPs_cytokines[,1:ncol(SNPs_cytokines)], function(x) as.integer(as.character(x)))
SNPs_cytokines$ID <- rownames(SNPs_cytokines)
SNPs_cytokines<- SNPs_cytokines[,c(ncol(SNPs_cytokines),1:ncol(SNPs_cytokines)-1)]
# Make sure that the SNPs and Cytokine files are in the same order
SNPs_cytokines<-SNPs_cytokines[,match(colnames(expression_cytokines), colnames(SNPs_cytokines))]

#####################################################################################
# MatrixEQTl creates a slicedData() from a file.
# We need to wrtie the filtered SNP and Expression files in a file

#Genotype
write.table(SNPs_cytokines, "DosageFile_temp_MatrixEQTL.txt", row.names = F, sep="\t")
SNP_file_name = paste(base.dir, "/DosageFile_temp_MatrixEQTL.txt", sep="");

#Expression
write.table(expression_cytokines, "Cytokines_temp.txt", row.names = F, sep="\t")
expression_file_name = paste(base.dir, "/Cytokines_temp.txt", sep="");

####################################################################################
# Make sure that the covariate files have the same columns as SNPs and Expression

# Covariates file name
# Set to character() for no covariates
#covariates_file_name = character()
covariates_file_name = paste(base.dir, "/Age_11_covariates.txt", sep="");

#We now need to filter out the samples not present in the sequencing data from the covariates file
covariates <- read_tsv(covariates_file_name)
colnames(covariates)[1] <- "ID"
covariates <- covariates[colnames(covariates) %in% colnames(SNPs_cytokines)]
SNPs_cytokines<-SNPs_cytokines[colnames(SNPs_cytokines) %in% colnames(covariates)]
write.table(covariates, "covariates_eqtl.txt", row.names = F, sep = "\t")
covariates<-read_tsv("covariates_eqtl.txt")

# Covariates file name
# Set to character() for no covariates
#covariates_file_name = character()
covariates_file_name = paste(base.dir, "/covariates_eqtl.txt", sep="");

#############################################################################
#SET UP THE PARAMETERS

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05;
pvOutputThreshold_tra = 0.005;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
#errorCovariance = read.table("Sample_Data/errorCovariance.txt");


# Distance for local gene-SNP pairs
cisDist = 1e6;

###########################################################################
## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

################### Filter for MAF  ######################################
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
#Filter for MAF > 5%
snps$RowReorder(maf>0.08);

#########################################################################

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the tab character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the comma character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
#If loading from file
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}



## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.csv(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
#useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelANOVA; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS


me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance,
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = T,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(me)

##########################
# BOXPLOT ######

df<-rbind(SNPs_cytokines, expression_cytokines)
df<-t(df)
df<-as.data.frame(df)
colnames(df)<-as.matrix(df[1,])
colnames(df)<-gsub("-","_", colnames(df))
colnames(df)<-gsub("\\(","_", colnames(df))
colnames(df)<-gsub("\\)","_", colnames(df))
df<-df[-1,]
df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))), row.names = rownames(df))

#Plot IL-15: significant at nominal p-value for RS6967330
subset(df,!is.na(df$chr2_191010309)) %>%
ggplot(aes(x=factor(chr2_191010309), y=IL_8_CXCL8__Fla, fill=as.factor(chr2_191010309))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha=0.3,position = position_jitter(width = 0.25), na.rm=T) +
  scale_y_continuous(breaks=round(seq(min(df$IL_8_CXCL8__Fla, na.rm = T), max(df$IL_8_CXCL8__Fla, na.rm = T), by = 1)))+
  labs(title="IL-8(CXCL8)_Fla vs chr2_191010309",
       x="chr2_191010309",
       y="IL-8(CXCL8)_Fla (Log2)") +
  guides(fill=guide_legend(title="Genotype")) +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"))
  


covar<-t(covariates)
covar<-as.data.frame(covar)
colnames(covar)<-as.matrix(covar[1,])
covar<-covar[-1,]
covar$cgender<-as.integer(as.character(covar$cgender))

# Check Lin. Reg
res.lin<-lm(df$IL_4 ~ df$chr7_106018005+covar$cgender+as.numeric(as.character(covar$viability_11)), na.action = na.exclude)
summary(res.lin)
# Check Anova
res.aov<-aov(df$IL_2 ~ as.factor(df$chr7_106018005)+covar$cgender+covar$viability_11)
summary(res.aov)

#####################################################################################

#Count genotype for the whole cohort
dosage<-read_tsv(paste(base.dir,"/CDHR3_Dosage_MatrixEQTL.txt",sep=""))
dosage<-t(dosage)
colnames(dosage)<-as.matrix(dosage[1,])
dosage<-as.data.frame(dosage)
table(dosage, useNA = "ifany")

#0              1              2  <NA> 
#662            268             25  7 




# Plot the cytokine expression by genotype for significant eQTLs
#Create the dataframe for plotting
df<-rbind(SNPs_cytokines, expression_cytokines)
df<-t(df)
df<-as.data.frame(df)
colnames(df)<-as.matrix(df[1,])
colnames(df)<-gsub("-","_", colnames(df))
df<-df[-1,]
df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))), row.names = rownames(df))
df$chr7_106018005<-as.factor(as.integer(as.character(df$chr7_106018005)))

#Count genotype for the cytokine subgroup
table(df$chr7_106018005, useNA = "ifany")
#   0    1     2  <NA> 
#  213   83    8    2 

#Plot IL-15: significant at nominal p-value for RS6967330
subset(df,!is.na(df$chr7_106018005)) %>% ggplot(aes(x=chr7_106018005, y=IL_15)) + 
  geom_boxplot(aes(fill=chr7_106018005)) + 
  geom_point() +
  labs(title="cQTL: IL_15 vs RS6967330 (Media)",
       x="CDHR3 (RS6967330)",
       y="IL-15 (Log2)")


#Plot IL-2: significant at nominal p-value for RS6967330
subset(df,!is.na(df$chr7_106018005)) %>% ggplot(aes(x=as.factor(chr7_106018005), y=IL_2_LPS)) + 
  geom_boxplot(aes(fill=as.factor(chr7_106018005))) + 
  geom_point() +
  labs(title="cQTL: IL_2 vs RS6967330 (Media)",
       x="CDHR3 (RS6967330)",
       y="IL-2 (Log2)")



Lin_reg<-lm(df$IL_15 ~ df$chr7_106018005, na.action = na.exclude)

####################################################################################

#############################################################
######### FILTER FOR ALLELE COUNT ###########################
dt <- df[,grep("chr", colnames(df))]
dt<-dt[,colSums(dt, na.rm = T) > 15]
SNPs_cytokines<-SNPs_cytokines[SNPs_cytokines$ID %in% colnames(dt),]

#write the sliced data in file
#CDHR3
write.table(SNPs_cytokines, "DosageFile_cdhr3_15_MatrixEQTL_sliced.txt", row.names = F, sep="\t")
SNP_file_name = paste(base.dir, "/DosageFile_cdhr3_15_MatrixEQTL_sliced.txt", sep="");
#IL-6
write.table(SNPs_cytokines, "DosageFile_il6_15_MatrixEQTL_sliced.txt", row.names = F, sep="\t")
SNP_file_name = paste(base.dir, "/DosageFile_il6_15_MatrixEQTL_sliced.txt", sep="");
#All SNPs
write.table(SNPs_cytokines, "DosageFile_AllSNPs_15_MatrixEQTL_sliced.txt", row.names = F, sep="\t")
SNP_file_name = paste(base.dir, "/DosageFile_AllSNPs_15_MatrixEQTL_sliced.txt", sep="");

#CDHR3
snps_location_file_name = paste(base.dir, "/SNP_location_cdhr3_MatrixEQTL.txt", sep="");
snploc<-read_tsv(snps_location_file_name)
snploc<-snploc[snploc$snp %in% SNPs_cytokines$ID,]
write.table(snploc, "SNP_location_cdhr3_15_MatrixEQTL.txt", row.names = F, sep="\t")
snps_location_file_name = paste(base.dir, "/SNP_location_cdhr3_15_MatrixEQTL.txt", sep="");
#IL-6
snps_location_file_name = paste(base.dir, "/SNP_location_il6_MatrixEQTL.txt", sep="");
snploc<-read_tsv(snps_location_file_name)
snploc<-snploc[snploc$snp %in% SNPs_cytokines$ID,]
write.table(snploc, "SNP_location_il6_15_MatrixEQTL.txt", row.names = F, sep="\t")
snps_location_file_name = paste(base.dir, "/SNP_location_il6_15_MatrixEQTL.txt", sep="");
#All SNPs
snps_location_file_name = paste(base.dir, "/SNP_location_AllSNPs_MatrixEQTL.txt", sep="");
snploc<-read_tsv(snps_location_file_name)
snploc<-snploc[snploc$snp %in% rownames(snps),]
write.table(snploc, "SNP_location_AllSNPs_15_MatrixEQTL.txt", row.names = F, sep="\t")
snps_location_file_name = paste(base.dir, "/SNP_location_AllSNPs_15_MatrixEQTL.txt", sep="");

##############################################################################
library(reshape2)

df_stim <- df
df_stim<- df_stim %>% select(chr7_106018005 ,IL_2_LPS)

df_med<-df
df_med<- df_med %>% select(chr7_106018005 ,IL_2)

df_merged<- cbind(df_med, df_stim)
df_merged<-df_merged[,-c(3)]

test<-melt(df_merged, id=c("chr7_106018005"))
rnames<-c(rownames(df_merged), rownames(df_merged))
test$IDs=rnames

#Plot the media vs stimulated values for IL-2
subset(test,!is.na(test$chr7_106018005)) %>% 
  ggplot(aes(x=as.factor(chr7_106018005), y=as.numeric(value), fill=variable)) + 
  geom_boxplot() +
 # geom_point(aes(group=variable), position=position_dodge(width = 0.75)) +
  geom_point(aes(group=variable), position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75)) +
  labs(title="cQTL: Baseline and LPS-stimulated IL-2 levels vs RS6967330",
       x="RS6967330 genotype",
       y="IL-2 (Log2)") +
  theme(plot.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) +
  scale_fill_discrete(labels = c("IL-2_Media", "IL-2_LPS")) 

#As previous but connect the dots
subset(test,!is.na(test$chr7_106018005) & !is.na(test$value)) %>% 
  ggplot(aes(x=variable, y=value)) + 
  facet_wrap(~as.factor(chr7_106018005)) +
  geom_boxplot(aes(fill=variable)) +
  geom_point(aes(group=variable), position=position_dodge(width = 0.75)) +
  geom_line(aes(group=IDs, colour=IDs)) +
  labs(title="cQTL: Baseline and LPS-stimulated IL-2 levels vs RS6967330",
       x="RS6967330 genotype",
       y="IL-2 (Log2)")


#Only Plot Hom. Alt.and connect the dots
subset(test,!is.na(test$chr7_106018005) & test$chr7_106018005 == 2) %>% 
  ggplot(aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  geom_point(aes(group=variable)) +
  geom_line(aes(group=IDs, colour=IDs)) +
  guides(colour=FALSE) +
  labs(title="Hom. Alt.for RS6967330: Baseline and LPS-stimulated IL-2 levels",
       x="RS6967330 genotype",
       y="IL-2 (Log2)") +
  theme(plot.title = element_text(size = 19, face = "bold"),
        axis.text = element_text(size=15, face = "bold"),
        axis.title = element_text(size = 25, face="bold"),
        legend.text = element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold")) #+
  scale_fill_discrete(labels = c("IL-2_Media", "IL-2_LPS"))



