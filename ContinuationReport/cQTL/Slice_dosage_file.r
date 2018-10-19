library(tidyverse)

setwd("B:/MauroTutino/VariantCalling")


# Genotype file name
#No HWE
SNP_file_name = "First_analysis_primary_assembly_only_only_SNP_HardFiltered_AC2_biallelic_db151_ann_noHLA_PASS_MissInd_CallR.recoded.renamed.dosage"
#AFTER HW Equilibrium filtering
SNP_file_name = "First_analysis_primary_assembly_only_only_SNP_HardFiltered_AC2_biallelic_db151_ann_noHLA_PASS_MissInd_CallR_HWE.recoded.dosage"

#Load the file in memory
dosage<-read_tsv(SNP_file_name)

#rename the #CHROM to CHROM to avoid problems with the "#" symbol
colnames(dosage)[1] <- "CHROM"

#Create column with SNP ID = CHROM_POSITION
dosage <- dosage %>% mutate(id=paste(CHROM,POS, sep="_"))
#Move the new column to the first position
dosage<- dosage[,c(ncol(dosage),1:ncol(dosage)-1)]



#Extract only the variants of interest
#Example: extarct only the variants on chromosome 4 and 5
#dosage_chr4_5dosage$id[grep("chr4|chr5",dosage$id)]

#Extract variants on chr4
#dosage_chr4 <- dosage[grep("chr4",dosage$id),]

#Then select only those variants between position 38772484 and 38773484
#dosage_sel <- dosage[dosage$POS >= 38772484 & dosage$POS < 38773484,]

#CDHR3
#chr7:105,963,239-106,033,772

#Extract all variants in CDHR3 +/- 5K bp
dosage_sel_cdhr3 <- dosage[dosage$CHROM == "chr7" & dosage$POS >= 105958239 & dosage$POS <= 106038772,]

#IL6
#chr7:22,725,884-22,732,002 +/- 5 kb
dosage_sel_il6 <- dosage[dosage$CHROM == "chr7" & dosage$POS >= 22720884 & dosage$POS <= 22738002,]


#Write SNP location for MatrixeQTL
setwd("C:/Users/mdxasmt5/Dropbox/Systems-Immunology/MauroPhD/Lijing_cytokines")
snplocation_cdhr3<-dosage_sel_cdhr3[,c(1,2,3)]
colnames(snplocation_cdhr3)<-c("snp","chr","pos")
write.table(snplocation_cdhr3, "SNP_location_cdhr3_MatrixEQTL.txt", row.names = F, sep="\t")

snplocation_il6<-dosage_sel_il6[,c(1,2,3)]
colnames(snplocation_il6)<-c("snp","chr","pos")
write.table(snplocation_il6, "SNP_location_il6_MatrixEQTL.txt", row.names = F, sep="\t")

#All SNPs
setwd("C:/Users/mdxasmt5/Dropbox/Systems-Immunology/MauroPhD/Lijing_cytokines")
dosage_sel <- dosage
snplocation<-dosage_sel[,c(1,2,3)]
colnames(snplocation)<-c("snp","chr","pos")
write.table(snplocation, "SNP_location_AllSNPs_MatrixEQTL.txt", row.names = F, sep="\t")



#Get rid of not useful columns for MatrixeQTL
dosage_sel_cdhr3<-dosage_sel_cdhr3[,-c(2,3,4,5)]
dosage_sel_il6<-dosage_sel_il6[,-c(2,3,4,5)]
dosage_sel<-dosage_sel[,-c(2,3,4,5)]

#Write dosage file
write.table(dosage_sel_cdhr3, "DosageFile_cdhr3_MatrixEQTL.txt", row.names = F, sep="\t")
write.table(dosage_sel_il6, "DosageFile_il6_MatrixEQTL.txt", row.names = F, sep="\t")
write.table(dosage_sel, "DosageFile_AllSNPs_MatrixEQTL.txt", row.names = F, sep="\t")


