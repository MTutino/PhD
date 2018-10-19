install.packages("ggrepel")
intall.packages("tidyverse")

library(tidyverse)
library(ggrepel)


setwd("//nask.man.ac.uk/home$/Desktop")

readcount<-read_tsv("readcount.txt",col_names = FALSE)
Cov_lt_10<-read_tsv("Cov_lt_10.txt",col_names = FALSE)
tot_mapped_reads_no_chr6<-read_tsv("Total_mapped_reads_no_chr6.txt",col_names = FALSE)

colnames(readcount)<-c("SampleID","readcount")
colnames(Cov_lt_10)<-c("SampleID","Coverage")
colnames(tot_mapped_reads_no_chr6)<-c("SampleID","mapped_reads_no_chr6")


readCount_cov <- inner_join(readcount, Cov_lt_10)
readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_no_chr6)



############################################################################################################

min_read_count<- min(readCount_cov$readcount, na.rm = T)
max_read_count<- max(readCount_cov$readcount, na.rm = T)

readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_ten=1 - Coverage)

readCount_cov <- readCount_cov %>% 
  mutate(low_cov_name = ifelse(perc_bases_cov_higher_ten <
                                0.8, readCount_cov$SampleID, NA))
num_samples_no_cov <- with(readCount_cov, table(perc_bases_cov_higher_ten < 0.8)["TRUE"])

pdf("number_of_reads_VS_coverage.pdf",width=20,height=14)

  ggplot(readCount_cov,aes(x=readcount, y=perc_bases_cov_higher_ten)) + 
    geom_point(size=2, colour="blue") +
 
     xlim(min_read_count,max_read_count) +
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 35, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold")) +
   
    # geom_label_repel(aes(label=low_cov_name), type = "closed", force = 5, colour = "red", segment.alpha = 0.5) +
   geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +
     
    annotate("text", x = 1.5e+7, y = 0.4, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in < 80% of\n targeted bases: ",num_samples_no_cov ), 
             colour="red",
             size=15) +
  labs(
    x="# Pair-end reads",
    y="% of bases with depth of coverage > 10X",
    title="Capture Targeted Exon Sequencing - no HLA region" )
     
   
  dev.off()
  
####################################################################################################################  
 
  
   min_read_count<- min(readCount_cov$readcount, na.rm = T)
  max_read_count<- max(readCount_cov$readcount, na.rm = T)
  
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_ten=1 - Coverage)
  
  readCount_cov <- readCount_cov %>% 
    mutate(low_cov_name = ifelse(perc_bases_cov_higher_ten <
                                   0.8, readCount_cov$SampleID, NA))
  num_samples_no_cov <- with(readCount_cov, table(perc_bases_cov_higher_ten < 0.8)["TRUE"])
  
  readCount_cov <- readCount_cov %>% mutate(perc_mapped_reads=mapped_reads_no_chr6/readcount*100)
  
  ggplot(readCount_cov,aes(x=readcount, y=perc_bases_cov_higher_ten)) + 
    geom_point(size=2, colour="blue") +
    
    xlim(min_read_count,max_read_count) +
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold")) +
    
    #geom_label_repel(aes(label=low_cov_name), type = "closed", force = 5, colour = "red", segment.alpha = 0.5) +
    
   geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +
   # geom_hline(aes(yintercept = mean(readCount_cov$perc_bases_cov_higher_ten), color="mean"),linetype="dashed", size=2 ) +
    
    annotate("text", x = 1.5e+7, y = 0.4, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in < 80% of\n targeted bases: ",num_samples_no_cov ), 
             colour="red",
             size=15) +
    labs(
      x="# Pair-end reads",
      y="% of bases with depth of coverage > 10X",
      title="Capture Targeted Exon Sequencing - no HLA region" )
  
  
  
  #############################################################################################################
  min_read_count<- min(readCount_cov$mapped_reads_no_chr6, na.rm = T)
  max_read_count<- max(readCount_cov$mapped_reads_no_chr6, na.rm = T)
  
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_ten=1 - Coverage)
  
  readCount_cov <- readCount_cov %>% 
    mutate(low_cov_name = ifelse(perc_bases_cov_higher_ten <
                                   0.8, readCount_cov$SampleID, NA))
  num_samples_no_cov <- with(readCount_cov, table(perc_bases_cov_higher_ten < 0.8)["TRUE"])
  
  readCount_cov <- readCount_cov %>% mutate(perc_mapped_reads=mapped_reads_no_chr6/(readcount*2)*100)
  
  ggplot(readCount_cov,aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_ten)) + 
    geom_point(size=2, colour="blue") +
    
  
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold")) +
    
    geom_label_repel(aes(label=low_cov_name), type = "closed", force = 5, colour = "red", segment.alpha = 0.5) +
    
    geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +
    # geom_hline(aes(yintercept = mean(readCount_cov$perc_bases_cov_higher_ten), color="mean"),linetype="dashed", size=2 ) +
    
    annotate("text", x = 2e+6, y = 0.5, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in < 80% of\n targeted bases: ",num_samples_no_cov ), 
             colour="red",
             size=6) +
    annotate("text", x = 2e+6, y = 0.1,
             label=paste("Mean % of reads mapped to\n chr other than chr6: ",mean(readCount_cov$perc_mapped_reads) ), 
             colour="red",
             size=6) +
    labs(
      x="# Mapped Pair-end reads",
      y="% of bases with depth of coverage > 10X",
      title="Capture Targeted Exon Sequencing - no HLA region" )

#########################################################################################
 #################################  FROM HERE ###########################################
  
  library(tidyverse)
  library(ggrepel)
  library(ggthemes)
  library(RColorBrewer)
  
  setwd("//nask.man.ac.uk/home$/Desktop")
  
  readcount<-read_tsv("readcount.txt",col_names = FALSE)
  #Change this to make plots at other coverage depths
  Cov_lt_10_no_chr6<-read_tsv("Cov_lt_10_no_chr6.txt",col_names = FALSE)
  Cov_lt_10_chr6<-read_tsv("Cov_lt_10_chr6.txt",col_names = FALSE)
  Cov_lt_15_no_chr6<-read_tsv("Cov_lt_15_no_chr6.txt",col_names = FALSE)
  Cov_lt_15_chr6<-read_tsv("Cov_lt_15_chr6.txt",col_names = FALSE)
  Cov_lt_20_no_chr6<-read_tsv("Cov_lt_20_no_chr6.txt",col_names = FALSE)
  Cov_lt_20_chr6<-read_tsv("Cov_lt_20_chr6.txt",col_names = FALSE)
  tot_mapped_reads_no_chr6<-read_tsv("Total_mapped_reads_no_chr6.txt",col_names = FALSE)
  tot_mapped_reads_chr6<-read_tsv("Total_mapped_reads_chr6.txt",col_names = FALSE)
  header_list <- read_tsv("Headers_list.txt", col_name=FALSE)
  
  
  
  colnames(readcount)<-c("SampleID","readcount")
  colnames(Cov_lt_10_no_chr6)<-c("SampleID","Coverage_no_chr6")
  colnames(Cov_lt_10_chr6)<-c("SampleID","Coverage_chr6")
  colnames(Cov_lt_15_no_chr6)<-c("SampleID","Coverage_15_no_chr6")
  colnames(Cov_lt_15_chr6)<-c("SampleID","Coverage_15_chr6")
  colnames(Cov_lt_20_no_chr6)<-c("SampleID","Coverage_20_no_chr6")
  colnames(Cov_lt_20_chr6)<-c("SampleID","Coverage_20_chr6")
  colnames(tot_mapped_reads_no_chr6)<-c("SampleID","mapped_reads_no_chr6")
  colnames(tot_mapped_reads_chr6)<-c("SampleID","mapped_reads_chr6")
  colnames(header_list)<-c("SampleID","header")
  
  header_list<- header_list %>%
    separate(SampleID, c("bin", "SampleID"), "_") %>%
    separate(header, c("machine", "run","flowcell","lane","bin2"), ":", extra="merge") %>%
    mutate(runlane=paste(run,flowcell,lane, sep=":")) %>%
    select(SampleID,runlane)
  
  readCount_cov <- inner_join(readcount, Cov_lt_10_no_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_10_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_15_no_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_15_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_20_no_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_20_chr6)
  readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_no_chr6)
  readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_chr6)
  readCount_cov <- inner_join(readCount_cov, header_list)
  
  ######################## Plot no-chr6 ###################
  
  min_read_count<- min(readCount_cov$mapped_reads_no_chr6, na.rm = T)
  max_read_count<- max(readCount_cov$mapped_reads_no_chr6, na.rm = T)
  
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_ten_no_chr6=1 - Coverage_no_chr6)
  
  readCount_cov <- readCount_cov %>% 
    mutate(low_cov_name_no_chr6 = ifelse(perc_bases_cov_higher_ten_no_chr6 <
                                   0.8, readCount_cov$SampleID, NA))
  num_samples_no_cov <- with(readCount_cov, table(perc_bases_cov_higher_ten_no_chr6 < 0.8)["FALSE"])
  
  readCount_cov <- readCount_cov %>% mutate(perc_mapped_reads_no_chr6=mapped_reads_no_chr6/(readcount*2)*100)
  
  colourCount <- length(unique(readCount_cov$runlane)) # number of levels
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  pdf("number_of_mapped_reads_VS_coverage_no_chr6.pdf",width=20,height=14)
  
  ggplot(readCount_cov,aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_ten_no_chr6)) + 
    geom_point(size=2, aes(colour=runlane)) +
    #scale_fill_manual(values = getPalette(colourCount)) +
   # scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(colourCount)) +
   # scale_color_manual(values = colorRampPalette(col_vector[15:])(colourCount)) +
    
  #xlim(min_read_count,max_read_count) +
    
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=15, face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
  
    
    #geom_label_repel(aes(label=low_cov_name), type = "closed", force = 5, colour = "red", segment.alpha = 0.5) +
    
    geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +
    # geom_hline(aes(yintercept = mean(readCount_cov$perc_bases_cov_higher_ten), color="mean"),linetype="dashed", size=2 ) +
    
    annotate("text", x = 2e+6, y = 0.5, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in > 80% of\n targeted bases: ",num_samples_no_cov, "(",round(num_samples_no_cov/length(readCount_cov$SampleID)*100), "% of samples )" ), 
             colour="red",
             size=8) +
    annotate("text", x = 2e+6, y = 0.1,
             label=paste("Mean % of reads mapped to\n chrs other than chr6: ",round(mean(readCount_cov$perc_mapped_reads_no_chr6), 2 ),"%"), 
             colour="red",
             size=8) +
    labs(
      x="# Mapped Pair-end reads",
      y="% of bases with depth of coverage > 10X",
      title="Capture Targeted Exon Sequencing - no chr6 region (no HLA)",
      colour="Run:Flowcell:Lane")
  
  dev.off()
  
  
  ################# plot chr6 ##################################
  
  min_read_count<- min(readCount_cov$mapped_reads_chr6, na.rm = T)
  max_read_count<- max(readCount_cov$mapped_reads_chr6, na.rm = T)
  
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_ten_chr6=1 - Coverage_chr6)
  
  readCount_cov <- readCount_cov %>% 
    mutate(low_cov_name = ifelse(perc_bases_cov_higher_ten_chr6 <
                                   0.8, readCount_cov$SampleID, NA))
  num_samples_no_cov <- with(readCount_cov, table(perc_bases_cov_higher_ten_chr6 < 0.8)["FALSE"])
  
  readCount_cov <- readCount_cov %>% mutate(perc_mapped_reads_chr6=mapped_reads_chr6/(readcount*2)*100)
  
  colourCount <- length(unique(readCount_cov$runlane)) # number of levels
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  
  pdf("number_of_mapped_reads_VS_coverage_chr6.pdf",width=20,height=14)
  
  ggplot(readCount_cov,aes(x=mapped_reads_chr6, y=perc_bases_cov_higher_ten_chr6)) + 
    geom_point(size=2, aes(colour=runlane)) +
    #scale_fill_manual(values = getPalette(colourCount)) +
    # scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(colourCount)) +
    # scale_color_manual(values = colorRampPalette(col_vector[15:])(colourCount)) +
    
    #xlim(min_read_count,max_read_count) +
    ylim(0.00,1) +
    
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=15, face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    
    
    #geom_label_repel(aes(label=low_cov_name), type = "closed", force = 5, colour = "red", segment.alpha = 0.5) +
    
    geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +
    # geom_hline(aes(yintercept = mean(readCount_cov$perc_bases_cov_higher_ten), color="mean"),linetype="dashed", size=2 ) +
    
    annotate("text", x = 2.0e+7, y = 0.5, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in > 80% of\n targeted bases: ",num_samples_no_cov, "(",round(num_samples_no_cov/length(readCount_cov$SampleID)*100), "% of samples )" ), 
             colour="red",
             size=8) +
    annotate("text", x = 2.0e+7, y = 0.1,
             label=paste("Mean % of reads mapped to\n chr6: ",round(mean(readCount_cov$perc_mapped_reads_chr6), 2 ),"%"), 
             colour="red",
             size=8) +
    labs(
      x="# Mapped Pair-end reads",
      y="% of bases with depth of coverage > 10X",
      title="Capture Targeted Exon Sequencing chr6 (HLA region)",
      colour="Run:Flowcell:Lane")
  
  dev.off()
  
  ##########################  Overlay different coverage depths #######################################
  
  min_read_count<- min(readCount_cov$mapped_reads_no_chr6, na.rm = T)
  max_read_count<- max(readCount_cov$mapped_reads_no_chr6, na.rm = T)
  
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_ten_no_chr6=1 - Coverage_no_chr6)
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_15_no_chr6=1 - Coverage_15_no_chr6)
  readCount_cov <- readCount_cov %>% mutate(perc_bases_cov_higher_20_no_chr6=1 - Coverage_20_no_chr6)
  
  pdf("Coverage_at_multiple_depths_no_chr6.pdf",width=20,height=14)
  

  num_samples_no_cov_10 <- with(readCount_cov, table(perc_bases_cov_higher_ten_no_chr6 < 0.8)["FALSE"])
  num_samples_no_cov_15 <- with(readCount_cov, table(perc_bases_cov_higher_15_no_chr6 < 0.8)["FALSE"])
  num_samples_no_cov_20 <- with(readCount_cov, table(perc_bases_cov_higher_20_no_chr6 < 0.8)["FALSE"])
  
  
  ggplot(readCount_cov,aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_ten_no_chr6)) + 
    geom_smooth(aes(colour="10X"), size=2) +
    geom_point(data=readCount_cov, aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_ten_no_chr6,colour="10X"),size=2, alpha=0.1) +
    geom_smooth(data=readCount_cov, aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_15_no_chr6, colour="15X"), size=2) +
    geom_point(data=readCount_cov, aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_15_no_chr6,colour="15X"),size=2, alpha=0.1) +
    geom_smooth(data=readCount_cov, aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_20_no_chr6, colour="20X"), size=2) +
    geom_point(data=readCount_cov, aes(x=mapped_reads_no_chr6, y=perc_bases_cov_higher_20_no_chr6,colour="20X"),alpha=0.1, size=2) +
   
    ylim(0.00,1) +
    
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold"),
          legend.text = element_text(size=18, face="bold"),
          legend.title = element_text(size=20, face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    
    geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +

    annotate("text", x = 2e+6, y = 0.6, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in > 80% of\n targeted bases: ",num_samples_no_cov_10, "(",round(num_samples_no_cov_10/length(readCount_cov$SampleID)*100), "% of samples )" ), 
             colour="coral2",
             size=8) +
    annotate("text", x = 2e+6, y = 0.4,
             label=paste("Number of samples with\ndepth of coverage > 15X\n in > 80% of\n targeted bases: ",num_samples_no_cov_15, "(",round(num_samples_no_cov_15/length(readCount_cov$SampleID)*100), "% of samples )" ), 
             colour="green3",
             size=8) +
    annotate("text", x = 2e+6, y = 0.2,
             label=paste("Number of samples with\ndepth of coverage > 20X\n in > 80% of\n targeted bases: ",num_samples_no_cov_20, "(",round(num_samples_no_cov_20/length(readCount_cov$SampleID)*100), "% of samples )" ), 
             colour="blue",
             size=8) +
    labs(
      x="# Mapped Pair-end reads",
      y="% of bases with depth of coverage > 10X/15X/20X",
      title="Capture Targeted Exon Sequencing - no chr6 region (no HLA)",
      colour="Coverage")
  
  dev.off()
  
  
  
  ########################### MAPPED READS ON CHR6 VS MAPPED READS NO-CHR6  ############################
  
  library(tidyverse)
  library(ggrepel)
  library(ggthemes)
  library(RColorBrewer)
  
  setwd("//nask.man.ac.uk/home$/Desktop")
  
  readcount<-read_tsv("readcount.txt",col_names = FALSE)
  Cov_lt_10_no_chr6<-read_tsv("Cov_lt_10_no_chr6.txt",col_names = FALSE)
  Cov_lt_10_chr6<-read_tsv("Cov_lt_10_chr6.txt",col_names = FALSE)
  tot_mapped_reads_no_chr6<-read_tsv("Total_mapped_reads_no_chr6.txt",col_names = FALSE)
  tot_mapped_reads_chr6<-read_tsv("Total_mapped_reads_chr6.txt",col_names = FALSE)
  header_list <- read_tsv("Headers_list.txt", col_name=FALSE)
  
  
  
  colnames(readcount)<-c("SampleID","readcount")
  colnames(Cov_lt_10_no_chr6)<-c("SampleID","Coverage_no_chr6")
  colnames(Cov_lt_10_chr6)<-c("SampleID","Coverage_chr6")
  colnames(tot_mapped_reads_no_chr6)<-c("SampleID","mapped_reads_no_chr6")
  colnames(tot_mapped_reads_chr6)<-c("SampleID","mapped_reads_chr6")
  colnames(header_list)<-c("SampleID","header")
  
  header_list<- header_list %>%
    separate(SampleID, c("bin", "SampleID"), "_") %>%
    separate(header, c("machine", "run","flowcell","lane","bin2"), ":", extra="merge") %>%
    mutate(runlane=paste(run,flowcell,lane, sep=":")) %>%
    select(SampleID,runlane)
  
  readCount_cov <- inner_join(readcount, Cov_lt_10_no_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_10_chr6)
  readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_no_chr6)
  readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_chr6)
  readCount_cov <- inner_join(readCount_cov, header_list)
  
  
  
  readCount_cov <- readCount_cov %>% mutate(norm_mapped_reads_no_chr6=mapped_reads_no_chr6/1.9)
  readCount_cov <- readCount_cov %>% mutate(norm_mapped_reads_chr6=mapped_reads_chr6/4.8)
  
  colourCount <- length(unique(readCount_cov$runlane)) # number of levels
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  
  
  ggplot(readCount_cov,aes(x=norm_mapped_reads_no_chr6, y=norm_mapped_reads_chr6)) + 
    geom_point(size=2, aes(colour=runlane)) +
    #scale_x_log10() +
    #scale_y_log10() +
    
    theme(panel.background = element_rect(fill="white", colour="white"),
          panel.grid.minor= element_line(colour = "grey", size = 0.1),
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size=15, face = "bold"),
          axis.title = element_text(size = 25, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=15, face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    
    
    geom_hline(yintercept=0.8,linetype="dashed", colour="red", size=2) +
    
    annotate("text", x = 1.5e+8, y = 0.5, 
             label=paste("Number of samples with\ndepth of coverage > 10X\n in < 80% of\n targeted bases: ",num_samples_no_cov, "(",round(num_samples_no_cov/length(readCount_cov$SampleID)*100), "% of samples )" ), 
             colour="red",
             size=8) +
    annotate("text", x = 1.5e+8, y = 0.1,
             label=paste("Mean % of reads mapped to\n chr6: ",round(mean(readCount_cov$perc_mapped_reads_chr6), 2 ),"%"), 
             colour="red",
             size=8) +
    labs(
      x="# mapped reads on chrms other than chr6/Mb",
      y="# mapped reads on chr6/Mb",
     # title="Capture Targeted Exon Sequencing chr6 (HLA region)",
      colour="Run:Flowcell:Lane")
  
  ########################## MAKE BOX PLOT INSTEAD  #########################################
  library(tidyverse)
  library(ggrepel)
  library(ggthemes)
  library(RColorBrewer)
  
  setwd("//nask.man.ac.uk/home$/Desktop")
  
  readcount<-read_tsv("readcount.txt",col_names = FALSE)
  Cov_lt_10_no_chr6<-read_tsv("Cov_lt_10_no_chr6.txt",col_names = FALSE)
  Cov_lt_10_chr6<-read_tsv("Cov_lt_10_chr6.txt",col_names = FALSE)
  tot_mapped_reads_no_chr6<-read_tsv("Total_mapped_reads_no_chr6.txt",col_names = FALSE)
  tot_mapped_reads_chr6<-read_tsv("Total_mapped_reads_chr6.txt",col_names = FALSE)
  header_list <- read_tsv("Headers_list.txt", col_name=FALSE)
  
  
  
  colnames(readcount)<-c("SampleID","readcount")
  colnames(Cov_lt_10_no_chr6)<-c("SampleID","Coverage_no_chr6")
  colnames(Cov_lt_10_chr6)<-c("SampleID","Coverage_chr6")
  colnames(tot_mapped_reads_no_chr6)<-c("SampleID","mapped_reads_no_chr6")
  colnames(tot_mapped_reads_chr6)<-c("SampleID","mapped_reads_chr6")
  colnames(header_list)<-c("SampleID","header")
  
  header_list<- header_list %>%
    separate(SampleID, c("bin", "SampleID"), "_") %>%
    separate(header, c("machine", "run","flowcell","lane","bin2"), ":", extra="merge") %>%
    mutate(runlane=paste(run,flowcell,lane, sep=":")) %>%
    select(SampleID,runlane)
  
  readCount_cov <- inner_join(readcount, Cov_lt_10_no_chr6)
  readCount_cov <- inner_join(readCount_cov, Cov_lt_10_chr6)
  readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_no_chr6)
  readCount_cov <- inner_join(readCount_cov, tot_mapped_reads_chr6)
  readCount_cov <- inner_join(readCount_cov, header_list)
  
  
  
  readCount_cov <- readCount_cov %>% mutate(norm_mapped_reads_no_chr6=mapped_reads_no_chr6/1.9)
  readCount_cov <- readCount_cov %>% mutate(norm_mapped_reads_chr6=mapped_reads_chr6/4.8)
  
  a = data.frame(group = "no_chr6", value=readCount_cov$norm_mapped_reads_no_chr6)
  b = data.frame(group = "chr6", value=readCount_cov$norm_mapped_reads_chr6)
    plot.data = rbind(a, b)
    
 pdf("number_of_mapped_reads_per_megabase_boxplot.pdf",width=20,height=14)
    
    ggplot(plot.data, aes(x=group, y=value, fill=group)) +  # This is the plot function
      geom_boxplot() +
      labs(
        x="Regions",
        y="Number of mapped reads/Mb",
        title="Number of reads mapped to targeted regions",
        fill="Regions") +
      
      theme(#panel.background = element_rect(fill="white", colour="white"),
          #  panel.grid.minor= element_line(colour = "grey", size = 0.1),
          #  panel.grid.major = element_line(colour = "grey", size = 0.1),
            plot.title = element_text(size = 30, face = "bold"),
            axis.text = element_text(size=20, face = "bold"),
            axis.title = element_text(size = 25, face="bold"),
            legend.text = element_text(size=18, face="bold"),
            legend.title = element_text(size=20, face="bold"))
  dev.off()
      
      
  
  