library(tidyverse)

#Set working directory
setwd("B:/MauroTutino/LD")

#Load the LD matrix computed using plink
LD_matrix<- read.table("plink_chr7.ld", header=T)

#Extract only the region and the columns of interest
LD_matrix_cdhr3 <- LD_matrix[as.numeric(as.character(LD_matrix$BP_A)) >= 105958239 & as.numeric(as.character(LD_matrix$BP_A)) <= 106038772 ,]

##########################################################
#Load the RS lookup
rs_lookup<-read.table("Chr_Pos_Rs_lookup.txt")
#Rename the columns to be meaninful
colnames(rs_lookup)<-c("CHR","POS","RS","REF","ALT")
#Create a column to match the one in he LD matrix
rs_lookup <- rs_lookup %>% mutate(CHR_POS=paste(CHR,POS, sep="_"))
#Subset to only the SNPs of interest
rs_lookup<-rs_lookup[rs_lookup$CHR_POS %in% LD_matrix_cdhr3$SNP_A | rs_lookup$CHR_POS %in% LD_matrix_cdhr3$SNP_B,]
#Rename the chr-pos to rs whenever there is a match in the lookup table
LD_matrix_cdhr3$SNP_A<-rs_lookup$RS[match(LD_matrix_cdhr3$SNP_A, rs_lookup$CHR_POS)]
#Change the factor levels to be in the same order of the elements in the column
#Before doing it, we need to reset the factor levels
#Do it for SNP_A ...
LD_matrix_cdhr3$SNP_A <- as.character(LD_matrix_cdhr3$SNP_A)
LD_matrix_cdhr3$SNP_A <- factor(LD_matrix_cdhr3$SNP_A, levels=unique(LD_matrix_cdhr3$SNP_A))
#And SNP_B
LD_matrix_cdhr3$SNP_B<-rs_lookup$RS[match(LD_matrix_cdhr3$SNP_B, rs_lookup$CHR_POS)]
LD_matrix_cdhr3$SNP_B <- as.character(LD_matrix_cdhr3$SNP_B)
LD_matrix_cdhr3$SNP_B <- factor(LD_matrix_cdhr3$SNP_B, levels=unique(LD_matrix_cdhr3$SNP_B))


#######################################################################

#Make the plot
ggheatmap <- ggplot(LD_matrix_cdhr3, aes(x=SNP_A, y=SNP_B, fill=R2)) + 
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
#Print PDF
pdf("LD_plot.pdf", width = 10, height = 10)
print(p)
dev.off()

#Print PNG
jpeg("LD_plot.jpeg", width = 1200, height = 1200, units = "px")
print(p)
dev.off()
