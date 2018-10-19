install.packages("devtools")
library(devtools);
install_github("sterding/powerEQTL")  
library('powerEQTL')


#################################################	CALCULATE POWER FOR DIFFERENT MAFs	############################################################
################################################# THIS ANALYSIS USES AN EFFECT SIZE OF 1	#########################################
# sample size
#N <- c(50,100,150,200,250,300)
N <-c(307)
nn <- length(N)

# MAF
MAF <- seq(0.5,20,0.1)/100
nq <- length(MAF)

# number of SNPs tested (10 SNPs per gene x 24 genes in total)
nSNP=240
#Unfiltered: 30500 SNPs X 6 viral-stimuli X 18 cytokines = 3.3M tests
nSNP=3300000
#AC2: 12953 SNPs X 6 viral-stimuli X 18 cytokines = 1.4 Milion tests
nSNP=1400000
#AC10:6812 SNPs X 6 viral-Stimuli X 18 cytokines = 736000 tests
nSNP=736000
#AF05: 3943 SNPs X 6 viral stimuli X 18 cytokines = 426000 tests
nSNP=425000
#GTEx trans-eqtl number of tests - it reproduces the paper graph
nSNP=100000000000

# significant level (FP)
a=0.05

# obtain power
power_unbalanced <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
  for (j in 1:nq){
    # unbalanced
    result <- powerEQTL.ANOVA(MAF=MAF[j], 
                              typeI=a, 
                              nTests=nSNP, 
                              myntotal=N[i], 
                              mystddev=0.13,
                              deltaVec = c(0.13, 0.13),
                              verbose = F)
    power_unbalanced[i,j] <-result;
  }
}

# set up graph
xrange <- range(MAF*100)
yrange <- c(0:1)
colors <- rainbow(length(N))
plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main=paste("Power Estimation for eQTL Studies\nSig=0.05, nSNP=",nSNP,"(one-way unbalanced ANOVA)",sep=" "))

abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
abline(h=0, v=c(1:10), lty=2,col="grey89")

# add power curves
for (i in 1:nn){
  lines(MAF*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
}

legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')


###################################################### 	Repeat the analysis at different effect sizes
ES<- seq(0.6,3,0.1)
ESl <- length(ES)

#AC10-filtered 6812 X 6 viral-stimuli X 18 cytokines = 735700 tests
nSNP=735700

nSNP=50000

power_unbalanced <- array(numeric(1*ESl), dim=c(1,ESl))
colnames(power_unbalanced) <-paste("EffectSize" ,ES,sep="_" )
rownames(power_unbalanced) <- "MAF"

#loop over the different effect sizes
for (i in 1:ESl){
result <- minMAFeQTL.ANOVA(effsize = ES[i],
                 typeI = 0.05,
                 nTests = nSNP,
                 myntotal = 307,
                 mypower = 0.8,
                 verbose = TRUE)
power_unbalanced[1,i] <-result;
}
#Print array
power_unbalanced

''
    EffectSize_0.6 EffectSize_0.7 EffectSize_0.8 EffectSize_0.9 EffectSize_1
MAF      0.2814119       0.181493      0.1308943     0.09984653   0.07905248
    EffectSize_1.1 EffectSize_1.2 EffectSize_1.3 EffectSize_1.4 EffectSize_1.5
MAF      0.0642792     0.05340753     0.04509721     0.03863782     0.03348637
    EffectSize_1.6 EffectSize_1.7 EffectSize_1.8 EffectSize_1.9 EffectSize_2
MAF     0.02926636     0.02586617     0.02299493     0.02056586   0.01853503
    EffectSize_2.1 EffectSize_2.2 EffectSize_2.3 EffectSize_2.4 EffectSize_2.5
MAF     0.01679551      0.0152795     0.01395514     0.01277915     0.01177674
    EffectSize_2.6 EffectSize_2.7 EffectSize_2.8 EffectSize_2.9 EffectSize_3
MAF     0.01088105     0.01007395    0.009402425    0.008743284  0.008177349
''

library(ggplot2)
df<-as.data.frame.table(power_unbalanced)
colnames(df)<-c("BIN","EffectSize", "MAF")

p<-ggplot(df,aes(x=EffectSize, y=MAF)) +
geom_bar(stat="identity",fill="blue") 
p + geom_text(aes(label=round(MAF,2), vjust=-0.3), col="red") +
  theme(axis.text =  element_text(angle = 90)) +
  labs(title="Minimum Detectable MAF for analysis power of 80% at different effect sizes")
  

xx<-barplot(power_unbalanced) + yrange(0:0.4)
text(x=xx, y=power_unbalanced[1,], label=round(power_unbalanced[1,],2),col="red", pos=3)

################################################	NOW ANALYSE THE EFFECT OF DIFFENT NUMBER OF TESTS	######################################################	
# MAF
MAF <- seq(0.5,20,0.1)/100
nq <- length(MAF)
#Unfiltered: 30414 SNPs X 6 viral-stimuli X 18 cytokines = 3.3M tests
nSNP=3300000
#AC2: 12953 SNPs X 6 viral-stimuli X 18 cytokines = 1.4 Milion tests
nSNP=1400000
#AC10:6812 SNPs X 6 viral-Stimuli X 18 cytokines = 736000 tests
nSNP=736000
#AF05: 3943 SNPs X 6 viral stimuli X 18 cytokines = 426000 tests
nSNP=425000
#Unfiltered: 30414 X 14 stimuli X 28 cytokines: 
11922288

nSNP <- c(425000, 736000, 1400000,3300000,11922288)
nS <- length(nSNP)


nSNP=500

# significant level (FP)
a=0.05

# obtain power
power_unbalanced <- array(numeric(nS*nq), dim=c(nS,nq))
for (i in 1:nS){
  for (j in 1:nq){
    # unbalanced
    result <- powerEQTL.ANOVA(MAF=MAF[j], 
                              typeI=a, 
                              nTests=nSNP[i], 
                              myntotal=307, 
                              mystddev=0.13,
                              deltaVec = c(0.13, 0.13),
                              verbose = F)
    power_unbalanced[i,j] <-result;
  }
}



# set up graph
xrange <- range(MAF*100)
yrange <- c(0:1)
colors <- rainbow(nS)
plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Power Estimation for eQTL Studies\nSig=0.05, SampleNumber=307 (one-way unbalanced ANOVA)")

abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
abline(h=0, v=c(1:10), lty=2,col="grey89")

# add power curves
for (i in 1:nS){
  lines(MAF*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
}

legend("topleft", title="Number Of Tests", as.character(nSNP),fill=colors, bty='n')

############################################################
###########################################################	SAME AS BEFORE BUT USING SIMPLE LINEAR REGRESSION INSTEAD OF ANOVA	##########################################
N <-c(307)
nn <- length(N)

# MAF
MAF <- seq(0.5,20,0.1)/100
nq <- length(MAF)

#Unfiltered: 30500 SNPs X 6 viral-stimuli X 18 cytokines = 3.3M tests
nSNP=3300000
#AC2: 12953 SNPs X 6 viral-stimuli X 18 cytokines = 1.4 Milion tests
nSNP=1400000
#AC10:6812 SNPs X 6 viral-Stimuli X 18 cytokines = 736000 tests
nSNP=736000
#AF05: 3943 SNPs X 6 viral stimuli X 18 cytokines = 426000 tests
nSNP=425000

power_unbalanced <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
  for (j in 1:nq){
     # unbalanced
   result<-  powerEQTL.SLR(
             MAF=MAF[j],
             typeI = 0.05,
             nTests = nSNP,
             slope = 0.13,
             myntotal = N[i],
             mystddev = 0.13,
             verbose = F)
             power_unbalanced[i,j] <-result;
  }
}

# set up graph
xrange <- range(MAF*100)
yrange <- c(0:1)
colors <- rainbow(length(nn))
plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main=paste("Power Estimation for eQTL Studies\nSig=0.05,nSNP=", nSNP," (Simple Linear Regression)"))

abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
abline(h=0, v=c(1:10), lty=2,col="grey89")

# add power curves
for (i in 1:nn){
  lines(MAF*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
}

legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')



