library(qtl)
library(viridis)
library(colourvalues)
library(fields)
#Plots 1A
setwd("/media/fredo/11a4cb61-5b90-4af8-b4fe-0296472c7fb1/home/fredo/Desktop_overflow/qtl_paper/upload_suppl_information/rqtl/")
phenotypes_1A<-read.table("data/1A_phenotypes.tsv", sep="\t", header=TRUE)
Seq_depth_1A<-read.table("data/1A_coverage.tsv", sep="\t", header=TRUE)
low_depth_1A<-Seq_depth_1A$sample[Seq_depth_1A$mean_cov<10]
#contaminants inferred from PCA
contaminants_1A<-c("1A_158","1A_187", "1A_67")
phenotypes_1A$symbol<-rep(16, length(phenotypes_1A$id))
phenotypes_1A$symbol[as.character(phenotypes_1A$id) %in% contaminants_1A]<-12
phenotypes_1A$symbol[phenotypes_1A$id %in% low_depth_1A]<-8

#Cross S1

phenotypes_S1<-read.table("data/S1_phenotypes.tsv", sep="\t", header=TRUE)

samples_S1_low_cov<-read.table("data/S1_low_cov.lst", header=FALSE)
samples_S1_contam<-read.table("data/S1_contam.lst", header=FALSE)

phenotypes_S1$symbol<-rep(16, length(phenotypes_S1$id))
phenotypes_S1$symbol[as.character(phenotypes_S1$id) %in% samples_S1_contam$V1]<-12
phenotypes_S1$symbol[as.character(phenotypes_S1$id) %in% samples_S1_low_cov$V1]<-8

#Cross 3A
phenotypes_3A<-read.table("data/3A_phenotypes.tsv", sep="\t", header=TRUE)
phenotypes_3A_low_cov<-read.table("data/3A_low_cov.lst", sep="\t", header=FALSE)
phenotypes_3A_contamn<-read.table("data/3A_contam.lst", sep="\t", header=FALSE)

phenotypes_3A$symbol<-rep(16, length(phenotypes_3A$id))
phenotypes_3A$symbol[phenotypes_3A$MRE>0.02]<-10

phenotypes_3A$symbol[as.character(phenotypes_3A$id) %in% phenotypes_3A_contamn$V1]<-12
phenotypes_3A$symbol[as.character(phenotypes_3A$id) %in% phenotypes_3A_low_cov$V1]<-8
#Cross 5A

phenotypes_5A<-read.table("data/5A_phenotypes.tsv", sep="\t", header=TRUE)
phenotypes_5A_low_cov<-read.table("data/5A_low_cov.lst", sep="\t", header=FALSE)
phenotypes_5A_contamn<-read.table("data/5A_contam.lst", sep="\t", header=FALSE)

phenotypes_5A$symbol<-rep(16, length(phenotypes_5A$id))
phenotypes_5A$symbol[as.character(phenotypes_5A$id) %in% phenotypes_5A_contamn$V1]<-12
phenotypes_5A$symbol[as.character(phenotypes_5A$id) %in% phenotypes_5A_low_cov$V1]<-8
phenotypes_5A$symbol[phenotypes_5A$FRE>0.1]<-10

FRE_max<-max(phenotypes_S1$FRE)
MRE_max<-max(phenotypes_1A$MRE)
  
#New colours from Ehouarn
x<-phenotypes_S1$rMA
ehouarn_color<-colorRamp(c("#64A395", "#E1B782","#E67C37","#8F2F0E"))


pdf("plots/phenotypes.pdf", width=8, height=6)
  #par(mfrow=c(1,4))
  layout(matrix(c(1,2,3,4), nrow=1), widths = c(1.3,1,1,1) )
  
  par(mar = c(5.1, 4.1, 4.1, 0))
  
  
  plot(phenotypes_1A$FRE, phenotypes_1A$MRE, pch=phenotypes_1A$symbol,  ylab="", col=rgb(ehouarn_color(phenotypes_1A$rMA), maxColorValue = 255), xlim=c(0, FRE_max), ylim=c(0, MRE_max), xlab="")
  title("Cross 1")
  text(0.07,0.096, paste("n = ",as.character(nrow(phenotypes_1A)),sep=""), cex=1.5)
  par(mar = c(5.1, 0, 4.1, 0))
  plot(phenotypes_S1$FRE, phenotypes_S1$MRE, pch=phenotypes_S1$symbol, col=rgb(ehouarn_color(phenotypes_S1$rMA), maxColorValue = 255), xlim=c(0, FRE_max), ylim=c(0, MRE_max), xlab="", ylab="", yaxt='n')
  text(0.07,0.096, paste("n = ",as.character(nrow(phenotypes_S1)),sep=""), cex=1.5)
  
  
  
  
  title("Cross 2")
  
  plot(phenotypes_3A$FRE, phenotypes_3A$MRE, pch=phenotypes_3A$symbol,  col=rgb(ehouarn_color(phenotypes_3A$rMA), maxColorValue = 255), xlim=c(0, FRE_max), ylim=c(0, MRE_max), xlab="", ylab="", yaxt='n')
  text(0.07,0.096, paste("n = ",as.character(nrow(phenotypes_3A)),sep=""), cex=1.5)
  
  title("Cross 3")
  
  plot(phenotypes_5A$FRE, phenotypes_5A$MRE, pch=phenotypes_5A$symbol,  col=rgb(ehouarn_color(phenotypes_5A$rMA), maxColorValue = 255), xlim=c(0, FRE_max), ylim=c(0, MRE_max), xlab="", ylab="", yaxt='n')
  title("Cross 4")
  legend("topright", legend=c("low coverage","outlier","contamination", "included"), pch=c(8,10,12,16))  
  mtext("Female reproductive effort",side=1,line=-2,outer=TRUE,cex=1)
  mtext("Male reproductive effort",side=2,line=-2,outer=TRUE,cex=1,las=0)
  text(0.12, 0.082, "relative male\nallocation")
  #colorbar.plot(0.15,0.07,(1:100)/100, col=colour_values((1:100)/100), horizontal=FALSE)
  colorbar.plot(0.15,0.07,(1:100)/100, col=rgb(ehouarn_color((1:100)/100), maxColorValue = 255), horizontal=FALSE)
  text(0.135,0.077, "1")
  text(0.135,0.063, "0")
  text(0.05,0.096, paste("n = ",as.character(nrow(phenotypes_5A)),sep=""), cex=1.5)
  
  
dev.off()
