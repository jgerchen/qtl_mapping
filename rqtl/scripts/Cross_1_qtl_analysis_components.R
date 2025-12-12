
setwd("/media/fredo/11a4cb61-5b90-4af8-b4fe-0296472c7fb1/home/fredo/Desktop_overflow/qtl_paper/upload_suppl_information/rqtl/")
library(qtl)
library(knitr)
library(dplyr)
cross_obj<-read.cross(format="csvr", file ="data/1A.csvr", alleles=c("A","B"))
cross_additonal_pheno<-read.table("data/cross_1A_components.tsv", header=T)
#multiply biomass of male flowers and female flowers for samples where only half the plant was phenotyped
cross_additonal_pheno$female_flowers_weight_adj<-cross_additonal_pheno$female_flowers_weight
cross_additonal_pheno$male_flowers_weight_adj<-cross_additonal_pheno$male_flowers_weight

for(add_pheno_i in 1:nrow(cross_additonal_pheno)){
  if(cross_additonal_pheno$phenotyped_at[add_pheno_i]=="half_plant"){
    cross_additonal_pheno$female_flowers_weight_adj[add_pheno_i]<-2*cross_additonal_pheno$female_flowers_weight[add_pheno_i]
    cross_additonal_pheno$male_flowers_weight_adj[add_pheno_i]<-2*cross_additonal_pheno$male_flowers_weight[add_pheno_i]
    
  }
  
}
cross_pheno_sub<-data.frame(id=cross_additonal_pheno$id, female_flowers_weight=cross_additonal_pheno$female_flowers_weight_adj, male_flowers_weight=cross_additonal_pheno$male_flowers_weight_adj)
cross_obj$pheno<-left_join(cross_obj$pheno, cross_pheno_sub, by="id")


cross_obj<-jittermap(cross_obj)
cross_obj<-calc.genoprob(cross_obj, step=1)

out_1A.hk_female_flowers <- scanone(cross_obj, pheno.col = 5, method="hk")
operm_1A.hk_female_flowers <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 5)

out_1A.hk_male_flowers <- scanone(cross_obj, pheno.col = 6, method="hk")
operm_1A.hk_male_flowers <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 6)


#############FINAL PLOTS!!!!!!!!!!!!!!!!!!!!!!!!####################################################

col_female_flowers="#E69F00"  
col_male_flowers="#56B4E9"


max(pull.map(cross_obj, chr="Lg1", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg2", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg3", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg4", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg5", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg6", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg7", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg8", as.table = TRUE)$pos)


pdf("plots/1A_single_qtl_components.pdf", width=13, height=6)
plot(out_1A.hk_female_flowers, out_1A.hk_male_flowers, lodcolumn = 1, col=c(col_female_flowers, col_male_flowers), ylab="LOD", xlab="")

abline(h=summary(operm_1A.hk_female_flowers)[1], col=col_female_flowers, lty=2, lwd=2)
abline(h=summary(operm_1A.hk_male_flowers)[1], col=col_male_flowers, lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_1A.hk_female_flowers, thechr = summary(out_1A.hk_female_flowers)[5,1], summary(out_1A.hk_female_flowers)[5,2]), col=col_female_flowers, lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_1A.hk_male_flowers, thechr = summary(out_1A.hk_male_flowers)[5,1], summary(out_1A.hk_male_flowers)[5,2]), col=col_male_flowers, lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_male_flowers, thechr = summary(out_1A.hk_male_flowers)[7,1], summary(out_1A.hk_male_flowers)[7,2]), col=col_male_flowers, lty=2, lwd=2)

dev.off()