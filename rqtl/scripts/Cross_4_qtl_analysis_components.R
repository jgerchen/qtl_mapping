setwd("/media/fredo/11a4cb61-5b90-4af8-b4fe-0296472c7fb1/home/fredo/Desktop_overflow/qtl_paper/upload_suppl_information/rqtl/")

library(qtl)
cross_obj<-read.cross(format="csvr", file ="data/5A.csvr", alleles=c("A","B"))
phenotypes_5A<-pull.pheno(cross_obj)
cross_obj_sub<-subset(cross_obj, ind=phenotypes_5A$FRE<0.1)

cross_obj_sub<-jittermap(cross_obj_sub)
cross_obj_sub<-calc.genoprob(cross_obj_sub, step=1)

qtl_table=data.frame(trait=c(),chr=c(), marker=c(), begin=c(), center=c(),end=c())



cross_additonal_pheno<-read.table("data/cross_5A_components.tsv", header=T)

cross_pheno_sub<-data.frame(id=cross_additonal_pheno$id, female_flowers_count=cross_additonal_pheno$female_flowers_count, male_flowers_weight=cross_additonal_pheno$male_flowers_weight)
cross_obj_sub$pheno<-left_join(cross_obj_sub$pheno, cross_pheno_sub, by="id")


out_5A.hk_female_flowers <- scanone(cross_obj_sub, pheno.col = 5, method="hk")
operm_5A.hk_female_flowers <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 5)
summary(operm_5A.hk_female_flowers)
summary(out_5A.hk_female_flowers)

out_5A.hk_male_flowers <- scanone(cross_obj_sub, pheno.col = 6, method="hk")
operm_5A.hk_male_flowers <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 6)
summary(operm_5A.hk_male_flowers)
summary(out_5A.hk_male_flowers)



col_female_flowers="#E69F00"  
col_male_flowers="#56B4E9"


max(pull.map(cross_obj_sub, chr="Lg1", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg2", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg3", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg4", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg5", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg6", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg7", as.table = TRUE)$pos)+
  max(pull.map(cross_obj_sub, chr="Lg8", as.table = TRUE)$pos)


pdf("plots/5A_single_qtl_components.pdf", width=13, height=6)
plot(out_5A.hk_female_flowers, out_5A.hk_male_flowers, lodcolumn = 1, col=c(col_female_flowers, col_male_flowers), ylab="LOD", xlab="")

abline(h=summary(operm_5A.hk_female_flowers)[1], col=col_female_flowers, lty=2, lwd=2)
abline(h=summary(operm_5A.hk_male_flowers)[1], col=col_male_flowers, lty=2, lwd=2)


abline(v=xaxisloc.scanone(out_5A.hk_male_flowers, thechr = summary(out_5A.hk_male_flowers)[1,1], summary(out_5A.hk_male_flowers)[1,2]), col=col_male_flowers, lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_male_flowers, thechr = summary(out_5A.hk_male_flowers)[4,1], summary(out_5A.hk_male_flowers)[4,2]), col=col_male_flowers, lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_male_flowers, thechr = summary(out_5A.hk_male_flowers)[5,1], summary(out_5A.hk_male_flowers)[5,2]), col=col_male_flowers, lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_male_flowers, thechr = summary(out_5A.hk_male_flowers)[7,1], summary(out_5A.hk_male_flowers)[7,2]), col=col_male_flowers, lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_male_flowers, thechr = summary(out_5A.hk_male_flowers)[8,1], summary(out_5A.hk_male_flowers)[8,2]), col=col_male_flowers, lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_5A.hk_female_flowers, thechr = summary(out_5A.hk_female_flowers)[4,1], summary(out_5A.hk_female_flowers)[4,2]), col=col_female_flowers, lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_female_flowers, thechr = summary(out_5A.hk_female_flowers)[8,1], summary(out_5A.hk_female_flowers)[8,2]), col=col_female_flowers, lty=2, lwd=2)


dev.off()


