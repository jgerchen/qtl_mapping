library(qtl)
library(knitr)
setwd("~/Desktop/qtl_paper/upload_suppl_information/rqtl")
cross_obj<-read.cross(format="csvr", file ="data/3A.csvr", alleles=c("A","B"))
phenotypes_3A<-pull.pheno(cross_obj)
#remove MRE outliers for QTL analysis!
cross_obj_sub<-subset(cross_obj, ind=phenotypes_3A$MRE<0.01)

cross_obj_sub<-jittermap(cross_obj_sub)
cross_obj_sub<-calc.genoprob(cross_obj_sub, step=1)
#save.image("data/cross_3A.Rdata")
qtl_table=data.frame(trait=c(),chr=c(), marker=c(), begin=c(), center=c(),end=c())



#run and plot FRE
out_3A.hk_FRE <- scanone(cross_obj_sub, pheno.col = 2, method="hk")
operm_3A.hk_FRE <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 2)
summary(operm_3A.hk_FRE)
summary(out_3A.hk_FRE)
#Here I am importing the results of the permutations for the two-QTL analysis, which I ran in chunks on the cluster using the following command: 
#scantwo_perm <- scantwopermhk(cross_obj_sub, pheno.col = 2, n.perm=100, verbose=TRUE, batchsize = 20)
#read FRE
x<-readRDS("data/cross_3A_out_FRE/1.rds")
for(i in 2:100){
  st<-paste("data/cross_3A_out_FRE/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}

#plot histograms of permutations
plot(x)


#out_3A.hk_FRE_two <- scantwo(cross_obj_sub, pheno.col = 2, method="hk")
#saveRDS(out_3A.hk_FRE_two, file = "out_3A_hk_FRE_two.rds")
out_3A.hk_FRE_two <-readRDS("data/out_3A_hk_FRE_two.rds")

fre_table_3A<-summary(out_3A.hk_FRE_two, perms=x, pvalues = TRUE)
#write.table(fre_table_3A,"data/3A_fre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

# model with Lg8 removed

qtl_3A_FRE_additive_nosd<-makeqtl(cross_obj_sub, chr=c("Lg4"), pos=c(40.5),  what="prob")
qtl_3A_FRE_additive_fit_nosd<-fitqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive_nosd, method="hk", pheno.col = 2, formula=y~Q1)
qtl_3A_FRE_additive_refine_nosd<-refineqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive_nosd, method="hk", pheno.col = 2, formula=y~Q1, keeplodprofile = TRUE)
qtl_3A_FRE_additive_fit_refine_nosd<-fitqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive_refine_nosd, method="hk", pheno.col = 2, formula=y~Q1)
qtl_3A_FRE_additive_fit_refine_ests_nosd<-fitqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive_refine_nosd, method="hk", pheno.col = 2, formula=y~Q1, dropone=FALSE, get.ests=TRUE)
bi_qtl_3A_FRE_additive_fit_refine_lg4_marker_nosd<-bayesint(results=qtl_3A_FRE_additive_refine_nosd, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)


# full model
qtl_3A_FRE_additive<-makeqtl(cross_obj_sub, chr=c("Lg4", "Lg8"), pos=c(40.5, 42.5),  what="prob")
qtl_3A_FRE_additive_fit<-fitqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive, method="hk", pheno.col = 2, formula=y~Q1+Q2)
qtl_3A_FRE_additive_refine<-refineqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive, method="hk", pheno.col = 2, formula=y~Q1+Q2, keeplodprofile = TRUE)
qtl_3A_FRE_additive_fit_refine<-fitqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive_refine, method="hk", pheno.col = 2, formula=y~Q1+Q2)
qtl_3A_FRE_additive_fit_refine_ests<-fitqtl(cross_obj_sub, qtl=qtl_3A_FRE_additive_refine, method="hk", pheno.col = 2, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)


summary(qtl_3A_FRE_additive_fit_refine)
summary(qtl_3A_FRE_additive_fit_refine_ests)
effectplot(cross = cross_obj_sub, mname1 = "Lg4@42.3", pheno.col = 2, main="Effect plot for FRE on Lg4")

effectplot(cross = cross_obj_sub, mname1 = "Lg8@6.2", pheno.col = 2, main="Effect plot for FRE on Lg8")

bi_qtl_3A_FRE_additive_fit_refine_lg4<-bayesint(results=qtl_3A_FRE_additive_refine, chr="Lg4", qtl.index = 1)
bi_qtl_3A_FRE_additive_fit_refine_lg4_marker<-bayesint(results=qtl_3A_FRE_additive_refine, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)
bi_qtl_3A_FRE_additive_fit_refine_lg8<-bayesint(results=qtl_3A_FRE_additive_refine, chr="Lg8", qtl.index = 2)
bi_qtl_3A_FRE_additive_fit_refine_lg8_marker<-bayesint(results=qtl_3A_FRE_additive_refine, chr="Lg8", qtl.index = 2, expandtomarkers = TRUE)

kable(bi_qtl_3A_FRE_additive_fit_refine_lg4, caption="Baysian interval FRE QTL on Lg4", digits = 3)
kable(bi_qtl_3A_FRE_additive_fit_refine_lg4_marker, caption="Baysian interval FRE QTL on Lg4 (physical loci)", digits = 3)


qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_3A_FRE_additive_fit_refine_lg4_marker$chr[1], begin=row.names(bi_qtl_3A_FRE_additive_fit_refine_lg4_marker)[1], 
                            center=row.names(bi_qtl_3A_FRE_additive_fit_refine_lg4_marker)[2], 
                            end=row.names(bi_qtl_3A_FRE_additive_fit_refine_lg4_marker)[3]))

qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_3A_FRE_additive_fit_refine_lg8_marker$chr[1], begin=row.names(bi_qtl_3A_FRE_additive_fit_refine_lg8_marker)[1], 
                            center=row.names(bi_qtl_3A_FRE_additive_fit_refine_lg8_marker)[2], 
                            end=row.names(bi_qtl_3A_FRE_additive_fit_refine_lg8_marker)[3]))


#MRE

#run and plot FRE
out_3A.hk_MRE <- scanone(cross_obj_sub, pheno.col = 3, method="hk")
operm_3A.hk_MRE <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 3)
summary(operm_3A.hk_MRE)

#Here I am importing the results of the permutations for the two-QTL analysis, which I ran in chunks on the cluster using the following command: 
#scantwo_perm <- scantwopermhk(cross_obj_sub, pheno.col = 2, n.perm=100, verbose=TRUE, batchsize = 20)
#read FRE
x<-readRDS("data/cross_3A_out_MRE/1.rds")
for(i in 2:100){
  st<-paste("data/cross_3A_out_MRE/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
#summary(x, alpha=0.05)

#plot histograms of permutations
plot(x)


#out_3A.hk_MRE_two <- scantwo(cross_obj_sub, pheno.col = 3, method="hk")
#saveRDS(out_3A.hk_FRE_two, file = "out_3A_hk_FRE_two.rds")
out_3A.hk_MRE_two <-readRDS("data/out_3A_hk_FRE_two.rds")
mre_table_3A<-summary(out_3A.hk_MRE_two, perms=x, pvalues = TRUE)
#write.table(mre_table_3A,"data/3A_mre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

###No signficance!

#rMA

out_3A.hk_rMA<- scanone(cross_obj_sub, pheno.col = 4, method="hk")
operm_3A.hk_rMA <- scanone(cross_obj_sub, method="hk", n.perm=10000, pheno.col = 4)
summary(operm_3A.hk_rMA)
x<-readRDS("data/cross_3A_out_rMA/1.rds")
for(i in 2:100){
  st<-paste("data/cross_3A_out_rMA/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
plot(x)

#out_3A.hk_rMA_two <- scantwo(cross_obj_sub, pheno.col = 4, method="hk")
#saveRDS(out_3A.hk_rMA_two, file = "out_3A_hk_rMA_two.rds")
out_3A.hk_rMA_two <-readRDS("data/out_3A_hk_rMA_two.rds")
rMA_table_3A<-summary(out_3A.hk_rMA_two, perms=x, pvalues = TRUE)
#write.table(rMA_table_3A,"data/3A_rma.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


#saveRDS(qtl_table, "cross_3A_qtls.Rdata")



pdf("plots/3A_single_qtl_final.pdf", width=13, height=6)
plot(out_3A.hk_FRE, out_3A.hk_MRE, out_3A.hk_rMA  ,lodcolumn = 1, col=c("blue", "red", "darkgreen"), ylab="LOD", xlab="", ylim=c(0,6))
abline(h=summary(operm_3A.hk_FRE)[1], col="blue", lty=2, lwd=2)
abline(h=summary(operm_3A.hk_MRE)[1], col="red", lty=2, lwd=2)
abline(h=summary(operm_3A.hk_rMA)[1], col="darkgreen", lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_3A.hk_FRE, thechr = summary(out_3A.hk_FRE)[4,1], summary(out_3A.hk_FRE)[4,2]), col="blue", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_3A.hk_FRE, thechr = summary(out_3A.hk_FRE)[8,1], summary(out_3A.hk_FRE)[8,2]), col="blue", lty=2, lwd=2)

dev.off()