
setwd("/home/fredo/Desktop/qtl_paper/upload_suppl_information/rqtl/")
library(qtl)
library(knitr)
cross_obj<-read.cross(format="csvr", file ="data/1A.csvr", alleles=c("A","B"))
cross_obj<-jittermap(cross_obj)
cross_obj<-calc.genoprob(cross_obj, step=1)
qtl_table=data.frame(trait=c(),chr=c(), marker=c(), begin=c(), center=c(),end=c())

#FRE
out_1A.hk_FRE <- scanone(cross_obj, pheno.col = 2, method="hk")
operm_1A.hk_FRE <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 2)
#summary(operm_1A.hk_FRE)
#summary(out_1A.hk_FRE)

#Here I am importing the results of the permutations for the two-QTL analysis, which I ran in chunks on the cluster using the following command: 
#scantwo_perm <- scantwopermhk(cross_obj, pheno.col = 2, n.perm=100, verbose=TRUE, batchsize = 20)
#read FRE
x<-readRDS("data/cross_1A_out_FRE/1.rds")
for(i in 2:100){
  st<-paste("data/cross_1A_out_FRE/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
#plot histograms of permutations
#plot(x)
#this is very slow, I just run it once, export the result as rds and reimport it for rerunning the whole script
#out_1A.hk_FRE_two <- scantwo(cross_obj, pheno.col = 2, method="hk")
#saveRDS(out_1A.hk_FRE_two, file = "data/out_1A_hk_FRE_two.rds")
out_1A.hk_FRE_two <-readRDS("data/out_1A_hk_FRE_two.rds")
fre_table_1A<-summary(out_1A.hk_FRE_two, perms=x, pvalues = TRUE)
#write.table(fre_table_1A,"data/1A_fre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


kable(fre_table_1A, caption="Output of the two-QTL analysis", digits = 3)
qtl_1a_FRE_additive<-makeqtl(cross_obj, chr=c("Lg3", "Lg5"), pos=c(9,38),  what="prob")
qtl_1a_FRE_additive_fit<-fitqtl(cross_obj, qtl=qtl_1a_FRE_additive, method="hk", pheno.col = 2, formula=y~Q1+Q2)
qtl_1a_FRE_additive_refine<-refineqtl(cross_obj, qtl=qtl_1a_FRE_additive, method="hk", pheno.col = 2, formula=y~Q1+Q2, keeplodprofile = TRUE)
qtl_1a_FRE_additive_fit_refine<-fitqtl(cross_obj, qtl=qtl_1a_FRE_additive_refine, method="hk", pheno.col = 2, formula=y~Q1+Q2)
qtl_1a_FRE_additive_fit_refine_ests<-fitqtl(cross_obj, qtl=qtl_1a_FRE_additive_refine, method="hk", pheno.col = 2, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)

summary(qtl_1a_FRE_additive_fit_refine)
summary(qtl_1a_FRE_additive_fit_refine_ests)
effectplot(cross = cross_obj, mname1 = "Lg3@9.6", pheno.col = 2, main="Effect plot for FRE on Lg3")
effectplot(cross = cross_obj, mname1 = "Lg5@38", pheno.col = 2, main="Effect plot for FRE on Lg5")

bi_qtl_1a_FRE_additive_fit_refine_lg3<-bayesint(results=qtl_1a_FRE_additive_refine, chr="Lg3", qtl.index = 1)
bi_qtl_1a_FRE_additive_fit_refine_lg3_marker<-bayesint(results=qtl_1a_FRE_additive_refine, chr="Lg3", qtl.index = 1, expandtomarkers = TRUE)

kable(bi_qtl_1a_FRE_additive_fit_refine_lg3, caption="Baysian interval FRE QTL on Lg3", digits = 3)
kable(bi_qtl_1a_FRE_additive_fit_refine_lg3_marker, caption="Baysian interval FRE QTL on Lg3 (physical loci)", digits = 3)

qtl_table<-rbind(qtl_table,
      data.frame(trait="FRE",chr=bi_qtl_1a_FRE_additive_fit_refine_lg3_marker$chr[1], begin=row.names(bi_qtl_1a_FRE_additive_fit_refine_lg3_marker)[1], 
                 center=row.names(bi_qtl_1a_FRE_additive_fit_refine_lg3_marker)[2], 
                 end=row.names(bi_qtl_1a_FRE_additive_fit_refine_lg3_marker)[3])) 


bi_qtl_1a_FRE_additive_fit_refine_lg5<-bayesint(results=qtl_1a_FRE_additive_refine, chr="Lg5", qtl.index = 2)
bi_qtl_1a_FRE_additive_fit_refine_lg5_marker<-bayesint(results=qtl_1a_FRE_additive_refine, chr="Lg5", qtl.index = 2, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_1a_FRE_additive_fit_refine_lg5_marker$chr[1], 
                            begin=row.names(bi_qtl_1a_FRE_additive_fit_refine_lg5_marker)[1], center=row.names(bi_qtl_1a_FRE_additive_fit_refine_lg5_marker)[2], 
                            end=row.names(bi_qtl_1a_FRE_additive_fit_refine_lg5_marker)[3])) 


kable(bi_qtl_1a_FRE_additive_fit_refine_lg5, caption="Baysian interval FRE QTL on Lg5", digits = 3)
kable(bi_qtl_1a_FRE_additive_fit_refine_lg5_marker, caption="Baysian interval FRE QTL on Lg5 (physical loci)", digits = 3)



#MRE
#run and plot FRE
out_1A.hk_MRE <- scanone(cross_obj, pheno.col = 3, method="hk")
operm_1A.hk_MRE <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 3)
summary(operm_1A.hk_MRE)
summary(out_1A.hk_MRE)
x<-readRDS("data/cross_1A_out_MRE/1.rds")
for(i in 2:100){
  st<-paste("data/cross_1A_out_MRE/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
plot(x)
out_1A.hk_MRE_two <-readRDS("data/out_1A_hk_MRE_two.rds")

mre_table_1A<-summary(out_1A.hk_MRE_two, perms=x, pvalues = TRUE)
#write.table(mre_table_1A,"data/1A_mre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

kable(mre_table_1A, caption="Output of the two-QTL analysis", digits = 3)



#QTL analysis MRE all loci
qtl_1a_MRE_additive<-makeqtl(cross_obj, chr=c("Lg4", "Lg5", "Lg7"), pos=c(9,38,64),  what="prob")
qtl_1a_MRE_additive_fit<-fitqtl(cross_obj, qtl=qtl_1a_MRE_additive, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3)
qtl_1a_MRE_additive_refine<-refineqtl(cross_obj, qtl=qtl_1a_MRE_additive, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3, keeplodprofile = TRUE)
qtl_1a_MRE_additive_fit_refine<-fitqtl(cross_obj, qtl=qtl_1a_MRE_additive_refine, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3)
qtl_1a_MRE_additive_fit_refine_ests<-fitqtl(cross_obj, qtl=qtl_1a_MRE_additive_refine, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3, dropone=FALSE, get.ests=TRUE)


summary(qtl_1a_MRE_additive_fit_refine)
summary(qtl_1a_MRE_additive_fit_refine_ests)
effectplot(cross = cross_obj, mname1 = "Lg4@49", pheno.col = 3, main="Effect plot for MRE on Lg4")
effectplot(cross = cross_obj, mname1 = "Lg5@40", pheno.col = 3, main="Effect plot for MRE on Lg5")
effectplot(cross = cross_obj, mname1 = "Lg7@59.1", pheno.col = 3, main="Effect plot for MRE on Lg7")
bi_qtl_1a_MRE_additive_fit_refine_lg4<-bayesint(results=qtl_1a_MRE_additive_refine, chr="Lg4", qtl.index = 1)
bi_qtl_1a_MRE_additive_fit_refine_lg4_marker<-bayesint(results=qtl_1a_MRE_additive_refine, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)


#MRE QTL analysis lg7 removed
qtl_1a_MRE_additive_nosd<-makeqtl(cross_obj, chr=c("Lg4", "Lg5"), pos=c(9,38),  what="prob")
qtl_1a_MRE_additive_fit_nosd<-fitqtl(cross_obj, qtl=qtl_1a_MRE_additive_nosd, method="hk", pheno.col = 3, formula=y~Q1+Q2)
qtl_1a_MRE_additive_refine_nosd<-refineqtl(cross_obj, qtl=qtl_1a_MRE_additive_nosd, method="hk", pheno.col = 3, formula=y~Q1+Q2, keeplodprofile = TRUE)
qtl_1a_MRE_additive_fit_refine_nosd<-fitqtl(cross_obj, qtl=qtl_1a_MRE_additive_refine_nosd, method="hk", pheno.col = 3, formula=y~Q1+Q2)
qtl_1a_MRE_additive_fit_refine_ests_nosd<-fitqtl(cross_obj, qtl=qtl_1a_MRE_additive_refine_nosd, method="hk", pheno.col = 3, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)

summary(qtl_1a_MRE_additive_fit_refine_nosd)
summary(qtl_1a_MRE_additive_fit_refine_ests_nosd)
bi_qtl_1a_MRE_additive_fit_refine_lg4_nosd<-bayesint(results=qtl_1a_MRE_additive_refine_nosd, chr="Lg4", qtl.index = 1)
bi_qtl_1a_MRE_additive_fit_refine_lg4_marker_nosd<-bayesint(results=qtl_1a_MRE_additive_refine_nosd, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)
bi_qtl_1a_MRE_additive_fit_refine_lg5_nosd<-bayesint(results=qtl_1a_MRE_additive_refine_nosd, chr="Lg5", qtl.index = 2)
bi_qtl_1a_MRE_additive_fit_refine_lg5_marker_nosd<-bayesint(results=qtl_1a_MRE_additive_refine_nosd, chr="Lg5", qtl.index = 2, expandtomarkers = TRUE)


#Find the closest physical marker to center: -> OW569317.1_36741192
MRE_lg4_map<-pull.map(cross_obj, as.table = TRUE, chr = "Lg4")

qtl_table<-rbind(qtl_table,
                      data.frame(trait="MRE",chr=bi_qtl_1a_MRE_additive_fit_refine_lg4_marker_nosd$chr[1], 
                            begin=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg4_marker_nosd)[1], 
                            center="OW569317.1_36583907", 
                            end=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg4_marker_nosd)[3])) 

kable(bi_qtl_1a_MRE_additive_fit_refine_lg4, caption="Baysian interval MRE QTL on Lg4", digits = 3)
kable(bi_qtl_1a_MRE_additive_fit_refine_lg4_marker, caption="Baysian interval MRE QTL on Lg4 (physical loci)", digits = 3)
bi_qtl_1a_MRE_additive_fit_refine_lg5<-bayesint(results=qtl_1a_MRE_additive_refine, chr="Lg5", qtl.index = 2)
bi_qtl_1a_MRE_additive_fit_refine_lg5_marker<-bayesint(results=qtl_1a_MRE_additive_refine, chr="Lg5", qtl.index = 2, expandtomarkers = TRUE)

#Find the closest physical marker to center: -> OW569314.1_15605703
MRE_lg5_map<-pull.map(cross_obj, as.table = TRUE, chr = "Lg5")

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_1a_MRE_additive_fit_refine_lg5_marker_nosd$chr[1], 
                            begin=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg5_marker_nosd)[1], 
                            center=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg5_marker_nosd)[2], 
                            end=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg5_marker_nosd)[3])) 



kable(bi_qtl_1a_MRE_additive_fit_refine_lg5, caption="Baysian interval MRE QTL on Lg5", digits = 3)
kable(bi_qtl_1a_MRE_additive_fit_refine_lg5_marker, caption="Baysian interval MRE QTL on Lg4 (physical loci)", digits = 3)
bi_qtl_1a_MRE_additive_fit_refine_lg7<-bayesint(results=qtl_1a_MRE_additive_refine, chr="Lg7", qtl.index = 3)
bi_qtl_1a_MRE_additive_fit_refine_lg7_marker<-bayesint(results=qtl_1a_MRE_additive_refine, chr="Lg7", qtl.index = 3, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_1a_MRE_additive_fit_refine_lg7_marker$chr[1], 
                            begin=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg7_marker)[1], 
                            center=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg7_marker)[2], 
                            end=row.names(bi_qtl_1a_MRE_additive_fit_refine_lg7_marker)[3])) 

kable(bi_qtl_1a_MRE_additive_fit_refine_lg7, caption="Baysian interval MRE QTL on Lg7", digits = 3)
kable(bi_qtl_1a_MRE_additive_fit_refine_lg7_marker, caption="Baysian interval MRE QTL on Lg4 (physical loci)", digits = 3)

#rMA
out_1A.hk_rMA<- scanone(cross_obj, pheno.col = 4, method="hk")
operm_1A.hk_rMA <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 4)
summary(operm_1A.hk_rMA)
summary(out_1A.hk_rMA)
x<-readRDS("data/cross_1A_out_rMA/1.rds")
for(i in 2:100){
  st<-paste("data/cross_1A_out_rMA/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
plot(x)

out_1A.hk_rMA_two<-readRDS("data/out_1A_hk_rMA_two.rds")
rMA_table_1A<-summary(out_1A.hk_rMA_two, perms=x, pvalues = TRUE)
#write.table(rMA_table_1A,"data/1A_rma.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


kable(rMA_table_1A, caption="Output of the two-QTL analysis", digits = 3)

#multi-QTL analysis Lg7 removed

qtl_1a_rMA_additive_nosd<-makeqtl(cross_obj, chr=c("Lg5", "Lg6"), pos=c(39.7,51.5),  what="prob")
qtl_1a_rMA_additive_fit_nosd<-fitqtl(cross_obj, qtl=qtl_1a_rMA_additive_nosd, method="hk", pheno.col = 4, formula=y~Q1+Q2)
qtl_1a_rMA_additive_refine_nosd<-refineqtl(cross_obj, qtl=qtl_1a_rMA_additive_nosd, method="hk", pheno.col = 4, formula=y~Q1+Q2, keeplodprofile = TRUE)
qtl_1a_rMA_additive_fit_refine_nosd<-fitqtl(cross_obj, qtl=qtl_1a_rMA_additive_refine_nosd, method="hk", pheno.col = 4, formula=y~Q1+Q2)
qtl_1a_rMA_additive_fit_refine_ests_nosd<-fitqtl(cross_obj, qtl=qtl_1a_rMA_additive_refine_nosd, method="hk", pheno.col = 4, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)

bi_qtl_1a_rMA_additive_fit_refine_lg5_marker_nosd<-bayesint(results=qtl_1a_rMA_additive_refine_nosd, chr="Lg5", qtl.index = 1, expandtomarkers = TRUE)
bi_qtl_1a_rMA_additive_fit_refine_lg6_marker_nosd<-bayesint(results=qtl_1a_rMA_additive_refine_nosd, chr="Lg6", qtl.index = 2, expandtomarkers = TRUE)



qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_1a_rMA_additive_fit_refine_lg5_marker_nosd$chr[1], 
                            begin=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg5_marker_nosd)[1], 
                            center=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg5_marker_nosd)[2], 
                            end=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg5_marker_nosd)[3])) 


#get locus for center

rMA_lg6_map<-pull.map(cross_obj, as.table = TRUE, chr = "Lg6")
qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_1a_rMA_additive_fit_refine_lg6_marker_nosd$chr[1], 
                            begin=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg6_marker_nosd)[1], 
                            center="OW569316.1_39321199", 
                            end=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg6_marker_nosd)[3])) 


#multi-QTL analysis all loci

qtl_1a_rMA_additive<-makeqtl(cross_obj, chr=c("Lg5", "Lg6", "Lg7"), pos=c(39.7,51.5,59.1),  what="prob")
qtl_1a_rMA_additive_fit<-fitqtl(cross_obj, qtl=qtl_1a_rMA_additive, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3)
qtl_1a_rMA_additive_refine<-refineqtl(cross_obj, qtl=qtl_1a_rMA_additive, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3, keeplodprofile = TRUE)
qtl_1a_rMA_additive_fit_refine<-fitqtl(cross_obj, qtl=qtl_1a_rMA_additive_refine, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3)
qtl_1a_rMA_additive_fit_refine_ests<-fitqtl(cross_obj, qtl=qtl_1a_rMA_additive_refine, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3, dropone=FALSE, get.ests=TRUE)

summary(qtl_1a_rMA_additive_fit_refine)
summary(qtl_1a_rMA_additive_fit_refine_ests)
bi_qtl_1a_rMA_additive_fit_refine_lg5<-bayesint(results=qtl_1a_rMA_additive_refine, chr="Lg5", qtl.index = 1)
bi_qtl_1a_rMA_additive_fit_refine_lg5_marker<-bayesint(results=qtl_1a_rMA_additive_refine, chr="Lg5", qtl.index = 1, expandtomarkers = TRUE)

#put QTLs

kable(bi_qtl_1a_rMA_additive_fit_refine_lg5, caption="Baysian interval rMA QTL on Lg5", digits = 3)
kable(bi_qtl_1a_rMA_additive_fit_refine_lg5_marker, caption="Baysian interval rMA QTL on Lg5 (physical loci)", digits = 3)
bi_qtl_1a_rMA_additive_fit_refine_lg6<-bayesint(results=qtl_1a_rMA_additive_refine, chr="Lg6", qtl.index = 2)
bi_qtl_1a_rMA_additive_fit_refine_lg6_marker<-bayesint(results=qtl_1a_rMA_additive_refine, chr="Lg6", qtl.index = 2, expandtomarkers = TRUE)




kable(bi_qtl_1a_rMA_additive_fit_refine_lg6, caption="Baysian interval rMA QTL on Lg6", digits = 3)
kable(bi_qtl_1a_rMA_additive_fit_refine_lg6_marker, caption="Baysian interval rMA QTL on Lg6 (physical loci)", digits = 3)
bi_qtl_1a_rMA_additive_fit_refine_lg7<-bayesint(results=qtl_1a_rMA_additive_refine, chr="Lg7", qtl.index = 3)
bi_qtl_1a_rMA_additive_fit_refine_lg7_marker<-bayesint(results=qtl_1a_rMA_additive_refine, chr="Lg7", qtl.index = 3, expandtomarkers = TRUE)


qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_1a_rMA_additive_fit_refine_lg7_marker$chr[1], 
                            begin=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg7_marker)[1], 
                            center=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg7_marker)[2], 
                            end=row.names(bi_qtl_1a_rMA_additive_fit_refine_lg7_marker)[3])) 




kable(bi_qtl_1a_rMA_additive_fit_refine_lg7, caption="Baysian interval rMA QTL on Lg7", digits = 3)
kable(bi_qtl_1a_rMA_additive_fit_refine_lg7_marker, caption="Baysian interval rMA QTL on Lg6 (physical loci)", digits = 3)
#export qtl table
#saveRDS(qtl_table, "cross_1A_qtls.Rdata")

#############FINAL PLOTS!!!!!!!!!!!!!!!!!!!!!!!!####################################################


max(pull.map(cross_obj, chr="Lg1", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg2", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg3", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg4", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg5", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg6", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg7", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg8", as.table = TRUE)$pos)


pdf("plots/1A_single_qtl_final.pdf", width=13, height=6)
plot(out_1A.hk_FRE, out_1A.hk_MRE, out_1A.hk_rMA  ,lodcolumn = 1, col=c("blue", "red", "darkgreen"), ylab="LOD", xlab="")
abline(h=summary(operm_1A.hk_FRE)[1], col="blue", lty=2, lwd=2)
#xaxisloc.scanone(out_1A.hk_FRE, thechr = summary(out_1A.hk_FRE)[3,1], summary(out_1A.hk_FRE)[3,2])
abline(v=xaxisloc.scanone(out_1A.hk_FRE, thechr = summary(out_1A.hk_FRE)[3,1], summary(out_1A.hk_FRE)[3,2]), col="blue", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_FRE, thechr = summary(out_1A.hk_FRE)[5,1], summary(out_1A.hk_FRE)[5,2]), col="blue", lty=2, lwd=2)

abline(h=summary(operm_1A.hk_MRE)[1], col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_MRE, thechr = summary(out_1A.hk_MRE)[4,1], summary(out_1A.hk_MRE)[4,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_MRE, thechr = summary(out_1A.hk_MRE)[5,1], summary(out_1A.hk_MRE)[5,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_MRE, thechr = summary(out_1A.hk_MRE)[7,1], summary(out_1A.hk_MRE)[7,2]), col="red", lty=2, lwd=2)

abline(h=summary(operm_1A.hk_rMA)[1], col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_MRE, thechr = summary(out_1A.hk_rMA)[5,1], summary(out_1A.hk_rMA)[5,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_1A.hk_MRE, thechr = summary(out_1A.hk_rMA)[6,1], summary(out_1A.hk_rMA)[6,2]), col="darkgreen", lty=2, lwd=2)

dev.off()