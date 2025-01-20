setwd("/home/fredo/Desktop/qtl_paper/upload_suppl_information/rqtl/")

library(qtl)
cross_obj<-read.cross(format="csvr", file ="data/5A.csvr", alleles=c("A","B"))
phenotypes_5A<-pull.pheno(cross_obj)
cross_obj_sub<-subset(cross_obj, ind=phenotypes_5A$FRE<0.1)

cross_obj_sub<-jittermap(cross_obj_sub)
cross_obj_sub<-calc.genoprob(cross_obj_sub, step=1)

qtl_table=data.frame(trait=c(),chr=c(), marker=c(), begin=c(), center=c(),end=c())



out_5A.hk_FRE_sub <- scanone(cross_obj_sub, pheno.col = 2, method="hk")
operm_5A.hk_FRE_sub <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 2)

out_5A.hk_MRE_sub <- scanone(cross_obj_sub, pheno.col = 3, method="hk")
operm_5A.hk_MRE_sub <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 3)

out_5A.hk_rMA_sub <- scanone(cross_obj_sub, pheno.col = 4, method="hk")
operm_5A.hk_rMA_sub <- scanone(cross_obj_sub, method="hk", n.perm=1000, pheno.col = 4)

# save.image("data/cross_5A_sub.Rdata")


x<-readRDS("data/cross_5A_out_FRE_sub/1.rds")
for(i in 2:100){
  st<-paste("data/cross_5A_out_FRE_sub/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
plot(x)
#runs out of memory...run on cluster!
#out_5A.hk_FRE_two <- scantwo(cross_obj, pheno.col = 2, method="hk")
out_5A.hk_FRE_two_sub<-readRDS("data/out_5A_hk_FRE_sub.rds")
fre_table_5A_sub<-summary(out_5A.hk_FRE_two_sub, perms=x, pvalues = TRUE)

#write.table(fre_table_5A_sub,"data/5A_fre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


x<-readRDS("data/cross_5A_out_MRE_sub/1.rds")
for(i in 2:100){
  st<-paste("data/cross_5A_out_MRE_sub/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
plot(x)
#runs out of memory...run on cluster!
#out_5A.hk_FRE_two <- scantwo(cross_obj, pheno.col = 2, method="hk")
out_5A.hk_MRE_two_sub<-readRDS("data/out_5A_hk_FRE_sub.rds")
mre_table_5A_sub<-summary(out_5A.hk_MRE_two_sub, perms=x, pvalues = TRUE)

#write.table(mre_table_5A_sub,"data/5A_mre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


x<-readRDS("data/cross_5A_out_rMA_sub/1.rds")
for(i in 2:100){
  st<-paste("data/cross_5A_out_rMA_sub/", i, ".rds", sep="")
  sobj<-readRDS(st)
  x<-c(x, sobj)
}
plot(x)
#runs out of memory...run on cluster!
#out_5A.hk_FRE_two <- scantwo(cross_obj, pheno.col = 2, method="hk")
out_5A.hk_rMA_two_sub<-readRDS("data/out_5A_hk_rMA_sub.rds")
rMA_table_5A_sub<-summary(out_5A.hk_rMA_two_sub, perms=x, pvalues = TRUE)

#multi-QTL model
#Lg8 removed
qtl_5a_FRE_additive_sub<-makeqtl(cross_obj_sub, chr=c("Lg4", "Lg7"), pos=c(12.7,72.4),  what="prob")
qtl_5a_FRE_additive_fit_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_sub, method="hk", pheno.col = 2, formula=y~Q1+Q2)
qtl_5a_FRE_additive_refine_sub<-refineqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_sub, method="hk", pheno.col = 2, formula=y~Q1+Q2, keeplodprofile = TRUE)
qtl_5a_FRE_additive_fit_refine_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_refine_sub, method="hk", pheno.col = 2, formula=y~Q1+Q2)
qtl_5a_FRE_additive_fit_refine_ests_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_refine_sub, method="hk", pheno.col = 2, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)

bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub<-bayesint(results=qtl_5a_FRE_additive_refine_sub, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub$chr[1], begin=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub)[3]))



bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub<-bayesint(results=qtl_5a_FRE_additive_refine_sub, chr="Lg7", qtl.index = 2, expandtomarkers = TRUE)


qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub$chr[1], begin=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub)[3]))



#multi-QTL model with two additonal QTLs on LG8
qtl_5a_FRE_additive_sub_sd<-makeqtl(cross_obj_sub, chr=c("Lg4", "Lg7","Lg8","Lg8"), pos=c(12.7,72.4,1,32),  what="prob")
qtl_5a_FRE_additive_fit_sub_sd<-fitqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_sub_sd, method="hk", pheno.col = 2, formula=y~Q1+Q2+Q3+Q4)
qtl_5a_FRE_additive_refine_sub_sd<-refineqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_sub_sd, method="hk", pheno.col = 2, formula=y~Q1+Q2+Q3+Q4, keeplodprofile = TRUE)
qtl_5a_FRE_additive_fit_refine_sub_sd<-fitqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_refine_sub_sd, method="hk", pheno.col = 2, formula=y~Q1+Q2+Q3+Q4)
qtl_5a_FRE_additive_fit_refine_ests_sub_sd<-fitqtl(cross_obj_sub, qtl=qtl_5a_FRE_additive_refine_sub_sd, method="hk", pheno.col = 2, formula=y~Q1+Q2+Q3+Q4, dropone=FALSE, get.ests=TRUE)

bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub_sd<-bayesint(results=qtl_5a_FRE_additive_refine_sub_sd, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)

#qtl_table<-rbind(qtl_table,
#                 data.frame(trait="FREsd",chr=bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub_sd$chr[1], begin=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub_sd)[1], 
#                            center=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub_sd)[2], 
#                            end=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg4_marke_sub_sd)[3]))

bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub_sd<-bayesint(results=qtl_5a_FRE_additive_refine_sub_sd, chr="Lg7", qtl.index = 2, expandtomarkers = TRUE)

#qtl_table<-rbind(qtl_table,
#                 data.frame(trait="FREsd",chr=bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub_sd$chr[1], begin=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub_sd)[1], 
#                            center=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub_sd)[2], 
#                            end=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg7_marke_sub_sd)[3]))

bi_qtl_5a_FRE_additive_fit_refine_lg8_1_marke_sub_sd<-bayesint(results=qtl_5a_FRE_additive_refine_sub_sd, chr="Lg8", qtl.index = 3, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_5a_FRE_additive_fit_refine_lg8_1_marke_sub_sd$chr[1], begin=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg8_1_marke_sub_sd)[1], 
                            center=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg8_1_marke_sub_sd)[2], 
                            end=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg8_1_marke_sub_sd)[3]))


bi_qtl_5a_FRE_additive_fit_refine_lg8_2_marke_sub_sd<-bayesint(results=qtl_5a_FRE_additive_refine_sub_sd, chr="Lg8", qtl.index = 4, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="FRE",chr=bi_qtl_5a_FRE_additive_fit_refine_lg8_2_marke_sub_sd$chr[1], begin=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg8_2_marke_sub_sd)[1], 
                            center="OW569318.1_39358793", 
                            end=row.names(bi_qtl_5a_FRE_additive_fit_refine_lg8_2_marke_sub_sd)[3]))


#MRE multi-QTL model

qtl_5a_MRE_additive_sub<-makeqtl(cross_obj_sub, chr=c("Lg1", "Lg3", "Lg4","Lg5", "Lg7"), pos=c(107.6,69.4,9.7,57.1,65.4),  what="prob")
qtl_5a_MRE_additive_fit_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_MRE_additive_sub, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3+Q4+Q5)
qtl_5a_MRE_additive_refine_sub<-refineqtl(cross_obj_sub, qtl=qtl_5a_MRE_additive_sub, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3+Q4+Q5, keeplodprofile = TRUE)
qtl_5a_MRE_additive_fit_refine_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_MRE_additive_refine_sub, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3+Q4+Q5)
qtl_5a_MRE_additive_fit_refine_ests_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_MRE_additive_refine_sub, method="hk", pheno.col = 3, formula=y~Q1+Q2+Q3+Q4+Q5, dropone=FALSE, get.ests=TRUE)

bi_qtl_5a_MRE_additive_fit_refine_lg1_marke_sub<-bayesint(results=qtl_5a_MRE_additive_refine_sub, chr="Lg1", qtl.index = 1, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_5a_MRE_additive_fit_refine_lg1_marke_sub$chr[1], begin=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg1_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg1_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg1_marke_sub)[3]))

bi_qtl_5a_MRE_additive_fit_refine_lg3_marke_sub<-bayesint(results=qtl_5a_MRE_additive_refine_sub, chr="Lg3", qtl.index = 2, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_5a_MRE_additive_fit_refine_lg3_marke_sub$chr[1], begin=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg3_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg3_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg3_marke_sub)[3]))

bi_qtl_5a_MRE_additive_fit_refine_lg4_marke_sub<-bayesint(results=qtl_5a_MRE_additive_refine_sub, chr="Lg4", qtl.index = 3, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_5a_MRE_additive_fit_refine_lg4_marke_sub$chr[1], begin=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg4_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg4_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg4_marke_sub)[3]))

bi_qtl_5a_MRE_additive_fit_refine_lg5_marke_sub<-bayesint(results=qtl_5a_MRE_additive_refine_sub, chr="Lg5", qtl.index = 4, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_5a_MRE_additive_fit_refine_lg5_marke_sub$chr[1], begin=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg5_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg5_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg5_marke_sub)[3]))


bi_qtl_5a_MRE_additive_fit_refine_lg7_marke_sub<-bayesint(results=qtl_5a_MRE_additive_refine_sub, chr="Lg7", qtl.index = 5, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="MRE",chr=bi_qtl_5a_MRE_additive_fit_refine_lg7_marke_sub$chr[1], begin=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg7_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg7_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_MRE_additive_fit_refine_lg7_marke_sub)[3]))


#rMA multi-QTL model QTL on Lg8 removed


qtl_5a_rMA_additive_sub<-makeqtl(cross_obj_sub, chr=c("Lg1", "Lg3", "Lg4","Lg5", "Lg7"), pos=c(107.41,68.97,11.73,61,62.71),  what="prob")
qtl_5a_rMA_additive_fit_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_sub, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5)
qtl_5a_rMA_additive_refine_sub<-refineqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_sub, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5, keeplodprofile = TRUE)
qtl_5a_rMA_additive_fit_refine_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_refine_sub, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5)
qtl_5a_rMA_additive_fit_refine_ests_sub<-fitqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_refine_sub, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5, dropone=FALSE, get.ests=TRUE)


bi_qtl_5a_rMA_additive_fit_refine_lg1_marke_sub<-bayesint(results=qtl_5a_rMA_additive_refine_sub, chr="Lg1", qtl.index = 1, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_5a_rMA_additive_fit_refine_lg1_marke_sub$chr[1], begin=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg1_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg1_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg1_marke_sub)[3]))

bi_qtl_5a_rMA_additive_fit_refine_lg3_marke_sub<-bayesint(results=qtl_5a_rMA_additive_refine_sub, chr="Lg3", qtl.index = 2, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_5a_rMA_additive_fit_refine_lg3_marke_sub$chr[1], begin=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg3_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg3_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg3_marke_sub)[3]))

bi_qtl_5a_rMA_additive_fit_refine_lg4_marke_sub<-bayesint(results=qtl_5a_rMA_additive_refine_sub, chr="Lg4", qtl.index = 3, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_5a_rMA_additive_fit_refine_lg4_marke_sub$chr[1], begin=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg4_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg4_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg4_marke_sub)[3]))

bi_qtl_5a_rMA_additive_fit_refine_lg5_marke_sub<-bayesint(results=qtl_5a_rMA_additive_refine_sub, chr="Lg5", qtl.index = 4, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_5a_rMA_additive_fit_refine_lg5_marke_sub$chr[1], begin=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg5_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg5_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg5_marke_sub)[3]))

bi_qtl_5a_rMA_additive_fit_refine_lg7_marke_sub<-bayesint(results=qtl_5a_rMA_additive_refine_sub, chr="Lg7", qtl.index = 5, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                 data.frame(trait="rMA",chr=bi_qtl_5a_rMA_additive_fit_refine_lg7_marke_sub$chr[1], begin=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg7_marke_sub)[1], 
                            center=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg7_marke_sub)[2], 
                            end=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg7_marke_sub)[3]))


#rMA with lg8 included



qtl_5a_rMA_additive_sub_sd<-makeqtl(cross_obj_sub, chr=c("Lg1", "Lg3", "Lg4","Lg5", "Lg7", "Lg8"), pos=c(107.41,68.97,11.73,61,62.71, 0),  what="prob")
qtl_5a_rMA_additive_fit_sub_sd<-fitqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_sub_sd, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5+Q6)
qtl_5a_rMA_additive_refine_sub_sd<-refineqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_sub_sd, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5+Q6, keeplodprofile = TRUE)
qtl_5a_rMA_additive_fit_refine_sub_sd<-fitqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_refine_sub_sd, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5+Q6)
qtl_5a_rMA_additive_fit_refine_ests_sub_sd<-fitqtl(cross_obj_sub, qtl=qtl_5a_rMA_additive_refine_sub_sd, method="hk", pheno.col = 4, formula=y~Q1+Q2+Q3+Q4+Q5+Q6, dropone=FALSE, get.ests=TRUE)


bi_qtl_5a_rMA_additive_fit_refine_lg1_marke_sub_sd<-bayesint(results=qtl_5a_rMA_additive_refine_sub_sd, chr="Lg1", qtl.index = 1, expandtomarkers = TRUE)
bi_qtl_5a_rMA_additive_fit_refine_lg3_marke_sub_sd<-bayesint(results=qtl_5a_rMA_additive_refine_sub_sd, chr="Lg3", qtl.index = 2, expandtomarkers = TRUE)
bi_qtl_5a_rMA_additive_fit_refine_lg4_marke_sub_sd<-bayesint(results=qtl_5a_rMA_additive_refine_sub_sd, chr="Lg4", qtl.index = 3, expandtomarkers = TRUE)
bi_qtl_5a_rMA_additive_fit_refine_lg5_marke_sub_sd<-bayesint(results=qtl_5a_rMA_additive_refine_sub_sd, chr="Lg5", qtl.index = 4, expandtomarkers = TRUE)
bi_qtl_5a_rMA_additive_fit_refine_lg7_marke_sub_sd<-bayesint(results=qtl_5a_rMA_additive_refine_sub_sd, chr="Lg7", qtl.index = 5, expandtomarkers = TRUE)
bi_qtl_5a_rMA_additive_fit_refine_lg8_marke_sub_sd<-bayesint(results=qtl_5a_rMA_additive_refine_sub_sd, chr="Lg8", qtl.index = 6, expandtomarkers = TRUE)

qtl_table<-rbind(qtl_table,
                  data.frame(trait="rMA",chr=bi_qtl_5a_rMA_additive_fit_refine_lg8_marke_sub_sd$chr[1], begin=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg8_marke_sub_sd)[1], 
                             center=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg8_marke_sub_sd)[2], 
                             end=row.names(bi_qtl_5a_rMA_additive_fit_refine_lg8_marke_sub_sd)[3]))
#saveRDS(qtl_table, "data/cross_5A_qtls.Rdata")





#write.table(rMA_table_5A_sub,"data/5A_rma.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



pdf("plots/5A_single_qtl_final.pdf", width=13, height=6)

plot(out_5A.hk_FRE_sub, out_5A.hk_MRE_sub, out_5A.hk_rMA_sub  ,lodcolumn = 1, col=c("blue", "red", "darkgreen"), ylab="LOD", xlab="")
abline(h=summary(operm_5A.hk_FRE_sub)[1], col="blue", lty=2, lwd=2)
abline(h=summary(operm_5A.hk_MRE_sub)[1], col="red", lty=2, lwd=2)
abline(h=summary(operm_5A.hk_rMA_sub)[1], col="darkgreen", lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_5A.hk_FRE_sub, thechr = summary(out_5A.hk_FRE_sub)[4,1], summary(out_5A.hk_FRE_sub)[4,2]), col="blue", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_FRE_sub, thechr = summary(out_5A.hk_FRE_sub)[7,1], summary(out_5A.hk_FRE_sub)[7,2]), col="blue", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_FRE_sub, thechr = summary(out_5A.hk_FRE_sub)[8,1], summary(out_5A.hk_FRE_sub)[8,2]), col="blue", lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_5A.hk_MRE_sub, thechr = summary(out_5A.hk_MRE_sub)[1,1], summary(out_5A.hk_MRE_sub)[1,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_MRE_sub, thechr = summary(out_5A.hk_MRE_sub)[3,1], summary(out_5A.hk_MRE_sub)[3,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_MRE_sub, thechr = summary(out_5A.hk_MRE_sub)[4,1], summary(out_5A.hk_MRE_sub)[4,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_MRE_sub, thechr = summary(out_5A.hk_MRE_sub)[5,1], summary(out_5A.hk_MRE_sub)[5,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_MRE_sub, thechr = summary(out_5A.hk_MRE_sub)[7,1], summary(out_5A.hk_MRE_sub)[7,2]), col="red", lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_5A.hk_rMA_sub, thechr = summary(out_5A.hk_rMA_sub)[1,1], summary(out_5A.hk_rMA_sub)[1,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_rMA_sub, thechr = summary(out_5A.hk_rMA_sub)[3,1], summary(out_5A.hk_rMA_sub)[3,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_rMA_sub, thechr = summary(out_5A.hk_rMA_sub)[4,1], summary(out_5A.hk_rMA_sub)[4,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_rMA_sub, thechr = summary(out_5A.hk_rMA_sub)[5,1], summary(out_5A.hk_rMA_sub)[5,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_rMA_sub, thechr = summary(out_5A.hk_rMA_sub)[7,1], summary(out_5A.hk_rMA_sub)[7,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_5A.hk_rMA_sub, thechr = summary(out_5A.hk_rMA_sub)[8,1], summary(out_5A.hk_rMA_sub)[8,2]), col="darkgreen", lty=2, lwd=2)

dev.off()