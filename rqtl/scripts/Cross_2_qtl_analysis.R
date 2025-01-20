  setwd("/home/fredo/Desktop/qtl_paper/upload_suppl_information/rqtl/")
  library(qtl)
  library(knitr)
  
  cross_obj<-read.cross(format="csvr", file ="data/S1.csvr", alleles=c("A","B"))
  cross_obj<-jittermap(cross_obj)
  cross_obj<-calc.genoprob(cross_obj, step=1)
  
  qtl_table=data.frame(trait=c(),chr=c(), marker=c(), begin=c(), center=c(),end=c())
  
  
  
  #pull.map(cross_obj, as.table = TRUE)
  #save.image("data/cross_S1.Rdata")
  out_S1.hk_FRE <- scanone(cross_obj, pheno.col = 2, method="hk")
  operm_S1.hk_FRE <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 2)
  summary(operm_S1.hk_FRE)
  summary(out_S1.hk_FRE)
  
  #Here I am importing the results of the permutations for the two-QTL analysis, which I ran in chunks on the cluster using the following command: 
  #scantwo_perm <- scantwopermhk(cross_obj, pheno.col = 2, n.perm=100, verbose=TRUE, batchsize = 20)
  x<-readRDS("data/cross_S1_out_FRE/1.rds")
  for(i in 2:100){
    st<-paste("data/cross_S1_out_FRE/", i, ".rds", sep="")
    sobj<-readRDS(st)
    x<-c(x, sobj)
  }
  #plot histograms of permutations
  plot(x)
  out_S1.hk_FRE_two <- scantwo(cross_obj, pheno.col = 2, method="hk")
  
  fre_table_S1<-summary(out_S1.hk_FRE_two, perms=x, pvalues = TRUE)
  #write.table(fre_table_S1,"data/S1_fre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
  

  kable(fre_table_S1, caption="Output of the two-QTL analysis", digits = 3)
  
  ### Multi-QTL analysis
  qtl_S1_FRE_additive<-makeqtl(cross_obj, chr="Lg7", pos=21.35,  what="prob")
  qtl_S1_FRE_additive_fit<-fitqtl(cross_obj, qtl=qtl_S1_FRE_additive, method="hk", pheno.col = 2, formula=y~Q1)
  qtl_S1_FRE_additive_refine<-refineqtl(cross_obj, qtl=qtl_S1_FRE_additive, method="hk", pheno.col = 2, formula=y~Q1, keeplodprofile = TRUE)
  qtl_S1_FRE_additive_fit_refine<-fitqtl(cross_obj, qtl=qtl_S1_FRE_additive_refine, method="hk", pheno.col = 2, formula=y~Q1)
  qtl_S1_FRE_additive_fit_refine_ests<-fitqtl(cross_obj, qtl=qtl_S1_FRE_additive_refine, method="hk", pheno.col = 2, formula=y~Q1, dropone=FALSE, get.ests=TRUE)
  
  
  summary(qtl_S1_FRE_additive_fit_refine)
  summary(qtl_S1_FRE_additive_fit_refine_ests)
  effectplot(cross = cross_obj, mname1 = "Lg7@24.4", pheno.col = 2, main="Effect plot for FRE on Lg7")
  
  
  bi_qtl_S1_FRE_additive_fit_refine_lg7<-bayesint(results=qtl_S1_FRE_additive_refine, chr="Lg7", qtl.index = 1)
  bi_qtl_S1_FRE_additive_fit_refine_lg7_marker<-bayesint(results=qtl_S1_FRE_additive_refine, chr="Lg7", qtl.index = 1, expandtomarkers = TRUE)
  
  kable(bi_qtl_S1_FRE_additive_fit_refine_lg7, caption="Baysian interval FRE QTL on Lg7", digits = 3)
  kable(bi_qtl_S1_FRE_additive_fit_refine_lg7_marker, caption="Baysian interval FRE QTL on Lg7 (physical loci)", digits = 3)
  
  
  qtl_table<-rbind(qtl_table,
                   data.frame(trait="FRE",chr=bi_qtl_S1_FRE_additive_fit_refine_lg7_marker$chr[1], begin=row.names(bi_qtl_S1_FRE_additive_fit_refine_lg7_marker)[1], 
                              center=row.names(bi_qtl_S1_FRE_additive_fit_refine_lg7_marker)[2], 
                              end=row.names(bi_qtl_S1_FRE_additive_fit_refine_lg7_marker)[3]))
  #MRE
  
  
  #run and plot FRE
  out_S1.hk_MRE <- scanone(cross_obj, pheno.col = 3, method="hk")
  operm_S1.hk_MRE <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 3)
  summary(operm_S1.hk_MRE)
  summary(out_S1.hk_MRE)
    
### Two-QTL analysis
    
    #Here I am importing the results of the permutations for the two-QTL analysis, which I ran in chunks on the cluster using the following command: 
    #scantwo_perm <- scantwopermhk(cross_obj, pheno.col = 3, n.perm=100, verbose=TRUE, batchsize = 20)
    x<-readRDS("data/cross_S1_out_MRE/1.rds")
    for(i in 2:100){
      st<-paste("data/cross_S1_out_MRE/", i, ".rds", sep="")
      sobj<-readRDS(st)
      x<-c(x, sobj)
    }
    plot(x)
    
    out_S1.hk_MRE_two <- scantwo(cross_obj, pheno.col = 3, method="hk")
    
    mre_table_S1<-summary(out_S1.hk_MRE_two, perms=x, pvalues = TRUE)
    #write.table(mre_table_S1,"data/S1_mre.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    
  kable(mre_table_S1, caption="Output of the two-QTL analysis", digits = 3)
  
  qtl_S1_MRE_additive<-makeqtl(cross_obj, chr=c("Lg4", "Lg7"), pos=c(42,8),  what="prob")
  qtl_S1_MRE_additive_fit<-fitqtl(cross_obj, qtl=qtl_S1_MRE_additive, method="hk", pheno.col = 3, formula=y~Q1+Q2)
  qtl_S1_MRE_additive_refine<-refineqtl(cross_obj, qtl=qtl_S1_MRE_additive, method="hk", pheno.col = 3, formula=y~Q1+Q2, keeplodprofile = TRUE)
  qtl_S1_MRE_additive_fit_refine<-fitqtl(cross_obj, qtl=qtl_S1_MRE_additive_refine, method="hk", pheno.col = 3, formula=y~Q1+Q2)
  qtl_S1_MRE_additive_fit_refine_ests<-fitqtl(cross_obj, qtl=qtl_S1_MRE_additive_refine, method="hk", pheno.col = 3, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)
  
  summary(qtl_S1_MRE_additive_fit_refine)
  summary(qtl_S1_MRE_additive_fit_refine_ests)
  effectplot(cross = cross_obj, mname1 = "Lg4@42.1", pheno.col = 3, main="Effect plot for MRE QTL on Lg4")
  effectplot(cross = cross_obj, mname1 = "Lg7@7.9", pheno.col = 3, main="Effect plot for MRE QTL on Lg7")
  
  bi_qtl_S1_MRE_additive_fit_refine_lg4<-bayesint(results=qtl_S1_MRE_additive_refine, chr="Lg4", qtl.index = 1)
  bi_qtl_S1_MRE_additive_fit_refine_lg4_marker<-bayesint(results=qtl_S1_MRE_additive_refine, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)
  kable(bi_qtl_S1_MRE_additive_fit_refine_lg4, caption="Baysian interval MRE QTL on Lg4", digits = 3)
  kable(bi_qtl_S1_MRE_additive_fit_refine_lg4_marker, caption="Baysian interval MRE QTL on Lg4 (physical loci)", digits = 3)
  bi_qtl_S1_MRE_additive_fit_refine_lg7<-bayesint(results=qtl_S1_MRE_additive_refine, chr="Lg7", qtl.index = 2)
  bi_qtl_S1_MRE_additive_fit_refine_lg7_marker<-bayesint(results=qtl_S1_MRE_additive_refine, chr="Lg7", qtl.index = 2, expandtomarkers = TRUE)
  kable(bi_qtl_S1_MRE_additive_fit_refine_lg7, caption="Baysian interval MRE QTL on Lg7", digits = 3)
  kable(bi_qtl_S1_MRE_additive_fit_refine_lg7_marker, caption="Baysian interval MRE QTL on Lg7 (physical loci)", digits = 3)
  
  
  qtl_table<-rbind(qtl_table,
                   data.frame(trait="MRE",chr=bi_qtl_S1_MRE_additive_fit_refine_lg4_marker$chr[1], begin=row.names(bi_qtl_S1_MRE_additive_fit_refine_lg4_marker)[1], 
                              center=row.names(bi_qtl_S1_MRE_additive_fit_refine_lg4_marker)[2], 
                              end=row.names(bi_qtl_S1_MRE_additive_fit_refine_lg4_marker)[3]))
  qtl_table<-rbind(qtl_table,
                   data.frame(trait="MRE",chr=bi_qtl_S1_MRE_additive_fit_refine_lg7_marker$chr[1], begin=row.names(bi_qtl_S1_MRE_additive_fit_refine_lg7_marker)[1], 
                              center=row.names(bi_qtl_S1_MRE_additive_fit_refine_lg7_marker)[2], 
                              end=row.names(bi_qtl_S1_MRE_additive_fit_refine_lg7_marker)[3]))
  
  
  
  ## Trait: relative male allocation rMA
  
  
  ### Single-QTL analysis
  
  out_S1.hk_rMA<- scanone(cross_obj, pheno.col = 4, method="hk")
  operm_S1.hk_rMA <- scanone(cross_obj, method="hk", n.perm=1000, pheno.col = 4)
  summary(operm_S1.hk_rMA)
  
  ### Two-QTL analysis
  
  x<-readRDS("data/cross_S1_out_rMA/1.rds")
  for(i in 2:100){
    st<-paste("data/cross_S1_out_rMA/", i, ".rds", sep="")
    sobj<-readRDS(st)
    x<-c(x, sobj)
  }
  plot(x)
  out_S1.hk_rMA_two <- scantwo(cross_obj, pheno.col = 4, method="hk")
  rMA_table_S1<-summary(out_S1.hk_rMA_two, perms=x, pvalues = TRUE)
  #write.table(rMA_table_S1,"data/S1_rma.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
  
  kable(rMA_table_S1, caption="Output of the two-QTL analysis", digits = 3)
  
  ### Multi-QTL analysis
  
  qtl_S1_rMA_additive<-makeqtl(cross_obj, chr=c("Lg4", "Lg7"), pos=c(41,7.92),  what="prob")
  qtl_S1_rMA_additive_fit<-fitqtl(cross_obj, qtl=qtl_S1_rMA_additive, method="hk", pheno.col = 4, formula=y~Q1+Q2)
  qtl_S1_rMA_additive_refine<-refineqtl(cross_obj, qtl=qtl_S1_rMA_additive, method="hk", pheno.col = 4, formula=y~Q1+Q2, keeplodprofile = TRUE)
  qtl_S1_rMA_additive_fit_refine<-fitqtl(cross_obj, qtl=qtl_S1_rMA_additive_refine, method="hk", pheno.col = 4, formula=y~Q1+Q2)
  qtl_S1_rMA_additive_fit_refine_ests<-fitqtl(cross_obj, qtl=qtl_S1_rMA_additive_refine, method="hk", pheno.col = 4, formula=y~Q1+Q2, dropone=FALSE, get.ests=TRUE)
  summary(qtl_S1_rMA_additive_fit_refine)
  summary(qtl_S1_rMA_additive_fit_refine_ests)
  effectplot(cross = cross_obj, mname1 = "Lg4@40", pheno.col = 4, main="Effect plot for rMA QTL on Lg4")
  effectplot(cross = cross_obj, mname1 = "Lg7@7.9", pheno.col = 4, main="Effect plot for rMA QTL on Lg7")
  
  bi_qtl_S1_rMA_additive_fit_refine_lg4<-bayesint(results=qtl_S1_rMA_additive_refine, chr="Lg4", qtl.index = 1)
  bi_qtl_S1_rMA_additive_fit_refine_lg4_marker<-bayesint(results=qtl_S1_rMA_additive_refine, chr="Lg4", qtl.index = 1, expandtomarkers = TRUE)
  kable(bi_qtl_S1_rMA_additive_fit_refine_lg4, caption="Baysian interval rMA QTL on Lg4", digits = 3)
  kable(bi_qtl_S1_rMA_additive_fit_refine_lg4_marker, caption="Baysian interval rMA QTL on Lg4 (physical loci)", digits = 3)
  bi_qtl_S1_rMA_additive_fit_refine_lg7<-bayesint(results=qtl_S1_rMA_additive_refine, chr="Lg7", qtl.index = 2)
  bi_qtl_S1_rMA_additive_fit_refine_lg7_marker<-bayesint(results=qtl_S1_rMA_additive_refine, chr="Lg7", qtl.index = 2, expandtomarkers = TRUE)
  kable(bi_qtl_S1_rMA_additive_fit_refine_lg7, caption="Baysian interval rMA QTL on Lg7", digits = 3)
  kable(bi_qtl_S1_rMA_additive_fit_refine_lg7_marker, caption="Baysian interval rMA QTL on Lg7 (physical loci)", digits = 3)
  
  qtl_table<-rbind(qtl_table,
                   data.frame(trait="rMA",chr=bi_qtl_S1_rMA_additive_fit_refine_lg4_marker$chr[1], begin=row.names(bi_qtl_S1_rMA_additive_fit_refine_lg4_marker)[1], 
                              center=row.names(bi_qtl_S1_rMA_additive_fit_refine_lg4_marker)[2], 
                              end=row.names(bi_qtl_S1_rMA_additive_fit_refine_lg4_marker)[3]))
  qtl_table<-rbind(qtl_table,
                   data.frame(trait="rMA",chr=bi_qtl_S1_rMA_additive_fit_refine_lg7_marker$chr[1], begin=row.names(bi_qtl_S1_rMA_additive_fit_refine_lg7_marker)[1], 
                              center=row.names(bi_qtl_S1_rMA_additive_fit_refine_lg7_marker)[2], 
                              end=row.names(bi_qtl_S1_rMA_additive_fit_refine_lg7_marker)[3]))
  
  
  
#  saveRDS(qtl_table, "cross_S1_qtls.Rdata")

#############FINAL PLOTS!!!!!!!!!!!!!!!!!!!!!!!!####################################################


max(pull.map(cross_obj, chr="Lg1", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg2", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg3", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg4", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg5", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg6", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg7", as.table = TRUE)$pos)+
  max(pull.map(cross_obj, chr="Lg8", as.table = TRUE)$pos)


pdf("plots/S1_single_qtl_final.pdf", width=13, height=6)
plot(out_S1.hk_FRE, out_S1.hk_MRE, out_S1.hk_rMA  ,lodcolumn = 1, col=c("blue", "red", "darkgreen"), ylab="LOD", xlab="")
abline(h=summary(operm_S1.hk_FRE)[1], col="blue", lty=2, lwd=2)
abline(h=summary(operm_S1.hk_MRE)[1], col="red", lty=2, lwd=2)
abline(h=summary(operm_S1.hk_rMA)[1], col="darkgreen", lty=2, lwd=2)

abline(v=xaxisloc.scanone(out_S1.hk_FRE, thechr = summary(out_S1.hk_FRE)[7,1], summary(out_S1.hk_FRE)[7,2]), col="blue", lty=2, lwd=2)


abline(v=xaxisloc.scanone(out_S1.hk_MRE, thechr = summary(out_S1.hk_MRE)[4,1], summary(out_S1.hk_MRE)[4,2]), col="red", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_S1.hk_MRE, thechr = summary(out_S1.hk_MRE)[7,1], summary(out_S1.hk_MRE)[7,2]), col="red", lty=2, lwd=2)



abline(v=xaxisloc.scanone(out_S1.hk_MRE, thechr = summary(out_S1.hk_rMA)[4,1], summary(out_S1.hk_rMA)[4,2]), col="darkgreen", lty=2, lwd=2)
abline(v=xaxisloc.scanone(out_S1.hk_MRE, thechr = summary(out_S1.hk_rMA)[7,1], summary(out_S1.hk_rMA)[7,2]), col="darkgreen", lty=2, lwd=2)



dev.off()



