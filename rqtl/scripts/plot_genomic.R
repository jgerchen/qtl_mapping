setwd("/home/fredo/Desktop/qtl_paper/upload_suppl_information/rqtl/")
library(qtl)
library(stringr)
#Genome information
#Chr1 OW569319.1: 74045953
#Chr2 OW569313.1: 58100813
#Chr3 OW569312.1: 76280018
#Chr4 OW569317.1: 41830924 
#Chr5 OW569314.1: 56051264 
#Chr6 OW569316.1: 47723244 
#Chr7 OW569315.1: 50031057 
#Chr8 OW569318.1: 43611580

genome_info_df<-data.frame(chr=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8"),
                           Lg=c("OW569319.1","OW569313.1","OW569312.1","OW569317.1","OW569314.1","OW569316.1","OW569315.1","OW569318.1"),
                           size=c(74045953, 58100813, 76280018, 41830924, 56051264, 47723244, 50031057, 43611580),
                           cum_size=1:8)
curr_size=0
for(chr_size_i in 1:8){
  genome_info_df$cum_size[chr_size_i]<-curr_size
  curr_size<-curr_size+genome_info_df$size[chr_size_i]
}

centromeric_repeats<-read.csv("data/Summary.of.repetitive.regions.GCA_937616625.1_ddMerAnnu1.1_genomic.fna.csv", header=TRUE)
most_common_centromeric_repeats<-centromeric_repeats[centromeric_repeats$most.freq.value.N==60 | centromeric_repeats$most.freq.value.N==99 ,]
most_common_centromeric_repeats$cum_start<-rep(NA, length(most_common_centromeric_repeats$name))
most_common_centromeric_repeats$cum_end<-rep(NA, length(most_common_centromeric_repeats$name))

for(cent_i in 1:length(most_common_centromeric_repeats$name)){
  most_common_centromeric_repeats$cum_start[cent_i]<-genome_info_df[genome_info_df$Lg==most_common_centromeric_repeats$name[cent_i],]$cum_size+most_common_centromeric_repeats$start[cent_i]
  most_common_centromeric_repeats$cum_end[cent_i]<-genome_info_df[genome_info_df$Lg==most_common_centromeric_repeats$name[cent_i],]$cum_size+most_common_centromeric_repeats$end[cent_i]
}

####Cross 1A
cross_obj_1A<-read.cross(format="csvr", file ="data/1A.csvr", alleles=c("A","B"))
#cross_obj_1A<-jittermap(cross_obj_1A)
cross_obj_1A<-calc.genoprob(cross_obj_1A, step=1)

genotypes_1A<-pull.geno(cross_obj_1A)
phenotypes_1A<-pull.pheno(cross_obj_1A)

cross_1A_df=data.frame(id=colnames(genotypes_1A), 
                       chrom=unlist(strsplit(colnames(genotypes_1A), split="_"))[ c(TRUE,FALSE) ], 
                       pos=as.numeric(unlist(strsplit(colnames(genotypes_1A), split="_"))[ c(FALSE,TRUE) ]),
                       FRE_mean_AA=rep(NA, length(colnames(genotypes_1A))),
                       FRE_mean_AB=rep(NA, length(colnames(genotypes_1A))), 
                       FRE_mean_BB=rep(NA, length(colnames(genotypes_1A))),
                       MRE_mean_AA=rep(NA, length(colnames(genotypes_1A))),
                       MRE_mean_AB=rep(NA, length(colnames(genotypes_1A))), 
                       MRE_mean_BB=rep(NA, length(colnames(genotypes_1A))),
                       rMA_mean_AA=rep(NA, length(colnames(genotypes_1A))),
                       rMA_mean_AB=rep(NA, length(colnames(genotypes_1A))), 
                       rMA_mean_BB=rep(NA, length(colnames(genotypes_1A))))

for(locus in 1:dim(genotypes_1A)[2]){
  cross_1A_df$FRE_mean_AA[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==1,2])
  cross_1A_df$FRE_mean_AB[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==2,2])
  cross_1A_df$FRE_mean_BB[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==3,2])
  cross_1A_df$MRE_mean_AA[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==1,3])
  cross_1A_df$MRE_mean_AB[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==2,3])
  cross_1A_df$MRE_mean_BB[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==3,3])
  cross_1A_df$rMA_mean_AA[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==1,4])
  cross_1A_df$rMA_mean_AB[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==2,4])
  cross_1A_df$rMA_mean_BB[locus]<-mean(phenotypes_1A[genotypes_1A[,locus]==3,4])
}

cross_1A_df$pos_cum<-cross_1A_df$pos
for(loc_i in 1:length(cross_1A_df$pos)){
  cross_1A_df$pos_cum[loc_i]<-genome_info_df[genome_info_df$Lg==cross_1A_df$chrom[loc_i],]$cum_size+cross_1A_df$pos[loc_i]
}
cross_1A_lm<-pull.map(cross_obj_1A, as.table = TRUE)
cross_1A_lm$pos_cum<-cross_1A_df$pos_cum
cross_1A_qtls<-readRDS("data/cross_1A_qtls.Rdata")

gt_1A<-geno.table(cross_obj_1A, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_1A<--log10(0.05/totmar(cross_obj_1A))


#plots here
#1A
pdf("plots/1A_plot.pdf", width = 8, height=12)
par(mfrow=c(5,1))
#bottom, left, top, right
#default_mar<-c(5.1, 4.1, 4.1, 2.1)
par(mar=c(0, 4.1,0.5, 2.1))

#FRE
plot(cross_1A_df$pos[cross_1A_df$chrom=="OW569319.1"], cross_1A_df$FRE_mean_AA[cross_1A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_1A_df$FRE_mean_AA,cross_1A_df$FRE_mean_AB,cross_1A_df$FRE_mean_BB)), max(c(cross_1A_df$FRE_mean_AA,cross_1A_df$FRE_mean_AB,cross_1A_df$FRE_mean_BB))), ylab="mean FRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?

#QTL1
cross_1A_FRE_qtl1_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[1]),]
cross_1A_FRE_qtl1_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[1]),]
cross_1A_FRE_qtl1_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[1]),]
cross_1A_FRE_qtl1_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_FRE_qtl1_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_FRE_qtl1_left_edge$pos_cum[1] & cross_1A_df$pos_cum<=cross_1A_FRE_qtl1_right_edge$pos_cum[1],]
cross_1A_FRE_qtl1_lines_AA<-lines(cross_1A_FRE_qtl1_lines_df$pos_cum, cross_1A_FRE_qtl1_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_1A_FRE_qtl1_lines_BB<-lines(cross_1A_FRE_qtl1_lines_df$pos_cum, cross_1A_FRE_qtl1_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_FRE_qtl1_lines_df$pos_cum, rev(cross_1A_FRE_qtl1_lines_df$pos_cum)), y=c(cross_1A_FRE_qtl1_lines_df$FRE_mean_BB, rev(cross_1A_FRE_qtl1_lines_df$FRE_mean_AA)),  col=rgb(0,0,1.0,alpha=0.5), lty=0)
abline(v=cross_1A_FRE_qtl1_center$pos_cum, col="blue", lty=2, lwd=2)

#QTL2
cross_1A_FRE_qtl2_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[2]),]
cross_1A_FRE_qtl2_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[2]),]
cross_1A_FRE_qtl2_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[2]),]
cross_1A_FRE_qtl2_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_FRE_qtl2_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_FRE_qtl2_left_edge$pos_cum[1] & cross_1A_df$pos_cum<=cross_1A_FRE_qtl2_right_edge$pos_cum[1],]
cross_1A_FRE_qtl2_lines_AA<-lines(cross_1A_FRE_qtl2_lines_df$pos_cum, cross_1A_FRE_qtl2_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_1A_FRE_qtl2_lines_BB<-lines(cross_1A_FRE_qtl2_lines_df$pos_cum, cross_1A_FRE_qtl2_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_FRE_qtl2_lines_df$pos_cum, rev(cross_1A_FRE_qtl2_lines_df$pos_cum)), y=c(cross_1A_FRE_qtl2_lines_df$FRE_mean_AA, rev(cross_1A_FRE_qtl2_lines_df$FRE_mean_BB)),  col=rgb(0,0,1.0,alpha=0.5), lty=0)
abline(v=cross_1A_FRE_qtl2_center$pos_cum, col="blue", lty=2, lwd=2)

abline(h=mean(phenotypes_1A$FRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$FRE_mean_AA[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$FRE_mean_AB[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$FRE_mean_BB[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}


#plot results
#MRE
plot(cross_1A_df$pos[cross_1A_df$chrom=="OW569319.1"], cross_1A_df$MRE_mean_AA[cross_1A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_1A_df$MRE_mean_AA,cross_1A_df$MRE_mean_AB,cross_1A_df$MRE_mean_BB)), max(c(cross_1A_df$MRE_mean_AA,cross_1A_df$MRE_mean_AB,cross_1A_df$MRE_mean_BB))), ylab="mean MRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)


#QTL3
cross_1A_MRE_qtl3_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[3]),]
cross_1A_MRE_qtl3_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[3]),]
cross_1A_MRE_qtl3_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[3]),]
cross_1A_MRE_qtl3_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_MRE_qtl3_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_MRE_qtl3_left_edge$pos_cum[1] 
                                        & cross_1A_df$pos_cum<=cross_1A_MRE_qtl3_right_edge$pos_cum[1],]
cross_1A_MRE_qtl3_lines_AA<-lines(cross_1A_MRE_qtl3_lines_df$pos_cum, cross_1A_MRE_qtl3_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_1A_MRE_qtl3_lines_BB<-lines(cross_1A_MRE_qtl3_lines_df$pos_cum, cross_1A_MRE_qtl3_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_MRE_qtl3_lines_df$pos_cum, rev(cross_1A_MRE_qtl3_lines_df$pos_cum)), 
        y=c(cross_1A_MRE_qtl3_lines_df$MRE_mean_BB, 
        rev(cross_1A_MRE_qtl3_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)
abline(v=cross_1A_MRE_qtl3_center$pos_cum, col="red", lty=2, lwd=2)

#QTL4

cross_1A_MRE_qtl4_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[4]),]
cross_1A_MRE_qtl4_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[4]),]
cross_1A_MRE_qtl4_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[4]),]
cross_1A_MRE_qtl4_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_MRE_qtl4_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_MRE_qtl4_left_edge$pos_cum[1] 
                                        & cross_1A_df$pos_cum<=cross_1A_MRE_qtl4_right_edge$pos_cum[1],]
cross_1A_MRE_qtl4_lines_AA<-lines(cross_1A_MRE_qtl4_lines_df$pos_cum, cross_1A_MRE_qtl4_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_1A_MRE_qtl4_lines_BB<-lines(cross_1A_MRE_qtl4_lines_df$pos_cum, cross_1A_MRE_qtl4_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_MRE_qtl4_lines_df$pos_cum, rev(cross_1A_MRE_qtl4_lines_df$pos_cum)), 
        y=c(cross_1A_MRE_qtl4_lines_df$MRE_mean_BB, 
            rev(cross_1A_MRE_qtl4_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)
abline(v=cross_1A_MRE_qtl4_center$pos_cum, col="red", lty=2, lwd=2)

#QTL5


cross_1A_MRE_qtl5_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[5]),]
cross_1A_MRE_qtl5_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[5]),]
cross_1A_MRE_qtl5_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[5]),]
cross_1A_MRE_qtl5_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_MRE_qtl5_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_MRE_qtl5_left_edge$pos_cum[1] 
                                        & cross_1A_df$pos_cum<=cross_1A_MRE_qtl5_right_edge$pos_cum[1],]
cross_1A_MRE_qtl5_lines_AA<-lines(cross_1A_MRE_qtl5_lines_df$pos_cum, cross_1A_MRE_qtl5_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_1A_MRE_qtl5_lines_BB<-lines(cross_1A_MRE_qtl5_lines_df$pos_cum, cross_1A_MRE_qtl5_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_MRE_qtl5_lines_df$pos_cum, rev(cross_1A_MRE_qtl5_lines_df$pos_cum)), 
        y=c(cross_1A_MRE_qtl5_lines_df$MRE_mean_AA, 
            rev(cross_1A_MRE_qtl5_lines_df$MRE_mean_BB)),  col=rgb(1,0,0.0,alpha=0.5), 
          lty=1,
        density=30,
        angle=45)
abline(v=cross_1A_MRE_qtl5_center$pos_cum, col="red", lty=2, lwd=2)

abline(h=mean(phenotypes_1A$MRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$MRE_mean_AA[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$MRE_mean_AB[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$MRE_mean_BB[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}





#rMA
plot(cross_1A_df$pos[cross_1A_df$chrom=="OW569319.1"], cross_1A_df$FRE_mean_AA[cross_1A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_1A_df$rMA_mean_AA,cross_1A_df$rMA_mean_AB,cross_1A_df$rMA_mean_BB)), max(c(cross_1A_df$rMA_mean_AA,cross_1A_df$rMA_mean_AB,cross_1A_df$rMA_mean_BB))), ylab="mean rMA", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)

#QTL6
cross_1A_rMA_qtl6_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[6]),]
cross_1A_rMA_qtl6_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[6]),]
cross_1A_rMA_qtl6_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[6]),]
cross_1A_rMA_qtl6_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_rMA_qtl6_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_rMA_qtl6_left_edge$pos_cum[1] 
                                        & cross_1A_df$pos_cum<=cross_1A_rMA_qtl6_right_edge$pos_cum[1],]
cross_1A_rMA_qtl6_lines_AA<-lines(cross_1A_rMA_qtl6_lines_df$pos_cum, cross_1A_rMA_qtl6_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_1A_MRE_qtl6_lines_BB<-lines(cross_1A_rMA_qtl6_lines_df$pos_cum, cross_1A_rMA_qtl6_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_rMA_qtl6_lines_df$pos_cum, rev(cross_1A_rMA_qtl6_lines_df$pos_cum)), 
        y=c(cross_1A_rMA_qtl6_lines_df$rMA_mean_AA, 
            rev(cross_1A_rMA_qtl6_lines_df$rMA_mean_BB)),  col=rgb(0,1,0.5,alpha=0.7), lty=0)
abline(v=cross_1A_rMA_qtl6_center$pos_cum, col="darkgreen", lty=2, lwd=2)


#QTL7
cross_1A_rMA_qtl7_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[7]),]
cross_1A_rMA_qtl7_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[7]),]
cross_1A_rMA_qtl7_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[7]),]
cross_1A_rMA_qtl7_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_rMA_qtl7_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_rMA_qtl7_left_edge$pos_cum[1] 
                                        & cross_1A_df$pos_cum<=cross_1A_rMA_qtl7_right_edge$pos_cum[1],]
cross_1A_rMA_qtl7_lines_AA<-lines(cross_1A_rMA_qtl7_lines_df$pos_cum, cross_1A_rMA_qtl7_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_1A_MRE_qtl7_lines_BB<-lines(cross_1A_rMA_qtl7_lines_df$pos_cum, cross_1A_rMA_qtl7_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_rMA_qtl7_lines_df$pos_cum, rev(cross_1A_rMA_qtl7_lines_df$pos_cum)), 
        y=c(cross_1A_rMA_qtl7_lines_df$rMA_mean_AA, 
            rev(cross_1A_rMA_qtl7_lines_df$rMA_mean_BB)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)
abline(v=cross_1A_rMA_qtl7_center$pos_cum, col="darkgreen", lty=2, lwd=2)

#QTL8
cross_1A_rMA_qtl8_left_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$begin[8]),]
cross_1A_rMA_qtl8_center<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$center[8]),]
cross_1A_rMA_qtl8_right_edge<-cross_1A_df[as.character(cross_1A_df$id)==as.character(cross_1A_qtls$end[8]),]
cross_1A_rMA_qtl8_lines_df<-cross_1A_df[as.character(cross_1A_df$chrom)==cross_1A_rMA_qtl8_left_edge$chrom[1] & cross_1A_df$pos_cum>=cross_1A_rMA_qtl8_left_edge$pos_cum[1] 
                                        & cross_1A_df$pos_cum<=cross_1A_rMA_qtl8_right_edge$pos_cum[1],]
cross_1A_rMA_qtl8_lines_AA<-lines(cross_1A_rMA_qtl8_lines_df$pos_cum, cross_1A_rMA_qtl8_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_1A_MRE_qtl8_lines_BB<-lines(cross_1A_rMA_qtl8_lines_df$pos_cum, cross_1A_rMA_qtl8_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_1A_rMA_qtl8_lines_df$pos_cum, rev(cross_1A_rMA_qtl8_lines_df$pos_cum)), 
        y=c(cross_1A_rMA_qtl8_lines_df$rMA_mean_AA, 
            rev(cross_1A_rMA_qtl8_lines_df$rMA_mean_BB)),  col=rgb(0,1,0.5),
        lty=1,
        density=30,
        angle=45)
abline(v=cross_1A_rMA_qtl8_center$pos_cum, col="darkgreen", lty=2, lwd=2)

abline(h=mean(phenotypes_1A$rMA), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$rMA_mean_AA[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$rMA_mean_AB[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_1A_df$pos_cum[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], cross_1A_df$rMA_mean_BB[cross_1A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}


#plot segregation distortion
plot(cross_1A_df$pos[cross_1A_df$chrom=="OW569319.1"], cross_1A_df$FRE_mean_AA[cross_1A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, max(gt_1A$neglog10P)), ylab="segregation distortion (-10*log p)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)

for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_1A_lm$pos_cum[cross_1A_lm$chr==chr_i], gt_1A$neglog10P[gt_1A$chr==chr_i], col="black", lwd=2)
}


abline(h=bonferroni_threshold_1A, col="red", lwd=2)



abline(v=genome_info_df$cum_size)



#plot linkage map
par(mar=c(5.1, 4.1,0.5, 2.1))
plot(cross_1A_df$pos[cross_1A_df$chrom=="OW569319.1"], cross_1A_df$FRE_mean_AA[cross_1A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, max(cross_1A_lm$pos)), ylab="genetic distance (cM)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)

rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end, 
     ytop = rep(max(cross_1A_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end)), col=rgb(1,0,0,alpha=0.3), border=NA)

rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end, 
     ytop = rep(max(cross_1A_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end)), col=rgb(0,0,1,alpha=0.3), border=NA)



for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_1A_lm$pos_cum[cross_1A_lm$chr==chr_i], cross_1A_lm$pos[cross_1A_lm$chr==chr_i], col="black", lwd=2)
}


abline(v=genome_info_df$cum_size)

axis(1, at = genome_info_df$cum_size+(genome_info_df$size/2),
     labels = genome_info_df$chr, cex.axis=1.4)



dev.off()



####Cross S1

cross_obj_S1<-read.cross(format="csvr", file ="data/S1.csvr", alleles=c("A","B"))
#cross_obj_S1<-jittermap(cross_obj_S1)
cross_obj_S1<-calc.genoprob(cross_obj_S1, step=1)

genotypes_S1<-pull.geno(cross_obj_S1)
phenotypes_S1<-pull.pheno(cross_obj_S1)





cross_S1_df=data.frame(id=colnames(genotypes_S1), 
                       chrom=unlist(strsplit(colnames(genotypes_S1), split="_"))[ c(TRUE,FALSE) ], 
                       pos=as.numeric(unlist(strsplit(colnames(genotypes_S1), split="_"))[ c(FALSE,TRUE) ]),
                       FRE_mean_AA=rep(NA, length(colnames(genotypes_S1))),
                       FRE_mean_AB=rep(NA, length(colnames(genotypes_S1))), 
                       FRE_mean_BB=rep(NA, length(colnames(genotypes_S1))),
                       MRE_mean_AA=rep(NA, length(colnames(genotypes_S1))),
                       MRE_mean_AB=rep(NA, length(colnames(genotypes_S1))), 
                       MRE_mean_BB=rep(NA, length(colnames(genotypes_S1))),
                       rMA_mean_AA=rep(NA, length(colnames(genotypes_S1))),
                       rMA_mean_AB=rep(NA, length(colnames(genotypes_S1))), 
                       rMA_mean_BB=rep(NA, length(colnames(genotypes_S1))))

for(locus in 1:dim(genotypes_S1)[2]){
  cross_S1_df$FRE_mean_AA[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==1,2])
  cross_S1_df$FRE_mean_AB[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==2,2])
  cross_S1_df$FRE_mean_BB[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==3,2])
  cross_S1_df$MRE_mean_AA[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==1,3])
  cross_S1_df$MRE_mean_AB[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==2,3])
  cross_S1_df$MRE_mean_BB[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==3,3])
  cross_S1_df$rMA_mean_AA[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==1,4])
  cross_S1_df$rMA_mean_AB[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==2,4])
  cross_S1_df$rMA_mean_BB[locus]<-mean(phenotypes_S1[genotypes_S1[,locus]==3,4])
}

cross_S1_df$pos_cum<-cross_S1_df$pos
for(loc_i in 1:length(cross_S1_df$pos)){
  cross_S1_df$pos_cum[loc_i]<-genome_info_df[genome_info_df$Lg==cross_S1_df$chrom[loc_i],]$cum_size+cross_S1_df$pos[loc_i]
}


cross_S1_lm<-pull.map(cross_obj_S1, as.table = TRUE)
cross_S1_lm$pos_cum<-cross_S1_df$pos_cum
cross_S1_qtls<-readRDS("data/cross_S1_qtls.Rdata")


gt_S1<-geno.table(cross_obj_S1, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_S1<--log10(0.05/totmar(cross_obj_S1))


#plot

pdf("plots/S1_plot.pdf", width = 8, height=12)
par(mfrow=c(5,1))
#bottom, left, top, right
#default_mar<-c(5.1, 4.1, 4.1, 2.1)
par(mar=c(0, 4.1,0.5, 2.1))

#FRE
plot(cross_S1_df$pos[cross_S1_df$chrom=="OW569319.1"], cross_S1_df$FRE_mean_AA[cross_S1_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_S1_df$FRE_mean_AA,cross_S1_df$FRE_mean_AB,cross_S1_df$FRE_mean_BB)), max(c(cross_S1_df$FRE_mean_AA,cross_S1_df$FRE_mean_AB,cross_S1_df$FRE_mean_BB))), ylab="mean FRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?


#QTL1
cross_S1_FRE_qtl1_left_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$begin[1]),]
cross_S1_FRE_qtl1_center<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$center[1]),]
cross_S1_FRE_qtl1_right_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$end[1]),]
cross_S1_FRE_qtl1_lines_df<-cross_S1_df[as.character(cross_S1_df$chrom)==cross_S1_FRE_qtl1_left_edge$chrom[1] & cross_S1_df$pos_cum>=cross_S1_FRE_qtl1_left_edge$pos_cum[1] & cross_S1_df$pos_cum<=cross_S1_FRE_qtl1_right_edge$pos_cum[1],]
cross_S1_FRE_qtl1_lines_AA<-lines(cross_S1_FRE_qtl1_lines_df$pos_cum, cross_S1_FRE_qtl1_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_S1_FRE_qtl1_lines_BB<-lines(cross_S1_FRE_qtl1_lines_df$pos_cum, cross_S1_FRE_qtl1_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_S1_FRE_qtl1_lines_df$pos_cum, rev(cross_S1_FRE_qtl1_lines_df$pos_cum)), y=c(cross_S1_FRE_qtl1_lines_df$FRE_mean_BB, rev(cross_S1_FRE_qtl1_lines_df$FRE_mean_AA)),  col=rgb(0,0,1.0,alpha=0.5), lty=0)

abline(v=cross_S1_FRE_qtl1_right_edge$pos_cum, col="blue", lty=2, lwd=2)

abline(h=mean(phenotypes_S1$FRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$FRE_mean_AA[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$FRE_mean_AB[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$FRE_mean_BB[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}


#MRE
plot(cross_S1_df$pos[cross_S1_df$chrom=="OW569319.1"], cross_S1_df$MRE_mean_AA[cross_S1_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_S1_df$MRE_mean_AA,cross_S1_df$MRE_mean_AB,cross_S1_df$MRE_mean_BB)), max(c(cross_S1_df$MRE_mean_AA,cross_S1_df$MRE_mean_AB,cross_S1_df$MRE_mean_BB))), ylab="mean MRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?


#QTL2
cross_S1_MRE_qtl2_left_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$begin[2]),]
cross_S1_MRE_qtl2_center<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$center[2]),]
cross_S1_MRE_qtl2_right_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$end[2]),]
cross_S1_MRE_qtl2_lines_df<-cross_S1_df[as.character(cross_S1_df$chrom)==cross_S1_MRE_qtl2_left_edge$chrom[1] & cross_S1_df$pos_cum>=cross_S1_MRE_qtl2_left_edge$pos_cum[1] & cross_S1_df$pos_cum<=cross_S1_MRE_qtl2_right_edge$pos_cum[1],]
cross_S1_MRE_qtl2_lines_AA<-lines(cross_S1_MRE_qtl2_lines_df$pos_cum, cross_S1_MRE_qtl2_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_S1_MRE_qtl2_lines_BB<-lines(cross_S1_MRE_qtl2_lines_df$pos_cum, cross_S1_MRE_qtl2_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_S1_MRE_qtl2_lines_df$pos_cum, rev(cross_S1_MRE_qtl2_lines_df$pos_cum)), y=c(cross_S1_MRE_qtl2_lines_df$MRE_mean_BB, rev(cross_S1_MRE_qtl2_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)
abline(v=cross_S1_MRE_qtl2_center$pos_cum, col="red", lty=2, lwd=2)



#QTL3
cross_S1_MRE_qtl3_left_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$begin[3]),]
cross_S1_MRE_qtl3_center<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$center[3]),]
cross_S1_MRE_qtl3_right_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$end[3]),]
cross_S1_MRE_qtl3_lines_df<-cross_S1_df[as.character(cross_S1_df$chrom)==cross_S1_MRE_qtl3_left_edge$chrom[1] & cross_S1_df$pos_cum>=cross_S1_MRE_qtl3_left_edge$pos_cum[1] & cross_S1_df$pos_cum<=cross_S1_MRE_qtl3_right_edge$pos_cum[1],]
cross_S1_MRE_qtl3_lines_AA<-lines(cross_S1_MRE_qtl3_lines_df$pos_cum, cross_S1_MRE_qtl3_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_S1_MRE_qtl3_lines_BB<-lines(cross_S1_MRE_qtl3_lines_df$pos_cum, cross_S1_MRE_qtl3_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_S1_MRE_qtl3_lines_df$pos_cum, rev(cross_S1_MRE_qtl3_lines_df$pos_cum)), y=c(cross_S1_MRE_qtl3_lines_df$MRE_mean_BB, rev(cross_S1_MRE_qtl3_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)
abline(v=cross_S1_MRE_qtl3_center$pos_cum, col="red", lty=2, lwd=2)
abline(h=mean(phenotypes_S1$MRE), lty=2)
abline(v=genome_info_df$cum_size)

for(chr_i in 1:8){
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$MRE_mean_AA[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$MRE_mean_AB[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$MRE_mean_BB[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}

#rMA
plot(cross_S1_df$pos[cross_S1_df$chrom=="OW569319.1"], cross_S1_df$rMA_mean_AA[cross_S1_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_S1_df$rMA_mean_AA,cross_S1_df$rMA_mean_AB,cross_S1_df$rMA_mean_BB)), max(c(cross_S1_df$rMA_mean_AA,cross_S1_df$rMA_mean_AB,cross_S1_df$rMA_mean_BB))), ylab="mean rMA", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?


#QTL4
cross_S1_rMA_qtl4_left_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$begin[4]),]
cross_S1_rMA_qtl4_center<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$center[4]),]
cross_S1_rMA_qtl4_right_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$end[4]),]
cross_S1_rMA_qtl4_lines_df<-cross_S1_df[as.character(cross_S1_df$chrom)==cross_S1_rMA_qtl4_left_edge$chrom[1] & cross_S1_df$pos_cum>=cross_S1_rMA_qtl4_left_edge$pos_cum[1] & cross_S1_df$pos_cum<=cross_S1_rMA_qtl4_right_edge$pos_cum[1],]
cross_S1_rMA_qtl4_lines_AA<-lines(cross_S1_rMA_qtl4_lines_df$pos_cum, cross_S1_rMA_qtl4_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_S1_rMA_qtl4_lines_BB<-lines(cross_S1_rMA_qtl4_lines_df$pos_cum, cross_S1_rMA_qtl4_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_S1_rMA_qtl4_lines_df$pos_cum, rev(cross_S1_rMA_qtl4_lines_df$pos_cum)), y=c(cross_S1_rMA_qtl4_lines_df$rMA_mean_BB, rev(cross_S1_rMA_qtl4_lines_df$rMA_mean_AA)),  
        col=rgb(0,1,0.5,alpha=0.5), lty=0)
abline(v=cross_S1_rMA_qtl4_center$pos_cum, col="darkgreen", lty=2, lwd=2)



#QTL5
cross_S1_rMA_qtl5_left_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$begin[5]),]
cross_S1_rMA_qtl5_center<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$center[5]),]
cross_S1_rMA_qtl5_right_edge<-cross_S1_df[as.character(cross_S1_df$id)==as.character(cross_S1_qtls$end[5]),]
cross_S1_rMA_qtl5_lines_df<-cross_S1_df[as.character(cross_S1_df$chrom)==cross_S1_rMA_qtl5_left_edge$chrom[1] & cross_S1_df$pos_cum>=cross_S1_rMA_qtl5_left_edge$pos_cum[1] & cross_S1_df$pos_cum<=cross_S1_rMA_qtl5_right_edge$pos_cum[1],]
cross_S1_rMA_qtl5_lines_AA<-lines(cross_S1_rMA_qtl5_lines_df$pos_cum, cross_S1_rMA_qtl5_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_S1_rMA_qtl5_lines_BB<-lines(cross_S1_rMA_qtl5_lines_df$pos_cum, cross_S1_rMA_qtl5_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_S1_rMA_qtl5_lines_df$pos_cum, rev(cross_S1_rMA_qtl5_lines_df$pos_cum)), y=c(cross_S1_rMA_qtl5_lines_df$rMA_mean_BB, rev(cross_S1_rMA_qtl5_lines_df$rMA_mean_AA)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)
abline(v=cross_S1_rMA_qtl5_center$pos_cum, col="darkgreen", lty=2, lwd=2)

abline(h=mean(phenotypes_S1$rMA), lty=2)
abline(v=genome_info_df$cum_size)


for(chr_i in 1:8){
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$rMA_mean_AA[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$rMA_mean_AB[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_S1_df$pos_cum[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], cross_S1_df$rMA_mean_BB[cross_S1_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}

#plot segregation distortion
plot(cross_S1_df$pos[cross_S1_df$chrom=="OW569319.1"], cross_S1_df$FRE_mean_AA[cross_S1_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, 6), ylab="segregation distortion (-10*log p)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)

for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_S1_lm$pos_cum[cross_S1_lm$chr==chr_i], gt_S1$neglog10P[gt_S1$chr==chr_i], col="black", lwd=2)
}


abline(h=bonferroni_threshold_S1, col="red", lwd=2)
abline(v=genome_info_df$cum_size)



#plot linkage map
par(mar=c(5.1, 4.1,0.5, 2.1))
plot(cross_S1_df$pos[cross_S1_df$chrom=="OW569319.1"], cross_S1_df$FRE_mean_AA[cross_S1_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, max(cross_S1_lm$pos)), ylab="genetic distance (cM)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)


#rect(xleft = most_common_centromeric_repeats$cum_start, ybottom = rep(0, length(most_common_centromeric_repeats$cum_start)), xright = most_common_centromeric_repeats$cum_end, ytop = rep(10, length(most_common_centromeric_repeats$cum_end))   )


rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end, 
     ytop = rep(max(cross_S1_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end)), col=rgb(1,0,0,alpha=0.5), border=NA)

rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end, 
     ytop = rep(max(cross_S1_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end)), col=rgb(0,0,1,alpha=0.5), border=NA)







for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_S1_lm$pos_cum[cross_S1_lm$chr==chr_i], cross_S1_lm$pos[cross_S1_lm$chr==chr_i], col="black")
}


abline(v=genome_info_df$cum_size)

axis(1, at = genome_info_df$cum_size+(genome_info_df$size/2),
     labels = genome_info_df$chr, cex.axis=1.4)
dev.off()
#cross 3A

#cross_obj_3A<-read.cross(format="csvr", file ="3A_rqtl/", alleles=c("A","B"))
load("data/cross_3A.Rdata")
cross_obj_3A<-cross_obj_sub

genotypes_3A<-pull.geno(cross_obj_3A)
phenotypes_3A<-pull.pheno(cross_obj_3A)



cross_3A_df=data.frame(id=colnames(genotypes_3A), 
                       chrom=unlist(strsplit(colnames(genotypes_3A), split="_"))[ c(TRUE,FALSE) ], 
                       pos=as.numeric(unlist(strsplit(colnames(genotypes_3A), split="_"))[ c(FALSE,TRUE) ]),
                       FRE_mean_AA=rep(NA, length(colnames(genotypes_3A))),
                       FRE_mean_AB=rep(NA, length(colnames(genotypes_3A))), 
                       FRE_mean_BB=rep(NA, length(colnames(genotypes_3A))),
                       MRE_mean_AA=rep(NA, length(colnames(genotypes_3A))),
                       MRE_mean_AB=rep(NA, length(colnames(genotypes_3A))), 
                       MRE_mean_BB=rep(NA, length(colnames(genotypes_3A))),
                       rMA_mean_AA=rep(NA, length(colnames(genotypes_3A))),
                       rMA_mean_AB=rep(NA, length(colnames(genotypes_3A))), 
                       rMA_mean_BB=rep(NA, length(colnames(genotypes_3A))))

for(locus in 1:dim(genotypes_3A)[2]){
  cross_3A_df$FRE_mean_AA[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==1,2])
  cross_3A_df$FRE_mean_AB[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==2,2])
  cross_3A_df$FRE_mean_BB[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==3,2])
  cross_3A_df$MRE_mean_AA[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==1,3])
  cross_3A_df$MRE_mean_AB[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==2,3])
  cross_3A_df$MRE_mean_BB[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==3,3])
  cross_3A_df$rMA_mean_AA[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==1,4])
  cross_3A_df$rMA_mean_AB[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==2,4])
  cross_3A_df$rMA_mean_BB[locus]<-mean(phenotypes_3A[genotypes_3A[,locus]==3,4])
}

cross_3A_df$pos_cum<-cross_3A_df$pos
for(loc_i in 1:length(cross_3A_df$pos)){
  cross_3A_df$pos_cum[loc_i]<-genome_info_df[genome_info_df$Lg==cross_3A_df$chrom[loc_i],]$cum_size+cross_3A_df$pos[loc_i]
}
cross_3A_lm<-pull.map(cross_obj_3A, as.table = TRUE)
cross_3A_lm$pos_cum<-cross_3A_df$pos_cum
cross_3A_qtls<-readRDS("data/cross_3A_qtls.Rdata")

gt_3A<-geno.table(cross_obj_3A, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_3A<--log10(0.05/totmar(cross_obj_3A))
pdf("plots/3A_plot.pdf", width = 8, height=12)
par(mfrow=c(5,1))
#bottom, left, top, right
#default_mar<-c(5.1, 4.1, 4.1, 2.1)
par(mar=c(0, 4.1,0.5, 2.1))

#FRE
plot(cross_3A_df$pos[cross_3A_df$chrom=="OW569319.1"], cross_3A_df$FRE_mean_AA[cross_3A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_3A_df$FRE_mean_AA,cross_3A_df$FRE_mean_AB,cross_3A_df$FRE_mean_BB)), max(c(cross_3A_df$FRE_mean_AA,cross_3A_df$FRE_mean_AB,cross_3A_df$FRE_mean_BB))), ylab="mean FRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?


#QTL1
cross_3A_FRE_qtl1_left_edge<-cross_3A_df[as.character(cross_3A_df$id)==as.character(cross_3A_qtls$begin[1]),]
cross_3A_FRE_qtl1_center<-cross_3A_df[as.character(cross_3A_df$id)==as.character(cross_3A_qtls$center[1]),]
cross_3A_FRE_qtl1_right_edge<-cross_3A_df[as.character(cross_3A_df$id)==as.character(cross_3A_qtls$end[1]),]
cross_3A_FRE_qtl1_lines_df<-cross_3A_df[as.character(cross_3A_df$chrom)==cross_3A_FRE_qtl1_left_edge$chrom[1] & cross_3A_df$pos_cum>=cross_3A_FRE_qtl1_left_edge$pos_cum[1] & cross_3A_df$pos_cum<=cross_3A_FRE_qtl1_right_edge$pos_cum[1],]
cross_3A_FRE_qtl1_lines_AA<-lines(cross_3A_FRE_qtl1_lines_df$pos_cum, cross_3A_FRE_qtl1_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_3A_FRE_qtl1_lines_BB<-lines(cross_3A_FRE_qtl1_lines_df$pos_cum, cross_3A_FRE_qtl1_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_3A_FRE_qtl1_lines_df$pos_cum, rev(cross_3A_FRE_qtl1_lines_df$pos_cum)), y=c(cross_3A_FRE_qtl1_lines_df$FRE_mean_BB, rev(cross_3A_FRE_qtl1_lines_df$FRE_mean_AA)),  col=rgb(0,0,1.0,alpha=0.5), lty=0)

abline(v=cross_3A_FRE_qtl1_center$pos_cum, col="blue", lty=2, lwd=2)


#QTL2
cross_3A_FRE_qtl2_left_edge<-cross_3A_df[as.character(cross_3A_df$id)==as.character(cross_3A_qtls$begin[2]),]
cross_3A_FRE_qtl2_center<-cross_3A_df[as.character(cross_3A_df$id)==as.character(cross_3A_qtls$center[2]),]
cross_3A_FRE_qtl2_right_edge<-cross_3A_df[as.character(cross_3A_df$id)==as.character(cross_3A_qtls$end[2]),]
cross_3A_FRE_qtl2_lines_df<-cross_3A_df[as.character(cross_3A_df$chrom)==cross_3A_FRE_qtl2_left_edge$chrom[1] & cross_3A_df$pos_cum>=cross_3A_FRE_qtl2_left_edge$pos_cum[1] & cross_3A_df$pos_cum<=cross_3A_FRE_qtl2_right_edge$pos_cum[1],]
cross_3A_FRE_qtl2_lines_AA<-lines(cross_3A_FRE_qtl2_lines_df$pos_cum, cross_3A_FRE_qtl2_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_3A_FRE_qtl2_lines_BB<-lines(cross_3A_FRE_qtl2_lines_df$pos_cum, cross_3A_FRE_qtl2_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_3A_FRE_qtl2_lines_df$pos_cum, rev(cross_3A_FRE_qtl2_lines_df$pos_cum)), y=c(cross_3A_FRE_qtl2_lines_df$FRE_mean_AB, rev(cross_3A_FRE_qtl2_lines_df$FRE_mean_AA)),  
        col=rgb(0,0,1.0,alpha=0.5),          
        lty=1,
        density=30,
        angle=45)

abline(v=cross_3A_FRE_qtl2_center$pos_cum, col="blue", lty=2, lwd=2)


abline(h=mean(phenotypes_3A$FRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$FRE_mean_AA[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$FRE_mean_AB[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$FRE_mean_BB[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}
#MRE
plot(cross_3A_df$pos[cross_3A_df$chrom=="OW569319.1"], cross_3A_df$MRE_mean_AA[cross_3A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
    ylim=c(min(c(cross_3A_df$MRE_mean_AA,cross_3A_df$MRE_mean_AB,cross_3A_df$MRE_mean_BB)), max(c(cross_3A_df$MRE_mean_AA,cross_3A_df$MRE_mean_AB,cross_3A_df$MRE_mean_BB))), ylab="mean MRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)

abline(h=mean(phenotypes_3A$MRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$MRE_mean_AA[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$MRE_mean_AB[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$MRE_mean_BB[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}
#rMA
plot(cross_3A_df$pos[cross_3A_df$chrom=="OW569319.1"], cross_3A_df$rMA_mean_AA[cross_3A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_3A_df$rMA_mean_AA,cross_3A_df$rMA_mean_AB,cross_3A_df$rMA_mean_BB)), max(c(cross_3A_df$rMA_mean_AA,cross_3A_df$rMA_mean_AB,cross_3A_df$rMA_mean_BB))), ylab="mean rMA", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)

abline(h=mean(phenotypes_3A$rMA), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$rMA_mean_AA[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$rMA_mean_AB[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_3A_df$pos_cum[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], cross_3A_df$rMA_mean_BB[cross_3A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}

#plot segregation distortion
plot(cross_3A_df$pos[cross_3A_df$chrom=="OW569319.1"], cross_3A_df$FRE_mean_AA[cross_3A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, max(gt_3A$neglog10P)), ylab="segregation distortion (-10*log p)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)

for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_3A_lm$pos_cum[cross_3A_lm$chr==chr_i], gt_3A$neglog10P[gt_3A$chr==chr_i], col="black", lwd=2)
}


abline(h=bonferroni_threshold_3A, col="red", lwd=2)

abline(v=genome_info_df$cum_size)


#plot linkage map
par(mar=c(5.1, 4.1,0.5, 2.1))
plot(cross_3A_df$pos[cross_3A_df$chrom=="OW569319.1"], cross_3A_df$FRE_mean_AA[cross_3A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(0, max(cross_3A_lm$pos)), ylab="genetic distance (cM)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)


#rect(xleft = most_common_centromeric_repeats$cum_start, ybottom = rep(0, length(most_common_centromeric_repeats$cum_start)), xright = most_common_centromeric_repeats$cum_end, ytop = rep(10, length(most_common_centromeric_repeats$cum_end))   )
rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end, 
     ytop = rep(max(cross_3A_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end)), col=rgb(1,0,0,alpha=0.5), border=NA)

rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end, 
     ytop = rep(max(cross_3A_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end)), col=rgb(0,0,1,alpha=0.5), border=NA)


for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_3A_lm$pos_cum[cross_3A_lm$chr==chr_i], cross_3A_lm$pos[cross_3A_lm$chr==chr_i], col="black")
}


abline(v=genome_info_df$cum_size)

axis(1, at = genome_info_df$cum_size+(genome_info_df$size/2),
     labels = genome_info_df$chr, cex.axis=1.4)
dev.off()
#Cross 5A

cross_obj_5A<-read.cross(format="csvr", file ="data/5A.csvr", alleles=c("A","B"))

genotypes_5A<-pull.geno(cross_obj_5A)
phenotypes_5A<-pull.pheno(cross_obj_5A)



cross_5A_df=data.frame(id=colnames(genotypes_5A), 
                       chrom=unlist(strsplit(colnames(genotypes_5A), split="_"))[ c(TRUE,FALSE) ], 
                       pos=as.numeric(unlist(strsplit(colnames(genotypes_5A), split="_"))[ c(FALSE,TRUE) ]),
                       FRE_mean_AA=rep(NA, length(colnames(genotypes_5A))),
                       FRE_mean_AB=rep(NA, length(colnames(genotypes_5A))), 
                       FRE_mean_BB=rep(NA, length(colnames(genotypes_5A))),
                       MRE_mean_AA=rep(NA, length(colnames(genotypes_5A))),
                       MRE_mean_AB=rep(NA, length(colnames(genotypes_5A))), 
                       MRE_mean_BB=rep(NA, length(colnames(genotypes_5A))),
                       rMA_mean_AA=rep(NA, length(colnames(genotypes_5A))),
                       rMA_mean_AB=rep(NA, length(colnames(genotypes_5A))), 
                       rMA_mean_BB=rep(NA, length(colnames(genotypes_5A))))

for(locus in 1:dim(genotypes_5A)[2]){
  cross_5A_df$FRE_mean_AA[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==1,2])
  cross_5A_df$FRE_mean_AB[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==2,2])
  cross_5A_df$FRE_mean_BB[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==3,2])
  cross_5A_df$MRE_mean_AA[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==1,3])
  cross_5A_df$MRE_mean_AB[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==2,3])
  cross_5A_df$MRE_mean_BB[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==3,3])
  cross_5A_df$rMA_mean_AA[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==1,4])
  cross_5A_df$rMA_mean_AB[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==2,4])
  cross_5A_df$rMA_mean_BB[locus]<-mean(phenotypes_5A[genotypes_5A[,locus]==3,4])
}

cross_5A_df$pos_cum<-cross_5A_df$pos
for(loc_i in 1:length(cross_5A_df$pos)){
  cross_5A_df$pos_cum[loc_i]<-genome_info_df[genome_info_df$Lg==cross_5A_df$chrom[loc_i],]$cum_size+cross_5A_df$pos[loc_i]
}
cross_5A_lm<-pull.map(cross_obj_5A, as.table = TRUE)
cross_5A_lm$pos_cum<-cross_5A_df$pos_cum
cross_5A_qtls<-readRDS("data/cross_5A_qtls.Rdata")
gt_5A<-geno.table(cross_obj_5A, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_5A<--log10(0.05/totmar(cross_obj_5A))
pdf("plots/5A_plot.pdf", width = 8, height=12)

par(mfrow=c(5,1))
#bottom, left, top, right
#default_mar<-c(5.1, 4.1, 4.1, 2.1)
par(mar=c(0, 4.1,0.5, 2.1))

#FRE
plot(cross_5A_df$pos[cross_5A_df$chrom=="OW569319.1"], cross_5A_df$FRE_mean_AA[cross_5A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(min(c(cross_5A_df$FRE_mean_AA,cross_5A_df$FRE_mean_AB,cross_5A_df$FRE_mean_BB)), max(c(cross_5A_df$FRE_mean_AA,cross_5A_df$FRE_mean_AB,cross_5A_df$FRE_mean_BB))), ylab="mean FRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?



#QTL1
cross_5A_FRE_qtl1_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[1]),]
cross_5A_FRE_qtl1_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[1]),]
cross_5A_FRE_qtl1_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[1]),]
cross_5A_FRE_qtl1_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_FRE_qtl1_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_FRE_qtl1_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_FRE_qtl1_right_edge$pos_cum[1],]
cross_5A_FRE_qtl1_lines_AA<-lines(cross_5A_FRE_qtl1_lines_df$pos_cum, cross_5A_FRE_qtl1_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_5A_FRE_qtl1_lines_BB<-lines(cross_5A_FRE_qtl1_lines_df$pos_cum, cross_5A_FRE_qtl1_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_FRE_qtl1_lines_df$pos_cum, rev(cross_5A_FRE_qtl1_lines_df$pos_cum)), y=c(cross_5A_FRE_qtl1_lines_df$FRE_mean_BB, rev(cross_5A_FRE_qtl1_lines_df$FRE_mean_AA)),  col=rgb(0,0,1.0,alpha=0.5), lty=0)

abline(v=cross_5A_FRE_qtl1_center$pos_cum, col="blue", lty=2, lwd=2)


#QTL2 -lg7
cross_5A_FRE_qtl2_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[2]),]
cross_5A_FRE_qtl2_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[2]),]
cross_5A_FRE_qtl2_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[2]),]
cross_5A_FRE_qtl2_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_FRE_qtl2_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_FRE_qtl2_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_FRE_qtl2_right_edge$pos_cum[1],]
cross_5A_FRE_qtl2_lines_AA<-lines(cross_5A_FRE_qtl2_lines_df$pos_cum, cross_5A_FRE_qtl2_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_5A_FRE_qtl2_lines_BB<-lines(cross_5A_FRE_qtl2_lines_df$pos_cum, cross_5A_FRE_qtl2_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_FRE_qtl2_lines_df$pos_cum, rev(cross_5A_FRE_qtl2_lines_df$pos_cum)), y=c(cross_5A_FRE_qtl2_lines_df$FRE_mean_BB, rev(cross_5A_FRE_qtl2_lines_df$FRE_mean_AA)),  col=rgb(0,0,1.0,alpha=0.5), lty=0)

abline(v=cross_5A_FRE_qtl2_center$pos_cum, col="blue", lty=2, lwd=2)




#QTL3 - lg8 two qtls
cross_5A_FRE_qtl3_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[3]),]
#Fix this!!!!
cross_5A_FRE_qtl3_center<-cross_5A_df[as.character(cross_5A_df$id)=="OW569318.1_282665",]
cross_5A_FRE_qtl3_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[3]),]
cross_5A_FRE_qtl3_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_FRE_qtl3_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_FRE_qtl3_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_FRE_qtl3_right_edge$pos_cum[1],]
cross_5A_FRE_qtl3_lines_AA<-lines(cross_5A_FRE_qtl3_lines_df$pos_cum, cross_5A_FRE_qtl3_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_5A_FRE_qtl3_lines_BB<-lines(cross_5A_FRE_qtl3_lines_df$pos_cum, cross_5A_FRE_qtl3_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_FRE_qtl3_lines_df$pos_cum, rev(cross_5A_FRE_qtl3_lines_df$pos_cum)), y=c(cross_5A_FRE_qtl3_lines_df$FRE_mean_AB, rev(cross_5A_FRE_qtl3_lines_df$FRE_mean_AA)),  
        col=rgb(0,0,1.0,alpha=0.5), 
        lty=1,
        density=30,
        angle=45)

abline(v=cross_5A_FRE_qtl3_center$pos_cum, col="blue", lty=2, lwd=2)


#QTL4 - lg 8 one qtl
cross_5A_FRE_qtl4_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[4]),]
cross_5A_FRE_qtl4_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[4]),]
cross_5A_FRE_qtl4_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[4]),]
cross_5A_FRE_qtl4_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_FRE_qtl4_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_FRE_qtl4_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_FRE_qtl4_right_edge$pos_cum[1],]
cross_5A_FRE_qtl4_lines_AA<-lines(cross_5A_FRE_qtl4_lines_df$pos_cum, cross_5A_FRE_qtl4_lines_df$FRE_mean_AA, col="black", lwd=0)
cross_5A_FRE_qtl4_lines_BB<-lines(cross_5A_FRE_qtl4_lines_df$pos_cum, cross_5A_FRE_qtl4_lines_df$FRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_FRE_qtl4_lines_df$pos_cum, rev(cross_5A_FRE_qtl4_lines_df$pos_cum)), y=c(cross_5A_FRE_qtl4_lines_df$FRE_mean_AB, rev(cross_5A_FRE_qtl4_lines_df$FRE_mean_AA)),  col=rgb(0,0,1.0,alpha=0.5),
        lty=1,
        density=30,
        angle=45)

abline(v=cross_5A_FRE_qtl4_center$pos_cum, col="blue", lty=2, lwd=2)

abline(h=mean(phenotypes_5A$FRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$FRE_mean_AA[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$FRE_mean_AB[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$FRE_mean_BB[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}



#MRE

plot(cross_5A_df$pos[cross_5A_df$chrom=="OW569319.1"], cross_5A_df$MRE_mean_AA[cross_5A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(0, max(c(cross_5A_df$MRE_mean_AA,cross_5A_df$MRE_mean_AB,cross_5A_df$MRE_mean_BB))), ylab="mean MRE", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?



#QTL1
cross_5A_MRE_qtl1_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[5]),]
cross_5A_MRE_qtl1_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[5]),]
cross_5A_MRE_qtl1_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[5]),]
cross_5A_MRE_qtl1_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_MRE_qtl1_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_MRE_qtl1_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_MRE_qtl1_right_edge$pos_cum[1],]
cross_5A_MRE_qtl1_lines_AA<-lines(cross_5A_MRE_qtl1_lines_df$pos_cum, cross_5A_MRE_qtl1_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_5A_MRE_qtl1_lines_BB<-lines(cross_5A_MRE_qtl1_lines_df$pos_cum, cross_5A_MRE_qtl1_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_MRE_qtl1_lines_df$pos_cum, rev(cross_5A_MRE_qtl1_lines_df$pos_cum)), y=c(cross_5A_MRE_qtl1_lines_df$MRE_mean_BB, rev(cross_5A_MRE_qtl1_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)

abline(v=cross_5A_MRE_qtl1_center$pos_cum, col="red", lty=2, lwd=2)


#qtl2
cross_5A_MRE_qtl2_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[6]),]
cross_5A_MRE_qtl2_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[6]),]
cross_5A_MRE_qtl2_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[6]),]
cross_5A_MRE_qtl2_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_MRE_qtl2_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_MRE_qtl2_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_MRE_qtl2_right_edge$pos_cum[1],]
cross_5A_MRE_qtl2_lines_AA<-lines(cross_5A_MRE_qtl2_lines_df$pos_cum, cross_5A_MRE_qtl2_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_5A_MRE_qtl2_lines_BB<-lines(cross_5A_MRE_qtl2_lines_df$pos_cum, cross_5A_MRE_qtl2_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_MRE_qtl2_lines_df$pos_cum, rev(cross_5A_MRE_qtl2_lines_df$pos_cum)), y=c(cross_5A_MRE_qtl2_lines_df$MRE_mean_BB, rev(cross_5A_MRE_qtl2_lines_df$MRE_mean_AB)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)

abline(v=cross_5A_MRE_qtl2_center$pos_cum, col="red", lty=2, lwd=2)


#qtl3
cross_5A_MRE_qtl3_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[7]),]
cross_5A_MRE_qtl3_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[7]),]
cross_5A_MRE_qtl3_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[7]),]
cross_5A_MRE_qtl3_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_MRE_qtl3_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_MRE_qtl3_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_MRE_qtl3_right_edge$pos_cum[1],]
cross_5A_MRE_qtl3_lines_AA<-lines(cross_5A_MRE_qtl3_lines_df$pos_cum, cross_5A_MRE_qtl3_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_5A_MRE_qtl3_lines_BB<-lines(cross_5A_MRE_qtl3_lines_df$pos_cum, cross_5A_MRE_qtl3_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_MRE_qtl3_lines_df$pos_cum, rev(cross_5A_MRE_qtl3_lines_df$pos_cum)), y=c(cross_5A_MRE_qtl3_lines_df$MRE_mean_BB, rev(cross_5A_MRE_qtl3_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)

abline(v=cross_5A_MRE_qtl3_center$pos_cum, col="red", lty=2, lwd=2)


#qtl4
cross_5A_MRE_qtl4_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[8]),]
cross_5A_MRE_qtl4_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[8]),]
cross_5A_MRE_qtl4_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[8]),]
cross_5A_MRE_qtl4_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_MRE_qtl4_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_MRE_qtl4_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_MRE_qtl4_right_edge$pos_cum[1],]
cross_5A_MRE_qtl4_lines_AA<-lines(cross_5A_MRE_qtl4_lines_df$pos_cum, cross_5A_MRE_qtl4_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_5A_MRE_qtl4_lines_BB<-lines(cross_5A_MRE_qtl4_lines_df$pos_cum, cross_5A_MRE_qtl4_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_MRE_qtl4_lines_df$pos_cum, rev(cross_5A_MRE_qtl4_lines_df$pos_cum)), y=c(cross_5A_MRE_qtl4_lines_df$MRE_mean_BB, rev(cross_5A_MRE_qtl4_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)

abline(v=cross_5A_MRE_qtl4_center$pos_cum, col="red", lty=2, lwd=2)


#qtl5
cross_5A_MRE_qtl5_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[9]),]
cross_5A_MRE_qtl5_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[9]),]
cross_5A_MRE_qtl5_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[9]),]
cross_5A_MRE_qtl5_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_MRE_qtl5_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_MRE_qtl5_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_MRE_qtl5_right_edge$pos_cum[1],]
cross_5A_MRE_qtl5_lines_AA<-lines(cross_5A_MRE_qtl5_lines_df$pos_cum, cross_5A_MRE_qtl5_lines_df$MRE_mean_AA, col="black", lwd=0)
cross_5A_MRE_qtl5_lines_BB<-lines(cross_5A_MRE_qtl5_lines_df$pos_cum, cross_5A_MRE_qtl5_lines_df$MRE_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_MRE_qtl5_lines_df$pos_cum, rev(cross_5A_MRE_qtl5_lines_df$pos_cum)), y=c(cross_5A_MRE_qtl5_lines_df$MRE_mean_BB, rev(cross_5A_MRE_qtl5_lines_df$MRE_mean_AA)),  col=rgb(1,0,0.0,alpha=0.5), lty=0)

abline(v=cross_5A_MRE_qtl5_center$pos_cum, col="red", lty=2, lwd=2)



abline(h=mean(phenotypes_5A$MRE), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$MRE_mean_AA[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$MRE_mean_AB[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$MRE_mean_BB[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}

#rMA


plot(cross_5A_df$pos[cross_5A_df$chrom=="OW569319.1"], cross_5A_df$rMA_mean_AA[cross_5A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), 
     ylim=c(0, max(c(cross_5A_df$rMA_mean_AA,cross_5A_df$rMA_mean_AB,cross_5A_df$rMA_mean_BB))), ylab="mean rMA", xlab ="", xaxt='n', cex.lab=1.4, cex.axis=1.3)
#plot QTLS first?




#qtl6
cross_5A_rMA_qtl6_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[10]),]
cross_5A_rMA_qtl6_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[10]),]
cross_5A_rMA_qtl6_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[10]),]
cross_5A_rMA_qtl6_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_rMA_qtl6_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_rMA_qtl6_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_rMA_qtl6_right_edge$pos_cum[1],]
cross_5A_rMA_qtl6_lines_AA<-lines(cross_5A_rMA_qtl6_lines_df$pos_cum, cross_5A_rMA_qtl6_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_5A_rMA_qtl6_lines_BB<-lines(cross_5A_rMA_qtl6_lines_df$pos_cum, cross_5A_rMA_qtl6_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_rMA_qtl6_lines_df$pos_cum, rev(cross_5A_rMA_qtl6_lines_df$pos_cum)), y=c(cross_5A_rMA_qtl6_lines_df$rMA_mean_AA, rev(cross_5A_rMA_qtl6_lines_df$rMA_mean_BB)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)

abline(v=cross_5A_rMA_qtl6_center$pos_cum, col="darkgreen", lty=2, lwd=2)





#QTL1
cross_5A_rMA_qtl1_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[11]),]
cross_5A_rMA_qtl1_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[11]),]
cross_5A_rMA_qtl1_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[11]),]
cross_5A_rMA_qtl1_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_rMA_qtl1_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_rMA_qtl1_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_rMA_qtl1_right_edge$pos_cum[1],]
cross_5A_rMA_qtl1_lines_AA<-lines(cross_5A_rMA_qtl1_lines_df$pos_cum, cross_5A_rMA_qtl1_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_5A_rMA_qtl1_lines_BB<-lines(cross_5A_rMA_qtl1_lines_df$pos_cum, cross_5A_rMA_qtl1_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_rMA_qtl1_lines_df$pos_cum, rev(cross_5A_rMA_qtl1_lines_df$pos_cum)), y=c(cross_5A_rMA_qtl1_lines_df$rMA_mean_BB, rev(cross_5A_rMA_qtl1_lines_df$rMA_mean_AA)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)

abline(v=cross_5A_rMA_qtl1_center$pos_cum, col="darkgreen", lty=2, lwd=2)


#QTL2
cross_5A_rMA_qtl2_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[12]),]
cross_5A_rMA_qtl2_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[12]),]
cross_5A_rMA_qtl2_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[12]),]
cross_5A_rMA_qtl2_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_rMA_qtl2_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_rMA_qtl2_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_rMA_qtl2_right_edge$pos_cum[1],]
cross_5A_rMA_qtl2_lines_AA<-lines(cross_5A_rMA_qtl2_lines_df$pos_cum, cross_5A_rMA_qtl2_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_5A_rMA_qtl2_lines_BB<-lines(cross_5A_rMA_qtl2_lines_df$pos_cum, cross_5A_rMA_qtl2_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_rMA_qtl2_lines_df$pos_cum, rev(cross_5A_rMA_qtl2_lines_df$pos_cum)), y=c(cross_5A_rMA_qtl2_lines_df$rMA_mean_BB, rev(cross_5A_rMA_qtl2_lines_df$rMA_mean_AA)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)

abline(v=cross_5A_rMA_qtl2_center$pos_cum, col="darkgreen", lty=2, lwd=2)

#QTL3
cross_5A_rMA_qtl3_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[13]),]
cross_5A_rMA_qtl3_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[13]),]
cross_5A_rMA_qtl3_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[13]),]
cross_5A_rMA_qtl3_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_rMA_qtl3_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_rMA_qtl3_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_rMA_qtl3_right_edge$pos_cum[1],]
cross_5A_rMA_qtl3_lines_AA<-lines(cross_5A_rMA_qtl3_lines_df$pos_cum, cross_5A_rMA_qtl3_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_5A_rMA_qtl3_lines_BB<-lines(cross_5A_rMA_qtl3_lines_df$pos_cum, cross_5A_rMA_qtl3_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_rMA_qtl3_lines_df$pos_cum, rev(cross_5A_rMA_qtl3_lines_df$pos_cum)), y=c(cross_5A_rMA_qtl3_lines_df$rMA_mean_BB, rev(cross_5A_rMA_qtl3_lines_df$rMA_mean_AA)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)

abline(v=cross_5A_rMA_qtl3_center$pos_cum, col="darkgreen", lty=2, lwd=2)


#QTL4
cross_5A_rMA_qtl4_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[14]),]
cross_5A_rMA_qtl4_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[14]),]
cross_5A_rMA_qtl4_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[14]),]
cross_5A_rMA_qtl4_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_rMA_qtl4_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_rMA_qtl4_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_rMA_qtl4_right_edge$pos_cum[1],]
cross_5A_rMA_qtl4_lines_AA<-lines(cross_5A_rMA_qtl4_lines_df$pos_cum, cross_5A_rMA_qtl4_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_5A_rMA_qtl4_lines_BB<-lines(cross_5A_rMA_qtl4_lines_df$pos_cum, cross_5A_rMA_qtl4_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_rMA_qtl4_lines_df$pos_cum, rev(cross_5A_rMA_qtl4_lines_df$pos_cum)), y=c(cross_5A_rMA_qtl4_lines_df$rMA_mean_BB, rev(cross_5A_rMA_qtl4_lines_df$rMA_mean_AA)),  col=rgb(0,1,0.5,alpha=0.5), lty=0)

abline(v=cross_5A_rMA_qtl4_center$pos_cum, col="darkgreen", lty=2, lwd=2)


#QTL5
cross_5A_rMA_qtl5_left_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$begin[15]),]
cross_5A_rMA_qtl5_center<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$center[15]),]
cross_5A_rMA_qtl5_right_edge<-cross_5A_df[as.character(cross_5A_df$id)==as.character(cross_5A_qtls$end[15]),]
cross_5A_rMA_qtl5_lines_df<-cross_5A_df[as.character(cross_5A_df$chrom)==cross_5A_rMA_qtl5_left_edge$chrom[1] & cross_5A_df$pos_cum>=cross_5A_rMA_qtl5_left_edge$pos_cum[1] & cross_5A_df$pos_cum<=cross_5A_rMA_qtl5_right_edge$pos_cum[1],]
cross_5A_rMA_qtl5_lines_AA<-lines(cross_5A_rMA_qtl5_lines_df$pos_cum, cross_5A_rMA_qtl5_lines_df$rMA_mean_AA, col="black", lwd=0)
cross_5A_rMA_qtl5_lines_BB<-lines(cross_5A_rMA_qtl5_lines_df$pos_cum, cross_5A_rMA_qtl5_lines_df$rMA_mean_BB, col="black", lwd=0)
polygon(x=c(cross_5A_rMA_qtl5_lines_df$pos_cum, rev(cross_5A_rMA_qtl5_lines_df$pos_cum)), y=c(cross_5A_rMA_qtl5_lines_df$rMA_mean_BB, rev(cross_5A_rMA_qtl5_lines_df$rMA_mean_AA)),
        col=rgb(0,1,0.5,alpha=0.5), 
        lty=1,
        density=30,
        angle=45)

abline(v=cross_5A_rMA_qtl5_center$pos_cum, col="darkgreen", lty=2, lwd=2)


abline(h=mean(phenotypes_5A$rMA), lty=2)
abline(v=genome_info_df$cum_size)
for(chr_i in 1:8){
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$rMA_mean_AA[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], col="red", lwd=2)
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$rMA_mean_AB[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], lwd=2)
  lines(cross_5A_df$pos_cum[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], cross_5A_df$rMA_mean_BB[cross_5A_df$chrom==genome_info_df$Lg[chr_i]], col="blue", lwd=2) 
}


#plot segregation distortion
plot(cross_5A_df$pos[cross_5A_df$chrom=="OW569319.1"], cross_5A_df$FRE_mean_AA[cross_5A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, max(gt_5A$neglog10P)), ylab="segregation distortion (-10*log p)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)

for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_5A_lm$pos_cum[cross_5A_lm$chr==chr_i], gt_5A$neglog10P[gt_5A$chr==chr_i], col="black", lwd=2)
}


abline(h=bonferroni_threshold_5A, col="red", lwd=2)
abline(v=genome_info_df$cum_size)

#plot linkage map
par(mar=c(5.1, 4.1,0.5, 2.1))
plot(cross_5A_df$pos[cross_5A_df$chrom=="OW569319.1"], cross_5A_df$FRE_mean_AA[cross_5A_df$chrom=="OW569319.1"], cex=0, xlim=c(0, sum(genome_info_df$size)), ylim=c(0, max(cross_5A_lm$pos)), ylab="genetic distance (cM)", xaxt="n", xlab="", cex.lab=1.4, cex.axis=1.3)


#rect(xleft = most_common_centromeric_repeats$cum_start, ybottom = rep(0, length(most_common_centromeric_repeats$cum_start)), xright = most_common_centromeric_repeats$cum_end, ytop = rep(10, length(most_common_centromeric_repeats$cum_end))   )
rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end, 
     ytop = rep(max(cross_5A_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==60,]$cum_end)), col=rgb(1,0,0,alpha=0.5), border=NA)

rect(xleft = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start, 
     ybottom = rep(0, length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_start)), 
     xright = most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end, 
     ytop = rep(max(cross_5A_lm$pos), length(most_common_centromeric_repeats[most_common_centromeric_repeats$most.freq.value.N==99,]$cum_end)), col=rgb(0,0,1,alpha=0.5), border=NA)





for(chr_i in c("Lg1","Lg2","Lg3","Lg4","Lg5","Lg6","Lg7","Lg8")){
  lines(cross_5A_lm$pos_cum[cross_5A_lm$chr==chr_i], cross_5A_lm$pos[cross_5A_lm$chr==chr_i], col="black")
}


abline(v=genome_info_df$cum_size)

axis(1, at = genome_info_df$cum_size+(genome_info_df$size/2),
     labels = genome_info_df$chr, cex.axis=1.4)
dev.off()
