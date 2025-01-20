options(scipen=999)
setwd("/home/fredo/Desktop/qtl_paper/upload_suppl_information/rqtl/")
library(qtl)
library(stringr)
library(shape)
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

####Cross 1A
cross_obj_1A<-read.cross(format="csvr", file ="data/1A.csvr", alleles=c("A","B"))
#cross_obj_1A<-jittermap(cross_obj_1A)
cross_obj_1A<-calc.genoprob(cross_obj_1A, step=1)


out_1A.hk_FRE <- scanone(cross_obj_1A, pheno.col = 2, method="hk")
operm_1A.hk_FRE <- scanone(cross_obj_1A, method="hk", n.perm=1000, pheno.col = 2)
out_1A.hk_MRE <- scanone(cross_obj_1A, pheno.col = 3, method="hk")
operm_1A.hk_MRE <- scanone(cross_obj_1A, method="hk", n.perm=1000, pheno.col = 3)
out_1A.hk_rMA<- scanone(cross_obj_1A, pheno.col = 4, method="hk")
operm_1A.hk_rMA <- scanone(cross_obj_1A, method="hk", n.perm=1000, pheno.col = 4)

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

cross_1A_lm<-pull.map(cross_obj_1A, as.table = TRUE)


cross_1A_lm$pos_cum<-cross_1A_df$pos_cum

for(loc_i in 1:length(cross_1A_df$pos)){
  cross_1A_df$pos_cum[loc_i]<-genome_info_df[genome_info_df$Lg==cross_1A_df$chrom[loc_i],]$cum_size+cross_1A_df$pos[loc_i]
}

#can we directly add the linkage map?
cross_1A_qtls<-readRDS("data/cross_1A_qtls.Rdata")
cross_1A_qtls$chr_number<-c(3,5,4,5,7,5,6,7)
cross_1A_qtls$direction<-c(2,1,1,2,1,2,2,1)
cross_1A_qtls$cross<-rep(1, nrow(cross_1A_qtls))
cross_1A_qtls$prop_variance<-c(7.9, 11.67, 10.53, 14.18,10.67, 24.09, 3.854, 10.04)


gt_1A<-geno.table(cross_obj_1A, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_1A<--log10(0.05/totmar(cross_obj_1A))



##########Cross S1

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
cross_S1_qtls$chr_number<-c(7,4,7,4,7)
cross_S1_qtls$direction<-c(1,2,2,2,2)
cross_S1_qtls$cross<-rep(2, nrow(cross_S1_qtls))
cross_S1_qtls$prop_variance<-c(20.78, 4.72, 30.23, 4.71, 31.99)

gt_S1<-geno.table(cross_obj_S1, scanone.output = TRUE)
#do bonferroni correction
bonferroni_threshold_S1<--log10(0.05/totmar(cross_obj_S1))

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
cross_3A_qtls$chr_number<-c(4,8)
cross_3A_qtls$direction<-c(2,1)
cross_3A_qtls$cross<-rep(3, nrow(cross_3A_qtls))
cross_3A_qtls$prop_variance<-c(5.23, 4.47)

gt_3A<-geno.table(cross_obj_3A, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_3A<--log10(0.05/totmar(cross_obj_3A))
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
cross_5A_qtls$chr_number<-c(4,7,8,8,1,3,4,5,7,1,3,4,5,7,8)
cross_5A_qtls$direction<-c(1,2,1,1,1,2,2,2,2,1,2,2,2,2,2)
cross_5A_qtls$cross<-rep(4, nrow(cross_5A_qtls))
cross_5A_qtls$prop_variance<-c(6.64, 8.34, 6.15, 7.27, 5.68, 3.61, 5.03, 12.01, 9.32, 8.06, 4.69, 8.03, 18.56, 4.83, 3.19)

gt_5A<-geno.table(cross_obj_5A, scanone.output = TRUE)

#do bonferroni correction
bonferroni_threshold_5A<--log10(0.05/totmar(cross_obj_5A))


qtl_table<-rbind(cross_1A_qtls, cross_S1_qtls, cross_3A_qtls, cross_5A_qtls)
qtl_table$begin_pos=as.numeric(unlist(strsplit(as.character(qtl_table$begin), split="_"))[ c(FALSE,TRUE) ])
qtl_table$center_pos=as.numeric(unlist(strsplit(as.character(qtl_table$center), split="_"))[ c(FALSE,TRUE) ])
qtl_table$end_pos=as.numeric(unlist(strsplit(as.character(qtl_table$end), split="_"))[ c(FALSE,TRUE) ])
  
  
  ##################Horizontal
    
    pdf("plots/main_figure.pdf", width=11,height=10)
    
    
    my_lay <- layout(mat = matrix(c(1,2,
                                    1,3), nrow = 2, byrow = TRUE))
    #layout.show(my_lay)                 
    
    #Top plot with qtl direction an dprop variation explained
    #botton, left, top, right
    #par(mar=c(1,1,1,1))
    par(mar=c(4,5,6.5,1))
    total_length=76
    #main plot
    tick_positions<-1:8
    par(xpd = T) 
    #plot(0, ylab="", xlab="", xaxt="n", yaxt="n")
    #text(par('usr')[1]-(0.1*par('usr')[1]),par('usr')[4]*1.06,labels=expression(bold("A")) , cex=1.5)
    plot(0, xlim=c(0,90000000), ylim= c(0,total_length), cex=0, xaxt="n", yaxt = "n", ylab="", xlab="Position (Mb)", frame.plot = FALSE, cex.lab=1.3)
    text(par('usr')[1]-(0.1*par('usr')[1]),par('usr')[4]*1.06,labels=expression(bold("A")) , cex=1.5)
    legend("topright", legend = c("Female reproductive effort","Male reproductive effort","Relative male allocation", "Physical chromosomal length","Regions of segregation distortion"), lty = c(1,1,1, NA, NA), col=c("blue","red","darkgreen",NA,NA), inset=c(0.06,-0.115), lwd=2, cex=1.2)
    rect(23500000,80.5 , 31000000,81, col="darkgray")
    rect(23500000,76.9 , 31000000,79.2,col="red", density=30, lty=1, angle=45)
    #SD chross 1
    sd_1A_Lg4_beg<-25262091
    sd_1A_Lg4_end<-29416081
    rect(sd_1A_Lg4_beg,(total_length-19),sd_1A_Lg4_end,(total_length-21), density=30, lty=1, angle=45, col="red")
    
    
    sd_1A_Lg7_beg<-11687577
    sd_1A_Lg7_end<-48863894
    rect(sd_1A_Lg7_beg,(total_length-51),sd_1A_Lg7_end,(total_length-55), col="red", density=30, lty=1, angle=45)
    #SD cross 3
    sd_3A_Lg8_beg<-153808
    sd_3A_Lg8_end<-40141015
    rect(sd_3A_Lg8_beg, (total_length-71), sd_3A_Lg8_end, (total_length-77), col="red", density=30, lty=1, angle=45)
    #SD cross 4
    sd_5A_Lg8_beg<-282665
    sd_5A_Lg8_end<-39620005
    rect(sd_5A_Lg8_beg,(total_length-69), sd_5A_Lg8_end, (total_length-71), col="red", density=30, lty=1, angle=45)
    
    x_position=total_length
    arrow_step=5
    for(chr_i in 1:8){
      #draw chrom
      rect(0, x_position, genome_info_df$size[chr_i], x_position-0.5, col="darkgray")
      
      
      segments(0,
               x_position-0.75,
               95000000,
               x_position-0.75, lty=2, col="darkgray")
      
      
      #abline(h=x_position-0.75, lty=2)
      # rect(x_position, min(most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df$Lg[chr_i],]$start), x_position+1,
      #       max(most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df$Lg[chr_i],]$end), col="gray")
      tick_positions[chr_i]<-x_position-(0.5)
      x_position<-x_position-2
      #subset table
      chr_table<-qtl_table[qtl_table$chr_number==chr_i,]
      prev_cross=0
      if(nrow(chr_table)>0){
        
        for(chr_qtl_i in 1:nrow(chr_table)){
          
          trait=as.character(chr_table$trait[chr_qtl_i])
          cross=as.numeric(chr_table$cross[chr_qtl_i])
          if(trait=="FRE"){
            arrows(
              chr_table$begin_pos[chr_qtl_i],
              x_position,
              chr_table$end_pos[chr_qtl_i],
              x_position,
              lwd=3, col="blue", angle=90, code=3, length=0.05
            )
            #Arrows(x_position, 
            #       chr_table$begin_pos[chr_qtl_i],
            #       x_position,
            #       chr_table$end_pos[chr_qtl_i],
            #       code = chr_table$direction[chr_qtl_i],
            #       lwd=4, arr.type="triangle", col="blue", arr.width = 0.2,arr.length = 0.2)
            
            #points(chr_table$begin_pos[chr_qtl_i],x_position, pch="|", cex=1, col="blue")
            #points(chr_table$end_pos[chr_qtl_i],x_position, pch="|", cex=1, col="blue")
            if(chr_table$direction[chr_qtl_i]==2){
              points(chr_table$center_pos[chr_qtl_i],x_position, pch=24,col="black", bg="blue")
              
            }
            else{
              points(chr_table$center_pos[chr_qtl_i],x_position, pch=25,col="black", bg="blue")
              
            }
            
          } 
          else if(trait=="MRE"){
            arrows(
                    chr_table$begin_pos[chr_qtl_i],
                    x_position,
                    chr_table$end_pos[chr_qtl_i],
                    x_position,
                    lwd=3, col="red", angle=90, code=3, length=0.05
            )
            #   Arrows(x_position, 
            #         chr_table$begin_pos[chr_qtl_i],
            #         x_position,
            #         chr_table$end_pos[chr_qtl_i],
            #         code = chr_table$direction[chr_qtl_i],
            #         lwd=4, arr.type="triangle", col="red", arr.width = 0.2 ,arr.length = 0.2)
            #points(chr_table$begin_pos[chr_qtl_i],x_position, pch="|", cex=1, col="red")
            #points(chr_table$end_pos[chr_qtl_i],x_position, pch="|", cex=1, col="red")
            if(chr_table$direction[chr_qtl_i]==2){
              points(chr_table$center_pos[chr_qtl_i],x_position, pch=24,col="black", bg="red")
              
            }
            else{
              points(chr_table$center_pos[chr_qtl_i],x_position, pch=25,col="black", bg="red")
              
            }
          }
          else if(trait=="rMA"){
            arrows(chr_table$begin_pos[chr_qtl_i],
                     x_position,
                     chr_table$end_pos[chr_qtl_i],
                     x_position,
                     lwd=3, col="darkgreen", angle=90, code=3, length=0.05
            )
            #Arrows(x_position, 
            #       chr_table$begin_pos[chr_qtl_i],
            #       x_position,
            #       chr_table$end_pos[chr_qtl_i],
            #       code = chr_table$direction[chr_qtl_i],
            #       lwd=4, arr.type="triangle", col="darkgreen", arr.width = 0.2 ,arr.length = 0.2)
            #points(chr_table$begin_pos[chr_qtl_i],x_position, pch="|", cex=1, col="darkgreen")
            #points(chr_table$end_pos[chr_qtl_i],x_position, pch="|", cex=1, col="darkgreen")
            if(chr_table$direction[chr_qtl_i]==2){
              points(chr_table$center_pos[chr_qtl_i],x_position, pch=24,col="black", bg="darkgreen")
              
            }
            else{
              points(chr_table$center_pos[chr_qtl_i],x_position, pch=25,col="black", bg="darkgreen")
              
            }
          }
          
          
          if(cross!=prev_cross && prev_cross!=0){
            #abline(h=x_position+1, lty=2)
            segments(0,
                     x_position+1,
                     95000000,
                     x_position+1, lty=2, col="darkgrey")
          }
          if(cross!=prev_cross){
            text( 85000000,x_position, labels = paste("Cross ", as.character(cross)), cex=1.2)
          }
          x_position<-x_position-2
          prev_cross=cross
        }
      }
      else{
        x_position<-x_position-2
      }
      
    }
    
    axis(2, at = tick_positions, labels = c("Chr1", "Chr2","Chr3", "Chr4","Chr5", "Chr6","Chr7", "Chr8" ), las=2, lwd=0, line = -1, cex.axis=1.3)
    axis(1, at = c(0, 20000000,40000000,60000000,80000000), labels = c(0, 20,40, 60,80), cex.axis=1.2)
    
    ############Second plot
    par(mar=c(4,5,2,1))
    
    out_1A_rMA_summary=summary(out_1A.hk_rMA)
    out_1A_rMA_qtl_pos=out_1A_rMA_summary$pos[out_1A_rMA_summary$chr=="Lg5"]
    out_1A_rMA_qtl_LOD=out_1A_rMA_summary$lod[out_1A_rMA_summary$chr=="Lg5"]
    
    out_1A_FRE_summary=summary(out_1A.hk_FRE)
    out_1A_FRE_qtl_pos=out_1A_FRE_summary$pos[out_1A_FRE_summary$chr=="Lg5"]
    out_1A_FRE_qtl_LOD=out_1A_FRE_summary$lod[out_1A_FRE_summary$chr=="Lg5"]
    
    out_1A_MRE_summary=summary(out_1A.hk_MRE)
    out_1A_MRE_qtl_pos=out_1A_MRE_summary$pos[out_1A_MRE_summary$chr=="Lg5"]
    out_1A_MRE_qtl_LOD=out_1A_MRE_summary$lod[out_1A_MRE_summary$chr=="Lg5"]
    
    plot(out_1A.hk_rMA, lodcolumn = 1, col="darkgreen", ylab="LOD", chr="Lg5", xlab="Map position (cM)", cex.lab=1.3, cex.axis=1.2)

    #limit lines
    #abline(v=out_1A_rMA_qtl_pos, col="darkgreen", lty=2, lwd=2)
    segments(out_1A_rMA_qtl_pos,out_1A_rMA_qtl_LOD ,out_1A_rMA_qtl_pos,out_1A_MRE_qtl_LOD,col="darkgreen", lty=2, lwd=2)
    #abline(h=summary(operm_1A.hk_FRE)[1], col="darkgray", lty=2, lwd=2)
    segments(0, summary(operm_1A.hk_FRE)[1], 80, summary(operm_1A.hk_FRE)[1], col="darkgray", lty=2, lwd=2)
    plot(out_1A.hk_MRE, lodcolumn = 1, col="red", ylab="LOD", chr="Lg5", add=TRUE)
    segments(out_1A_MRE_qtl_pos,out_1A_MRE_qtl_LOD ,out_1A_MRE_qtl_pos,out_1A_FRE_qtl_LOD,col="red", lty=2, lwd=2)
    #abline(h=summary(operm_1A.hk_FRE)[1], col="red", lty=2, lwd=2)
    plot(out_1A.hk_FRE, lodcolumn = 1, col="blue", ylab="LOD", chr="Lg5", add=TRUE)
    segments(out_1A_FRE_qtl_pos,out_1A_FRE_qtl_LOD ,out_1A_FRE_qtl_pos,0,col="blue", lwd=2, lty=2)
    #abline(h=summary(operm_1A.hk_FRE)[1], col="blue", lty=2, lwd=2)
    #xaxisloc.scanone(out_1A.hk_FRE, thechr = summary(out_1A.hk_FRE)[3,1], summary(out_1A.hk_FRE)[3,2])
    #abline(v=xaxisloc.scanone(out_1A.hk_FRE, thechr = summary(out_1A.hk_FRE)[3,1], summary(out_1A.hk_FRE)[3,2]), col="darkgreen", lty=2, lwd=2)
    abline(v=xaxisloc.scanone(out_1A.hk_FRE, thechr = summary(out_1A.hk_FRE)[5,1], summary(out_1A.hk_FRE)[5,2]), col="blue", lty=2, lwd=2)
    title("Cross 1, Chr5")
    
    text(par('usr')[1]-(0.1*par('usr')[1]),par('usr')[4]*1.035,labels=expression(bold("B")) , cex=1.5)
    
    #cross_1A_df
    #Chr 5 center qtl OW569314.1_15482270
    genotypes_1A_chr5_rMA_qtl<-genotypes_1A[,"OW569314.1_15482270"]
    df_gt_1A_chr5_rMA_qtl_FRE<-data.frame(genotypes=as.factor( genotypes_1A_chr5_rMA_qtl), phenotypes=phenotypes_1A$FRE)
    df_gt_1A_chr5_rMA_qtl_MRE<-data.frame(genotypes=as.factor( genotypes_1A_chr5_rMA_qtl), phenotypes=phenotypes_1A$MRE)
    df_gt_1A_chr5_rMA_qtl_colours<-genotypes_1A_chr5_rMA_qtl
    #set colours manually
    df_gt_1A_chr5_rMA_qtl_colours[df_gt_1A_chr5_rMA_qtl_colours==1]<-"#884ea0"
    #purple?
    df_gt_1A_chr5_rMA_qtl_colours[df_gt_1A_chr5_rMA_qtl_colours==2]<-"darkgrey"
    #orange?
    df_gt_1A_chr5_rMA_qtl_colours[df_gt_1A_chr5_rMA_qtl_colours==3]<-"#e67e22"
    #set symbols manually  
    df_gt_1A_chr5_rMA_qtl_pchs<-genotypes_1A_chr5_rMA_qtl
    df_gt_1A_chr5_rMA_qtl_pchs[df_gt_1A_chr5_rMA_qtl_pchs==1]<-17
    df_gt_1A_chr5_rMA_qtl_pchs[df_gt_1A_chr5_rMA_qtl_pchs==2]<-18
    df_gt_1A_chr5_rMA_qtl_pchs[df_gt_1A_chr5_rMA_qtl_pchs==3]<-19
    df_gt_1A_chr5_rMA_qtl_FRE_MRE<-data.frame(genotypes=as.factor( genotypes_1A_chr5_rMA_qtl), MRE=phenotypes_1A$MRE, FRE=phenotypes_1A$FRE)
    
    #Pretty nice...add histograms to sides?
    plot(df_gt_1A_chr5_rMA_qtl_FRE_MRE$MRE, df_gt_1A_chr5_rMA_qtl_FRE_MRE$FRE, col=df_gt_1A_chr5_rMA_qtl_colours, pch=df_gt_1A_chr5_rMA_qtl_pchs, xlab="Male reproductive effort", ylab="Female reproductive effort", cex.lab=1.3, cex.axis=1.2)
    legend("topright", legend = c("AA", "AB", "BB"), pch = c(17,18,19), col=c("#884ea0","darkgrey","#e67e22"), title = "Genotypes", cex=1.3)
    text(par('usr')[1]-(0.1*par('usr')[1]),par('usr')[4]*1.035,labels=expression(bold("C")) , cex=1.5)
    
    
    dev.off()
    
    
    
    
    
    