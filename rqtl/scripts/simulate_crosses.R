options(scipen=999)
setwd("/home/fredo/Stuff/qtl_paper/upload_suppl_information/rqtl/")
genome_info_df<-data.frame(chr=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8"),
                           Lg=c("OW569319.1","OW569313.1","OW569312.1","OW569317.1","OW569314.1","OW569316.1","OW569315.1","OW569318.1"),
                           size=c(74045953, 58100813, 76280018, 41830924, 56051264, 47723244, 50031057, 43611580))
plot_color<-function(trait){
  if(trait=="FRE"){
    return("#E69F00")
  }
  if(trait=="MRE"){
    return("#56B4E9")  
  }
  if(trait=="rMA"){
    return("#009E73")
  }
  #col_FRE="#E69F00"
  #rgb_FRE=rgb(0.9,0.62,0, alpha=0.5)
  #col_MRE="#56B4E9"
  #rgb_MRE=rgb(0.34,0.71,0.91, alpha=0.5)
  #col_rMA="#009E73"
  #rgb_rMA=rgb(0,0.62,0.45, alpha=0.5)
  
}

plot_color_alpha<-function(trait){
  if(trait=="FRE"){
    return(rgb(0.9,0.62,0, alpha=0.5))
  }
  if(trait=="MRE"){
    return(rgb(0.34,0.71,0.91, alpha=0.5))  
  }
  if(trait=="rMA"){
    return(rgb(0,0.62,0.45, alpha=0.5))
  }
}


centromeric_repeats<-read.csv("data/Summary.of.repetitive.regions.GCA_937616625.1_ddMerAnnu1.1_genomic.fna.csv", header=TRUE)
most_common_centromeric_repeats<-centromeric_repeats[centromeric_repeats$most.freq.value.N==60 | centromeric_repeats$most.freq.value.N==99 ,]


col_FRE="#E69F00"
#rgb_FRE=rgb(0.9,0.62,0, alpha=0.5)
col_MRE="#56B4E9"
#rgb_MRE=rgb(0.34,0.71,0.91, alpha=0.5)
col_rMA="#009E73"
#rgb_rMA=rgb(0,0.62,0.45, alpha=0.5)



library(qtl)
n_its<-100

cross_obj_1A<-read.cross(format="csvr", file ="data/1A.csvr", alleles=c("A","B"))
cross_obj_1A<-jittermap(cross_obj_1A)
cross_1A_lm<-pull.map(cross_obj_1A, as.table = TRUE)
cross_1A_LGs<-unique(cross_1A_lm$chr)
cross_1A_lm_list<-list()
cross_1A_lm_list$`1`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg1"]
names(cross_1A_lm_list$`1`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg1",])

cross_1A_lm_list$`2`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg2"]
names(cross_1A_lm_list$`2`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg2",])

cross_1A_lm_list$`3`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg3"]
names(cross_1A_lm_list$`3`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg3",])

cross_1A_lm_list$`4`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg4"]
names(cross_1A_lm_list$`4`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg4",])

cross_1A_lm_list$`5`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg5"]
names(cross_1A_lm_list$`5`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg5",])

cross_1A_lm_list$`6`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg6"]
names(cross_1A_lm_list$`6`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg6",])

cross_1A_lm_list$`7`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg7"]
names(cross_1A_lm_list$`7`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg7",])

cross_1A_lm_list$`8`<-cross_1A_lm$pos[cross_1A_lm$chr=="Lg8"]
names(cross_1A_lm_list$`8`)<-rownames(cross_1A_lm[cross_1A_lm$chr=="Lg8",])

cross_1_qtl<-rbind(c(3,9.6,0.079,21.67,6,2386203,5708265,1432222,0),
                   c(5,37.97,0.1167,63,29.99,13897875,52855638,8332668,1),
                   c(4,48,0.1053,52.9,40.55,36583907,38902564,34430945,0),
                   c(5,39.92,0.1418,45.57,24.32,15605703,42893618,6318978,1),
                   c(5,39.73,0.2409,42.07,36.81,15482270,16308679,13172098,0),
                   c(6,51,0.0385,70.56,36.66,39321199,45632260,11609975,1))
rownames(cross_1_qtl)<-c("FRE","FRE","MRE","MRE","rMA","rMA")
colnames(cross_1_qtl)<-c("Lg","Pos","Var", "cm_end","cm_start", "bp_center","bp_end","bp_start", "cent_span")


cross_1_qtl_intervals_start<-vector("list", nrow(cross_1_qtl))
cross_1_qtl_intervals_end<-vector("list", nrow(cross_1_qtl))
cross_1_qtl_intervals_bp_start<-vector("list", nrow(cross_1_qtl))
cross_1_qtl_intervals_bp_end<-vector("list", nrow(cross_1_qtl))

#iterate over QTL
for(cross_1_qtl_i in 1:nrow(cross_1_qtl)){
  #iterate over replicate simulations
  qtl_intervals_start<-1:n_its
  qtl_intervals_end<-1:n_its
  qtl_intervals_bp_start<-1:n_its
  qtl_intervals_bp_end<-1:n_its
  
  qtl_alpha<-sqrt(2*cross_1_qtl[cross_1_qtl_i,3]/(1-cross_1_qtl[cross_1_qtl_i,3]))
  for(it_i in 1:n_its){
    print(it_i)
    sim_cross<-sim.cross(cross_1A_lm_list, n.ind=247, type="f2", c(cross_1_qtl[cross_1_qtl_i,1],cross_1_qtl[cross_1_qtl_i,2],qtl_alpha,0)  )
    sim_cross<-calc.genoprob(sim_cross, step=0.5)
    sim_hk<-scanone(sim_cross, method="hk")
    sim_bi<-bayesint(results=sim_hk, chr=cross_1_qtl[cross_1_qtl_i,1],expandtomarkers = T)
    qtl_length_cm_end<-sim_bi$pos[3]
    qtl_length_cm_start<-sim_bi$pos[1]
    qtl_intervals_start[it_i]<-round(qtl_length_cm_start, digits = 2)
    qtl_intervals_end[it_i]<-round(qtl_length_cm_end, digits=2)
    
    qtl_length_bp_end<-as.numeric(unlist(strsplit(row.names(sim_bi)[3],split = "_"))[2])
    qtl_length_bp_start<-as.numeric(unlist(strsplit(row.names(sim_bi)[1],split = "_"))[2])
    qtl_intervals_bp_start[it_i]<-qtl_length_bp_start
    qtl_intervals_bp_end[it_i]<-qtl_length_bp_end
    
  }
  cross_1_qtl_intervals_start[[cross_1_qtl_i]]<-qtl_intervals_start
  cross_1_qtl_intervals_end[[cross_1_qtl_i]]<-qtl_intervals_end
  
  cross_1_qtl_intervals_bp_start[[cross_1_qtl_i]]<-qtl_intervals_bp_start
  cross_1_qtl_intervals_bp_end[[cross_1_qtl_i]]<-qtl_intervals_bp_end
}


pdf("plots/Cross1_simulations_2.pdf", width=8, height=12)

#try layout with titles instead...
#par(mfrow=c(nrow(cross_1_qtl),2))
#par(mar = c(2, 2, 2, 2))

layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,13,14,15,16,16,17,18), byrow=T,ncol=2),heights=c(0.75,3,0.75,3,0.75,3,0.75,3,0.75,3,0.75,3))
sub_ids<-c("A","B","C","D","E","F")

for(cross_1_plot_i in 1:nrow(cross_1_qtl)){
  #cM
  #First title
  par(mar=c(0,0.5,0,0.5))
  plot.new()
  text(0.5,0.5,paste("Chromosome ",cross_1_qtl[cross_1_plot_i,1],", ",rownames(cross_1_qtl)[cross_1_plot_i],sep=""),cex=1.5,font=1.5)
  text(-0.01, 0.5, sub_ids[cross_1_plot_i], cex=1.5,font=1.5)
  if(cross_1_plot_i<nrow(cross_1_qtl)){
    par(mar=c(2,0.5,0,0.5))
  }
  else{
    par(mar=c(5,0.5,0,0.5))
  }
  
  qtl_lengths<-cross_1_qtl_intervals_end[[cross_1_plot_i]]-cross_1_qtl_intervals_start[[cross_1_plot_i]]
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, max(unlist(cross_1A_lm_list[cross_1_qtl[cross_1_plot_i,1]]))), ylab="", yaxt='n', frame.plot = FALSE, xlab="", xaxt="n")
  chrom_ticks_cm<-seq(0,max(unlist(cross_1A_lm_list[cross_1_qtl[cross_1_plot_i,1]])),10)
  axis(1, at = c(chrom_ticks_cm,max(unlist(cross_1A_lm_list[cross_1_qtl[cross_1_plot_i,1]]))), labels = c(chrom_ticks_cm,""), cex.axis=1.2)
  
  Plot_dists=seq(0.1,10,0.1)
  segment_colors<-rep("black", n_its)
  segment_colors[cross_1_qtl_intervals_start[[cross_1_plot_i]]<=cross_1_qtl[cross_1_plot_i,5] & cross_1_qtl_intervals_end[[cross_1_plot_i]]>=cross_1_qtl[cross_1_plot_i,4]]<-plot_color(rownames(cross_1_qtl)[cross_1_plot_i])
  segments(cross_1_qtl_intervals_start[[cross_1_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists,cross_1_qtl_intervals_end[[cross_1_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors[order(qtl_lengths, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_1_qtl[cross_1_plot_i,2], col=plot_color(rownames(cross_1_qtl)[cross_1_plot_i]), lwd=1)
  abline(v=cross_1_qtl[cross_1_plot_i,4], col=plot_color(rownames(cross_1_qtl)[cross_1_plot_i]), lwd=1, lty=2)
  abline(v=cross_1_qtl[cross_1_plot_i,5], col=plot_color(rownames(cross_1_qtl)[cross_1_plot_i]), lwd=1, lty=2)
  if(cross_1_plot_i==nrow(cross_1_qtl)){
    title(xlab="genetic distance (cM)", cex.lab=1.5)
  }
  
  #bp
  qtl_lengths_bp<-cross_1_qtl_intervals_bp_end[[cross_1_plot_i]]-cross_1_qtl_intervals_bp_start[[cross_1_plot_i]]
  chrom_ticks<-seq(0,ceiling(genome_info_df$size[cross_1_qtl[cross_1_plot_i,1]]/1000000), 10)
  
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, genome_info_df$size[cross_1_qtl[cross_1_plot_i,1]]), ylab="", yaxt='n', frame.plot = FALSE, xlab="",xaxt="n")
  axis(1, at = c(chrom_ticks*1000000,genome_info_df$size[cross_1_qtl[cross_1_plot_i,1]]), labels = c(chrom_ticks,""), cex.axis=1.2)
  
  segment_colors_bp<-rep("black", n_its)
  segment_colors_bp[cross_1_qtl_intervals_bp_start[[cross_1_plot_i]]<=cross_1_qtl[cross_1_plot_i,8] & cross_1_qtl_intervals_bp_end[[cross_1_plot_i]]>=cross_1_qtl[cross_1_plot_i,7]]<-plot_color(rownames(cross_1_qtl)[cross_1_plot_i])
  
  segments(cross_1_qtl_intervals_bp_start[[cross_1_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists,cross_1_qtl_intervals_bp_end[[cross_1_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors_bp[order(qtl_lengths_bp, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_1_qtl[cross_1_plot_i,6], col=plot_color(rownames(cross_1_qtl)[cross_1_plot_i]), lwd=1)
  abline(v=cross_1_qtl[cross_1_plot_i,7], col=plot_color(rownames(cross_1_qtl)[cross_1_plot_i]), lwd=1, lty=2)
  abline(v=cross_1_qtl[cross_1_plot_i,8], col=plot_color(rownames(cross_1_qtl)[cross_1_plot_i]), lwd=1, lty=2)
  #get centromeres
  chrom_centromeres<-most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df[cross_1_qtl[cross_1_plot_i,1],2],]
  rect(chrom_centromeres$start, -0.2,chrom_centromeres$start,-0.1, col="purple", border = "purple")
}
title(xlab="physical distance (Mb)", cex.lab=1.5)

dev.off()



#Cross 2
cross_obj_S1<-read.cross(format="csvr", file ="data/S1.csvr", alleles=c("A","B"))
cross_obj_S1<-jittermap(cross_obj_S1)
cross_S1_lm<-pull.map(cross_obj_S1, as.table = TRUE)
cross_S1_LGs<-unique(cross_S1_lm$chr)
cross_S1_lm_list<-list()
cross_S1_lm_list$`1`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg1"]
names(cross_S1_lm_list$`1`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg1",])

cross_S1_lm_list$`2`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg2"]
names(cross_S1_lm_list$`2`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg2",])

cross_S1_lm_list$`3`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg3"]
names(cross_S1_lm_list$`3`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg3",])

cross_S1_lm_list$`4`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg4"]
names(cross_S1_lm_list$`4`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg4",])

cross_S1_lm_list$`5`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg5"]
names(cross_S1_lm_list$`5`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg5",])

cross_S1_lm_list$`6`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg6"]
names(cross_S1_lm_list$`6`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg6",])

cross_S1_lm_list$`7`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg7"]
names(cross_S1_lm_list$`7`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg7",])

cross_S1_lm_list$`8`<-cross_S1_lm$pos[cross_S1_lm$chr=="Lg8"]
names(cross_S1_lm_list$`8`)<-rownames(cross_S1_lm[cross_S1_lm$chr=="Lg8",])

cross_2_qtl<-rbind(c(7,24.24,0.2078,24.24,14.48,8706833, 8706833,2916894,0),
                   c(4,39.37,0.0471,48.42,22.3, 30877784,33657137,25565755,0),
                   c(7,7.26,0.3023,13.13,1.89,1010464,2550366,849466,0),
                   c(4,37.34,0.0471,44.24,23.95,30465469,32646736,26044219,0),
                   c(7,7.26,0.3199,16.26,1.89,1010464,4184838,849466,0))

rownames(cross_2_qtl)<-c("FRE","MRE","MRE","rMA","rMA")
colnames(cross_2_qtl)<-c("Lg","Pos","Var", "cm_end","cm_start", "bp_center","bp_end","bp_start", "cent_span")


cross_2_qtl_intervals_start<-vector("list", nrow(cross_2_qtl))
cross_2_qtl_intervals_end<-vector("list", nrow(cross_2_qtl))
cross_2_qtl_intervals_bp_start<-vector("list", nrow(cross_2_qtl))
cross_2_qtl_intervals_bp_end<-vector("list", nrow(cross_2_qtl))

#iterate over QTL
for(cross_2_qtl_i in 1:nrow(cross_2_qtl)){
  #iterate over replicate simulations
  qtl_intervals_start<-1:n_its
  qtl_intervals_end<-1:n_its
  qtl_intervals_bp_start<-1:n_its
  qtl_intervals_bp_end<-1:n_its
  
  qtl_alpha<-sqrt(2*cross_2_qtl[cross_2_qtl_i,3]/(1-cross_2_qtl[cross_2_qtl_i,3]))
  for(it_i in 1:n_its){
    print(it_i)
    sim_cross<-sim.cross(cross_S1_lm_list, n.ind=368, type="f2", c(cross_2_qtl[cross_2_qtl_i,1],cross_2_qtl[cross_2_qtl_i,2],qtl_alpha,0)  )
    sim_cross<-calc.genoprob(sim_cross, step=0.5)
    sim_hk<-scanone(sim_cross, method="hk")
    sim_bi<-bayesint(results=sim_hk, chr=cross_2_qtl[cross_2_qtl_i,1],expandtomarkers = T)
    qtl_length_cm_end<-sim_bi$pos[3]
    qtl_length_cm_start<-sim_bi$pos[1]
    qtl_intervals_start[it_i]<-qtl_length_cm_start
    qtl_intervals_end[it_i]<-qtl_length_cm_end
    
    qtl_length_bp_end<-as.numeric(unlist(strsplit(row.names(sim_bi)[3],split = "_"))[2])
    qtl_length_bp_start<-as.numeric(unlist(strsplit(row.names(sim_bi)[1],split = "_"))[2])
    qtl_intervals_bp_start[it_i]<-qtl_length_bp_start
    qtl_intervals_bp_end[it_i]<-qtl_length_bp_end
    
  }
  cross_2_qtl_intervals_start[[cross_2_qtl_i]]<-round(qtl_intervals_start, digits = 2)
  cross_2_qtl_intervals_end[[cross_2_qtl_i]]<-round(qtl_intervals_end, digits = 2)
  
  cross_2_qtl_intervals_bp_start[[cross_2_qtl_i]]<-qtl_intervals_bp_start
  cross_2_qtl_intervals_bp_end[[cross_2_qtl_i]]<-qtl_intervals_bp_end
}

pdf("plots/Cross2_simulations_2.pdf", width=8, height=10)
#par(mfrow=c(nrow(cross_2_qtl),2))
#par(mar = c(2, 2, 2, 2))

layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,13,14,15), byrow=T,ncol=2),heights=c(0.75,3,0.75,3,0.75,3,0.75,3,0.75,3))
sub_ids<-c("A","B","C","D","E")

for(cross_2_plot_i in 1:nrow(cross_2_qtl)){
  #First title
  par(mar=c(0,0.5,0,0.5))
  plot.new()
  text(0.5,0.5,paste("Chromosome ",cross_2_qtl[cross_2_plot_i,1],", ",rownames(cross_2_qtl)[cross_2_plot_i],sep=""),cex=1.5,font=1.5)
  text(-0.01, 0.5, sub_ids[cross_2_plot_i], cex=1.5,font=1.5)
  if(cross_2_plot_i<nrow(cross_2_qtl)){
    par(mar=c(2,0.5,0,0.5))
  }
  else{
    par(mar=c(5,0.5,0,0.5))
  }
  
  #cM
  qtl_lengths<-cross_2_qtl_intervals_end[[cross_2_plot_i]]-cross_2_qtl_intervals_start[[cross_2_plot_i]]
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, max(unlist(cross_S1_lm_list[cross_2_qtl[cross_2_plot_i,1]]))), ylab="", yaxt='n',, xlab="", frame.plot = FALSE,xaxt="n")
  chrom_ticks_cm<-seq(0,max(unlist(cross_S1_lm_list[cross_2_qtl[cross_2_plot_i,1]])),10)
  axis(1, at = c(chrom_ticks_cm,max(unlist(cross_S1_lm_list[cross_2_qtl[cross_2_plot_i,1]]))), labels = c(chrom_ticks_cm,""), cex.axis=1.2)
  
  Plot_dists=seq(0.1,10,0.1)
  segment_colors<-rep("black", n_its)
  segment_colors[cross_2_qtl_intervals_start[[cross_2_plot_i]]<=cross_2_qtl[cross_2_plot_i,5] & cross_2_qtl_intervals_end[[cross_2_plot_i]]>=cross_2_qtl[cross_2_plot_i,4]]<-plot_color(rownames(cross_2_qtl)[cross_2_plot_i])
  segments(cross_2_qtl_intervals_start[[cross_2_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists,cross_2_qtl_intervals_end[[cross_2_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors[order(qtl_lengths, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_2_qtl[cross_2_plot_i,2], col=plot_color(rownames(cross_2_qtl)[cross_2_plot_i]), lwd=1)
  abline(v=cross_2_qtl[cross_2_plot_i,4], col=plot_color(rownames(cross_2_qtl)[cross_2_plot_i]), lwd=1, lty=2)
  abline(v=cross_2_qtl[cross_2_plot_i,5], col=plot_color(rownames(cross_2_qtl)[cross_2_plot_i]), lwd=1, lty=2)
  if(cross_2_plot_i==nrow(cross_2_qtl)){
    title(xlab="genetic distance (cM)", cex.lab=1.5)
  }
  
  #bp
  qtl_lengths_bp<-cross_2_qtl_intervals_bp_end[[cross_2_plot_i]]-cross_2_qtl_intervals_bp_start[[cross_2_plot_i]]
  chrom_ticks<-seq(0,ceiling(genome_info_df$size[cross_2_qtl[cross_2_plot_i,1]]/1000000), 10)
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, genome_info_df$size[cross_2_qtl[cross_2_plot_i,1]]), ylab="", yaxt='n', xlab="", frame.plot = FALSE, xaxt="n")
  axis(1, at = c(chrom_ticks*1000000,genome_info_df$size[cross_2_qtl[cross_2_plot_i,1]]), labels = c(chrom_ticks,""), cex.axis=1.2)
  segment_colors_bp<-rep("black", n_its)
  segment_colors_bp[cross_2_qtl_intervals_bp_start[[cross_2_plot_i]]<=cross_2_qtl[cross_2_plot_i,8] & cross_2_qtl_intervals_bp_end[[cross_2_plot_i]]>=cross_2_qtl[cross_2_plot_i,7]]<-plot_color(rownames(cross_2_qtl)[cross_2_plot_i])
  segments(cross_2_qtl_intervals_bp_start[[cross_2_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists,cross_2_qtl_intervals_bp_end[[cross_2_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors_bp[order(qtl_lengths_bp, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_2_qtl[cross_2_plot_i,6], col=plot_color(rownames(cross_2_qtl)[cross_2_plot_i]), lwd=1)
  abline(v=cross_2_qtl[cross_2_plot_i,7], col=plot_color(rownames(cross_2_qtl)[cross_2_plot_i]), lwd=1, lty=2)
  abline(v=cross_2_qtl[cross_2_plot_i,8], col=plot_color(rownames(cross_2_qtl)[cross_2_plot_i]), lwd=1, lty=2)
  chrom_centromeres<-most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df[cross_2_qtl[cross_2_plot_i,1],2],]
  rect(chrom_centromeres$start, -0.2,chrom_centromeres$start,-0.1, col="purple", border = "purple")
  
}

title(xlab="physical distance (Mb)", cex.lab=1.5)

dev.off()

#Cross 3
cross_obj_3A<-read.cross(format="csvr", file ="data/3A.csvr", alleles=c("A","B"))
phenotypes_3A<-pull.pheno(cross_obj_3A)
#remove MRE outliers for QTL analysis!
cross_obj_sub_3A<-subset(cross_obj_3A, ind=phenotypes_3A$MRE<0.01)
cross_obj_sub_3A<-jittermap(cross_obj_sub_3A)
cross_3A_lm<-pull.map(cross_obj_sub_3A, as.table = TRUE)

cross_3A_LGs<-unique(cross_3A_lm$chr)
cross_3A_lm_list<-list()
cross_3A_lm_list$`1`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg1"]
names(cross_3A_lm_list$`1`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg1",])

cross_3A_lm_list$`2`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg2"]
names(cross_3A_lm_list$`2`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg2",])

cross_3A_lm_list$`3`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg3"]
names(cross_3A_lm_list$`3`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg3",])

cross_3A_lm_list$`4`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg4"]
names(cross_3A_lm_list$`4`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg4",])

cross_3A_lm_list$`5`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg5"]
names(cross_3A_lm_list$`5`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg5",])

cross_3A_lm_list$`6`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg6"]
names(cross_3A_lm_list$`6`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg6",])

cross_3A_lm_list$`7`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg7"]
names(cross_3A_lm_list$`7`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg7",])

cross_3A_lm_list$`8`<-cross_3A_lm$pos[cross_3A_lm$chr=="Lg8"]
names(cross_3A_lm_list$`8`)<-rownames(cross_3A_lm[cross_3A_lm$chr=="Lg8",])


cross_3_qtl<-rbind(c(4,42.33,0.0523,50.05,2.85, 31616069 ,34602851,4340717,1))


rownames(cross_3_qtl)<-c("FRE")
colnames(cross_3_qtl)<-c("Lg","Pos","Var", "cm_end","cm_start", "bp_center","bp_end","bp_start", "cent_span")
cross_3_qtl_intervals_start<-1:n_its
cross_3_qtl_intervals_end<-1:n_its
cross_3_qtl_intervals_bp_start<-1:n_its
cross_3_qtl_intervals_bp_end<-1:n_its
qtl_alpha<-sqrt(2*cross_3_qtl[1,3]/(1-cross_3_qtl[1,3]))
for(it_i in 1:n_its){
  print(it_i)
  sim_cross<-sim.cross(cross_3A_lm_list, n.ind=349, type="f2", c(cross_3_qtl[1,1],cross_3_qtl[1,2],qtl_alpha,0)  )
  sim_cross<-calc.genoprob(sim_cross, step=0.5)
  sim_hk<-scanone(sim_cross, method="hk")
  sim_bi<-bayesint(results=sim_hk, chr=cross_3_qtl[1,1],expandtomarkers = T)
  qtl_length_cm_end<-sim_bi$pos[3]
  qtl_length_cm_start<-sim_bi$pos[1]
  cross_3_qtl_intervals_start[it_i]<-round(qtl_length_cm_start, digits=2)
  cross_3_qtl_intervals_end[it_i]<-round(qtl_length_cm_end, digits=2)
  
  qtl_length_bp_end<-as.numeric(unlist(strsplit(row.names(sim_bi)[3],split = "_"))[2])
  qtl_length_bp_start<-as.numeric(unlist(strsplit(row.names(sim_bi)[1],split = "_"))[2])
  cross_3_qtl_intervals_bp_start[it_i]<-qtl_length_bp_start
  cross_3_qtl_intervals_bp_end[it_i]<-qtl_length_bp_end
  
}


pdf("plots/Cross3_simulations_2.pdf", width=8, height=2)
#par(mfrow=c(nrow(cross_2_qtl),2))
#par(mar = c(2, 2, 2, 2))

layout(matrix(c(1,1,2,3), byrow=T,ncol=2),heights=c(0.75,3))

par(mar=c(0,0.5,0,0.5))
plot.new()
text(0.5,0.5,paste("Chromosome ",cross_3_qtl[1,1],", ",rownames(cross_3_qtl)[1],sep=""),cex=1.5,font=1.5)
par(mar=c(5,0.5,0,0.5))
qtl_lengths<-cross_3_qtl_intervals_end-cross_3_qtl_intervals_start
plot(0, cex=0, ylim=c(0,10), xlim=c(0, max(unlist(cross_3A_lm_list[cross_3_qtl[1,1]]))), ylab="", yaxt='n',xlab="", frame.plot = FALSE, xaxt="n")
chrom_ticks_cm<-seq(0,max(unlist(cross_3A_lm_list[cross_3_qtl[1,1]])),10)
axis(1, at = c(chrom_ticks_cm,max(unlist(cross_3A_lm_list[cross_3_qtl[1,1]]))), labels = c(chrom_ticks_cm,""), cex.axis=1.2)

Plot_dists=seq(0.1,10,0.1)
segment_colors<-rep("black", n_its)
segment_colors[cross_3_qtl_intervals_start<=cross_3_qtl[1,5] & cross_3_qtl_intervals_end>=cross_3_qtl[1,4]]<-plot_color(rownames(cross_3_qtl)[1])
segments(cross_3_qtl_intervals_start[order(qtl_lengths, decreasing = T)],Plot_dists,cross_3_qtl_intervals_end[order(qtl_lengths, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors[order(qtl_lengths, decreasing = T)])
#plot empirical position and qtl location
abline(v=cross_3_qtl[1,2], col=plot_color(rownames(cross_3_qtl)[1]), lwd=1)
abline(v=cross_3_qtl[1,4], col=plot_color(rownames(cross_3_qtl)[1]), lwd=1, lty=2)
abline(v=cross_3_qtl[1,5], col=plot_color(rownames(cross_3_qtl)[1]), lwd=1, lty=2)
title(xlab="genetic distance (cM)", cex.lab=1.5)

#bp
qtl_lengths_bp<-cross_3_qtl_intervals_bp_end-cross_3_qtl_intervals_bp_start
plot(-0.01, cex=0, ylim=c(0,10), xlim=c(0, genome_info_df$size[cross_3_qtl[1,1]]), ylab="", yaxt='n', xlab="", frame.plot = FALSE, xaxt="n")
chrom_ticks<-seq(0,ceiling(genome_info_df$size[cross_3_qtl[1,1]]/1000000), 10)
axis(1, at = c(chrom_ticks*1000000,genome_info_df$size[cross_3_qtl[1,1]]), labels = c(chrom_ticks,""), cex.axis=1.2)



segment_colors_bp<-rep("black", n_its)
segment_colors_bp[cross_3_qtl_intervals_bp_start<=cross_3_qtl[1,8] & cross_3_qtl_intervals_bp_end>=cross_3_qtl[1,7]]<-plot_color(rownames(cross_3_qtl)[1])
segments(cross_3_qtl_intervals_bp_start[order(qtl_lengths_bp, decreasing = T)],Plot_dists,cross_3_qtl_intervals_bp_end[order(qtl_lengths_bp, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors_bp[order(qtl_lengths_bp, decreasing = T)])
#plot empirical position and qtl location
abline(v=cross_3_qtl[1,6], col=plot_color(rownames(cross_3_qtl)[1]), lwd=1)
abline(v=cross_3_qtl[1,7], col=plot_color(rownames(cross_3_qtl)[1]), lwd=1, lty=2)
abline(v=cross_3_qtl[1,8], col=plot_color(rownames(cross_3_qtl)[1]), lwd=1, lty=2)
chrom_centromeres<-most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df[cross_3_qtl[1,1],2],]
rect(chrom_centromeres$start, -0.2,chrom_centromeres$start,-0.1, col="purple", border = "purple")
title(xlab="physical distance (Mb)", cex.lab=1.5)

dev.off()
#Cross 4

cross_obj_5A<-read.cross(format="csvr", file ="data/5A.csvr", alleles=c("A","B"))
phenotypes_5A<-pull.pheno(cross_obj_5A)
cross_obj_5A_sub<-subset(cross_obj_5A, ind=phenotypes_5A$FRE<0.1)
cross_obj_5A_sub<-jittermap(cross_obj_5A_sub)

cross_5A_lm<-pull.map(cross_obj_5A_sub, as.table = TRUE)

cross_5A_LGs<-unique(cross_5A_lm$chr)
cross_5A_lm_list<-list()
cross_5A_lm_list$`1`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg1"]
names(cross_5A_lm_list$`1`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg1",])

cross_5A_lm_list$`2`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg2"]
names(cross_5A_lm_list$`2`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg2",])

cross_5A_lm_list$`3`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg3"]
names(cross_5A_lm_list$`3`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg3",])

cross_5A_lm_list$`4`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg4"]
names(cross_5A_lm_list$`4`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg4",])

cross_5A_lm_list$`5`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg5"]
names(cross_5A_lm_list$`5`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg5",])

cross_5A_lm_list$`6`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg6"]
names(cross_5A_lm_list$`6`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg6",])

cross_5A_lm_list$`7`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg7"]
names(cross_5A_lm_list$`7`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg7",])

cross_5A_lm_list$`8`<-cross_5A_lm$pos[cross_5A_lm$chr=="Lg8"]
names(cross_5A_lm_list$`8`)<-rownames(cross_5A_lm[cross_5A_lm$chr=="Lg8",])

cross_4_qtl<-rbind(c(4,12.67,0.0664,39.4,7.2,18962716,30528464,15400265,0),
                   c(7,72.42,0.0834,76.03,68.51,48389478,49795787,46913936,0),
                   c(1,107.41,0.0568,108.19,26.93,73153124,73999319,7174544,1),
                   c(3,69.44,0.0361,82.57,10.01,66670375,72420091,3720630,1),
                   c(4,46.43,0.0503,50.02,6.11,33482941,34415133,15294180,0),
                   c(5,57.85,0.1201,61.13,55.98,47058661,51295682,42860012,0),
                   c(7,65.37,0.0932,70.23,61.45,45234516,47610546,45234516,0),
                   c(1,107.41,0.0806,108.19,104.92,73153124,73999319,71824179,0),
                   c(3,73.2,0.0469,77.88,61.28,67773366,70015234,64377429,0),
                   c(4,11.73,0.0803,43.14,5.95,18685092,32496038,15074052,0),
                   c(5,58.16,0.1856,60.03,55.8,48799986,50965214,42389979,0),
                   c(7,34,0.0483,67.73,2.19,14106314,46685197,912681,1))

rownames(cross_4_qtl)<-c("FRE","FRE","MRE","MRE","MRE","MRE","MRE","rMA","rMA","rMA","rMA", "rMA")
colnames(cross_4_qtl)<-c("Lg","Pos","Var", "cm_end","cm_start", "bp_center","bp_end","bp_start", "cent_span")

cross_4_qtl_intervals_start<-vector("list", nrow(cross_4_qtl))
cross_4_qtl_intervals_end<-vector("list", nrow(cross_4_qtl))
cross_4_qtl_intervals_bp_start<-vector("list", nrow(cross_4_qtl))
cross_4_qtl_intervals_bp_end<-vector("list", nrow(cross_4_qtl))

#iterate over QTL
for(cross_4_qtl_i in 1:nrow(cross_4_qtl)){
  #iterate over replicate simulations
  qtl_intervals_start<-1:n_its
  qtl_intervals_end<-1:n_its
  qtl_intervals_bp_start<-1:n_its
  qtl_intervals_bp_end<-1:n_its
  
  qtl_alpha<-sqrt(2*cross_4_qtl[cross_4_qtl_i,3]/(1-cross_4_qtl[cross_4_qtl_i,3]))
  for(it_i in 1:n_its){
    print(it_i)
    sim_cross<-sim.cross(cross_5A_lm_list, n.ind=321, type="f2", c(cross_4_qtl[cross_4_qtl_i,1],cross_4_qtl[cross_4_qtl_i,2],qtl_alpha,0)  )
    sim_cross<-calc.genoprob(sim_cross, step=0.5)
    sim_hk<-scanone(sim_cross, method="hk")
    sim_bi<-bayesint(results=sim_hk, chr=cross_4_qtl[cross_4_qtl_i,1],expandtomarkers = T)
    qtl_length_cm_end<-sim_bi$pos[3]
    qtl_length_cm_start<-sim_bi$pos[1]
    qtl_intervals_start[it_i]<-qtl_length_cm_start
    qtl_intervals_end[it_i]<-qtl_length_cm_end
    
    qtl_length_bp_end<-as.numeric(unlist(strsplit(row.names(sim_bi)[3],split = "_"))[2])
    qtl_length_bp_start<-as.numeric(unlist(strsplit(row.names(sim_bi)[1],split = "_"))[2])
    qtl_intervals_bp_start[it_i]<-qtl_length_bp_start
    qtl_intervals_bp_end[it_i]<-qtl_length_bp_end
    
  }
  cross_4_qtl_intervals_start[[cross_4_qtl_i]]<-round(qtl_intervals_start, digits = 2)
  cross_4_qtl_intervals_end[[cross_4_qtl_i]]<-round(qtl_intervals_end, digits = 2)
  
  cross_4_qtl_intervals_bp_start[[cross_4_qtl_i]]<-qtl_intervals_bp_start
  cross_4_qtl_intervals_bp_end[[cross_4_qtl_i]]<-qtl_intervals_bp_end
}


pdf("plots/Cross4_simulations_2_1.pdf", width=8, height=12)
#par(mfrow=c(nrow(cross_2_qtl),2))
#par(mar = c(2, 2, 2, 2))

layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,13,14,15,16,16,17,18), byrow=T,ncol=2),heights=c(0.75,3,0.75,3,0.75,3,0.75,3,0.75,3,0.75,3))
sub_ids<-c("A","B","C","D","E","F")

for(cross_4_plot_i in 1:6){
  #First title
  par(mar=c(0,0.5,0,0.5))
  plot.new()
  text(0.5,0.5,paste("Chromosome ",cross_4_qtl[cross_4_plot_i,1],", ",rownames(cross_4_qtl)[cross_4_plot_i],sep=""),cex=1.5,font=1.5)
  text(-0.01, 0.5, sub_ids[cross_4_plot_i], cex=1.5,font=1.5)
  if(cross_4_plot_i<6){
    par(mar=c(2,0.5,0,0.5))
  }
  else{
    par(mar=c(5,0.5,0,0.5))
  }

  #cM
  qtl_lengths<-cross_4_qtl_intervals_end[[cross_4_plot_i]]-cross_4_qtl_intervals_start[[cross_4_plot_i]]
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, max(unlist(cross_5A_lm_list[cross_4_qtl[cross_4_plot_i,1]]))), ylab="", yaxt='n', xlab="",xaxt="n", frame.plot = FALSE)
  chrom_ticks_cm<-seq(0,max(unlist(cross_5A_lm_list[cross_4_qtl[cross_4_plot_i,1]])),10)
  axis(1, at = c(chrom_ticks_cm,max(unlist(cross_5A_lm_list[cross_4_qtl[cross_4_plot_i,1]]))), labels = c(chrom_ticks_cm,""), cex.axis=1.2)
  Plot_dists=seq(0.1,10,0.1)
  segment_colors<-rep("black", n_its)
  segment_colors[cross_4_qtl_intervals_start[[cross_4_plot_i]]<=cross_4_qtl[cross_4_plot_i,5] & cross_4_qtl_intervals_end[[cross_4_plot_i]]>=cross_4_qtl[cross_4_plot_i,4]]<-plot_color(rownames(cross_4_qtl)[cross_4_plot_i])
  segments(cross_4_qtl_intervals_start[[cross_4_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists,cross_4_qtl_intervals_end[[cross_4_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors[order(qtl_lengths, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_4_qtl[cross_4_plot_i,2], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1)
  abline(v=cross_4_qtl[cross_4_plot_i,4], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  abline(v=cross_4_qtl[cross_4_plot_i,5], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  if(cross_4_plot_i==6){
    title(xlab="genetic distance (cM)", cex.lab=1.5)
  }
  
  
  #bp
  qtl_lengths_bp<-cross_4_qtl_intervals_bp_end[[cross_4_plot_i]]-cross_4_qtl_intervals_bp_start[[cross_4_plot_i]]
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, genome_info_df$size[cross_4_qtl[cross_4_plot_i,1]]), ylab="", yaxt='n', xlab="",xaxt="n", frame.plot = FALSE)
  chrom_ticks<-seq(0,ceiling(genome_info_df$size[cross_4_qtl[cross_4_plot_i,1]]/1000000), 10)
  axis(1, at = c(chrom_ticks*1000000,genome_info_df$size[cross_4_qtl[cross_4_plot_i,1]]), labels = c(chrom_ticks,""), cex.axis=1.2)
  
  segment_colors_bp<-rep("black", n_its)
  segment_colors_bp[cross_4_qtl_intervals_bp_start[[cross_4_plot_i]]<=cross_4_qtl[cross_4_plot_i,8] & cross_4_qtl_intervals_bp_end[[cross_4_plot_i]]>=cross_4_qtl[cross_4_plot_i,7]]<-plot_color(rownames(cross_4_qtl)[cross_4_plot_i])
  segments(cross_4_qtl_intervals_bp_start[[cross_4_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists,cross_4_qtl_intervals_bp_end[[cross_4_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors_bp[order(qtl_lengths_bp, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_4_qtl[cross_4_plot_i,6], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1)
  abline(v=cross_4_qtl[cross_4_plot_i,7], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  abline(v=cross_4_qtl[cross_4_plot_i,8], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  chrom_centromeres<-most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df[cross_4_qtl[cross_4_plot_i,1],2],]
  rect(chrom_centromeres$start, -0.2,chrom_centromeres$start,-0.1, col="purple", border = "purple")
  
}
title(xlab="physical distance (Mb)", cex.lab=1.5)

dev.off()


pdf("plots/Cross4_simulations_2_2.pdf", width=8, height=12)
#par(mfrow=c(nrow(cross_2_qtl),2))
#par(mar = c(2, 2, 2, 2))

layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,13,14,15,16,16,17,18), byrow=T,ncol=2),heights=c(0.75,3,0.75,3,0.75,3,0.75,3,0.75,3,0.75,3))
sub_ids<-c("A","B","C","D","E","F")

for(cross_4_plot_i in 7:12){
  #First title
  par(mar=c(0,0.5,0,0.5))
  plot.new()
  text(0.5,0.5,paste("Chromosome ",cross_4_qtl[cross_4_plot_i,1],", ",rownames(cross_4_qtl)[cross_4_plot_i],sep=""),cex=1.5,font=1.5)
  text(-0.01, 0.5, sub_ids[cross_4_plot_i-6], cex=1.5,font=1.5)
  if(cross_4_plot_i<12){
    par(mar=c(2,0.5,0,0.5))
  }
  else{
    par(mar=c(5,0.5,0,0.5))
  }
  
  #cM
  qtl_lengths<-cross_4_qtl_intervals_end[[cross_4_plot_i]]-cross_4_qtl_intervals_start[[cross_4_plot_i]]
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, max(unlist(cross_5A_lm_list[cross_4_qtl[cross_4_plot_i,1]]))), ylab="", yaxt='n', xlab="",xaxt="n", frame.plot = FALSE)
  chrom_ticks_cm<-seq(0,max(unlist(cross_5A_lm_list[cross_4_qtl[cross_4_plot_i,1]])),10)
  axis(1, at = c(chrom_ticks_cm,max(unlist(cross_5A_lm_list[cross_4_qtl[cross_4_plot_i,1]]))), labels = c(chrom_ticks_cm,""), cex.axis=1.2)
  Plot_dists=seq(0.1,10,0.1)
  segment_colors<-rep("black", n_its)
  segment_colors[cross_4_qtl_intervals_start[[cross_4_plot_i]]<=cross_4_qtl[cross_4_plot_i,5] & cross_4_qtl_intervals_end[[cross_4_plot_i]]>=cross_4_qtl[cross_4_plot_i,4]]<-plot_color(rownames(cross_4_qtl)[cross_4_plot_i])
  segments(cross_4_qtl_intervals_start[[cross_4_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists,cross_4_qtl_intervals_end[[cross_4_plot_i]][order(qtl_lengths, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors[order(qtl_lengths, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_4_qtl[cross_4_plot_i,2], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1)
  abline(v=cross_4_qtl[cross_4_plot_i,4], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  abline(v=cross_4_qtl[cross_4_plot_i,5], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  if(cross_4_plot_i==12){
    title(xlab="genetic distance (cM)", cex.lab=1.5)
  }
  
  
  
  #bp
  qtl_lengths_bp<-cross_4_qtl_intervals_bp_end[[cross_4_plot_i]]-cross_4_qtl_intervals_bp_start[[cross_4_plot_i]]
  plot(0, cex=0, ylim=c(0,10), xlim=c(0, genome_info_df$size[cross_4_qtl[cross_4_plot_i,1]]), ylab="", yaxt='n', xlab="",xaxt="n", frame.plot = FALSE)
  chrom_ticks<-seq(0,ceiling(genome_info_df$size[cross_4_qtl[cross_4_plot_i,1]]/1000000), 10)
  axis(1, at = c(chrom_ticks*1000000,genome_info_df$size[cross_4_qtl[cross_4_plot_i,1]]), labels = c(chrom_ticks,""), cex.axis=1.2)
  segment_colors_bp<-rep("black", n_its)
  segment_colors_bp[cross_4_qtl_intervals_bp_start[[cross_4_plot_i]]<=cross_4_qtl[cross_4_plot_i,8] & cross_4_qtl_intervals_bp_end[[cross_4_plot_i]]>=cross_4_qtl[cross_4_plot_i,7]]<-plot_color(rownames(cross_4_qtl)[cross_4_plot_i])
  segments(cross_4_qtl_intervals_bp_start[[cross_4_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists,cross_4_qtl_intervals_bp_end[[cross_4_plot_i]][order(qtl_lengths_bp, decreasing = T)],Plot_dists, lwd=0.5, col=segment_colors_bp[order(qtl_lengths_bp, decreasing = T)])
  #plot empirical position and qtl location
  abline(v=cross_4_qtl[cross_4_plot_i,6], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1)
  abline(v=cross_4_qtl[cross_4_plot_i,7], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  abline(v=cross_4_qtl[cross_4_plot_i,8], col=plot_color(rownames(cross_4_qtl)[cross_4_plot_i]), lwd=1, lty=2)
  chrom_centromeres<-most_common_centromeric_repeats[most_common_centromeric_repeats$name==genome_info_df[cross_4_qtl[cross_4_plot_i,1],2],]
  rect(chrom_centromeres$start, -0.2,chrom_centromeres$start,-0.1, col="purple", border = "purple")
  
}
title(xlab="physical distance (Mb)", cex.lab=1.5)
dev.off()



