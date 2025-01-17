require(adegenet)
require(pegas)
require(vcfR)
args = commandArgs(trailingOnly=TRUE)

cross_name<-args[1]
#cross_name<-"3A"
vcf_input<-args[2]
#vcf_input<-"~/Desktop/qtl_paper/vcfs/3A-round1_low_cov_filtered.recode.vcf"
pca_output_pdf<-args[3]
#pca_output_pdf<-"~/Desktop/qtl_paper/test/3A_test_pca.pdf"
pca_output_titles_pdf<-args[4]
#pca_output_titles_pdf<-"~/Desktop/qtl_paper/test/3A_test_pca_titles.pdf"
pca_output_table<-args[5]
#pca_output_table<-"~/Desktop/qtl_paper/test/3A_test_pca_table.tsv"

A1_vcfR<-read.vcfR(vcf_input)
A1_vcfR_genind<-vcfR2genind(A1_vcfR)
pca_obj<-scaleGen(A1_vcfR_genind, NA.method="mean")
pca1 <- dudi.pca(pca_obj,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
pca1_variation_explained_axis1<-pca1$eig[1]/sum(pca1$eig)
pca1_variation_explained_axis2<-pca1$eig[2]/sum(pca1$eig)
plot_data_pca1<-data.frame(samples=rownames(pca1$li),x1=pca1$li$Axis1, x2=pca1$li$Axis2, colors=rep("black", length(pca1$li$Axis1)), symbols=rep(19, length(pca1$li$Axis1)), stringsAsFactors = FALSE)

#make plot of PCA
pdf(pca_output_pdf, width=7, height=5)
plot(plot_data_pca1$x1, plot_data_pca1$x2, xlab=paste("PC1 [",round(pca1_variation_explained_axis1,2),"]", sep=""), ylab=paste("PC2 [",round(pca1_variation_explained_axis2,2),"]", sep=""), col=plot_data_pca1$colors, pch=plot_data_pca1$symbols)
title(cross_name)
dev.off()

#make plot with sample names
pdf(pca_output_titles_pdf, width=7, height=5)
plot(plot_data_pca1$x1, plot_data_pca1$x2, xlab=paste("PC1 [",round(pca1_variation_explained_axis1,2),"]", sep=""), ylab=paste("PC2 [",round(pca1_variation_explained_axis2,2),"]", sep=""), cex=0)
text(plot_data_pca1$x1, plot_data_pca1$x2, labels = plot_data_pca1$samples)
title(cross_name)
dev.off()

#write table with PCA values
output_table<-data.frame(samples=plot_data_pca1$samples, PC1=plot_data_pca1$x1, PC2=plot_data_pca1$x2)
write.table(output_table, pca_output_table, sep="\t", quote=FALSE, row.names=FALSE)

