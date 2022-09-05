library(ChIPseeker)
library(ggplot2)
library(GenomicFeatures)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(intervals)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

#Sys.setenv("/ddn/LiB/softwares/Homer/bin")
annotate_peak <- function(bed_file,result_dir,txdb,prefix,tss_start=-2500,tss_end=2500,flank_Distance=50000){
	print(prefix)
	
	#step1.visualize the peak distribution on chromosomes
	covplot(bed_file, weightCol="V5")
	ggsave(filename = paste(prefix,'_coverage_plot.png',sep=''),device = "png",path = result_dir,width = 12,height = 10)
	
	#step2:annotate peaks gene	
	annotate_peaks <- annotatePeak(bed_file, tssRegion=c(tss_start, tss_end), TxDb=txdb, 
                  addFlankGeneInfo=TRUE, flankDistance=flank_Distance, annoDb="org.Mm.eg.db")
	annotate_peaks <- as.data.frame(annotate_peaks)
	save_annotate_peak_file <- paste(result_dir,'/',prefix,'_annotate_peak_gene.txt',sep='')
	write.table(annotate_peaks,save_annotate_peak_file,sep="\t",quote=FALSE)

	#step3:VISUALIZE peaks gene region
	pdf(paste(result_dir,'/',prefix,'_peak_region_bar.pdf',sep=''))
	peakAnno <- annotatePeak(bed_file, 
                         tssRegion=c(tss_start, tss_end),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
	plotAnnoBar(peakAnno)
	dev.off()

	#step4: visualize the peak distribution neighbor TSSs
	#1.heatmap
	png(paste(result_dir,'/',prefix,'_peak_region_heatmap.png',sep=''))
	peakHeatmap(bed_file, weightCol="V5", TxDb=txdb, 
            upstream=abs(tss_start), downstream=abs(tss_end), 
            color=rainbow(length(bed_file)))
	dev.off()
	
	
	#2.line
	#添加置信区间并分面
	promoter <- getPromoters(TxDb=txdb, 
                  upstream=abs(tss_start), downstream=abs(tss_end))
	tagMatrix <- getTagMatrix(bed_file, 
	                        windows=promoter)
	plotAvgProf(tagMatrix, xlim=c(tss_start, tss_end), 
		            conf=0.95,resample=500, facet="row")
	ggsave(filename = paste(prefix,'_peak_region_wave_line.png',sep=''),device = "png",path = result_dir,width = 10,height = 10)
}


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
bed_file <- "/ddn/LiB/project/chip-seq/data/analysis/examples/callpeak/SRR13292279/SRR13292279_peaks.narrowPeak"
result_dir <- "/ddn/LiB/project/chip-seq/data/analysis/examples/peaks_annotation"
tss_start=-2500
tss_end=2500
flank_Distance = 500000
annotate_peak(bed_file,result_dir,txdb,"SRR13292279",tss_start,tss_end,flank_Distance)

