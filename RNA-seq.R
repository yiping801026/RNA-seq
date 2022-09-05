install.packages('DESeq2')
BiocManager::install('DESeq2')
BiocManager::install('clusterProfiler')
BiocManager::install('stringr')
BiocManager::install('enrichplot')
BiocManager::install('biomaRt')
BiocManager::install('KEGGREST')
install.packages("Cairo")
install.packages('devtools')
install.packages('ggvenn')
install.packages("ggVennDiagram")

library(edgeR)
library(DESeq2)
library(ggrepel)
library(clusterProfiler)
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(KEGGREST)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Cairo)
library(dplyr)
library(usethis)
library(devtools)
library(ggvenn)
library(ggVennDiagram)

#指定富集分析的物种库
GO_database <- 'org.Mm.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
#GO_database <- 'org.Rn.eg.db'
KEGG_database <- 'mmu' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#####
##TPM
input_data_raw <- read.table("/Users/pingyi/Desktop/GWAS/RNA_seq_basic/final_featureCounts.TPM.txt",header = TRUE,row.names = 1)
colnames(input_data_raw)
input_data_raw <- input_data_raw[,6:17]



input_data_raw <- read.table("/Users/pingyi/Desktop/GWAS/RNA_seq_basic/count_matrix_name.txt",header = TRUE,row.names = 1)
#info_data_raw <- read.table('/Users/pingyi/Desktop/GWAS/RNA_seq_basic/sample_info.txt',header = TRUE,row.names = 1)
### use PT01 as case 
info_data_raw <- read.table('/Users/pingyi/Desktop/GWAS/RNA_seq_basic/sample_info_PT01.txt',header = TRUE,row.names = 1)

###### 2 vs 2
##R848 7,8
input_data <- input_data_raw[,c(2,3,10,11)]
info_data <- info_data_raw[c(2,3,10,11),]
a <- Run_DESeq2(info_data,input_data,'R848')
##ADU 9,10
input_data <- input_data_raw[,c(2,3,12,1)]
info_data <- info_data_raw[c(2,3,12,1),]
b <- Run_DESeq2(info_data,input_data,'ADU')
##ADP 3,4
input_data <- input_data_raw[,c(2,3,6,7)]
info_data <- info_data_raw[c(2,3,6,7),]
c <- Run_DESeq2(info_data,input_data,'ADP')
##PT01_0010
input_data <- input_data_raw[,c(2,3,4,5)]
info_data <- info_data_raw[c(2,3,4,5),]
d <- Run_DESeq2(info_data,input_data,'PT01_0010')
##PT01_0079
input_data <- input_data_raw[,c(2,3,8,9)]
info_data <- info_data_raw[c(2,3,8,9),]
e <- Run_DESeq2(info_data,input_data,'PT01_0079')

###################################################
###################################################

##PT01_0010 vs R848
input_data <- input_data_raw[,c(4,5,10,11)]
info_data <- info_data_raw[c(4,5,10,11),]
f <- Run_DESeq2(info_data, input_data, med_name = 'PT01_0010', control_name = 'R848')

##PT01_0010 vs ADU
input_data <- input_data_raw[,c(4,5,1,12)]
info_data <- info_data_raw[c(4,5,1,12),]
g <- Run_DESeq2(info_data, input_data, med_name = 'PT01_0010', control_name = 'ADU')

##PT01_0010 vs ADP
input_data <- input_data_raw[,c(4,5,6,7)]
info_data <- info_data_raw[c(4,5,6,7),]
h <- Run_DESeq2(info_data, input_data, med_name = 'PT01_0010', control_name = 'ADP')

##PT01_0079 vs R848
input_data <- input_data_raw[,c(8,9,10,11)]
info_data <- info_data_raw[c(8,9,10,11),]
i <- Run_DESeq2(info_data, input_data, med_name = 'PT01_0079', control_name = 'R848')

##PT01_0079 vs ADU
input_data <- input_data_raw[,c(8,9,1,12)]
info_data <- info_data_raw[c(8,9,1,12),]
j <- Run_DESeq2(info_data, input_data, med_name = 'PT01_0079', control_name = 'ADU')

##PT01_0079 vs ADP
input_data <- input_data_raw[,c(8,9,6,7)]
info_data <- info_data_raw[c(8,9,6,7),]
k <- Run_DESeq2(info_data, input_data, med_name = 'PT01_0079', control_name = 'ADP')




#####remove smaller than 10
#input_data <- input_data [which(rowSums(input_data) > 9),]
#input_data <- round(input_data,digits=0)



#input_data <- as.matrix(input_data)
#input_data
#write.table(input_data,'input_data.csv')
#info_data <- as.matrix(info_data)
#info_data

############################################################## get data prepared and write table
#### DESeq2 use the TPM inside (no need to do any other steps)
Run_DESeq2<- function(info_data,input_data,med_name,control_name){
  #med_name = 'PT01_0010'
  #control_name = 'R848'
  ### DESeq2 if use control as control
  info_data$label <- factor(info_data$label ,levels =c('CON','CASE') )
  dds <- DESeqDataSetFromMatrix(countData = input_data,colData = info_data,design = ~label)
  dds <- DESeq(dds)
  result <- results(dds,alpha = 0.05)
  result <- result[order(-abs(result$log2FoldChange)),]
  ### get differencial expressed genes
  diff_gene <- subset(result, pvalue < 0.05 & (log2FoldChange > 0 | log2FoldChange < 0))
  ### get heatmap_plot
  heatmap_plot(dds,diff_gene,med_name,control_name)
  ### write table into file
  result_data <- merge(as.data.frame(result),as.data.frame(counts(dds,normalized=TRUE)),by='row.names',sort=FALSE)
  names(result_data)[1] <- 'Gene'
  write.table(result_data,file= paste('gene_result_file_',med_name,'_vs_',control_name,'.txt',sep = ''),sep = '\t',quote=F,row.names = F)
  ### get VOL plot
  VOL_plot(result,med_name,result_data,control_name)
  ### get GO enrichment plot 
  GO_run(result,med_name,control_name)
  
  
  ### 也可以输出差异基因
  #write.table(x=as.data.frame(diff_gene), file="results_gene_annotated_significant_ADP.txt", sep="\t", quote=F, col.names=NA)
  return(result)
}

############################# heatmap function 

heatmap_plot <- function(dds,diff_gene,med_name,control_name){
  dds_rlog <- rlog(dds, blind=FALSE)
  colData(dds_rlog)$sampleid <- rownames(colData(dds_rlog))
  # select top 40 genes to make the heatmap plot
  mat <- assay(dds_rlog[row.names(diff_gene)])[1:40, ]
  # choose the coloumns for heatmap 
  annotation_col <- data.frame(
    Group=factor(colData(dds_rlog)$label),
    medcine=factor(colData(dds_rlog)$med),
    row.names=colData(dds_rlog)$sampleid
  )
  # change heatmap 'med' colors
  ann_colors <- list(
    Group=c(CASE="darkorange", CON="lightblue"),
    if(med_name == 'R848'){ medcine=c(NC="forestgreen",R848 = 'yellow')}
    else if (med_name == 'PT01') { medcine=c(NC="forestgreen",PT01='blue')}
    else if (med_name == 'ADU') {medcine=c(NC="forestgreen",ADU="darkred")}
    else {medcine=c(NC="forestgreen",ADP = 'red')}
  )
  # make heatmap plot
  pheatmap(mat=mat, 
           color=colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
           scale="row", # Scale genes to Z-score (how many standard deviations)
           annotation_col=annotation_col, # Add multiple annotations to the samples
           annotation_colors=ann_colors,# Change the default colors of the annotations
           fontsize=6.5,
           filename = paste(med_name,'_vs_',control_name,'_Heatmap_top40genes.pdf',sep = '')) # Make fonts smaller
}


plotPCA(dds_rlog, intgroup='label', ntop=500) +
  #theme_bw() + # 修改主体
  #geom_point(size=5) + # 增加点大小
  #scale_y_continuous(limits=c(-5, 5)) +
  ggtitle(label="Principal Component Analysis (PCA)",
          subtitle="Top 500 most variable genes")

###################################

#### write table 
vol_data_out <- data.frame(gene=row.names(result), pval=result$pvalue, log2FoldChange=result$log2FoldChange)
vol_data_out <- mutate(vol_data_out, regulate=case_when(
  vol_data_out$log2FoldChange > 0 & vol_data_out$pval < 0.05 ~ "up-regulated",
  vol_data_out$log2FoldChange < 0 & vol_data_out$pval < 0.05 ~ "down-regulated",
  vol_data_out$pval > 0.05 ~ "nonsignificant"))
write.table(vol_data_out,file='gene_file_R848.txt',sep = '\t',quote=F,col.names = T,row.names = F)


############################# VOL plot function
VOL_plot <- function(result,med_name,result_data,control_name){
  ### make vol_data as input for plot making use padj
  #vol_data <- data.frame(gene=row.names(result), pval=-log10(result$padj), lfc=result$log2FoldChange)
  ### make vol_data as input for plot making use p-value
  vol_data <- data.frame(gene=row.names(result), pval=-log10(result$pvalue), lfc=result$log2FoldChange)
  # remove zreo
  vol_data <- na.omit(vol_data)
  # set colors for up_regulated and down_regulated
  vol_data <- mutate(vol_data, color=case_when(
    vol_data$lfc >= 1 & vol_data$pval > -log10(0.05) ~ "UP",
    vol_data$lfc <= -1 & vol_data$pval > -log10(0.05) ~ "DOWN",
    vol_data$lfc > -1 & vol_data$lfc < 1 ~ "nonsignificant",
    vol_data$pval <= -log10(0.05) ~ "nonsignificant"))
  # select genes we interested + top 40 genes
  diff_gene <- subset(vol_data, vol_data$pval > -log10(0.05) )
  top_gene <- diff_gene[1:40,]$gene
  vol_data$label = ifelse(vol_data$gene%in%top_gene,
                          #| vol_data$gene %in% c('Alpk1','Tifa','Traf6'), 
                          as.character(vol_data$gene),"")
  ### make VOL plot
  vol <- ggplot(vol_data, 
                aes(x=lfc, y=pval, color=color , label = vol_data$label ) )
  # change bands for VOL plot
  vol + ggtitle(label=paste(med_name,'_vs_',control_name," Volcano Plot",sep = ''), subtitle="Colored by fold-change direction") +
    geom_point(size=2.5, alpha=0.8, na.rm=T) +
    scale_color_manual(name="Regulated",
                       values=c(DOWN="#008B00",UP="#CD4F39", 
                                nonsignificant="darkgray")) +
    theme_bw(base_size=14) +
    theme(legend.position="right") +
    #xlab(expression(log[2]("Case" / "Control"))) +
    xlab(label = paste('log2(',med_name,'/',control_name,')',sep = '')) +
    ylab(expression(-log[10]("p-value"))) +
    geom_hline(yintercept= -log10(0.05), colour="darkgrey",show.legend = TRUE) +
    geom_vline(xintercept=1, colour="darkgrey",show.legend = TRUE) +
    geom_vline(xintercept=-1, colour="darkgrey",show.legend = TRUE) +
    scale_y_continuous(trans="log1p",breaks = c(0,1.3,100,200,300,400))+
    scale_x_continuous(breaks = c(-5,-1,0,1,5,10))+
    geom_text_repel(data = vol_data, aes(x = lfc, y = pval, label = label),size = 4)+
    ggsave(paste(med_name,'_vs_',control_name,'_VOL.pdf'),width = 20,height =10 )
}



############################### GO_all_in_one
GO_all_in_one <- function(gene_all,regulate_name,med_name,control_name){
  gene_all_go <-as.character(gene_all$gene)
  g_list <- gene_all$log2FoldChange
  names(g_list) <- gene_all$gene
  GO_en<-enrichGO( gene_all_go, #GO富集分析
                   OrgDb = GO_database,
                   keyType = "SYMBOL",# using the gene name as input
                   ont = "ALL",# including BP,CC,MF (Biological Process,Cellular Component,Mollecular Function)
                   pvalueCutoff = 0.05,# p cutoff
                   qvalueCutoff = 0.05,# q cutoff
                   readable = F)
  ### make matrix & output
  GO_en@result <- GO_en@result[order(GO_en@result$ONTOLOGY,-GO_en@result$Count),]
  GO_out <- as.data.frame(GO_en)
  write.table(GO_out,paste(med_name,'_vs_',control_name,regulate_name,'_GO_out',sep = ''))
  ### make dotplot for GO
  dotplot(GO_en, split="ONTOLOGY") + 
    facet_grid(ONTOLOGY~., scale="free")+
    ggtitle(label=paste(med_name,'_vs_',control_name," DotPlot",sep = ''))+ 
              ggsave(paste(med_name,'_vs_',control_name,regulate_name,'_dotplot.pdf',sep = ''),width = 10,height =10 )
  ### make net plot for GO 
  cnetplot(GO_en, foldChange=g_list,categorySize='geneNum',showCategory = 5, colorEdge = TRUE,size=15) + 
    ggtitle(label=paste(med_name,'_vs_',control_name," NetPlot",sep = ''))+ 
    ggsave(paste(med_name,'_vs_',control_name,regulate_name,'_net.pdf',sep = ''),width = 20,height =20 )
  ### make circle net plot for GO
  #cnetplot(GO_en, foldChange=g_list,circular = TRUE, colorEdge = TRUE) + ggsave(paste(med_name,'_',regulate_name,'_cir.pdf',sep = ''),width = 20,height =20 )
  #dotplot(GO_down, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")+ggsave('t.pdf',width = 10,height =10 )
  #cnetplot(GO_down, foldChange=g_down)+ggsave('t2.pdf',width = 20,height =20 )
  #heatplot(GO_down, foldChange=g_down)+ggsave('t3.pdf',width = 20,height =20 )
  
}



############################### GO_run function 

GO_run <- function(result,med_name,control_name){
  ### prepare data for GO
  #info <- data.frame(gene=row.names(result), pval=result$padj, log2FoldChange=result$log2FoldChange)
  info <- data.frame(gene=row.names(result), pval=result$pvalue, log2FoldChange=result$log2FoldChange)
  info <- mutate(info, regulated=case_when(
    info$log2FoldChange > 1 & info$pval < 1e-5 ~ "up-regulated",
    info$log2FoldChange < -1 & info$pval < 1e-5 ~ "down-regulated",
    info$pval > 1e-5 ~ "nonsignificant"))
  ### GO_up 
  gene_up_all <- filter(info,regulated == 'up-regulated')
  GO_all_in_one(gene_up_all,'up',med_name,control_name )
  ### GO_down 
  gene_down_all <- filter(info,regulated == 'down-regulated')
  GO_all_in_one(gene_down_all, 'down',med_name ,control_name)

}

################################### Venn_plot function 

a1 <-  data.frame(a)
a_up <- filter(a1,log2FoldChange > 1 & a1$pvalue< 1e-5)
a_down <- filter(a1,log2FoldChange < -1 & a1$pvalue< 1e-5)
b1 <-  data.frame(b)
b_up <- filter(b1,log2FoldChange > 1 & b1$pvalue< 1e-5)
b_down <- filter(b1,log2FoldChange < -1 & b1$pvalue< 1e-5)
c1 <-  data.frame(c)
c_up <- filter(c1,log2FoldChange > 1 & c1$pvalue< 1e-5)
c_down <- filter(c1,log2FoldChange < -1 & c1$pvalue< 1e-5)
d1 <-  data.frame(d)
d_up <- filter(d1,log2FoldChange > 1 & d1$pvalue< 1e-5)
d_down <- filter(d1,log2FoldChange < -1 & d1$pvalue< 1e-5)
e1 <-  data.frame(e)
e_up <- filter(e1,log2FoldChange > 1 & e1$pvalue< 1e-5)
e_down <- filter(e1,log2FoldChange < -1 & e1$pvalue< 1e-5)

venn_list_up <- list(
  rownames(a_up),
  rownames(b_up),
  rownames(c_up),
  rownames(d_up),
  rownames(e_up)
) 
names(venn_list_up) <- c('R848 (485 genes)','ADU (988 genes)','ADP (358 genes)','PT01_0010 (378 genes)','PT01_0079 (360 genes)')

dev.new()

ggVennDiagram(venn_list_up[1:5],label_alpha=0,label_color = "white",
              label_size = 6,
              set_size = 9 ) + scale_fill_gradient(low="blue",high = "red")+ggsave('vv.pdf',width = 30,height =20 )



















gene_down <- as.character(gene_down_all$gene)
g_down <- gene_down_all$log2FoldChange
names(g_down) <- gene_down_all$gene
#de <- names(g_down)[abs(geneList) > 2]
GO_down<-enrichGO( gene_down,#GO富集分析
                   OrgDb = GO_database,
                   keyType = "SYMBOL",#设定读取的gene ID类型
                   ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                   pvalueCutoff = 0.05,#设定p值阈值
                   qvalueCutoff = 0.05,#设定q值阈值
                   readable = F)
a <- as.data.frame(GO_down)
a[order(a,by)]

GO_down@result <- GO_down@result[order(GO_down@result$ONTOLOGY,-GO_down@result$Count),]


cnetplot(GO_down, foldChange=g_down,colorEdge = TRUE,label_size = 20,set_size = 20) + theme(legend.text=element_text(size=10)) 
dotplot(GO_down, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

mtcars[order(mtcars$cyl,mtcars$disp),] 
a1 <- a[order(a$ONTOLOGY,-a$Count),]








GO_try <- enrichGO(gene_down,OrgDb = GO_database,
                   keyType = "SYMBOL",#设定读取的gene ID类型
                   ont = "BP",
                   readable = F
                   )
GO_try2 <- simplify(GO_try)
data(geneList, package="DOSE")
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
## remove redundent GO terms
ego2 <- simplify(ego)
cnetplot(ego2, foldChange=geneList,categorySize='geneNum',)+ggsave('t2.pdf',width = 10,height =10 )
cnetplot(ego, foldChange=geneList)+ggsave('t.pdf',width = 10,height =10 )


#+ggsave('t.pdf',width = 10,height =10 )

info2 <- info[1:40,]
info3 <- info2[,c(1,3)]
ego2 <- as.character(info3$gene)
geneList <- info$log2FoldChange



#enrichKK = DOSE::setReadable(enrichKK,OrgDb = 'org.Mm.eg.db',keyType = 'ENTREZID')
#enrichKK


barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,title = 'KEGG Pathway')

dotplot(KEGG)
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
GO2 <- pairwise_termsim(GO)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图


write.table(KEGG$ID, file = "/Users/ZYP/Downloads/KEGG_GO/KEGG_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


################################ try this TPM 
################################
data <- read.delim("~/Desktop/GWAS/GWAS/matrix2.txt", comment.char="#", row.names=1)
countdata <- data[ , 6:ncol(data)]
## remove the genes which have fewer than 10 hits  
countdata <- countdata [which(rowSums(countdata) > 9),]
name <- rownames(countdata)
data <- data[which(rownames(data)%in%name),]
## calculate TPM
KB <- data$Length / 1000
RPK <- countdata / KB
TPM <- t(t(RPK) / colSums(RPK) * 1000000)
TPM <- merge(data[,1:5], as.data.frame(TPM), by="row.names", sort=FALSE)
#TPM <- TPM[, 1:ncol(TPM)]
write.table(TPM, "final_featureCounts.TPM.txt", sep="\t", quote=FALSE, row.names=FALSE)

##
info_data <- read.table('/Users/pingyi/Desktop/GWAS/RNA_seq_basic/sample_info.txt',header = TRUE)
colnames(info_data) <- c('sampleid','med','Group')
rownames(info_data) <- info_data$sampleid
dds@metadata <- info_data

#############################
#############################







