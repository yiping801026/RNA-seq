
######## install library


######## load library
library(Rsubread)
library(edgeR)
suppressMessages(library(org.Mm.eg.db))
library("AnnotationDbi")

########

#rownames(count_matrix$counts) <- mapIds(org.Mm.eg.db, keys = rownames(count_matrix$counts), keytype = "ENSEMBL", column="SYMBOL")
	
	# symbol_org <- keys(org.Mm.eg.db, keytype = "SYMBOL")
	
	# save_matrix_file <- paste(save_dir,'/gene_sample_count_matrix_symbol.txt',sep='')
	# write.table(count_matrix$counts,save_matrix_file,sep="\t",quote=FALSE)

calculate_counts <- function(sample_info_file,mapping_dir,gtf_file,save_dir,save_symbol_file){
	gene_symbol_matrix <- read.delim(save_symbol_file,header=T,row.names=1)

	bam_file_list <- c()
	sample_info <- read.delim(sample_info_file,header=FALSE)
	sample_names <- as.character(sample_info[,1])
	sample_id_list <- c()
	for(i in 1:length(sample_names)){
		sample_id <- paste(unlist(strsplit(sample_names[i],split='_'))[1:3],collapse='_')
		c_bam_file <- paste(mapping_dir,'/',sample_id,'/',sample_id,'Aligned.sortedByCoord.out.bam',sep='')	
		bam_file_list <- c(bam_file_list,c_bam_file)
		sample_id_list <- c(sample_id_list,sample_id)
	}
	
	# count_matrix <- featureCounts(files=bam_file_list,annot.ext=gtf_file,
	# 	isPairedEnd=FALSE,nthreads=10,isGTFAnnotationFile=TRUE)	
	# colnames(count_matrix$counts) <- sample_id_list
	save_matrix_file <- paste(save_dir,'/gene_sample_count_matrix.txt',sep='')
	#write.table(count_matrix$counts,save_matrix_file,sep="\t",quote=FALSE)

	save_matrix_annot_file <- paste(save_dir,'/gene_sample_annotation_matrix.txt',sep='')
	#rownames(count_matrix$annotation) <- count_matrix$annotation$GeneID
	#write.table(count_matrix$annotation,save_matrix_annot_file,sep="\t",quote=FALSE,row.names=FALSE)

	count_matrix_counts <- read.delim(save_matrix_file,header=T,row.names=1)

	#filter low-expression genesin all samples
	print(paste(nrow(count_matrix_counts),'genes'))
	count_gene_exp <- count_matrix_counts[rowMax(as.matrix(count_matrix_counts))>=10,]
	save_filter_matrix_file <- paste(save_dir,'/gene_sample_count_matrix_min10.txt',sep='')
	write.table(count_gene_exp,save_filter_matrix_file,sep="\t",quote=FALSE)
	print(paste('after removing low-expression genes',nrow(count_gene_exp),'genes'))

	#rownames(count_matrix_counts) <- as.character(gene_symbol_matrix[rownames(count_matrix_counts),'gene_symbol'])
	count_matrix_length <- read.delim(save_matrix_annot_file,header=T,row.names=1)
	count_matrix_length_filter <- count_matrix_length[rownames(count_gene_exp),]
	count_matrix_length_filter$symbol <- as.character(gene_symbol_matrix[rownames(count_gene_exp),]$gene_symbol)
	save_filter_matrix_annot_file <- paste(save_dir,'/gene_sample_annotation_matrix_min10.txt',sep='')
	#colnames(count_matrix_length_filter) <- c("symbol","Chr","Start","End","Strand","Length")

	#count_matrix_length_filter$Length <- as.numeric(as.character(count_matrix_length_filter$Length))

	write.table(count_matrix_length_filter,save_filter_matrix_annot_file,sep="\t",quote=FALSE)
	group_list <- apply(as.matrix(colnames(count_gene_exp)),1,function(x){
		return(paste(unlist(strsplit(x[1],split='_'))[1:2],collapse='_'))
	})	
	value_list <- calculate_rpkm(as.matrix(count_gene_exp),group_list,count_matrix_length_filter)	
	save_rpkm_matrix_file <- paste(save_dir,'/gene_sample_rpkm_matrix_min10.txt',sep='')
	write.table(value_list[['rpkm']],save_rpkm_matrix_file,sep="\t",quote=FALSE)
	save_cpm_matrix_file <- paste(save_dir,'/gene_sample_cpm_matrix_min10.txt',sep='')
	write.table(value_list[['cpm']],save_cpm_matrix_file,sep="\t",quote=FALSE)
	
}

calculate_rpkm <- function(counts,group_list,genes_annot){
	y <- DGEList(counts=counts,group=group_list,genes=genes_annot)
	y <- calcNormFactors(y)
	RPKM <- rpkm(y)
	CPM <- cpm(y)
	return(list(rpkm=RPKM,cpm=CPM))
}


get_gene_symbol_map <- function(gtf_file,save_file) {
  if (is.character(gtf_file)) {
    if(!file.exists(gtf_file)) stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(gtf_file, header = FALSE)
  } else {
    data.table::setDT(gtf_file)
  }
  
  input = input[input[[3]] == "transcript", ]
  
  pattern_id = ".*gene_id \"(ENSM[0-9]+)\";.*"
  pattern_name = ".*gene_id \"([^;]+)\";.*"
  pattern_symbol = ".*gene_name \"([^;]+)\";.*"
  
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  gene_symbol = sub(pattern_symbol, "\\1", input[[9]])
  gene_start <- input[[4]]
  gene_end <- input[[5]]
  gene_strand <- as.character(input[[7]])
  gene_chr <- apply(as.matrix(as.character(input[[1]])),1,function(x){
  	if(as.character(x[1]) %in% as.character(1:23)){
  		return(paste('chr',as.character(x[1]),sep=''))
  	}else{
  		return(as.character(x[1]))
  	}
  })
  gene_table <- data.frame(
             gene_name = gene_name,
             gene_id = gene_id,
             gene_symbol = gene_symbol,
             gene_start = gene_start,
             gene_end = gene_end,
             gene_strand = gene_strand,
             gene_chr = gene_chr,
             stringsAsFactors = FALSE)

 write.table(gene_table,file=save_file,quote=F,row.names=FALSE,sep="\t")
}

main <- function(){
	sample_info_file = "/ddn/LiB/pingyi/RNA_seq/RNA_lists.txt"
    sample_dir = "/ddn/LiB/pingyi/RNA_seq/data_RNA_seq"
    result_dir = "/ddn/LiB/pingyi/RNA_seq/results"
    mouse_gtf_file = "/ddn/LiB/data/ref_genome/mouse/genome/Mus_musculus.GRCm38.102.gtf"
    mapping_dir <- paste(result_dir,'/STAR',sep='')
    save_count_dir <- paste(result_dir,'/gene_count',sep='')
    dir.create(save_count_dir)
    save_symbol_file <- "/ddn/LiB/data/ref_genome/mouse/genome/mm10_knowngene_ENSEMBL_TO_SYMBOL.txt"
	
	get_gene_symbol_map(mouse_gtf_file,save_symbol_file)

	#calculate_counts(sample_info_file,mapping_dir,mouse_gtf_file,save_count_dir,save_symbol_file)
	


}
main()