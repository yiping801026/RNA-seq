import os

def fastqc(script,input_fastq,outdir):
	script = "/ddn/LiB/softwares/fastqc_v0.11.9/FastQC/fastqc"
	cmd = "%s -o %s -t 20 %s"%(script,outdir,input_fastq)
	os.system(cmd)

def multiqc(fastqc_dir,outdir,out_name):
	cmd = "multiqc %s -o %s -n %s"%(fastqc_dir,outdir,out_name)
	os.system(cmd)

def trimmomatic(script,input_fastq,output_fastq):
	cmd = "java -jar %s SE -phred33 %s %s \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"%(script,input_fastq,output_fastq)
	os.system(cmd)

def trim_galore(input_fastq,outdir,pair_flag,adapter_seq):
	script = "/ddn/LiB/softwares/TrimGalore-0.6.6/trim_galore"
	if adapter_seq!='':
		if pair_flag=='paired':
			cmd = "%s -q 25 --phred33 -e 0.1 --stringency 4 --adapter %s --paired --retain_unpaired --cores 20 --gzip --fastqc -o %s %s"%(script,adapter_seq,outdir,input_fastq)
		else:
			cmd = "%s -q 25 --phred33 -e 0.1 --stringency 4 --adapter %s --gzip --fastqc -o %s %s"%(script,adapter_seq,outdir,input_fastq)
		os.system(cmd)
	else:
		if pair_flag=='paired':
			cmd = "%s -q 25 --phred33 -e 0.1 --stringency 4 --paired --retain_unpaired --cores 20 --gzip --fastqc -o %s %s"%(script,outdir,input_fastq)
		else:
			cmd = "%s -q 25 --phred33 -e 0.1 --stringency 4 --gzip --fastqc -o %s %s"%(script,outdir,input_fastq)
		os.system(cmd)

def run_fastqc(sample_info_file,sample_dir,fastqc_dir,prefix):
	# with open(sample_info_file) as f:
	#     contents = f.readlines()
	# for line in contents[1:]:
	#     line = line.strsip().split('\t')
	#     sample_name = line[0]
	#     c_sample_file = os.path.join(sample_dir.sample_name)
	fastqc_script = "/ddn/LiB/softwares/fastqc_v0.11.9/FastQC/fastqc"
	input_fastq = os.path.join(sample_dir,'*.fastq*')
	fastqc(fastqc_script,input_fastq,fastqc_dir)
	multiqc_dir = os.path.join(fastqc_dir,'multiqc')
	mkdir(multiqc_dir)
	multiqc(fastqc_dir,multiqc_dir,prefix)

def get_sample_info(sample_info_file,save_sample_info_file):
	f_save = open(save_sample_info_file,'w')
	f_save.write('sample_name\tsample_id\tsample_type\tphenotype\n')
	with open(sample_info_file) as f:
		contents = f.readlines()
	for line in contents:
		sample_name = line.strip().split('\t')[0]
		if 'CHIP' in sample_name:
			sample_id = '_'.join(sample_name.split('_')[1:3])
			sample_type = sample_name.split('_')[1]
			phenotype = sample_name.split('_')[2]
		else:
			sample_id = sample_name.split('_')[0]
			sample_type = 'input'
			phenotype = sample_name.split('_')[0]
		f_save.write(sample_name+'\t'+sample_id+'\t'+sample_type+'\t'+phenotype+'\n')
		f_save.flush()
	f_save.close()

def mkdir(dirname):
	cmd = "mkdir -p %s"%(dirname)
	os.system(cmd)

def bowtie2_build(fasta_file,index_file):
	cmd = "bowtie2-build %s %s"%(fasta_file,index_file)
	os.system(cmd)

def sam2bam(sam_file,bam_file):
    cmd_convert = 'samtools view -@ 25 -bS %s > %s'%(sam_file,bam_file)
    os.system(cmd_convert)

def bowtie2_mapping(read1_file,read2_file,sam_file,ref_db,bowtie2_log_file,pair_flag):
	if pair_flag =='single':
		cmd_bowtie2 = 'bowtie2 -q -p 20 -x %s -U %s -S %s > %s 2>&1'%(ref_db,read1_file,sam_file,bowtie2_log_file)
		print(cmd_bowtie2)
		os.system(cmd_bowtie2)
	else:
		cmd_bowtie2 = 'bowtie2 -q -p 20 -x %s -1 %s -2 %s -S %s > %s 2>&1'%(ref_db,read1_file,read2_file,sam_file,bowtie2_log_file)
		print(cmd_bowtie2)
		os.system(cmd_bowtie2)
	sam2bam(sam_file,sam_file+'.bam')
	samtools_sort_index(sam_file+'.bam',sam_file+'.bam'+'_sort')
	get_mapped_bam_file(sam_file+'.bam'+'_sort',sam_file+'.bam'+'_sort_mapped')

def samtools_sort_index(bam_file,sort_bam_file):
	command = "samtools sort -@ 20 %s -o %s"%(bam_file,sort_bam_file)
	os.system(command)
	#print(command)
	command = "samtools index %s"%(sort_bam_file)
	os.system(command)
	print(command)

def cutadapt(output_fastq,input_fastq,adapter_seq):
	if adapter_seq=='':
		cmd = "cutadapt --cores=20 -o %s %s"%(output_fastq,input_fastq)
	else:
		cmd = "cutadapt --cores=20 -a %s -o %s %s"%(adapter_seq,output_fastq,input_fastq)
	os.system(cmd)

def run_cutadapt(save_sample_info_file,sample_dir,cutadapt_dir,default_falg):

	with open(save_sample_info_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		sample_name = line[header.index('sample_name')]
		sample_id = line[header.index('sample_id')]
		c_fastq_file = os.path.join(sample_dir,sample_name)
		
		c_outdir = os.path.join(cutadapt_dir,'trim_galore',sample_id)
		rm_adapter_fastq = os.path.join(cutadapt_dir,sample_name)
		mkdir(c_outdir)
		trim_galore(c_fastq_file,c_outdir,'single','')		
		c_trim_file = os.path.join(cutadapt_dir,'trim_galore',sample_id,sample_name.strip('.fastq.gz')+'_trimmed.fq.gz')
		os.rename(c_trim_file,rm_adapter_fastq)

def run_bowtie2(sample_info_file,sample_dir,bowtie2_dir,ref_db,thread_num):
	with open(sample_info_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		sample_name = line[0]
		sample_id = line[header.index('sample_id')]
		c_sample_file = os.path.join(sample_dir,sample_name)
		save_dir = os.path.join(bowtie2_dir,sample_id)
		mkdir(save_dir)
		save_prefix = os.path.join(save_dir,sample_id)
		sam_file = os.path.join(save_dir,sample_id+'_bowtie2_mm10.sam')
		bowtie2_log_file = os.path.join(save_dir,sample_id+'_bowtie2_mm10.log')
		bowtie2_mapping(c_sample_file,"",sam_file,ref_db,bowtie2_log_file,'single')	

def macs(query_bam_file,control_bam_file,outdir,prefix,species='mouse'):
	if control_bam_file!='na':
		if species=='mouse':
			cmd = "macs3 callpeak -t %s -c %s -f BAM -g mm -n %s -B --outdir %s"%(query_bam_file,control_bam_file,prefix,outdir)
			os.system(cmd)
		if species == 'human':
			cmd = "macs3 callpeak -t %s -c %s -f BAM -g hs -n %s -B --outdir %s"%(query_bam_file,control_bam_file,prefix,outdir)
			os.system(cmd)
	else:
		if species=='mouse':
			cmd = "macs3 callpeak -t %s -f BAM -g mm -n %s -B --outdir %s"%(query_bam_file,prefix,outdir)
			os.system(cmd)
		if species == 'human':
			cmd = "macs3 callpeak -t %s -f BAM -g hs -n %s -B --outdir %s"%(query_bam_file,prefix,outdir)
			os.system(cmd)

def sambamba(bam_file,unique_bam_file):
	cmd = "sambamba view -h -t 10 -f bam -F '[XS] == null and not unmapped  and not duplicate' %s > %s"%(bam_file,unique_bam_file)
	os.system(cmd)

def run_callpeak(sample_info_file,bowtie2_dir,callpeak_dir,control_bam_file="na"):
	with open(sample_info_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		sample_name = line[header.index('sample_name')]
		sample_id = line[header.index('sample_id')]
		sample_type = line[header.index('sample_type')]
		c_bowtie2_dir = os.path.join(bowtie2_dir,sample_id)
		save_prefix = os.path.join(c_bowtie2_dir,sample_id)
		bam_file = os.path.join(c_bowtie2_dir,sample_id+'_bowtie2_mm10.sam.bam_sort')
		c_callpeak_dir = os.path.join(callpeak_dir,sample_id)
		mkdir(c_callpeak_dir)
		if sample_type != 'input':
			macs(bam_file,control_bam_file,c_callpeak_dir,sample_id,'mouse')

def bamcoverage(bam_file,bigwig_file):
	cmd = "bamCoverage --bam %s -o %s \
	    --binSize 20 \
	    --normalizeUsing CPM \
	    --smoothLength 60 \
		--extendReads 150 \
		--centerReads"%(bam_file,bigwig_file)
	os.system(cmd)
 
def bamCompare(bam_file,input_bam_file,bigwig_file):
	cmd = "bamCompare -b1 %s -b2 %s -o %s \
	    --binSize 20 \
	    --normalizeUsing CPM \
	    --smoothLength 60 \
		--extendReads 150 \
		--centerReads"%(bam_file,input_bam_file,bigwig_file)
	os.system(cmd)

def getfasta_from_bed(input_fasta_file,bed_file,out_file):
	cmd = "bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>"
	os.system(cmd)

def calc_cov(sort_bam_file,cov_file):
    cmd_bedtools = 'bedtools genomecov -ibam %s > %s'%(sort_bam_file,cov_file)
    print(cmd_bedtools)
    os.system(cmd_bedtools)

def calc_cov_read(sort_bam_file,cov_file):
    cmd_bedtools = 'bedtools genomecov -dz -ibam %s > %s'%(sort_bam_file,cov_file)
    print(cmd_bedtools)
    os.system(cmd_bedtools)

def get_mapped_bam_file(bam_file,map_bam_file):
    cmd_samtools = 'samtools view -@ 25 -bF 12 %s > %s'%(bam_file,map_bam_file)
    os.system(cmd_samtools)

def homer_motif(input_file,outdir,genome_file,size):
	homer_path = "/ddn/LiB/softwares/Homer/bin"
	cmd = "%s/findMotifsGenome.pl %s %s %s -size %s"%(homer_path,input_file,genome_file,outdir,str(size))
	os.system(cmd)
	#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]

def run_homer(sample_info_file,callpeak_dir,homer_dir,genome_file):
	with open(sample_info_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		sample_name = line[header.index('sample_name')]
		sample_id = line[header.index('sample_id')]
		sample_type = line[header.index('sample_type')]
		c_callpeak_dir = os.path.join(callpeak_dir,sample_id)
		c_callpeak_file = os.path.join(c_callpeak_dir,sample_id+'_peaks.narrowPeak')
		c_homer_dir = os.path.join(homer_dir,sample_id)
		mkdir(c_homer_dir)
		homer_motif(c_callpeak_file,c_homer_dir,genome_file,200)

if __name__ == '__main__':
	#step1.create sample information

	save_sample_info_file = "/ddn/LiB/project/chip-seq/data/examples/chipseq_sample.txt"
	result_dir = "/ddn/LiB/project/chip-seq/data/analysis/examples"
	sample_dir = "/ddn/LiB/project/chip-seq/data/examples"	
	prefix = "chipseq"
	mouse_bowtie2_index_file = "/ddn/LiB/data/ref_genome/mouse/bowtie2_index/mm10"

	#step2: run fastqc
	fastqc_dir = os.path.join(result_dir,'FastQC')
	mkdir(fastqc_dir)
	#run_fastqc(save_sample_info_file,sample_dir,fastqc_dir,'CHIP-seq_fastqc')

	#step3: remove adapter
	cutadapt_dir = os.path.join(result_dir,'cutadapt')
	mkdir(cutadapt_dir)
	#run_cutadapt(save_sample_info_file,sample_dir,cutadapt_dir,'no')

	#step3:run bowtie2
	#bowtie2_build(mouse_genome_file,mouse_bowtie2_index_file)
	bowtie2_dir = os.path.join(result_dir,'bowtie2')
	mkdir(bowtie2_dir)
	#run_bowtie2(save_sample_info_file,cutadapt_dir,bowtie2_dir,mouse_bowtie2_index_file,20)

	#step4: call peaks
	callpeak_dir = os.path.join(result_dir,'callpeak')
	mkdir(callpeak_dir)
	run_callpeak(save_sample_info_file,bowtie2_dir,callpeak_dir,'na')

	#step5: motif enrichment analysis
	homer_dir = os.path.join(result_dir,'homer')
	run_homer(save_sample_info_file,callpeak_dir,homer_dir,'mm10')
