
#!/usr/bin
#-*-coding:utf-8-*-

import os
import threading
import multiprocessing


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

def trim_galore(input_fastq,outdir,pair_flag):
    if pair_flag=='paired':
        cmd = "trim_galore --length 35 --adapter %s --paired --retain_unpaired --cores 20 --gzip --fastqc -o %s %s"%(adapter_seq,outdir,input_fastq)
    else:
        cmd = "trim_galore --length 35 --adapter %s --cores 20 --gzip --fastqc -o %s %s"%(adapter_seq,outdir,input_fastq)
    os.system(cmd)
    
def star_index(script,ref_star_path,fasta_file,gtf_file,thread_num=20):
    command = "%s --runThreadN %s --runMode genomeGenerate\
     --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --limitGenomeGenerateRAM 214748364800"%(script,str(thread_num),ref_star_path,fasta_file,gtf_file)
    os.system(command)
    print(command)

def star_align(script,fastq_file,ref_fasta_path,ref_gtf_path,save_prefix,thread_num=30,gz_flag='no'):
    command = "%s --runThreadN %s --genomeDir %s --readFilesIn %s --outFileNamePrefix %s --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 314748364800 "%(script,str(thread_num),ref_fasta_path,fastq_file,save_prefix)
    #if gz_flag=='yes':
    #    command = command +" --readFilesCommand zcat"
    os.system(command)
    bam_file = save_prefix+'Aligned.sortedByCoord.out.bam'
    command = "samtools index %s"%(bam_file)
    os.system(command)
    print(command)

def featureCounts(script,ref_gtf_file,save_prefix,star_align_bam_file):
    command = "%s -a %s -o %s -O -R BAM %s -T 4 -t CDS"%(script,ref_gtf_file,save_prefix,star_align_bam_file)
    os.system(command)
    print(command)

def samtools_sort_index(featureCounts_bam_file,align_gene_sort_index_bam_file):
    command = "samtools sort -@ 20 %s -o %s"%(featureCounts_bam_file,align_gene_sort_index_bam_file)
    os.system(command)
    print(command)
    command = "samtools index %s"%(align_gene_sort_index_bam_file)
    os.system(command)
    print(command)

def get_gene_count(bam_file,ref_gtf_database,save_dir,prefix):
    featureCounts_script = "/ddn/LiB/softwares/subread-2.0.2-Linux-x86_64/bin/featureCounts"
    save_map_gene_prefix = os.path.join(save_dir,prefix+'_map_gene')
    featureCounts(featureCounts_script,ref_gtf_database,save_map_gene_prefix,bam_file)

    sort_index_bam_file = os.path.join(save_dir,prefix+'_map_gene_sort')
    map_gene_bam_file = os.path.join(save_dir,prefix+'_map_virusAligned.sortedByCoord.out.bam.featureCounts.bam')
    samtools_sort_index(map_gene_bam_file,sort_index_bam_file)
    
    #get gene counts from mapped gene results
    gene_count_file = os.path.join(save_dir,prefix+'_gene_count.tsv.gz')
    umi_count(sort_index_bam_file,gene_count_file)

def parse_gene_count(gene_count_gz_file,gtf_file,save_prefix):
    un_gz_file = save_prefix+'_gene_count.txt'
    gzip(gene_count_gz_file,un_gz_file)
    gtf_info_dict = pasre_gtf(gtf_file)
    save_gene_count_file = save_prefix+'_gene_cell_count_matrix.txt'
    with open(un_gz_file) as f:
        contents = f.readlines()    
    gene_cell_count_dict = {}
    gene_id_list = []
    for line in contents[1:]:
        line = line.strip().split('\t')
        gene_ids = line[0].strip(',').split(',')
        cell_id = line[1]
        read_count = line[2]
        if cell_id not in gene_cell_count_dict.keys():
            gene_cell_count_dict.update({cell_id:{}})
        for gene_id in gene_ids:
            gene_id = gtf_info_dict[gene_id]
            if gene_id not in gene_cell_count_dict[cell_id].keys():
                gene_cell_count_dict[cell_id].update({gene_id:read_count})
            if gene_id not in gene_id_list:
                gene_id_list.append(gene_id)
    f_save = open(save_gene_count_file,'w')
    f_save.write('cell_id\t'+'\t'.join(gene_id_list)+'\n')
    for cell_id in gene_cell_count_dict.keys(): 
        gene_count_list = []
        for gene_id in gene_id_list:
            if gene_id in gene_cell_count_dict[cell_id].keys():
                gene_count_list.append(gene_cell_count_dict[cell_id][gene_id])
            else:
                gene_count_list.append('0')
        f_save.write(cell_id+'\t'+'\t'.join(gene_count_list)+'\n')
        f_save.flush()
    f_save.close()

def bowtie2_mapping_trans(read1_file,read2_file,sam_file,ref_db,bowtie2_log_file,pair_flag):
    if pair_flag.upper()=='SINGLE':
        cmd_bowtie2 = 'bowtie2 -q -p 10 -sensitive -X 2000 --non-deterministic -x %s -U %s -S %s > %s 2>&1'%(ref_db,read1_file,sam_file,bowtie2_log_file)
        print(cmd_bowtie2)
        os.system(cmd_bowtie2)
    else:
        cmd_bowtie2 = 'bowtie2 -q -p 10 -sensitive -X 2000 --non-deterministic -x %s -1 %s -2 %s -S %s > %s 2>&1'%(ref_db,read1_file,read2_file,sam_file,bowtie2_log_file)
        print(cmd_bowtie2)
        os.system(cmd_bowtie2)

def bowtie2_mapping(read1_file,read2_file,sam_file,ref_db,bowtie2_log_file,pair_flag):
    if pair_flag.upper()=='SINGLE':
        cmd_bowtie2 = 'bowtie2 -q -p 10 -sensitive -X 2000 --non-deterministic -x %s -U %s -S %s > %s 2>&1'%(ref_db,read1_file,sam_file,bowtie2_log_file)
        print(cmd_bowtie2)
        os.system(cmd_bowtie2)
    else:
        cmd_bowtie2 = 'bowtie2 -q -p 10 -sensitive -X 2000 --non-deterministic -x %s -1 %s -2 %s -S %s > %s 2>&1'%(ref_db,read1_file,read2_file,sam_file,bowtie2_log_file)
        print(cmd_bowtie2)
        os.system(cmd_bowtie2)

def sam2bam(sam_file,bam_file):
    cmd_convert = 'samtools view -@ 25 -bS %s > %s'%(sam_file,bam_file)
    os.system(cmd_convert)

def rm_file(file_path):
    if os.path.exists(file_path):
        cmd_rm = 'rm -rf %s'%file_path
        os.system(cmd_rm)

def check_mapping(bam_file,mapping_quality_file):
    cmd_flgstat = 'samtools flagstat %s > %s'%(bam_file,mapping_quality_file)
    os.system(cmd_flgstat)

def get_mapped_bam_file(bam_file,map_bam_file):
    cmd_samtools = 'samtools view -@ 25 -bF 12 %s > %s'%(bam_file,map_bam_file)
    os.system(cmd_samtools)

def bam2sam(bam_file,sam_file):
    cmd_bam2sam = 'samtools view -h %s > %s'%(bam_file,sam_file)
    os.system(cmd_bam2sam)

def sort_bam_file(bam_file,sort_bam_file_prefix):
    cmd_sort = 'samtools sort -@ 25 -o %s %s'%(sort_bam_file_prefix,bam_file)
    os.system(cmd_sort)
    print(cmd_sort)

def index_bam_file(bam_sort_file):
    cmd_index = 'samtools index %s'%(bam_sort_file)
    os.system(cmd_index)
    print(cmd_index)

def calc_cov(sort_bam_file,cov_file):
    cmd_bedtools = 'bedtools genomecov -ibam %s > %s'%(sort_bam_file,cov_file)
    print(cmd_bedtools)
    os.system(cmd_bedtools)

def get_mapped_bam_sam_file(bam_file, mapped_bam_file,mapped_sam_file):
    get_mapped_bam_file(bam_file, mapped_bam_file)
    bam2sam(mapped_bam_file, mapped_sam_file)

def count_map_bam_forach(bam_file,count_dir,count_file):
    cmd = "samtools idxstats %s > %s"%(bam_file,count_file)
    os.system(cmd)
    with open(count_file) as f:
        contents = f.readlines()
    for line in contents[1:]:
        line = line.strip().split('\t')
        match_region = line[0]
        c_sort_bam_file = os.path.join(count_dir,match_region+'.bam')
        cmd = "samtools view -b %s %s > %s"%(bam_file,match_region,c_sort_bam_file)
        os.system(cmd)
        c_cov_file = os.path.join(count_dir,match_region+'_cov.txt')
        calc_cov(c_sort_bam_file,c_cov_file)

def sub_job_only_one(read1_file, read2_file, sam_file, bam_file, mapping_quality_file, mapped_bam_file,
                     mapped_sam_file, bam_sort_file_prefix, bam_sort_file, cov_file,mapping_result_file,bowtie2_log_file,pair_flag):
    cmd_srun = 'srun -N 1 -p cpuall.q python bowtie2_mapping_analysis.py %s %s %s %s %s %s %s %s %s %s %s %s'%(read1_file, read2_file, sam_file, bam_file, mapping_quality_file, mapped_bam_file,
                     mapped_sam_file, bam_sort_file_prefix, bam_sort_file, cov_file,mapping_result_file,bowtie2_log_file,pair_flag)
    os.system(cmd_srun)

def mkdir(dirPath):
    cmd_mkdir = 'mkdir -p %s'%dirPath
    os.system(cmd_mkdir)

def rm_dir(dirPath):
    if os.path.exists(dirPath):
        cmd_rm = 'rm -rf %s'%dirPath
        os.system(cmd_rm)


def fastq_dump(sra_file,outdir):
    cmd = "/home/chencg/projects/sc_metabolism/multiorgan_tcell/sratoolkit.2.11.0-ubuntu64/bin/fastq-dump --split-3 -O %s %s"%(outdir,sra_file)
    os.system(cmd)
    print(cmd)

def calc_cov(sort_bam_file,cov_file):
    cmd_bedtools = 'bedtools genomecov -dz -ibam %s > %s'%(sort_bam_file,cov_file)
    print(cmd_bedtools)
    os.system(cmd_bedtools)

def count_map_bam_forach(bam_file,count_dir,count_file):
    cmd = "samtools idxstats %s > %s"%(bam_file,count_file)
    os.system(cmd)
    with open(count_file) as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip().split('\t')
        match_region = line[0]
        c_sort_bam_file = os.path.join(count_dir,match_region+'.bam')
        cmd = "samtools view -b %s %s > %s"%(bam_file,match_region,c_sort_bam_file)
        os.system(cmd)
        c_cov_file = os.path.join(count_dir,match_region+'_cov.txt')
        calc_cov(c_sort_bam_file,c_cov_file)

def sort_bam_file(bam_file,sort_bam_file_prefix):
    cmd_sort = 'samtools sort -@ 25 -o %s %s'%(sort_bam_file_prefix,bam_file)
    os.system(cmd_sort)
    print(cmd_sort)

def merge_results(raw_data_metadata,result_save_root_dir,save_dir):
    #step1: merge all bam files
    mapped_bam_file = os.path.join(result_save_root_dir,'*','*_bowtie2.bam.mapped.bam')
    merge_bam_file = os.path.join(save_dir,'merged_bowtie2.bam.mapped.bam')
    cmd = "samtools merge -f %s %s"%(merge_bam_file,mapped_bam_file)
    os.system(cmd)
    merge_bam_file_sort = os.path.join(save_dir,'merged_bowtie2.bam.mapped_sort.bam')
    sort_bam_file(merge_bam_file,merge_bam_file_sort)
    cmd = "samtools index %s"%(merge_bam_file_sort)
    os.system(cmd)
    count_dir = os.path.join(result_save_root_dir,'count_each_item')
    mkdir(count_dir)
    count_file = os.path.join(result_save_root_dir,'count_mapped_number.txt')
    count_map_bam_forach(merge_bam_file_sort,count_dir,count_file)
    cov_file = os.path.join(save_dir,'merged_bowtie2_mapped_cov.txt')
    calc_cov(merge_bam_file_sort,cov_file)
    save_cov_add_flag_file = os.path.join(save_dir,'merged_bowtie2_mapped_cov_add_flag.txt')
    get_plot_R_map(cov_file,save_cov_add_flag_file)

def get_plot_R_map(cov_file,save_file):
    f_save = open(save_file,'w')
    f_save.write('crispr_cas_id\thit_position\thit_read_number\tc2c10\thit_true_position\tregion_class\n')
    with open(cov_file) as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip().split('\t')
        array_id = line[0]
        hit_position = line[1]
        hit_read_number = line[2]
        hit_position = int(hit_position)+int(array_id.split('-')[0].split('_')[-1])-1
        if 'up' in array_id:
            pro_expand_region = array_id.split('_')[-2]
            pro_region = array_id.split('_')[-3]
            middle_region = array_id.split('_')[-4]
            crispr_region = array_id.split('_')[-5]
            crispr_expand_region = array_id.split('_')[-6]
        else:
            pro_expand_region = array_id.split('_')[-6]
            pro_region = array_id.split('_')[-5]
            middle_region = array_id.split('_')[-4]
            crispr_region = array_id.split('_')[-3]
            crispr_expand_region = array_id.split('_')[-2]

        if int(hit_position)>=int(pro_expand_region.split('-')[0]) and int(hit_position)<=int(pro_expand_region.split('-')[1]):
            region_class = 'c2c10_expand_region'
        else:
            if int(hit_position)>=int(pro_region.split('-')[0]) and int(hit_position)<=int(pro_region.split('-')[1]):
                region_class = 'c2c10'
            else:
                if int(hit_position)>=int(middle_region.split('-')[0]) and int(hit_position)<=int(middle_region.split('-')[1]):
                    region_class = 'c2c10_crispr_non_coding_region'
                else:
                    if int(hit_position)>=int(crispr_region.split('-')[0]) and int(hit_position)<=int(crispr_region.split('-')[1]):
                        region_class = 'crispr'
                    else:
                        region_class = 'crispr_expand_region'
        f_save.write('\t'.join(line).strip()+'\t'+array_id+'\t'+str(hit_position)+'\t'+region_class+'\n')
        f_save.flush()
    f_save.close()

def before():
    raw_data_metadata = '/100T/ganr/bgc_metatranscriptomics/attention_bgc_IBD_metatranscriptome/run_metatranscriptome42_list.txt'
    result_save_root_dir = '/100T/ganr/bgc_metatranscriptomics/attention_bgc_IBD_metatranscriptome/bowtie2_mapping'
    sample_save_root_dir = "/100T/ganr/metabolic/huada/restart_bgc/metranscriptomics/run_diamond_mapping/diamond_blastx_metatranscriptomics"
    
    raw_data_metadata = '/home/data/metatranscriptomics/run_list_merge.txt'
    result_save_root_dir = '/home/ganr/project/crispr/metagenome_analysis/analysis/human_gut/analysis_result/round2/metatranscriptome_analysis/bowtie2_mapping'
    sample_save_root_dir = "/home/data/metatranscriptomics/raw_data"
    sra_dir = "/home/data/metatranscriptomics/data"

    # pool = multiprocessing.Pool(processes=int(5))
    # with open(raw_data_metadata,'r') as fin:
    #     lines = fin.readlines()        
    # for line in lines[1:]:
    #     cur_sra_id = line.strip().split('\t')[0]
    #     pair_flag = line.strip().split('\t')[1]
    #     download_path = line.strip().split('\t')[2]
    #     c_sra_file = os.path.join(sra_dir,cur_sra_id,download_path.split('/')[-1])
    #     c_outdir = os.path.join(sample_save_root_dir,cur_sra_id)
    #     mkdir(c_outdir)
    #     if os.path.exists(c_outdir):
    #         if len(os.listdir(c_outdir))>=1:
    #             continue
    #     pool.apply_async(fastq_dump,(c_sra_file,c_outdir))
    # pool.close()
    # pool.join()

    tsk = []
    with open(raw_data_metadata,'r') as fin:
        lines = fin.readlines()        
        for line in lines[816:]:
            while True:
                if len(threading.enumerate()) < 6:
                    break
            cur_sra_id = line.strip().split('\t')[0]
            pair_flag = line.strip().split('\t')[1]
            download_path = line.strip().split('\t')[2]
            print(cur_sra_id)
            name_id = download_path.split('/')[-1]
            cur_result_save_dir = '%s/%s'%(result_save_root_dir,cur_sra_id)
            mkdir(cur_result_save_dir)
            if pair_flag.upper()=='PAIRED':
                cur_read1_fastq = os.path.join(sample_save_root_dir,cur_sra_id,name_id+'_1.fastq')
                cur_read2_fastq = os.path.join(sample_save_root_dir,cur_sra_id,name_id+'_2.fastq')
            else:
                cur_read1_fastq = os.path.join(sample_save_root_dir,cur_sra_id,name_id+'.fastq')
                cur_read2_fastq = 'NA'
            cur_sam_file = '%s/%s_bowtie2.sam'%(cur_result_save_dir,cur_sra_id)
            cur_bam_file = '%s/%s_bowtie2.bam'%(cur_result_save_dir,cur_sra_id)
            cur_mapping_quality_file = '%s/%s_bowtie2.bam.quality'%(cur_result_save_dir,cur_sra_id)
            cur_mapped_bam_file = '%s/%s_bowtie2.bam.mapped.bam'%(cur_result_save_dir,cur_sra_id)
            cur_mapped_sam_file = '%s/%s_bowtie2.bam.mapped.sam'%(cur_result_save_dir,cur_sra_id)
            cur_bam_sort_file_prefix = '%s/%s_bowtie2.bam.sort'%(cur_result_save_dir,cur_sra_id)
            cur_bam_sort_file = '%s/%s_bowtie2.bam.sort.bam'%(cur_result_save_dir,cur_sra_id)
            cur_cov_file = '%s/%s_bowtie2_cov_result.txt'%(cur_result_save_dir,cur_sra_id)
            cur_mapping_result_file = '%s/%s_bowtie2_mapping_result.txt'%(cur_result_save_dir,cur_sra_id)
            cur_bowtie2_log_file =  '%s/%s_bowtie2_log.txt'%(cur_result_save_dir,cur_sra_id)
            cur_mappint_cov_result_file = '%s/%s_bowtie2_mapping_coverage_result.txt'%(cur_result_save_dir,cur_sra_id)
            cur_sample_id = cur_sra_id
            # if os.path.exists(cur_mapping_quality_file):
            #     continue
            thread_name = MyThread(cur_read1_fastq, cur_read2_fastq, cur_sam_file, cur_bam_file, cur_mapping_quality_file,cur_mapped_bam_file,
                     cur_mapped_sam_file, cur_bam_sort_file_prefix, cur_bam_sort_file, cur_cov_file,cur_mapping_result_file,
                     cur_bowtie2_log_file,cur_mappint_cov_result_file,cur_sample_id,pair_flag)
            tsk.append(thread_name)
            thread_name.start()
    for t in tsk:
        t.join()

    save_dir = "/home/ganr/project/crispr/metagenome_analysis/analysis/human_gut/analysis_result/round2/metatranscriptome_analysis/analysis1"
    crispr_cas_fa_file = '/home/ganr/project/crispr/metagenome_analysis/analysis/human_gut/analysis_result/round2/result/specific_analysis/c2c10/human_gut_round12_c2c10_filter_distance_crispr_homo_c2c10_filter_add_location_nucl_no_spacer_flag_expand_protein.fa'
    merge_results(raw_data_metadata,result_save_root_dir,save_dir)

def run_fastqc(sample_info_file,sample_dir,fastqc_dir,prefix):
    # with open(sample_info_file) as f:
    #     contents = f.readlines()
    # for line in contents[1:]:plot
    #     line = line.strsip().split('\t')
    #     sample_name = line[0]
    #     c_sample_file = os.path.join(sample_dir.sample_name)
    fastqc_script = "/ddn/LiB/softwares/fastqc_v0.11.9/FastQC/fastqc"
    input_fastq = os.path.join(sample_dir,'*.fastq.gz')
    #fastqc(fastqc_script,input_fastq,fastqc_dir)
    multiqc_dir = os.path.join(fastqc_dir,'multiqc')
    mkdir(multiqc_dir)
    multiqc(fastqc_dir,multiqc_dir,prefix)

def run_star_counts(star_script,sample_info_file,sample_dir,star_dir,ref_star_path,ref_gtf_file,thread_num):
    with open(sample_info_file) as f:
        contents = f.readlines()
    for line in contents[0:]:
        line = line.strip().split('\t')
        sample_name = line[0]
        sample_id = '_'.join(sample_name.split('_')[0:3])
        c_sample_file = os.path.join(sample_dir,sample_name)
        save_dir = os.path.join(star_dir,sample_id)
        mkdir(save_dir)
        save_prefix = os.path.join(save_dir,sample_id)
        star_align(star_script,c_sample_file,ref_star_path,ref_gtf_file,save_prefix,thread_num,'yes')
        
        star_align_bam_file = os.path.join(save_dir,sample_id+'Aligned.sortedByCoord.out.bam')
        counts_dir = os.path.join(save_dir,'featureCounts')
        mkdir(counts_dir)
        save_counts_prefix = os.path.join(counts_dir,sample_id+'_featureCounts_')
        featureCounts_script = "/ddn/LiB/softwares/subread-2.0.2-Linux-x86_64/bin/featureCounts"
        featureCounts(featureCounts_script,ref_gtf_file,save_counts_prefix,star_align_bam_file)

def get_sample_info(sample_info_file,save_sample_info_file):
    f_save = open(save_sample_info_file,'w')
    f_save.write('sample_name\ttissue\tpatient\tphenotype\tphenotype1\n')
    with open(sample_info_file) as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip().split('\t')
        sample_name = '_'.join(line[0].split('_')[0:3])
        tissue_name = line[0].split('_')[1]
        patient = '_'.join(line[0].split('_')[0:2])
        phenotype = line[0].split('_')[0]
        f_save.write(sample_name+'\t'+tissue_name+'\t'+patient+'\t'+patient+'\t'+phenotype+'\n')
        f_save.flush()
    f_save.close()



if __name__ == '__main__':
    sample_info_file = "/ddn/LiB/pingyi/RNA_seq/RNA_lists.txt"
    #save_sample_info_file = "/ddn/LiB/ganrui/project/RNA-seq/sleep/liuqinghua/analysis/liuqinghua_rnaseq_sample_info.txt"
    sample_dir = "/ddn/LiB/pingyi/RNA_seq/data_RNA_seq"
    result_dir = "/ddn/LiB/pingyi/RNA_seq/results"
    star_script = "/ddn/LiB/softwares/STAR-2.7.9a/bin/Linux_x86_64/STAR"
    ref_star_path = "/ddn/LiB/data/ref_genome/mouse/star_index/Mus_musculus.GRCm38_star"
    mouse_fasta_file = "/ddn/LiB/data/ref_genome/mouse/genome/Mus_musculus.GRCm38.dna.toplevel.fa"
    mouse_gtf_file = "/ddn/LiB/data/ref_genome/mouse/genome/Mus_musculus.GRCm38.102.gtf"
    
    #get_sample_info(sample_info_file,save_sample_info_file)

    #step1:run fastqc
    #fastqc_dir = os.path.join(result_dir,'FastQC')
    #mkdir(fastqc_dir)
    #run_fastqc(sample_info_file,sample_dir,fastqc_dir,'RNA-seq_fastqc')

    #step2:build star index
    #star_index(star_script,ref_star_path,mouse_fasta_file,mouse_gtf_file,thread_num=20)

    #step3:star mapping and featureCounts
    star_dir = os.path.join(result_dir,'STAR')
    mkdir(star_dir)
    run_star_counts(star_script,sample_info_file,sample_dir,star_dir,ref_star_path,mouse_gtf_file,thread_num=20)
    