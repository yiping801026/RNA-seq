#!/bin/bash
#
#
nohup python /home/lib/data/softwares/SurVirus-master/surveyor.py /ddn/LiB/pingyi/HBV/try/B993-2_L4_R1.fastq,/ddn/LiB/pingyi/HBV/try/B993-2_L4_R2.fastq /home/lib/data/pingyi/HBV/try/B993  /home/lib/data/softwares/genome_fasta/hg38.fa /home/lib/data/softwares/genome_fasta/U95551.1.HBV.fa /home/lib/data/softwares/genome_fasta/hg38+HBV.fa --wgs --fq --dust /ddn/LiB/softwares/sdust-master/sdust --threads 32 > B993report.log 2>&1 &
