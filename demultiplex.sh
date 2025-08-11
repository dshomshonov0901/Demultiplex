#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=16                 #optional: number of cpus, default is 1
#SBATCH --mem=8GB                         #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=dansho@uoregon.edu   #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demultiplex                   #optional: job name
#SBATCH --output=demultiplex_%j.out              #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demultiplex_%j.err               #optional: file to store stderr from job, %j adds the assigned jobID

/usr/bin/time -v sbatch Assignment-the-third/demultiplex.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -idx Assignment-the-third/index.txt -o Assignment-the-third/output_final
