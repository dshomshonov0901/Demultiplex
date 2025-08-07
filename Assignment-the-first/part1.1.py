#!/usr/bin/env python3
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH -o slurm.out
#SBATCH -e slurm.err
import gzip
import matplotlib.pyplot as plt


#file1 = "test1.fq"

file1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
file2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"   
file3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"    
file4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"     

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33

def populate_list(file_path: str, length):
    my_list_r1 = [0] * length
    line_count = 0
    with gzip.open(file_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            if i % 4 == 0:  # quality line
                line = line.strip()
                for n, score in enumerate(line):
                    my_list_r1[n] += ord(score) - 33
            line_count += 1
    num_reads = line_count // 4
    return my_list_r1, num_reads


def plot_mean_quality(mean_scores, title):
    plt.figure()
    plt.plot(range(len(mean_scores)), mean_scores)
    plt.title(title)
    plt.xlabel("Base Position")
    plt.ylabel("Mean Quality Score")
    plt.savefig(f"{title.replace(' ', '_')}.png")
    plt.close()

my_list_r1, total_reads = populate_list(file4, 101)
mean_scores = [round(val / total_reads, 2) for val in my_list_r1]
plot_mean_quality(mean_scores, "Read 2 Quality Scores")

