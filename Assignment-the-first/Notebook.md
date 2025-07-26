To find sequence length:
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
101 for reads, 8 for indexes

To find phred encoding:
zcat 1294_S1_L008_R1_001.fastq.gz | head -4 | tail -1 (gets a score line)
result: A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ
no lower case letters, capital J indicates phred +33

Psuedocode:

We are demultiplexing because our library prep has some unknown sequences/indexes and some index hopping. To separate these from our properly matched indexes we are going to open
4 fastq files which contain the forward read, index for forward read, index for reverse read, and reverse read in 1294_S1_L008_R1_001.fastq.gz, 1294_S1_L008_R2_001.fastq.gz,
1294_S1_L008_R3_001.fastq.gz, and 1294_S1_L008_R4_001.fastq.gz respectively. To start we should open these 4 fastqs as readable and open 6 new fastqs to write into.
There will be 2 new fastqs ((one for forward read, and one for reverse) for each of the following categories:
1. read-pairs with properly matched indexes
2. read pairs with index-hopping observed
3. read-pairs with unknown index(es)

We are demultiplexing to separate properly matched index pairs from index-hopped or unknown index reads to ensure accuracy in downstream analysis.
Our sequencing data consists of four parallel FASTQ files: Read 1 (R1), Read 2 (R2), Index 1 (I1), and Index 2 (I2). Each read has a pair of index sequences
that should match the original sample assignment. We will compare the index pairs in each read to a list of 24 valid index sequences. Based on the comparison:
1. If both indexes are valid and match → write to matched FASTQ files
2. If both are valid but do not match → write to index-hopped files
3. If either index is invalid (contains 'N') or has low quality → write to unknown files

All read headers will be modified to include the index pair (e.g., @read1 ACGTACGT-TGCATGCA). Output files will be separated by category and used for later quality filtering and analysis.

Once all these files are open we should loop through the sequence line which is the 2nd line in the two index fastq records of 4 lines. We can use modulo to get this line.
We then need to check if the indexes to the 24 available ones in the text file 

1. open the 4 fq files to read and 6 files to write
2. loop through the fastq files at the same time
3. to get the sequence/index lines loop through the 2nd line in every 4 by using modulo
4. specifically for R2 and R3 files, search for whether the index line matches a index from the index text file
5. If there are no matches add the info to unkown fq
6.If there is a match and R2 and R3 also are equal to each other add to prop.matched fq
7. If there is a match but R2 and R3 are not equal add to index hopping.




STRATEGY:
We are demultiplexing to separate properly matched index pairs from index-hopped or unknown index reads to ensure accuracy in downstream analysis.
Our sequencing data consists of four parallel FASTQ files: Read 1 (R1), Read 2 (R2), Index 1 (I1), and Index 2 (I2). Each read has a pair of index sequences
that should match the original sample assignment. We will compare the index pairs in each read to a list of 24 valid index sequences. Based on the comparison:
1. If both indexes are valid and match → write to matched FASTQ files
2. If both are valid but do not match → write to index-hopped files
3. If either index is invalid (contains 'N') or has low quality → write to unknown files

All read headers will be modified to include the index pair (e.g., @read1 ACGTACGT-TGCATGCA). Output files will be separated by category and used for later quality filtering and analysis.

PSEUDOCODE:
DEFINE FUNCTION demultiplex_reads:
    INPUT:
        - r1_file, r2_file, i1_file, i2_file: paths to the 4 input FASTQ files
        - valid_indexes:  file of known good index sequences
        - quality_cutoff: minimum average Phred quality for indexes

    OPEN all input files for reading

    FOR EACH index IN valid_indexes:
        OPEN matched_r1_<index>.fq and matched_r2_<index>.fq for writing

    OPEN:
        - index_hop_r1.fq and index_hop_r2.fq
        - unknown_r1.fq and unknown_r2.fq

    WHILE not at end of input files:
        READ every 4 lines from each file (R1, R2, I1, I2) as one read block

        EXTRACT:
            - i1_seq = line 2 of I1
            - i2_seq = line 2 of I2, reverse complement
            - i1_quality = line 4 of I1
            - i2_quality = line 4 of I2

        IF i1_seq or i2_seq contains 'N' OR average quality < cutoff:
            WRITE R1/R2 read pair to unknown files
            CONTINUE

        IF i1_seq == i2_seq AND i1_seq in valid_indexes:
            WRITE R1/R2 to matched_r1_<i1_seq>.fq and matched_r2_<i1_seq>.fq
        ELSE IF i1_seq and i2_seq are both in valid_indexes:
            WRITE R1/R2 to index_hop_r1.fq and index_hop_r2.fq
        ELSE:
            WRITE R1/R2 to unknown_r1.fq and unknown_r2.fq

        APPEND "i1_seq-i2_seq" to header lines of R1 and R2

    CLOSE all files

    RETURN counts of: matched, hopped, unknown reads

HEADER:
   
```python
        
        def(demultiplex_reads(r1_file: str, r2_file: str, i1_file: str, i2_file: str, valid_indexes: set[str], quality_cutoff: float):
    
        Demultiplexes paired-end FASTQ reads into matched, index-hopped, and unknown categories.
    
        Parameters:
            r1_file, r2_file: Paths to biological read FASTQ files
            i1_file, i2_file: Paths to index read FASTQ files
            valid_indexes: Set of 24 valid index sequences
            quality_cutoff: Minimum average Phred score for index sequences

        Returns:
            A dictionary with counts of matched, hopped, and unknown read pairs
```

