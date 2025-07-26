# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. A good quality score cutoff for index reads is 30. I chose such a high score because as index reads are used to demultiplex reads to the correct sample, misassigning even a single base can result in incorrect sample assignment. A threshold of 30 ensures a 99.9% accuracy. For biological read pairs, I think a 20 cutoff is sufficient, since downstream tools like aligners can handle base errors. A quality score of 20 corresponds to a 99% base call accuracy and is a common threshold.
    3. i1 N count = 3976613 , i2 N count - 3328051
    
## Part 2
1. Define the problem  / Describe Ouput

```
We are demultiplexing to separate properly matched index pairs from index-hopped or unknown index reads to ensure accuracy in downstream analysis. Our sequencing data consists of four parallel FASTQ files: Read 1 (R1), Read 2 (R4), Index 1 (R2), and Index 2 (R3). Each read has a pair of index sequences that should match the original sample assignment. We will compare the index pairs in each read to a list of 24 valid index sequences. Based on these comparison:

1. If both indexes are valid and match we write to matched FASTQ files
2. If both are valid but do not match  we write to index-hopped files
3. If either index is invalid (contains 'N') or has low quality we write to unknown files

All read headers will be modified to include the index pair (e.g., @read1 ACGTACGT-TGCATGCA). Output files will be separated by category and used later.
```


2. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
3. Pseudocode
 ```
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
```
6. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

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

