# Hi-C Pre-processing Steps

## Step 1: Download SRA File (SRA Toolkit)

Use Docker in WSL to avoid installation issues:

```bash
docker pull ncbi/sra-tools
docker images
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp --name yy-sra ncbi/sra-tools
docker ps -a
docker attach <container_id>
```

In the container, use the SRA Toolkit CLI:

```bash
# Check bin
which prefetch
/usr/local/bin/prefetch

prefetch SRR5579177
2024-12-26T13:38:50 prefetch.3.1.0:  HTTPS download succeed
2024-12-26T13:41:08 prefetch.3.1.0:  'SRR5579177' is valid
2024-12-26T13:41:08 prefetch.3.1.0: 1) 'SRR5579177' was downloaded successfully
```

## Step 2: Convert SRA to FASTQ (SRA Toolkit)

```bash
# With progress bars
fasterq-dump SRR5579177 -p
join   :|-------------------------------------------------- 100%   
concat :|-------------------------------------------------- 100%   
spots read      : 301,260,192
reads read      : 602,520,384
reads written   : 602,520,384
```

Default is `--split-3`, which means:
The spots are split into (biological) reads, for each read: 4 lines of FASTQ or 2 lines of FASTA are written. For spots having 2 reads, the reads are written into the *_1.fastq and *_2.fastq files. Unmated reads are placed in *.fastq. If the accession has no spots with one single read, the *.fastq-file will not be created.

## Step 3: Quality Control (FastQC)

### Quality Control with FastQC

Review FastQC Reports: Open the generated HTML reports to review the quality of your reads. Look for issues such as low-quality scores, adapter contamination, or overrepresented sequences.

Install FastQC:

```bash
# Clean environment
docker pull ubuntu
docker images
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp/ --name yy-biofp ubuntu
docker attach yy-biofp

# In container
apt update
apt install fastqc
```

Run FastQC:

```bash
fastqc reads_1.fastq reads_2.fastq
null
Started analysis of SRR5579177_1.fastq
null
Approx 5% complete for SRR5579177_1.fastq
Approx 10% complete for SRR5579177_1.fastq
Approx 15% complete for SRR5579177_1.fastq
Approx 20% complete for SRR5579177_1.fastq
Approx 25% complete for SRR5579177_1.fastq
Approx 30% complete for SRR5579177_1.fastq
Approx 35% complete for SRR5579177_1.fastq
Approx 40% complete for SRR5579177_1.fastq
Approx 45% complete for SRR5579177_1.fastq
Approx 50% complete for SRR5579177_1.fastq
Approx 55% complete for SRR5579177_1.fastq
Approx 60% complete for SRR5579177_1.fastq
Approx 65% complete for SRR5579177_1.fastq
Approx 70% complete for SRR5579177_1.fastq
Approx 75% complete for SRR5579177_1.fastq
Approx 80% complete for SRR5579177_1.fastq
Approx 85% complete for SRR5579177_1.fastq
Approx 90% complete for SRR5579177_1.fastq
Approx 95% complete for SRR5579177_1.fastq
Analysis complete for SRR5579177_1.fastq
Started analysis of SRR5579177_2.fastq
Approx 5% complete for SRR5579177_2.fastq
Approx 10% complete for SRR5579177_2.fastq
Approx 15% complete for SRR5579177_2.fastq
Approx 20% complete for SRR5579177_2.fastq
Approx 25% complete for SRR5579177_2.fastq
Approx 30% complete for SRR5579177_2.fastq
Approx 35% complete for SRR5579177_2.fastq
Approx 40% complete for SRR5579177_2.fastq
Approx 45% complete for SRR5579177_2.fastq
Approx 50% complete for SRR5579177_2.fastq
Approx 55% complete for SRR5579177_2.fastq
Approx 60% complete for SRR5579177_2.fastq
Approx 65% complete for SRR5579177_2.fastq
Approx 70% complete for SRR5579177_2.fastq
Approx 75% complete for SRR5579177_2.fastq
Approx 80% complete for SRR5579177_2.fastq
Approx 85% complete for SRR5579177_2.fastq
Approx 90% complete for SRR5579177_2.fastq
Approx 95% complete for SRR5579177_2.fastq
Analysis complete for SRR5579177_2.fastq
```

the outputs are 2 reports, I placed it in: fastqc-report dir

### Trimming and Filtering with Cutadapt

Install Cutadapt:

```bash
apt update
apt install cutadapt

# Check
which cutadapt
/usr/bin/cutadap
```

Run Cutadapt:

```bash
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20,20 -m 36 \
    -o trimmed_reads_SRR5579177_1.fastq -p trimmed_readsSRR5579177_2.fastq \
    SRR5579177_1.fastq SRR5579177_2.fastq

Processing paired-end reads on 1 core ...
[---------=8 ] 01:22:36    24,900,000 reads @ 180.0 µs/read;   0.33 M reads/minute
```

* `-a ADAPTER_FWD`: Adapter sequence for the forward read.
* `-A ADAPTER_REV`: Adapter sequence for the reverse read.
* `-q 20,20`: Trim low-quality bases from the ends of each read before adapter removal (Phred score < 20).
* `-m 36`: Discard reads shorter than 36 bases after trimming.
* `-o trimmed_reads_1.fastq`: Output file for trimmed forward reads.
* `-p trimmed_reads_2.fastq`: Output file for trimmed reverse reads.
* `reads_1.fastq reads_2.fastq`: Input FASTQ files.

This step may take a long time. If you accidentally close the terminal, you can use `reptyr`:

```bash
apt install reptyr
reptyr <pid>
```

After trimming, we got 2 trimed fastq files (trimmed_reads_<SRR_id>_<forward/reverse>.fastq):

```bash
[>8          ] 23:54:58   184,350,000 reads @ 179.2 µs/read;   0.33 M reads/minute
Done           29:18:49   301,260,192 reads @ 350.3 µs/read;   0.17 M reads/minute
Finished in 105529.934 s (350.295 µs/read; 0.17 M reads/minute).

=== Summary ===

Total read pairs processed:        301,260,192
  Read 1 with adapter:               5,912,381 (2.0%)
  Read 2 with adapter:               5,993,054 (2.0%)

== Read fate breakdown ==
Pairs that were too short:           6,402,578 (2.1%)
Pairs written (passing filters):   294,857,614 (97.9%)

Total basepairs processed: 30,126,019,200 bp
  Read 1: 15,063,009,600 bp
  Read 2: 15,063,009,600 bp
Quality-trimmed:             449,006,638 bp (1.5%)
  Read 1:   148,885,410 bp
  Read 2:   300,121,228 bp
Total written (filtered):  29,248,747,190 bp (97.1%)
  Read 1: 14,648,224,791 bp
  Read 2: 14,600,522,399 bp

=== First read: Adapter 1 ===

Sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC; Type: regular 3'; Length: 34; Trimmed: 5912381 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-34 bp: 3

Bases preceding removed adapters:
  A: 32.5%
  C: 24.6%
  G: 19.8%
  T: 22.7%
  none/other: 0.4%

Overview of removed sequences
length  count   expect  max.err error counts
3       4486107 4707190.5       0       4486107
4       1084796 1176797.6       0       1084796
5       163324  294199.4        0       163324
6       134437  73549.9 0       134437
7       2987    18387.5 0       2987
8       1234    4596.9  0       1234
9       3861    1149.2  0       575 3286
10      5126    287.3   1       185 4941
11      2278    71.8    1       173 2105
12      1334    18.0    1       145 1189
13      259     4.5     1       133 126
14      171     1.1     1       139 32
15      156     0.3     1       141 15
16      118     0.1     1       106 12
17      142     0.0     1       128 14
18      134     0.0     1       122 12
19      136     0.0     1       114 20 2
20      134     0.0     2       123 8 3
21      146     0.0     2       113 22 11
22      153     0.0     2       124 18 11
23      129     0.0     2       111 10 8
24      131     0.0     2       111 13 7
25      124     0.0     2       108 9 7
26      152     0.0     2       128 18 6
27      110     0.0     2       92 8 10
28      105     0.0     2       86 9 10
29      105     0.0     2       83 11 9 2
30      94      0.0     3       68 12 6 8
31      108     0.0     3       85 13 9 1
32      93      0.0     3       73 10 6 4
33      88      0.0     3       66 13 5 4
34      74      0.0     3       47 15 10 2
35      72      0.0     3       43 14 13 2
36      73      0.0     3       45 15 10 3
37      58      0.0     3       41 3 10 4
38      40      0.0     3       29 5 4 2
39      58      0.0     3       34 9 6 9
40      53      0.0     3       21 15 12 5
41      69      0.0     3       35 9 22 3
42      23      0.0     3       7 8 5 3
43      38      0.0     3       8 11 12 7
44      46      0.0     3       15 9 13 9
45      77      0.0     3       15 20 20 22
46      189     0.0     3       30 65 52 42
47      330     0.0     3       1 120 101 108
48      973     0.0     3       9 415 282 267
49      2642    0.0     3       3 1476 945 218
50      19294   0.0     3       12 17703 1385 194


=== Second read: Adapter 2 ===

Sequence: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT; Type: regular 3'; Length: 33; Trimmed: 5993054 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 33.0%
  C: 24.2%
  G: 19.9%
  T: 22.5%
  none/other: 0.4%

Overview of removed sequences
length  count   expect  max.err error counts
3       4539125 4707190.5       0       4539125
4       1080699 1176797.6       0       1080699
5       179833  294199.4        0       179833
6       150051  73549.9 0       150051
7       3614    18387.5 0       3614
8       1386    4596.9  0       1386
9       3857    1149.2  0       655 3202
10      5152    287.3   1       204 4948
11      2271    71.8    1       195 2076
12      1355    18.0    1       189 1166
13      282     4.5     1       151 131
14      212     1.1     1       153 59
15      194     0.3     1       167 27
16      123     0.1     1       109 14
17      162     0.0     1       145 17
18      147     0.0     1       127 18 2
19      143     0.0     1       122 19 2
20      167     0.0     2       134 26 7
21      131     0.0     2       113 11 7
22      160     0.0     2       128 18 14
23      155     0.0     2       127 19 9
24      150     0.0     2       121 20 9
25      126     0.0     2       104 17 5
26      135     0.0     2       109 17 9
27      114     0.0     2       92 16 6
28      123     0.0     2       99 14 10
29      99      0.0     2       82 13 2 2
30      113     0.0     3       84 16 7 6
31      99      0.0     3       74 14 7 4
32      109     0.0     3       77 16 9 7
33      99      0.0     3       73 14 7 5
34      71      0.0     3       52 9 8 2
35      78      0.0     3       45 16 6 11
36      74      0.0     3       47 11 10 6
37      59      0.0     3       36 13 7 3
38      57      0.0     3       21 22 7 7
39      65      0.0     3       25 22 9 9
40      51      0.0     3       11 17 10 13
41      81      0.0     3       18 30 19 14
42      70      0.0     3       15 23 15 17
43      91      0.0     3       8 46 23 14
44      166     0.0     3       9 76 49 32
45      128     0.0     3       9 53 35 31
46      142     0.0     3       22 45 44 31
47      250     0.0     3       0 133 70 47
48      454     0.0     3       6 128 142 178
49      1608    0.0     3       4 572 778 254
50      19223   0.0     3       7 17395 1400 421
```

### Run FastQC Again on Trimmed Reads to Verify Quality

```bash
fastqc trimmed_reads_1_paired.fastq trimmed_reads_2_paired.fastq
```

## Step 4: Align Reads (BWA-MEM, HiC-Pro)

### Prepare the referance genome (dm3)

according to the paper, we need to use the dm3 for the reference genome which can find in the UCSC Genome Browser:
https://hgdownload.soe.ucsc.edu/downloads.html

and the dm3 is here:
https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/

```bash
# download through http
wget https://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz

# unzip
gunzip dm3.fa.gz

ls -al 
-rwxrwxrwx 1 jasonyy jasonyy   172111398 Jan 23  2020 dm3.fa

head dm3.fa
>chr2L
Cgacaatgcacgacagaggaagcagaacagatatttagattgcctctcat
tttctctcccatattatagggagaaatatgatcgcgtatgcgagagtagt
gccaacatattgtgctctttgattttttggcaacccaaaatggtggcgga
tgaaCGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCA
TTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTG
CCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCG
CAAACCCAAAAAgacaatacacgacagagagagagagcagcggagatatt
tagattgcctattaaatatgatcgcgtatgcgagagtagtgccaacatat
tgtgctctCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAA
```

then the fastq file with a bit over 17 MB is prepared.

### Index the Reference Genome

```bash
sudo apt-get install bwa

# index the reference genome
bwa index dm3.fa
[bwa_index] Pack FASTA... 0.86 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=337473074, availableWord=35745488
[BWTIncConstructFromPacked] 10 iterations done. 58964098 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 108932194 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 153339938 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 192805618 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 227878786 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 259047874 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 286747026 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 311362066 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 333235922 characters processed.
[bwt_gen] Finished constructing BWT in 93 iterations.
[bwa_index] 88.56 seconds elapse.
[bwa_index] Update BWT... 0.85 sec
[bwa_index] Pack forward-only FASTA... 0.64 sec
[bwa_index] Construct SA from BWT and Occ... 43.81 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index dm3.fa
[main] Real time: 145.710 sec; CPU: 134.736 sec

#  align the reads to the reference genome, for paired-end reads
bwa mem -SP5M reference_genome.fa reads_1.fastq reads_2.fastq > aligned_reads.sam

bwa mem -SP5M dm3.fa SRR5579177_1.fastq SRR5579177_2.fastq > dm3-aligned_reads.sam
```

* `-S`: Output in SAM format.
* `-P`: Assume the input consists of interleaved paired-end sequences.
* `5M`: Set the mismatch penalty to 5.
* `reference_genome.fa`: The reference genome file.
* `reads_1.fastq reads_2.fastq`: Input FASTQ files.

you can find the flags here:
https://bio-bwa.sourceforge.net/bwa.shtml

### make .pairs file

find the `.sizes` file in: 
https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/

```bash
# Find ligation pairs in .sam data, make .pairs. SAM_PATH : an input .sam/.bam file with paired-end sequence alignments of Hi-C molecules. If the path ends with .bam, the input is decompressed from bam with samtools.
pairtools parse -c dm3.chrom.sizes -o output.pairsam input.sam

# Sort pairs in the lexicographic order along chrom1 and chrom2, in the numeric order along pos1 and pos2 and in the lexicographic order along pair_type.
pairtools sort -o sorted.pairsam output.pairsam

# Find and remove PCR/optical duplicates.
pairtools dedup -o dedup.pairsam sorted.pairsam

# Select pairs with "unique-unique"
pairtools select '(pair_type == "UU")' -o output.pairs dedup.pairsam

```





```bash
# build bowtie index of dm3
bowtie-build ../dm3.fa dm3_index

# you will get 6 index files
root@08e290bd86fc:/home/yy/bio-fp/bowtie# ls
dm3_index.1.ebwt  dm3_index.2.ebwt  dm3_index.3.ebwt  dm3_index.4.ebwt  dm3_index.rev.1.ebwt  dm3_index.rev.2.ebwt

# check the summary of the dm3 bowtie index
bowtie-inspect dm3_index -s
SA-Sample       1 in 32
FTab-Chars      10
Sequence-1      chr2L   23011544
Sequence-2      chr2LHet        368872
Sequence-3      chr2R   21146708
Sequence-4      chr2RHet        3288761
Sequence-5      chr3L   24543557
Sequence-6      chr3LHet        2555491
Sequence-7      chr3R   27905053
Sequence-8      chr3RHet        2517507
Sequence-9      chr4    1351857
Sequence-10     chrU    10049037
Sequence-11     chrUextra       29004656
Sequence-12     chrX    22422827
Sequence-13     chrXHet 204112
Sequence-14     chrYHet 347038
Sequence-15     chrM    19517


bowtie dm3_index -1 ../SRR5579177_1.fastq -2 ../SRR5579177_2.fastq  -t -a -m 1 --best -S dm3_bowtie_align_output.sam

# Align paired-end reads with Bowtie1
bowtie reference_index \ # use the prefix is ok
    -1 reads_1.fastq \
    -2 reads_2.fastq \
    -S \              # SAM format output
    -m 1 \           # Report only unique alignments
    --best \         # Report best alignment
    -t \            # Print timing info
    -a \
    output.sam


```

```bash
# sam to bam
samtools view -bS dm3_bowtie_align_output.sam > output.bam

# sort bam
samtools sort output.bam -o sorted.bam

# bam to pairs
pairtools parse sorted.bam \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --output output.pairs
```


alternatively:

```bash
# Build the Bowtie2 index
bowtie2-build reference.fasta reference_index


# Align paired-end reads (FASTQ) to the reference
bowtie2 -x reference_index \
    -1 reads_1.fastq \
    -2 reads_2.fastq \
    -S output.sam
```

bowtie manual: https://bowtie-bio.sourceforge.net/manual.shtml

| flags  |                                                          |
|--------|----------------------------------------------------------|
| -t     | Print the amount of wall-clock time taken by each phase. |
| -a     | Full path the bowtie2 installation directory             |
| -m 1   | Full path to the samtools installation directory         |
| --best | Full path to the R installation directory                |

### Process Aligned Reads with HiC-Pro

#### Installation

```bash
docker pull nservant/hicpro:latest
docker images
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp --name yy-hicpro nservant/hicpro
```

#### Annotation Files

##### BED File

BED file of the restriction fragments after digestion. This file depends both on the restriction enzyme and the reference genome.

##### Table File

A table file of chromosomes' size. This file can be easily found on the UCSC genome browser. Note: Pay attention to the contigs or scaffolds, and be aware that HiC-Pro will generate a map per chromosome pair.

##### Bowtie2 Indexes

Mentioned in the previous section.

#### Config File

```bash
# config-hicpro.txt
# Paths
OUTPUT_DIR = output_data
INPUT_DIR = input_data
GENOME_FASTA = reference_genome.fa
BOWTIE2_IDX_PATH = /path/to/bowtie2/index
# Other parameters
MIN_MAPQ = 30
```

Run HiC-Pro:

```bash
hic-pro -i input_data -o output_data -c config-hicpro.txt
```

|               | SYSTEM CONFIGURATION                                                       |
|---------------|----------------------------------------------------------------------------|
| PREFIX        | Path to installation folder                                                |
| BOWTIE2_PATH  | Full path the bowtie2 installation directory                               |
| SAMTOOLS_PATH | Full path to the samtools installation directory                           |
| R_PATH        | Full path to the R installation directory                                  |
| PYTHON_PATH   | Full path to the python installation directory                             |
| CLUSTER_SYS   | Scheduler to use for cluster submission. Must be TORQUE, SGE, SLURM or LSF |

## Step 5: Filter and Deduplicate (HiC-Pro)

```bash
hic-pro -i aligned_reads.sam -o output_data -c config_file
```

## Step 6: Generate Contact Matrix (HiC-Pro, Juicer, HiCExplorer)

```bash
hic-pro -i filtered_data -o contact_matrix -c config_file
```

## Step 7: Normalize Contact Matrix (HiC-Pro, Juicer, HiCExplorer)

```bash
hic-pro -i contact_matrix -o normalized_matrix -c config_file
```

## Step 8: Visualize Contact Map (Juicebox, HiCPlotter, HiCExplorer)

```bash
java -jar Juicebox.jar -g reference_genome.fa -n normalized_matrix
```

######################################################################

## from .cool

```mermaid

```

```bash
cd /code

# build .cool from bin & contact
# output: output.cool
python3 build-cooler.py

# Generate a multi-resolution cooler file by coarsening
# output: output.mcool
cooler zoomify output.cool

# build a contact map with .mcool
python3 mcool_map.py
```

cooler doc:
https://cooler.readthedocs.io/en/latest/cli.html#cooler-zoomify


## Reference

### SRA-Tools Docker

* [SRA-Tools GitHub](https://github.com/ncbi/sra-tools)
* [SRA-Tools Docker Hub](https://hub.docker.com/r/ncbi/sra-tools)
* [SRA-Tools Docker Wiki](https://github.com/ncbi/sra-tools/wiki/SRA-tools-docker)
* [Prefetch and Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
* [HowTo: Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)

