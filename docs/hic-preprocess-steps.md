# Hi-C pre-processing steps

## Step 1: Download SRA file (SRA Toolkit)

I use docker in wsl to avoid any instalation problem:

```
docker pull ncbi/sra-tools
docker images
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp --name yy-sra ncbi/sra-tools
docker ps -a
docker attach <container_id>
```

In the container, use the sra-toolkit cli:

```
# check bin
which prefetch
/usr/local/bin/prefetch

prefetch SRR5579177
2024-12-26T13:38:50 prefetch.3.1.0:  HTTPS download succeed
2024-12-26T13:41:08 prefetch.3.1.0:  'SRR5579177' is valid
2024-12-26T13:41:08 prefetch.3.1.0: 1) 'SRR5579177' was downloaded successfully
```

## Step 2: Convert SRA to FASTQ (SRA Toolkit)

```
# with progress bars
fasterq-dump SRR5579177 -p
join   :|-------------------------------------------------- 100%   
concat :|-------------------------------------------------- 100%   
spots read      : 301,260,192
reads read      : 602,520,384
reads written   : 602,520,384
```

default is `--split-3`, which means:
The spots are split into ( biological ) reads, for each read : 4 lines of FASTQ or 2 lines of FASTA are written. For spots having 2 reads, the reads are written into the *_1.fastq and *_2.fastq files. Unmated reads are placed in *.fastq. If the accession has no spots with one single read, the *.fastq-file will not be created.

## Step 3: Quality control (FastQC)

### Quality Control with FastQC

Review FastQC Reports: Open the generated HTML reports to review the quality of your reads. Look for issues such as low-quality scores, adapter contamination, or overrepresented sequences.

install fastqc:

```
# clean env
docker pull ubuntu
docker images
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp/ --name yy-biofp ubuntu
docker attach yy-biofp

# in container
apt update
apt install fastqc
```

```
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

the outputs are 2 reports, I put it in results: fastqc-report

### Trimming and Filtering with Cutadapt 

cutadapt instalation:

```
apt update
apt install cutadapt

# check
which cutadapt
/usr/bin/cutadap
```

```
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20,20 -m 36 \
    -o trimmed_reads_SRR5579177_1.fastq -p trimmed_readsSRR5579177_2.fastq \
    SRR5579177_1.fastq SRR5579177_2.fastq
```

* `-a ADAPTER_FWD`: Adapter sequence for the forward read.
* `-A ADAPTER_REV`: Adapter sequence for the reverse read.
* -`q 20,20`: Trim low-quality bases from the ends of each read before adapter removal (Phred score < 20).
* `-m 36`: Discard reads shorter than 36 bases after trimming.
* `-o trimmed_reads_1.fastq`: Output file for trimmed forward reads.
* `-p trimmed_reads_2.fastq`: Output file for trimmed reverse reads.
* `reads_1.fastq reads_2.fastq`: Input FASTQ files.

### Run FastQC again on trimmed reads to verify quality

```
fastqc trimmed_reads_1_paired.fastq trimmed_reads_2_paired.fastq
```

## Step 4: Align reads (BWA-MEM, HiC-Pro)

### Index the Reference Genome

```
sudo apt-get install bwa
bwa index reference_genome.fa


bwa mem -SP5M reference_genome.fa reads_1.fastq reads_2.fastq > aligned_reads.sam
```
-S: Output in SAM format.
-P: Assume the input consists of interleaved paired-end sequences.
5M: Set the mismatch penalty to 5.
reference_genome.fa: The reference genome file.
reads_1.fastq reads_2.fastq: Input FASTQ files.

### Process Aligned Reads with HiC-Pro

```
# config-hicpro.txt
# Paths
OUTPUT_DIR = output_data
INPUT_DIR = input_data
GENOME_FASTA = reference_genome.fa
BOWTIE2_IDX_PATH = /path/to/bowtie2/index
# Other parameters
MIN_MAPQ = 30
```

```
hic-pro -i input_data -o output_data -c config-hicpro.txt
```

## Step 5: Filter and deduplicate (HiC-Pro)

```
hic-pro -i aligned_reads.sam -o output_data -c config_file
```

## Step 6: Generate contact matrix (HiC-Pro, Juicer, HiCExplorer)

```
hic-pro -i filtered_data -o contact_matrix -c config_file
```

## Step 7: Normalize contact matrix (HiC-Pro, Juicer, HiCExplorer)

```
hic-pro -i contact_matrix -o normalized_matrix -c config_file
```

## Step 8: Visualize contact map (Juicebox, HiCPlotter, HiCExplorer)

```
java -jar Juicebox.jar -g reference_genome.fa -n normalized_matrix
```

## Reference

### sra-tools docker

* https://github.com/ncbi/sra-tools
* https://hub.docker.com/r/ncbi/sra-tools
* https://github.com/ncbi/sra-tools/wiki/SRA-tools-dockerhttps://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
* https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
* 

###

