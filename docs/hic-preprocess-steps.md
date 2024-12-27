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

### Run FastQC Again on Trimmed Reads to Verify Quality

```bash
fastqc trimmed_reads_1_paired.fastq trimmed_reads_2_paired.fastq
```

## Step 4: Align Reads (BWA-MEM, HiC-Pro)

### Index the Reference Genome

```bash
sudo apt-get install bwa
bwa index reference_genome.fa

bwa mem -SP5M reference_genome.fa reads_1.fastq reads_2.fastq > aligned_reads.sam
```

* `-S`: Output in SAM format.
* `-P`: Assume the input consists of interleaved paired-end sequences.
* `5M`: Set the mismatch penalty to 5.
* `reference_genome.fa`: The reference genome file.
* `reads_1.fastq reads_2.fastq`: Input FASTQ files.

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
| ------------- | -------------------------------------------------------------------------- |
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

## Reference

### SRA-Tools Docker

* [SRA-Tools GitHub](https://github.com/ncbi/sra-tools)
* [SRA-Tools Docker Hub](https://hub.docker.com/r/ncbi/sra-tools)
* [SRA-Tools Docker Wiki](https://github.com/ncbi/sra-tools/wiki/SRA-tools-docker)
* [Prefetch and Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
* [HowTo: Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)

