# Hi-C pre-processing steps

## Step 1: Download SRA file (SRA Toolkit)

```
prefetch SRR5579177
```

## Step 2: Convert SRA to FASTQ (SRA Toolkit)

```
fastq-dump --split-3 SRR5579177.sra
```

## Step 3: Quality control (FastQC)

### Quality Control with FastQC

```
fastqc reads_1.fastq reads_2.fastq
```

Review FastQC Reports: Open the generated HTML reports to review the quality of your reads. Look for issues such as low-quality scores, adapter contamination, or overrepresented sequences.

### Trimming and Filtering with Cutadapt 

```
cutadapt -a ADAPTER_FWD -A ADAPTER_REV \
    -q 20,20 -m 36 \
    -o trimmed_reads_1.fastq -p trimmed_reads_2.fastq \
    reads_1.fastq reads_2.fastq
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
