# TADs are 3D structural units of higher-order chromosome organization in Drosophila

### Contribute Members
* 李柏漢
* 林穎彥
* 黃宇秀
* 邱淦均

### Demo 

** For the full demo version, please check: docs/hic-preprocess-steps.md

#### Data Proccess

```bash
# sra-tool env
docker pull ncbi/sra-tools
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp --name yy-sra ncbi/sra-tools
docker attach yy-sra 

prefetch SRR5579177
fasterq-dump SRR5579177 -p

# other process env
docker pull ubuntu
docker images
docker run -itd -v /mnt/e/workspace/bio-fp/:/home/yy/bio-fp/ --name yy-biofp ubuntu
docker attach yy-biofp

apt install fastqc
fastqc reads_1.fastq reads_2.fastq
apt install cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20,20 -m 36 \
    -o trimmed_reads_SRR5579177_1.fastq -p trimmed_readsSRR5579177_2.fastq \
    SRR5579177_1.fastq SRR5579177_2.fastq

wget https://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz
gunzip dm3.fa.gz
bowtie-build ../dm3.fa dm3_index
bowtie dm3_index -1 ../SRR5579177_1.fastq -2 ../SRR5579177_2.fastq  -t -a -m 1 --best -S dm3_bowtie_align_output.sam

samtools view -bS dm3_bowtie_align_output.sam > output.bam
samtools sort output.bam -o sorted.bam
    
pairtools parse -c dm3.chrom.sizes -o output.pairsam dm3_bowtie_align_output.sam
pairtools sort -o sorted.pairsam output.pairsam
pairtools dedup -o dedup.pairsam sorted.pairsam
pairtools select '(pair_type == "UU")' -o output.pairs dedup.pairsam
```

#### Visulaize: R

```bash
cd code
Rscript contact_file_generate.R

Rscript contact_map_generate.R 
```

#### Visulaize: Python with .cool

```bash
python3 build-cooler.py

cooler zoomify output.cool

python3 mcool_map.py
```

## Folder organization and its related information

idea by Noble WS (2009) [A Quick Guide to Organizing Computational Biology Projects.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424) PLoS Comput Biol 5(7): e1000424.

### docs

* Presentation:, 1131_bioinformatics_FP_group1.pdf
  
* Related Document
  * docs/hic-preprocess-steps.md
  * data/data-src.md

### data (do not upload fastq file)

* Source
  * data/data-src.md (experiment data source details)
* Format
  * sra
  * fastq
  * fa (fasta)
  * sam
  * bam
  * txt
  * ebwt
* Size
  * SRA: ~ 15.3 GB
  * FASTQ: ~ 68 GB
  * FASTA (dm3): ~ 164 MB\
  * EBWT Bowtie Index: ~ 1 KB ~ 161 MB
  * SAM: ~ 115 GB
  * Sizes (dm3): ~ 1 KB
  * PairSAM: ~ 60.8 ~ 133 GB
  * Pairs: ~ 60.8 GB
  * BINS: ~ 332 KB
  

### code

* Which packages do you use? 
  * original packages in the paper
    * bowtie
    * 
  * additional packages you found
    * sra-tools 
    * wget
    * fastqc
    * cutadapt
    * reptyr
    * samtools
    * pairtools
    * R Lib: ggplot2
    * R Lib: reshape2
    * Cooler
    * Python Lib: cooler
    
* Analysis steps
  * Download SRA/Reference Genome Files
  * Convert SRA to FASTQ
  * FASTQ Quality Control
  * Build Bowtie Index
  * Trimming
  * Alignment
  * Build Pairs
  * Store SAM
  * Create Contact File
  * Build Contact Matrix
  * Visualize Contact Map
  
### results

* Reproduce Part
  * Figure 1A Hi-C Contact Map
  * result/res_10000_contact_heatmap.png
* QC Reports:
  * results/fastqc-report/SRR5579177_1_fastqc.html
  * results/fastqc-report/SRR5579177_2_fastqc.html
  
|                          | Information                                   |
|--------------------------|-----------------------------------------------|
| Organism                 | Drosophila melanogaster                       |
| Instrument Model         | Illumina HiSeq 2500                           |
| Mapping Genome           | Drosophila melanogaster genome (assembly dm3) |
| Data Processing Software | bowtie (params: -t -a -m 1 --best)            |

* Any improvement or change in your package?

## References

### Packages you use

#### SRA-Tools Docker

* [SRA-Tools GitHub](https://github.com/ncbi/sra-tools)
* [SRA-Tools Docker Hub](https://hub.docker.com/r/ncbi/sra-tools)
* [SRA-Tools Docker Wiki](https://github.com/ncbi/sra-tools/wiki/SRA-tools-docker)
* [Prefetch and Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
* [HowTo: Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)

#### Reference Genome
* [dm3 related](https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/)

#### Illumina

* [Illumina sequencing libraries](https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html)

#### Bowtie
* [bowtie1 comparison](https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html#differences-between-bowtie-and-bowtie2)
* [installation](https://www.metagenomics.wiki/tools/bowtie2/install)
* [bowtie2 docs](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner)
* [bowtie 1 docs](https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-inspect-index-inspector)

#### bwa

* [github](https://github.com/lh3/bwa#multihit)

#### Pairs

* [pairix](https://github.com/4dn-dcic/pairix)
* [pypairix](https://pypi.org/project/pypairix/0.1.0/)
* [pairtools](https://pairtools.readthedocs.io/en/latest/cli_tools.html)
* [micro-c](https://micro-c.readthedocs.io/en/latest/contact_map.html#from-pairs-to-cooler-contact-matrix)
* [pairs file format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)

#### Cool
* [cooler doc](https://cooler.readthedocs.io/en/stable/cli.html)
* [cooler github](https://github.com/open2c/cooler/tree/master/docs)

#### HIC-PRO
* [hic-pro docs](https://github.com/nservant/HiC-Pro/blob/master/doc)
* [github](https://github.com/nservant/HiC-Pro)

#### GenPipes
* [docs](https://genpipes.readthedocs.io/en/genpipes-v-3.6.2/user_guide/pipelines/gp_hicseq.html)
* Related publications
