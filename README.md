[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/7YvIlKos)
# TADs are 3D structural units of higher-order chromosome organization in Drosophila

### Members
* 李柏漢, 113753218
* 林穎彥, 113971012
* name, student ID2
* name, student ID3
* ...

### Demo 

You might provide an example command or a few commands to reproduce your analysis, i.e., the following R script
```R
Rscript code/your_script.R --input data/training --output results/performance.tsv
```

## Folder organization and its related information

idea by Noble WS (2009) [A Quick Guide to Organizing Computational Biology Projects.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424) PLoS Comput Biol 5(7): e1000424.

### docs

* Your presentation, 1131_bioinformatics_FP_groupID.ppt/pptx/pdf (i.e.,1131_bioinformatics_FP_group1.ppt), by **01.02 9:20 am**
  
* Any related document for the project
  * i.e., software user guide
  * hic-preprocess-steps.md

### data (do not upload fastq file)

* Source
  * data/data-src.md (experiment source details)
* Format
  * sra
  * fastq
  * fa (fasta)
  * sam
  * txt
* Size
  * sra: approx. 15.4 GB
  * fastq: approx. 68 GB
  * 

### code

* Which packages do you use? 
  * original packages in the paper
    * sra-tools
    * bowtie2 (bwa)
    * 
  * additional packages you found
    * fastqc
    * cutadapt
    * reptyr
    * hic-pro
    * 
* Analysis steps
  * 

### results

* Which part of the paper do you reproduce?
  * <experiment title>

|                          | Information                                   |
|--------------------------|-----------------------------------------------|
| Organism                 | Drosophila melanogaster                       |
| Instrument Model         | Illumina HiSeq 2500                           |
| Mapping Genome           | Drosophila melanogaster genome (assembly dm3) |
| Data Processing Software | bowtie (params: -t -a -m 1 --best)            |

* Any improvement or change in your package?

## References

* Packages you use
  * SRA-Tools Docker
    * [SRA-Tools GitHub](https://github.com/ncbi/sra-tools)
    * [SRA-Tools Docker Hub](https://hub.docker.com/r/ncbi/sra-tools)
    * [SRA-Tools Docker Wiki](https://github.com/ncbi/sra-tools/wiki/SRA-tools-docker)
    * [Prefetch and Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
    * [HowTo: Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
  * 
* Related publications
