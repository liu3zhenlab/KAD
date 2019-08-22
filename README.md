# KAD
Assessing genome assemblies by comparing k-mer copies in assemblies and K-mer abundance in reads

### Introduction
KAD is designed for evaluating the accuracy of nucleotide base quality of genome assemblies. Briefly, abundance of k-mers are quantified for both sequencing reads and assembly sequences. Comparison of the two values results in a single value per k-mer, K-mer Abundance Difference (KAD), which indicates how well the assembly matches read data for each k-mer.


<img src="https://latex.codecogs.com/svg.latex?\Large&space;KAD=log2\frac{c+m}{m*(n+1)}" />

where, _c_ is the count of a k-mer from reads, _m_ the mode of counts of read k-mers, _n_ is the copy of the k-mer in the assembly. 

### Requirements
The script was written with Perl and R is invoked. Both Perl and R are generally installed. If needed, please refer to [Perl](https://www.perl.org/) and [R](https://www.r-project.org/) for installation guides. To generate reports, the R packages [knitr](https://github.com/yihui/knitr) and [rmarkdown](https://rmarkdown.rstudio.com) are needed to be installed.

[Jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) was used to generated k-mers for either FASTA or FASTQ data. The binary executable is in the _bin_ directory of the KAD package.

### Installation
git clone https://github.com/liu3zhenlab/KAD.git  
perl ./KAD/seqKADprofile.pl

### Running guide
1. **seqKADprofile.pl**: producing KAD profiles for input assemblies.
Usage: perl seqKADprofile.pl [options]
[Options]
--read <file>:	FASTQ/A read file for k-mer generation; the parameter can be used multiple times; zip files with the suffix of .gz are allowed;  
--minc <num>:	minimal number of counts per k-mer from reads; k-mers with counts smaller than <num> are not output. default=5.  
--asm <file>:	FASTA sequence file for k-mer generation; the parameter can be used multiple times to allow using multiple FASTA file;
					each file is considered an indepedent assembly.  
--rid <str>:	ID used in the header of the k-mer table generated from reads.  
--aid <str>:	ID used in the header of the asm k-mer table to be generated from each assembly; the parameter can be used multiple times to match --asm input. By default, a header ID is generated from the file name of each assembly by removing PATH and the suffix of .fa, .fas, or .fasta.  
--prefix <str>:the output directory and the prefix for output files; default=kad.  
--klen <num>:  length of k-mers; default=25.  
--readdepth <num>: estimated depth of reads; not required; if specified, it will be compared to the mode of read k-mers.  
--kadcutoff <str of nums>: a set of numbers to define k-mer categories; default="-2 -0.5, 0.5, 0.75, 2".  
--binlen <num>:		bin length to count KAD; default=0.05.  
--threads <num>:		number of cpus; default=1.  
--version:		version  
--help:			help information

2. **seqKADprofile.pl**: producing KAD profiles for input assemblies.

### Walk-through example
Let us say you have three assembly versions:
1. asm1.fas
2. asm2.fas
3. asm3.fas

You also have a read set:
1. read1.fq
2. read2.fq

Assuming the Perl script was in the directory of _scriptpath_, run the following script to generate KAD profiles for all three assemblies.
```
perl scriptpath/seqKADprofile.pl --read read1.fq --read read2.fq \
                                 --asm asm1.fas --asm asm2.fas --asm asm3.fas
```

You might want to assign names all three assemblies with a1, a2, and a3.
```
perl scriptpath/seqKADprofile.pl --read read1.fq --read read2.fq \
                                 --asm asm1.fas --asm asm2.fas --asm asm3.fas \
                                 --aid a1 --aid a2 --a3
```
You need to be carefully use --aid, which must matches with --asm order.



