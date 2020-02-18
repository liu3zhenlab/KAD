### case I: Assessing an E. coli genome assembly
An example to show how to use KAD scripts and what results are expected

#### step 1: data preparation
1. genome assemblies
For this example, we used two genome assemblies of *Escherichia coli* K-12 MG1655, which were from Genbank accession U00096.1 and U00096.3.
2. Illumina read data
Illumina data for the same bacterial strain were downloaded from NCBI with [sra-tools](https://github.com/ncbi/sra-tools].
```
fasterq-dump --skip-technical --split-3 ERR022075
```
3. trimming of read data with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```
#!/bin/bash
########### subject to change #############
trimmomaticDir=<trimmomatic directory>
trimmomatic=$trimmomaticDir/trimmomatic-0.38.jar
adaptor=$trimmomaticDir/adapters/TruSeq3-PE.fa
nthreads=4
read1=ERR022075_1.fastq
read2=ERR022075_2.fastq
prefix=MG1655
minlen=60
###########################################
java -jar $trimmomatic PE \
    -threads $nthreads \
    $read1 $read2 \
    -baseout $prefix \
    ILLUMINACLIP:$adaptor:3:20:10:1:true \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:13 \
    MINLEN:$minlen
```
4. Subsampling of read data with [seqtk](https://github.com/lh3/seqtk) to have ~45x sequencing depth
```
nreads=1000000
seqtk sample -s 100 MG1655_1P $nreads | gzip > MG1655_1.fq.gz
seqtk sample -s 100 MG1655_2P $nreads | gzip > MG1655_2.fq.gz
```


#### step 2: KAD profiling with [KADprofile.pl](KADprofile.pl)


#### step 3: KAD distribution with [KADdist.pl](KADdist.pl)  


#### step 4: KAD comparison with [KADcompare.pl](KADcompare.pl)

