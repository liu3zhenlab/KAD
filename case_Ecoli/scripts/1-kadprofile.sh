#!/bin/bash

fq1=MG1655_1.fq.gz
fq2=MG1655_2.fq.gz

ref1=U00096.1.fasta
ref3=U00096.3.fasta

perl ~/scripts2/KAD/KADprofile.pl \
	--read $fq1 --read $fq2 \
	--asm $ref1 --asm $ref3 \
	--aid U00096.1  --aid U00096.3  \
	--prefix MG1655

