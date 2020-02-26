#!/bin/bash

ddir=../data
fq1=$ddir/read1.fq.gz
fq2=$ddir/read2.fq.gz

asm0=$ddir/asm0.fas
asm1=$ddir/asm1.fas
asm2=$ddir/asm2.fas

perl ../KADprofile.pl \
	--read $fq1 --read $fq2 \
	--asm $asm0 --asm $asm1 --asm $asm2 \
	--aid a0  --aid a1 --aid a2 \
	--prefix 1-KADprofile

