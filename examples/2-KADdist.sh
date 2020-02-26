#!/bin/bash
perl /homes/liu3zhen/scripts2/KAD/KADdist.pl \
	--kad 1-KADprofile/1-KADprofile_4_kad.txt \
	--aid a0 \
	--winsize 500 \
	--asm ../data/asm0.fas \
	--prefix a0_dist

perl /homes/liu3zhen/scripts2/KAD/KADdist.pl \
	--kad 1-KADprofile/1-KADprofile_4_kad.txt \
	--aid a1 \
	--winsize 500 \
	--asm ../data/asm1.fas \
	--prefix a1_dist

perl /homes/liu3zhen/scripts2/KAD/KADdist.pl \
	--kad 1-KADprofile/1-KADprofile_4_kad.txt \
	--aid a2 \
	--winsize 500 \
	--asm ../data/asm2.fas \
	--prefix a2_dist

