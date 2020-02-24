#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-23:00:00
perl ~/scripts2/KAD/KADdist.pl --prefix U00096.1dist --kad ./MG1655/MG1655_4_kad.txt --aid U00096.1 --asm U00096.1.fasta
perl ~/scripts2/KAD/KADdist.pl --prefix U00096.3dist --kad ./MG1655/MG1655_4_kad.txt --aid U00096.3 --asm U00096.3.fasta
