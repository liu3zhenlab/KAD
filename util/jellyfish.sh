#!/bin/bash
data=fasta
out=
$binPath/jellyfish count -m 25 -t 8 -o $out.jf -L 1 -s 5000M -C $data 
echo 
$binPath/jellyfish dump -c -T $out.jf > $out.txt
rm $out.jf

