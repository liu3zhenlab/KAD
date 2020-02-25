#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($fas, $wig, $chr, $start, $end, $help);

sub prompt {
	print <<EOF;
	Usage: perl wig.region.pl --wig <wig file> --chr <str> --start <num> --end <num>
	Output: BED formatted regional data
	--wig <file>:  wig file, required
	--chr <str>:   chromosome or contig ID, required
	--start <num>: start position (1-based), required
	--end <num>:   end position (1-based), required
	--fas <file>:  fasta file, required if sequence bases will be added
	--help
EOF
exit;
}

# read the parameters:
&GetOptions("fas|f=s" => \$fas,
            "wig|w=s" => \$wig,
			"chr|c=s" => \$chr,
            "start|s=i" => \$start,
			"end|e=i" => \$end,
			"help" => \$help);

&prompt if ($help);
if (!defined $wig or !defined $chr or !defined $start or !defined $end) {
	&prompt;
}

# sequencing data
my ($seq_name, $seq, $extract);
my $extract_done = 0;

if (defined $fas) {
	open(IN, "<", $fas) || die;
	while (<IN>) {
		$_ =~ s/\R//g;
		chomp;
		if (/^>(\S+)/) {
			if (defined $seq_name and $seq_name eq $chr) {
				$extract = substr($seq, $start - 1, $end - $start + 1);
				$extract_done = 1;
				last;
			}
			$seq_name = $1;
			$seq = "";
		} else {
    		$_ =~ s/\s//g; # remove spaces
		  	$seq .= $_;
		}	
	}
	# last element:
	if (!$extract_done and defined $seq_name and $seq_name eq $chr) {
		$extract = substr($seq, $start - 1, $end - $start + 1);    
	}
	close IN;
}
# wig region
my $wig_seqname;
my $start2read = 0;
open(WIG, "<", $wig) || die;
while (<WIG>) {
	chomp;
	#variableStep  chrom=chr1
	if (/^variableStep.*chrom=(\S+)/) {	
		if (!$start2read) {
			$wig_seqname = $1;
			if ($wig_seqname eq $chr) {
				$start2read = 1;
			}
		} else {
			last;	
		}
	} elsif ($start2read) {
		my @line = split;
		my $pos = $line[0];
		my $value = $line[1];
		if ($pos >= $start and $pos<= $end) {
			my $zero_based_pos = $pos - 1;
			print "$chr\t$zero_based_pos\t$pos";
			if (defined $extract) {
				my $base = substr($extract, $pos - $start, 1);
				print "\t$base";
			} else {
				print "\t.";
			}
			print "\t$value\n";
		} elsif ($line[0] > $end) {
			last;
		}
	}
}
close WIG;

