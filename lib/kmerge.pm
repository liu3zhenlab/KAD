#!/usr/bin/perl -w
# ===============================================================
# kmerge.pm
# Sanzhen Liu
# 8/9/2019
# ===============================================================

package kmerge;

sub kmerge {
	use strict;
	use warnings;
	use File::Temp;
	use Getopt::Long;

	my ($kmer01, $kmer02, $merge) = @_;
	
	# ouput file
	open(OUT, ">$merge") || die;

	### headers:
	my $header01;
	chomp($header01 = `head $kmer01 -n 1`);
	my $header02;
	chomp($header02 = `head $kmer02 -n 1`);
	my @header01 = split(/\t/, $header01);
	my $ncol01 = $#header01;
	my @header02 = split(/\t/, $header02);
	my $ncol02 = $#header02;

	# query header
	my $header01_add = join("\t", @header01[0..$ncol01]);
	print OUT "$header01_add";

	# target header
	my $header02_add = join("\t", @header02[1..$ncol02]);
	print OUT "\t$header02_add\n";
	
	close OUT;

	# missing values
	my $mvalue = 0;
	my $missing_value01 = $mvalue;
	if ($ncol01 > 1) { 
		for (my $i = 2; $i <= $ncol01; $i++) {
			$missing_value01 .= "\t$mvalue";	
		}
	}

	my $missing_value02 = $mvalue;
	if ($ncol02 > 1) {
		for (my $i = 2; $i <= $ncol02; $i++) {
			$missing_value02 .= "\t$mvalue";
		}	
	}

	###################################################
	### divide a large file into small pieces for efficient merging
	###################################################
	my @nt = ("A", "T", "C", "G");
	my @indexnt = @nt;
	my (@newindex, %kmerfiles, $fh01in, $fh02in);
	for (my $i = 1; $i < 3; $i++) {
		foreach my $nt1 (@indexnt) {
			if (!exists $kmerfiles{1}{$nt1}) {
				$fh01in = $kmer01;
			} else {
				$fh01in = $kmerfiles{1}{$nt1};
			}

			if (!exists $kmerfiles{2}{$nt1}) {
				$fh02in = $kmer02;
			} else {
				$fh02in = $kmerfiles{2}{$nt1};
			}

			my $fh01out = File::Temp->new(SUFFIX => '.kmertmp');
			`grep "^$nt1" $fh01in > $fh01out`;
			my $fh02out = File::Temp->new(SUFFIX => '.kmertmp');
			`grep "^$nt1" $fh02in > $fh02out`;

			foreach my $nt2 (@nt) {
				my $comb = $nt1.$nt2;
				my $fh01out2 = File::Temp->new(SUFFIX => '.kmertmp');
				`grep "^$comb" $fh01out > $fh01out2`;
				$kmerfiles{1}{$comb} = $fh01out2;
			
				my $fh02out2 = File::Temp->new(SUFFIX => '.kmertmp');
				`grep "^$comb" $fh02out > $fh02out2`;
				$kmerfiles{2}{$comb} = $fh02out2;
			
				push(@newindex, $comb);
			}
		}
		@indexnt = @newindex;
		@newindex = ();
	}

	foreach (@indexnt) {
	# merge two sets of kmer counts
		my $infile01 = $kmerfiles{1}{$_};
		my @k1 = `grep $_ $infile01`;
		my $infile02 = $kmerfiles{2}{$_};
		my @k2 = `grep $_ $infile02`;
		&merge(\@k1, \@k2, $merge, $missing_value01, $missing_value02);
	}
}


sub merge {
# merge each subset:
	my ($inkm01, $inkm02, $inmerge, $inmvalue1, $inmvalue2) = @_;
	my @km01 = @{$inkm01};
	my @km02 = @{$inkm02};
	my %kmer_hash01 = ();
	my %kmer_hash02 = ();

	foreach (@km01) {
		chomp;
		my @kma = split(/\t/, $_);
		my $km = $kma[0];
		my $kmac = join("\t", @kma[1..$#kma]);
		$kmer_hash01{$km} = $kmac;
	}

	### adding data:
	my $kmbc;
	foreach (@km02) {
		chomp;
		my @kmb = split(/\t/, $_);
		my $km = $kmb[0];
		$kmbc = join("\t", @kmb[1..$#kmb]);

		if (exists $kmer_hash01{$km}) {
			$kmer_hash02{$km} = $kmbc;
		} else {
			$kmer_hash01{$km} = $inmvalue1;
			$kmer_hash02{$km} = $kmbc;
		}
	}
	
	### output:
	open(OUT, ">>$inmerge") || die;
	foreach my $eachkmer (sort {$a cmp $b} keys %kmer_hash01) {
		print OUT "$eachkmer";
		print OUT "\t$kmer_hash01{$eachkmer}";
		if (exists $kmer_hash02{$eachkmer}) {
			print OUT "\t$kmer_hash02{$eachkmer}";
		} else {
			print OUT "\t$inmvalue2";
		}
		print OUT "\n";
	}
	close OUT;
}
###################################################

1;

