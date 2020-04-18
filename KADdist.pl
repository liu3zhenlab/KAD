#!/usr/bin/perl -w
#
#=========================================================
# KADdist.pl
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 8/31/2019
#
# ========================================================
use strict;
use warnings;
use Getopt::Long;
use FindBin;

my $version = 0.1.6;

# 4/17/2020: fix the running failure due to shared assembly names on the KAD file

sub prompt {
	print <<EOF;
	Usage KADdist.pl [options]
	[Options]
	--kad|k <file>:      KAD output file from KADprofile.pl (required)
	--aid|i <str>:       assembly ID in the header of KAD file (required)
	--asm|a <file>:      assembly FASTA file, including path (required)
	--mincopy|m <num>: k-mers  with at least --mincopy in the assembly will be aligned to the assembly (1)
	--maxcopy|i <num>: k-mers  with at most --maxcopy in the assembly will be aligned to the assembly (100)
	--winsize|w <num>: window size on which the number of each KAD type is counted (50000)
	--kadcutoff|s <str>: same to --kadcutoff in the KADprofile.pl ("-0.8 -0.5 0.5 0.75 2")
	--prefix|p <str>:  the output directory and the prefix for output files (KADdist)
	--minwin4plot|n <num>: contigs or chromosomes with minimum window number (--minwin4plot) will be plotted (10)
	--pdfoutdir|o <str>: the subdirectory under --prefix directory for PDF outputs (pdf)
	--help
EOF
exit;
}

###############################################
# parameters
###############################################
my ($kad_file, $aid, $asm_file, $mincopy, $maxcopy,
    $prefix, $winsize, $kadcutoff, 
	$minwin4plot, $pdfoutdir, $help);
&GetOptions("kad|k=s" => \$kad_file,
			"aid|i=s" => \$aid,
            "asm|a=s" => \$asm_file,
			"mincopy|m=i" => \$mincopy,
			"maxcopy|x=i" => \$maxcopy,
			"winsize|w=i" => \$winsize,
			"kadcutoff=s" => \$kadcutoff,
			"prefix|p=s" => \$prefix,
			"minwin4plot|n=i" => \$minwin4plot,
			"pdfoutdir|o=s" => \$pdfoutdir,
			"help|h" => \$help) || &prompt;
&prompt if ($help or !defined $kad_file or !defined $aid or !defined $asm_file);
$mincopy = 1 if !defined $mincopy;
$maxcopy = 100 if !defined $maxcopy;
$kadcutoff = "-0.8 -0.5 0.5 0.75 2" if !defined $kadcutoff;
$winsize = 50000 if !defined $winsize;
$prefix = "KADdist" if !defined $prefix;
$minwin4plot = 10 if !defined $minwin4plot;
$pdfoutdir = "pdf" if !defined $pdfoutdir;

###############################################
# preparation
###############################################
if (-d $prefix) {
	print STDERR "error: the directory $prefix exists.\n";
	#exit;
} else {
	`mkdir $prefix`;
}

# script path:
my $scriptPath = $FindBin::Bin;
my $binPath = $scriptPath."/bin/";
#
# # log file
my $logfile = $prefix."/".$prefix."_KADdist.log";
open(LOG, ">$logfile") || die;

###############################################
# k-mers to fasta
###############################################
my $fas_file = $prefix."/".$prefix."_1_kmer.fasta";
open(FAS, ">$fas_file") || die;

open(KAD, $kad_file) || die;
# header and determine select cols:
$_ = <KAD>;
chomp;
my @head = split("\t", $_);
my @selcols = ();
my $kad_id = $aid.".KAD";
for (my $i=0; $i<=$#head; $i++) {
	my $head_value = $head[$i];
	$head_value =~ s/\.KAD$//g;
	if ($head_value eq $aid) {
		push(@selcols, $i);
	}
}

if (!@selcols) {
	print STDERR "error: no columns were identified to match --aid\n";
	exit;
} elsif ($#selcols == 0) {
	print STDERR "error: only 1 column was identified to match --aid\n";
	print STDERR "       there should be one column for --aid and one for its KAD values\n";
	exit;
} elsif ($#selcols > 1) {
	print STDERR "error: >2 columns were identified to match --aid\n";
	exit;
}

# k-mers to fasta
my $rnum = 0;
while(<KAD>) {
	chomp;
	$rnum++;
	my @line = split;
	my $kmer = $line[0];
	my $kinfo = join("_", @line[@selcols]);
	my $asmcopy =  @line[$selcols[0]];
	my $asmkad =  @line[$selcols[1]];
	if ($asmcopy >= $mincopy and $asmcopy <= $maxcopy) {
		print FAS ">k$rnum";
		print FAS "_$kinfo\n";
		print FAS "$kmer\n";
	}
}

close KAD;
close FAS;

###############################################
# aln k-mers to ref
###############################################
my $kmeraln_file = $prefix."/".$prefix."_2_kmer.aln";
my $bowtie_dbidx = $prefix."/bowtie";
if (! -d $bowtie_dbidx) {
	`mkdir $bowtie_dbidx`;
}

# bowtie path:
my $bowtiebuild = `which bowtie-build || echo $binPath/bowtie/bowtie-build`;
chomp $bowtiebuild;
my $bowtie = `which bowtie || echo $binPath/bowtie/bowtie`;
chomp $bowtie;

# index asm
`$bowtiebuild $asm_file $bowtie_dbidx/$aid`;
# aln
`$bowtie -n 0 -v 0 -k $maxcopy --quiet --no-unal -a -B 1 --sam --sam-nohead -f $bowtie_dbidx/$aid $fas_file | cut -f 1,3,4,6 > $kmeraln_file`;


###############################################
# generate bed from bowtie alignments
###############################################
my $kmerkad_bedfile = $prefix."/".$prefix."_3_kmer.kad.bed";
&aln2bed($kmeraln_file, $kmerkad_bedfile);

###############################################
# error count (KAD=-1)
###############################################
my $errpos_file = $prefix."/".$prefix."_4_error.pos.bed";
my $errpos_merge_file = $prefix."/".$prefix."_4_error.pos.merge.bed";
`awk '{ if (\$5 == -1) print }' $kmerkad_bedfile | cut -f 1-3 > $errpos_file`;
`bedtools merge -i $errpos_file > $errpos_merge_file`;

###############################################
# generate a kad-wig file
###############################################
my $chrlen_file = $prefix."/".$prefix."_5_asm.lengths";
my $wig_file = $prefix."/".$prefix."_5_kad.wig";
my $bigwig_file = $prefix."/".$prefix."_5_kad.bigwig";

&bed2kadwig($kmerkad_bedfile, $wig_file);

# chrlen
`$scriptPath/util/fastaSize.pl $asm_file > $chrlen_file`;
`$binPath/wigToBigWig $wig_file $chrlen_file $bigwig_file`;

close LOG;

###############################################
# aln to dist
###############################################
my $dist_file = $prefix."/".$prefix."_6_KADtype.dist.txt";
my @kadtype = qw(good error overrep lounrep hiunrep others);
open(DIST, ">$dist_file") || die;
print DIST "chr\twindow\tstart\tend\t";
print DIST join("\t", @kadtype);
print DIST "\n";
my ($dist, $outchrsize) = &aln2dist($chrlen_file, $winsize, $kmeraln_file, $kadcutoff);
my %dist = %{$dist};
my %outchrsize = %{$outchrsize};
foreach my $echr (sort {$a cmp $b} keys %dist) {
	my %echrdist = %{$dist{$echr}};
	foreach my $ewin (sort {$a <=> $b} keys %echrdist) {
		print DIST "$echr\t$ewin\t";
		my $winstart = $ewin * $winsize + 1;
		my $winend = ($ewin + 1) * $winsize;
		if ($winend > $outchrsize{$echr}) {
			$winend = $outchrsize{$echr};
		}
		print DIST "$winstart\t$winend";
		for (my $i=0; $i<=$#kadtype; $i++) {
			print DIST "\t$echrdist{$ewin}{$kadtype[$i]}";
		}
		print DIST "\n";
	}
}
close DIST;

###############################################
# plot KAD distribution on contigs or chromosomes
###############################################
my $pdfdir = $prefix."/".$pdfoutdir;
if (-d $pdfdir) {
	`rm $pdfdir -rf`;
}
`mkdir $pdfdir`;

`Rscript $scriptPath/util/kaddistPlot.R $dist_file $minwin4plot $pdfdir`;

###############################################
## module to bed
################################################
sub aln2bed {
	my ($in_aln, $out_bed) = @_;
	my %pos_kad;

	# read aln and analyze KAD at each pos
	open(INALN, "<", $in_aln) || die;
	while(<INALN>) {
		chomp;
		my ($inkinfo, $inchr, $inpos, $nummatch) = split;
		my @kmer_kad = split("_", $inkinfo);
		my $inkad = $kmer_kad[2];
		$nummatch =~ s/M//g;
		my $instart = $inpos - 1; # 0-based
		$pos_kad{$inchr}{$instart} = $nummatch."_".$inkad;
	}
	close INALN;

	# output
	open(BED, ">", $out_bed) || die; # sorted BED
	foreach my $inec (sort {$a cmp $b} keys %pos_kad) {
		my %ec_pos_kad = %{$pos_kad{$inec}};
		foreach my $estart (sort {$a <=> $b} keys %ec_pos_kad) {
			my $match_kad = $ec_pos_kad{$estart};
			my ($nummatch_value, $kad_value) = split("_", $match_kad);
			my $eend = $estart + $nummatch_value;
			print BED "$inec\t$estart\t$eend\t\.\t$kad_value\n";
		}
	}
	close BED;
}

###############################################
## module to generate a wig file and error data
################################################
sub bed2kadwig {
	my ($in_bed, $out_wig) = @_;
	my %pos_nkad;
	my %pos_kadsum;
	my %inctghash;
	my $curctg;
	
	# WIG output
	open(WIG, ">$out_wig") || die;

	# read aln and analyze KAD at each pos
	open(BED, $in_bed) || die;
	while(<BED>) {
		chomp;
		my @inline = split;
		my ($inctg, $instart, $inend, $inkad) = @inline[0,1,2,4];
		$instart += 1; # adjust to 1-based

		if (defined $curctg and !exists $inctghash{$inctg}) {
			print WIG "variableStep  chrom=$curctg\n";
			foreach my $epos (sort {$a <=> $b} keys %pos_kadsum) {
				my $kad_mean = $pos_kadsum{$epos} / $pos_nkad{$epos};
				print WIG "$epos\t$kad_mean\n";
			}
			%pos_nkad = ();
			%pos_kadsum = ();
			$curctg = "";
		}

		$curctg = $inctg;
		$inctghash{$curctg}++;
		for (my $cpos = $instart; $cpos <= $inend; $cpos++) {
			$pos_nkad{$cpos}++;
			if (exists $pos_kadsum{$cpos}) {
				$pos_kadsum{$cpos} += $inkad;
			} else {
				$pos_kadsum{$cpos} = $inkad;
			}
		}
	}
	close BED;
	# last one
	print WIG "variableStep  chrom=$curctg\n";
	foreach my $epos (sort {$a <=> $b} keys %pos_kadsum) {
		my $kad_mean = $pos_kadsum{$epos} / $pos_nkad{$epos};
		print WIG "$epos\t$kad_mean\n";
	}
	%pos_nkad = ();
	%pos_kadsum = ();
	close WIG;
}

###############################################
## module to generate distribution of 4 categories
################################################
sub aln2dist {
	my ($in_len, $in_win, $in_aln, $inkad_cutoff) = @_;
	my @inkad_cutoff = split(/[ ,]+/, $inkad_cutoff);
	my %wins;
	my %inchrsize;
	my @kadtype = qw(good error overrep lounrep hiunrep others);
	# chrlen to windows
	open(INLEN, $in_len) || die;
	while (<INLEN>) {
		chomp;
		my ($inchr, $inlen) = split;
		$inchrsize{$inchr} = $inlen;
		my $max_win = int(($inlen - 1) / $in_win);
		for (my $i = 0; $i <= $max_win; $i++) {
			foreach (@kadtype) {
				$wins{$inchr}{$i}{$_} = 0;
			}
		}
	}
	close INLEN;

	# read aln and analyze KAD at each pos
	open(INALN, $in_aln) || die;
	while(<INALN>) {
		chomp;
		my ($inkinfo, $inchr, $inpos, $nummatch) = split;
		my @kmer_kad = split("_", $inkinfo);
		my $inkad = $kmer_kad[2];
		my $cwin = int(($inpos - 1) / $in_win);
		my $type = "others";
		if ($inkad <= $inkad_cutoff[0] and $inkad != -1) {
			$type = "overrep";
		} elsif ($inkad >= $inkad_cutoff[1] and $inkad <= $inkad_cutoff[2]) {
			$type = "good";
		} elsif ($inkad == -1) {
			$type = "error";
		} elsif ($inkad >= $inkad_cutoff[3] and $inkad < $inkad_cutoff[4]) {
			$type = "lounrep";
		} elsif ($inkad >= $inkad_cutoff[4]) {
			$type = "hiuprep";
		}
		$wins{$inchr}{$cwin}{$type}++;
	}
	close INALN;
	return(\%wins, \%inchrsize);
}

# cleanup
`rm $bowtie_dbidx -rf`;
#`rm $chrlen_file`;

