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

my $version = 0.10;

sub prompt {
	print <<EOF;
	Usage KADdist.pl [options]
	[Options]
	--kad|k <file>:      KAD output file from seqKADprofile.pl; required
	--aid|i <str>:       assembly ID in the header of KAD file; required
	--asm|a <file>:      assembly FASTA file, including path; required
	--mincopy|m <num>: k-mers  with at least --mincopy in the assembly will be aligned to the assembly; default=1
	--maxcopy|i <num>: k-mers  with at most --maxcopy in the assembly will be aligned to the assembly; default=100
	--winsize|w <num>: window size on which the number of each KAD type is counted; default=50000
	--kadcutoff <str>: same to --kadcutoff in the seqKADprofile.pl
	                   default="-2 -0.5 0.5 0.75 2"
	--prefix|p <str>:  the output directory and the prefix for output files; default=KADdist
	--minwin4plot|n <num>: contigs or chromosomes with minimum window number (--minwin4plot) will be plotted; default=10
	--pdfoutdir|o <str>: the subdirectory under --prefix directory for PDF outputs
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
$kadcutoff = "-2 -0.5 0.5 0.75 2" if !defined $kadcutoff;
$winsize = 50000 if !defined $winsize;
$prefix = "KADdist" if !defined $prefix;
$minwin4plot = 10 if !defined $minwin4plot;
$pdfoutdir = "pdf" if !defined $pdfoutdir;

###############################################
# preparation
###############################################
if (-d $prefix) {
	print STDERR "error: the directory $prefix exists.\n";
	exit;
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
my @head = split("\t", $_);
my @selcols = ();
for (my $i=0; $i<=$#head; $i++) {
	if ($head[$i] =~ /$aid/) {
		push(@selcols, $i);
	}
}

if (!@selcols) {
	print STDERR "error: no columns were identified to match --aid\n";
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
`mkdir $bowtie_dbidx`;

# index asm
`$binPath/bowtie/bowtie-build $asm_file $bowtie_dbidx/$aid`;
# aln
`$binPath/bowtie/bowtie -n 0 -v 0 -k $maxcopy --quiet --no-unal -a -B 1 --sam --sam-nohead -f $bowtie_dbidx/$aid $fas_file | cut -f 1,3,4,6 > $kmeraln_file`;


###############################################
# generate a kad-wig file
###############################################
my $chrlen_file = $prefix."/".$prefix."_3_asm.lengths";
my $wig_file = $prefix."/".$prefix."_4_kad.wig";
my $bigwig_file = $prefix."/".$prefix."_4_kad.bigwig";

my $errors = &aln2kadwig_error($kmeraln_file, $wig_file);

# wig to bigwig
`$scriptPath/util/fastaSize.pl $asm_file > $chrlen_file`;
`$binPath/wigToBigWig $wig_file $chrlen_file $bigwig_file`;

close LOG;

###############################################
# error count (KAD=-1)
###############################################
my $errpos_file = $prefix."/".$prefix."_5_error.pos.txt";
open(ERRPOS, ">$errpos_file") || die;
my %errors = %{$errors};
foreach my $echr (sort {$a cmp $b} keys %errors) {
	my %echrpos = %{$errors{$echr}};
	foreach (sort {$a <=> $b} keys %echrpos) {
		print ERRPOS "$echr\t$_\n";
	}
}
close ERRPOS;

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
## module to generate a wig file and error data
################################################
sub aln2kadwig_error {
	my ($in_aln, $out_wig) = @_;
	my %pos_nkad;
	my %pos_kadsum;
	my %error_pos;

	# read aln and analyze KAD at each pos
	open(INALN, $in_aln) || die;
	while(<INALN>) {
		chomp;
		my ($inkinfo, $inchr, $inpos, $nummatch) = split;
		my @kmer_kad = split("_", $inkinfo);
		my $inkad = $kmer_kad[2];
		$nummatch =~ s/M//g;
		for (my $cpos = $inpos; $cpos <= $inpos + $nummatch - 1; $cpos++) {
			$pos_nkad{$inchr}{$cpos}++;
			if (exists $pos_kadsum{$inchr."_".$cpos}) {
				$pos_kadsum{$inchr."_".$cpos} += $inkad;
			} else {
				$pos_kadsum{$inchr."_".$cpos} = $inkad;
			}
		}
		
		# errors
		if ($inkad == -1) {
			$error_pos{$inchr}{$inpos}++;
		}

	}
	close INALN;

	# output a WIG file
	open(WIG, ">$out_wig") || die;
	foreach my $echr (sort { $a cmp $b } keys %pos_nkad) {
		print WIG "variableStep  chrom=$echr\n";
		my %echrPos = %{$pos_nkad{$echr}};
		foreach my $epos (sort {$a <=> $b} keys %echrPos) {
			my $kad_mean = $pos_kadsum{$echr."_".$epos} / $echrPos{$epos};
			print WIG "$epos\t$kad_mean\n";
		}
	}
	close WIG;
	return(\%error_pos);
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
		if ($inkad <= $inkad_cutoff[0]) {
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


