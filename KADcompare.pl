#!/use/bin/perl -w
#
#======================================================================
# KADcompare.pl
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 8/17/2019
#
# The script is to compare KADs from two assemblies
#======================================================================

use strict;
use warnings;
use Getopt::Long;
use FindBin;

my $version = "0.10";

sub prompt {
	print <<EOF;
	Usage: perl KADcompare.pl [options] <kad>
	[Options]
	--set1:		*assembly ID 1 in the <kad> file; required
	--set2:		*assembly ID 2 in the <kad> file; required
	--prefix:	prefix for all outputs; including output directory and files; default=KC
	--binlen:	bin length to count KAD; default=0.05;
	--force:	overwrite the existing directory if specified; if not, quit if the output directory exists.
	--help: 	help information
EOF
exit;
}

###############################################
# parameters:
###############################################
my %opts = ();
my ($set1, $set2, $prefix, $binlen, $force, $help);
&GetOptions("set1=s" => \$set1,
            "set2=s" => \$set2,
			"prefix=s" => \$prefix,
			"binlen=i" => \$binlen,
			"force" => \$force,
			"help" => \$help);

&prompt if $help or @ARGV == 0;
$prefix = "KC" if (!defined $prefix); 
$binlen = 0.05 if (!defined $binlen);

# log file:
my $logfile = $prefix."_KADcompare.log";
open(LOG, ">$logfile") || die;

###############################################
# preparation
###############################################

# create a directory for outputs
if (-d $prefix) {
	if (!$force) {
		print LOG "ERR: Run stopped due to existing output directory: $prefix\n";
		exit;
	} else {
		print LOG "WARNING: the directory $prefix exists.\n";
	}
} else {
	`mkdir $prefix`;
}

# script path:
my $scriptPath = $FindBin::Bin;
#my $binPath = $scriptPath."/bin/";
#my $utilPath = $scriptPath."/util/";


###############################################
# diff data
###############################################
my $kad_file = $ARGV[0];
open(IN, $ARGV[0]) || die;
$_ = <IN>; # header
chomp;
my @header = split;
my $kad_count = 0;

foreach my $ehead (@header) {
	if ($ehead =~ /KAD$/) {
		$kad_count++;
	}
}

my $ncol1 = 0;
my $ncol2 = 0;
my $kadcol1 = 0;
my $kadcol2 = 0;

if ($kad_count < 2) {
	print STDERR "Less than 2 KAD columns. Check the input $ARGV[0]\n"; 
} elsif ($kad_count >= 2) {
	if (!defined $set1 or !defined $set2) {
		print STDERR "Both --set1 and --set2 must be provided\n";
	} else {
		for (my $i=2; $i<=$#header; $i++) {
			if ($header[$i] eq $set1) {
				$ncol1 = $i;
			} elsif ($header[$i] eq $set2) {
				$ncol2 = $i;
			} elsif ($header[$i] eq $set1.".KAD") {
				$kadcol1 = $i;
			} elsif ($header[$i] eq $set2.".KAD") {
				$kadcol2 = $i;
			}
		}
	}
}

# diff KAD output file
my $kaddiff_file = $prefix."/".$prefix."_1_".$set1."-".$set2.".KADdiff.txt";
open(KD, ">$kaddiff_file") || die;
print KD join("\t", @header[0,1,$ncol1,$ncol2,$kadcol1,$kadcol2]);
print KD "\n";

# binhash
my %allbinnum;

# from the row after the header
while(<IN>) {
	chomp;
	my @line = split;

	# change NA to 0
	$line[$kadcol1] = 0 if $line[$kadcol1] eq "NA";
	$line[$kadcol2] = 0 if $line[$kadcol2] eq "NA";

	# output:
	if ($line[$ncol1] != $line[$ncol2]) {
		my $kadbinnum1 = &num2bin($line[$kadcol1], $binlen);
		my $kadbinnum2 = &num2bin($line[$kadcol2], $binlen);
		$allbinnum{$kadbinnum1}{1}++;
		$allbinnum{$kadbinnum2}{2}++;
		print KD join("\t", @line[0,1,$ncol1,$ncol2,$kadcol1,$kadcol2]);
		print KD "\n";
	}

}
close IN;
close KD;


# bincount output file
my $bincount_file = $prefix."/".$prefix."_2_".$set1."-".$set2.".bincount.txt";
open(BIN, ">$bincount_file") || die;
print BIN "BIN\t$set1\t$set2\n";

foreach (sort {$a <=> $b} keys %allbinnum) {
	my $binc1 = 0;
	my $binc2 = 0;
	$binc1 = $allbinnum{$_}{1} if exists $allbinnum{$_}{1};
	$binc2 = $allbinnum{$_}{2} if exists $allbinnum{$_}{2};
	print BIN "$_\t$binc1\t$binc2\n";
}
close BIN;

###############################################
# KAD diff profiles
###############################################
`Rscript $scriptPath/util/KADcompare.plot.R $bincount_file $prefix`;

# Rmd
my $rmd = $scriptPath."\/util\/KADcompare.report.Rmd";

# report
my $reportdir= $prefix."/report";
if ( ! -d $reportdir) {
	`mkdir $prefix/report`;
}

my $jellyfish_version = `$scriptPath/bin/jellyfish  --version`;
chomp($jellyfish_version);

# html report:
my $htmlout = $prefix."_".$set1."-".$set2.".report.html";
my $tmpRscript = $prefix."/".$prefix."_X_Rmd.render.R";
open(TMPR, ">$tmpRscript") || die;
print TMPR "library\(rmarkdown\)\n";
#print TMPR "library\(knitr\)\n";
print TMPR "\n";
print TMPR "render\(\'$rmd\',\n";
print TMPR "  params = list(prefix=\"$prefix\",\n";
print TMPR "    jfversion=\"$jellyfish_version\",\n";
print TMPR "    version=\"$version\",\n";
print TMPR "    kad=\"$kad_file\",\n";
print TMPR "    kdbin=\"$bincount_file\",\n";
print TMPR "    set1=\"$set1\",\n";
print TMPR "    set2=\"$set2\",\n";
print TMPR "    binlen=\"$binlen\"\),\n";
print TMPR "  knit_root_dir=getwd\(\),\n";
print TMPR "  output_dir=paste0\(getwd\(\), \"\/$reportdir\"\),\n";
print TMPR "  output_format=\"html_document\",\n";
print TMPR "  output_file=\"$htmlout\"\)\n";
close TMPR;

# run
`Rscript $tmpRscript`;
#`rm $tmpRscript`;

close LOG;

###############################################
# module 1: num2bin
###############################################
sub num2bin {
# round a number using to binsize;
# # 1.54 and binsize 0.5; return 1.5
# # 1.04 and binsize 0.5; return 1.0
	my ($innum, $inbinsize) = @_; 
	my $binned_num = int($innum / $inbinsize + 0.5) * $inbinsize;
	return $binned_num;
}


