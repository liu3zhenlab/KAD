#!/use/bin/perl -w
#
#======================================================================
# KADprofile.pl
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 7/27/2019
#
# The script is to generate KAD profiles for 1 or multiple assemblies
#======================================================================

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use kmerge;
use kad;

use constant VER => "0.01";

sub prompt {
	print <<EOF;
	Usage: perl KADprofile [options]
	[Options]
	--read <file>:	FASTQ/A read file for k-mer generation; overridden by --rkmer;
					the parameter can be used multiple times;
					Jellyfish is used to generate counts of k-mers.
	--minc <num>:	minimal number of counts per k-mer from reads;
					k-mers with counts smaller than <num> are not output. default=5
	--asm <file>:	FASTA sequence file for k-mer generation; overridden by --akmer
					the parameter can be used multiple times to allow using multiple FASTA file;
					each file is considered an indepedent assembly.
	--rkmer <file>:	a k-mer table file of reads; the file contains two columns and the header (e.g., Kmer	SampleA)
					the 1st column contains k-mer sequence; the 2nd is for counts per k-mer.
	--akmer <file>:	a k-mer table file of one or multiple assemblies;
					each file at least has two columns with k-mer sequences at the 1st column;
					the rest of columns are counts for each k-mer in each assembly;
					more than two columns to represent multiple assemblies are allowed;
					the header is required to contain information for each column.	
	--rid <str>:	ID used in the header of the k-mer table generated from reads;
					invalid when --rkmer is specified;
	--aid <str>:	ID used in the header of the asm k-mer table to be generated from each assembly;
					the parameter can be used multiple times to match --asm input.
	                By default, a header ID is generated from the file name of each assembly
					by removing PATH and the suffix of .fa, .fas, or .fasta.
					IMPORTANT: If --aid is specified, it must match the --asm parameter. One --aid is used
					for one k-mer table generated from one assembly file specified by --asm. To match --asm and --aid input,
					the order of two inputs should have the matching order;
					invalid when --akmer is specified.
	--prefix		the prefix for all output files; also it is the output directory; default=kad.
	--klen			length of k-mers; default=25.
	--readdepth		estimated depth of reads; not required; if specified, it will be compared to the mode of read k-mers.
	--kadcutoff		a set of numbers to define k-mer categories:
					1. Good: k-mers basically containing no errors, for which, by default, KADs are between -0.5 and 0.5;
					         Note: some k-mers with low counts but not presented in the assembly are in this category;
					2. SingleError: k-mers showing a single copy in the assembly but with no reads supported, for which
					                KADs equals to -1; this value is fixed; 
					3. MultiError: k-mers showing multiple locations in the assembly but read depths indicate lower copies, for
					               which, by default, KADs are smaller than -2;
					4. LowMiss: k-mers showing less copies in the assembly compared to copies indicated by read depths, for
					            which, by default, KADs are between 0.75 and 2;
					5. HighMiss: k-mers showing less copies at a high degree in the assembly as compared to copies
					             indicated by read depths, for which, by default, KADs are higher than 2;
					default=(-2, -0.5, 0.5, 0.75, 2)
	--threads		number of cpus; default=1.
	--help:			help information

	o example 1: to generate KAD profiles with read and assemblies
		perl KADprofile.pl --read data1.fq --read data2.fq --rid PE250 --minc 5 \
		                   --asm asm1.fasta --asm asm2.fasta --asm asm3.fasta \
						   --aid asm1 --aid asm2 --aid asm3

	o example 2: to generate KAD profiles wtih counts of k-mers
		perl KADprofile.pl --rkmer read.kmer.txt --akmer asms.kmer.txt \

EOF
exit;
}

###############################################
# parameters:
###############################################
my %opts = ();
my (@read, $minc, @asm, $rkmer, $akmer, $rid, @aid, $readdepth, $threads,
    $kadcutoff, $klen);
&GetOptions(\%opts, "read=s@", "minc=i", "asm=s@", "rkmer=s",
                    "akmer=s", "rid=s", "aid=s@", "prefix=s",
					"readdepth=i", "threads=i", "klen=i",
					"kadcutoff=s", "help");

&prompt if exists $opts{help};
@read = @{$opts{read}} if exists $opts{read};
$minc = exists $opts{minc} ? $opts{minc} : 5;
$rkmer = $opts{rkmer} if exists $opts{rkmer};
@asm = @{$opts{asm}} if exists $opts{asm};
$akmer = $opts{akmer} if exists $opts{akmer};
$rid = $opts{rid} if exists $opts{rid};
@aid = @{$opts{aid}} if exists $opts{aid};
my $prefix = "kad"; # default
$prefix = $opts{prefix} if exists $opts{prefix};
$readdepth = $opts{readdepth} if exists $opts{readdepth};
$threads = exists $opts{threads} ? $opts{threads} : 1;
$klen = exists $opts{klen} ? $opts{klen} : 25;
$kadcutoff = exists $opts{kadcutoff} ? $opts{kadcutoff} : "-2 -0.5 0.5 0.75 2";
my @kadcutoff = split(" ", $kadcutoff);

###############################################
# preparation
###############################################

# create a directory for outputs
if (-d $prefix) {
	print LOG "the directory $prefix exists. Run stopped!\n";
	exit;
} else {
	`mkdir $prefix`;
}

# script path:
my $scriptPath = $FindBin::Bin;
my $binPath = $scriptPath."/bin/";
#my $utilPath = $scriptPath."/util/";


# log file
my $logfile = $prefix."_KADrun.log";
open(LOG, ">$logfile") || die;

###############################################
# read to kmers:
###############################################
my $rkout; # read kmer table file name
if (!exists $opts{rkmer}) {
	$rkout = $prefix."/".$prefix."_0_reads.kmer.table.txt"; # reads k-mer table output
	if (exists $opts{read}) {
		print LOG "o generate counts of read k-mers\n";
		foreach (@read) {
			print LOG "$_\n";
		}
		my $reads_files = join(" ", @read);
		if (!exists $opts{rid}) {
			$rid = "reads";
		}
		&kgen($reads_files, $rid, $minc, $rkout); # input: files, id, minimal count, outfile #output: k-mer table file
	} else {
		print LOG "ERROR: read files or the read k-mer file must be provided. Run stopped\n"; 
		exit;
	}
} else {
	$rkout = $rkmer;
}

sub kgen {
# fasta/q to kmer table
	my ($infiles, $in_rid, $in_minc, $in_out) = @_;
	my $jf_out = $in_out.".jf"; 
	
	# generate jf file
	`$binPath/jellyfish count -C -s 1000M -m $klen -t $threads -L  $in_minc -o $jf_out $infiles`;
	
	# add header
	open(OUT, ">$in_out") || die;
	print OUT "Kmer\t$in_rid\n";
	close OUT;
	
	# jf to tab-separated text format
	`$binPath/jellyfish dump -c -t $jf_out >> $in_out`;
	`rm $jf_out`;
}

###############################################
# asm to kmers:
###############################################
my @akout; # assembly k-mer table outputs
if (!exists $opts{akmer}) {
	if (exists $opts{asm}) {
		if (exists $opts{aid}) {
			if ($#aid != $#asm) {
				print LOG "ERROR: --aid and --asm should be invoked equal times. Run stopped\n";
				exit;
			} else {
				foreach (@aid) {  # set asm kmer output files
					my $akout = $prefix."/".$prefix."_1_".$_.".kmer.table.txt";
					push(@akout, $akout);
				}
			}
		} else { # by default, using asm file names as IDs
			foreach (@asm) {
				my $asm_fn = $_;
				$asm_fn =~ s/.*\///g;
				$asm_fn =~ s/\.fa$|\.fas$|\.fasta$//g;
				push(@aid, $asm_fn);
				my $akout = $prefix."/".$prefix."_1_".$asm_fn.".kmer.table.txt";
				push(@akout, $akout);
			}
		}

		# asm2kmer
		print LOG "o generate counts of k-mers of assemblies\n";
		
		foreach (@asm) {
			print LOG "$_\n";
		}
		
		# generate k-mer talbe for each assembly:
		foreach (my $i=0; $i<=$#asm; $i++) {
			&kgen($asm[$i], $aid[$i], 1, $akout[$i]); # input: files, id, minimal count, outfile #output: k-mer table file
		}
	}
} else {
	push(@akout, $akmer);
}

###############################################
# sequentially merge k-mer tables (taking time)
###############################################
my $mergekmer = $prefix."/".$prefix."_2_merge.kmer.table.";
my $mergekmer_out = $mergekmer."txt";
`cp $rkout $mergekmer_out`;
print LOG "o merge k-mer tables\n";
for (my $i=0; $i<=$#akout; $i++) {
	my $tmpout = $mergekmer."$i";
	kmerge::kmerge($mergekmer_out, $akout[$i], $tmpout);
	`mv $tmpout $mergekmer_out`;
}


# read count mode:
print LOG "o determine mode of counts of read k-mers\n";
my $cmodev = kad::cmode($mergekmer_out, 2);
print LOG "the mode of counts of read k-mers is $cmodev\n";
if (defined $readdepth and abs($readdepth - $cmodev) / $cmodev > 0.5) {
	print LOG "WARNING: the estimated mode is quite different from the estimated depth that is $readdepth.\n";
}

###############################################
# KAD calculation
###############################################
my $kadout = $prefix."/".$prefix."_3_kad.txt";
open(KAD_OUT, ">$kadout") || die;

open(MERGE_IN, $mergekmer_out) || die;
$_ = <MERGE_IN>;
chomp;
my @header = split;

print LOG "o calculate KAD\n";

# print header
my @kad_cols; # cols of KAD values in the new output
my @kad_colnames; # col names
print KAD_OUT "$_";
for (my $i=2; $i<=$#header; $i++) {
	my $new_colname = $header[$i].".KAD";
	print KAD_OUT "\t$new_colname";
	push(@kad_cols, $#header + $i);
	push(@kad_colnames, $new_colname);
}
print KAD_OUT "\n";

while(<MERGE_IN>) {
	chomp;
	my @line = split;
	my $read_cok = $line[1];
	print KAD_OUT $_;
	for (my $i=2; $i<=$#line; $i++) {
		my $asm_cok = $line[$i];
		my $kad_value = "NA";
		if ($read_cok > 0 or $asm_cok > 0) {
			$kad_value = kad::kad_cal($read_cok, $asm_cok, $cmodev);
		}
		print KAD_OUT "\t$kad_value";
	}
	print KAD_OUT "\n";
}

close MERGE_IN;
close KAD_OUT;

###############################################
# KAD statistics
###############################################
my %plusHighNum;
my %minusHighAbsNum;
my %negative1count;
my %around0; # -.5 to .5
my @around0_range = @kadcutoff[1..2];
my %around1;
my @around1_range = @kadcutoff[3..4];
my %total;

print LOG "o Summarize KAD values\n";

#default: @kadcutoff = (-2 -0.5 0.5 0.75 2)

open(KAD_IN, $kadout) || die;
$_ = <KAD_IN>; # skip header
while(<KAD_IN>) {
	chomp;
	my @line = split;
	for (my $i=0; $i<=$#kad_cols; $i++) {
		my $colname = $kad_colnames[$i];
		my $count = $line[$kad_cols[$i] - 1];
		if ($count ne "NA") {
			# total
			$total{$colname}++;
			# values close to certain numbers
			if ($count >= $around0_range[0] and $count <= $around0_range[1]) {
				$around0{$colname}++;
			} elsif ($count >= $around1_range[0] and $count <= $around1_range[1]) {
				$around1{$colname}++;
			} elsif ($count == -1) {
				$negative1count{$colname}++;
			}
			
			# plus and minus
			if ($count >= $kadcutoff[4]) {
				$plusHighNum{$colname}++;;
			} elsif ($count <= $kadcutoff[0]) {
				$minusHighAbsNum{$colname}++;
			}
		}
	}
}
close KAD_IN;

# summary output
my $kadstat = $prefix."/".$prefix."_4_kad.stat.txt";
open(KAD_STAT, ">$kadstat") || die;
print KAD_STAT "Data\tTotal\tGood\tSingleError\tMultiError\tLowMiss\tHighMiss\n";
foreach my $ecolname (@kad_colnames) {
	exists $total{$ecolname} ? print KAD_STAT "$ecolname\t$total{$ecolname}\t" : print KAD_STAT "$ecolname\t0\t";
	exists $around0{$ecolname} ? print KAD_STAT "$around0{$ecolname}\t" : print KAD_STAT "0\t";
	exists $negative1count{$ecolname} ? print KAD_STAT "$negative1count{$ecolname}\t" : print KAD_STAT "0\t";
	exists $minusHighAbsNum{$ecolname} ? print KAD_STAT "$minusHighAbsNum{$ecolname}\t" : print KAD_STAT "0\t";
	exists $around1{$ecolname} ? print KAD_STAT "$around1{$ecolname}\t" : print KAD_STAT "0\t";
	exists $plusHighNum{$ecolname} ? print KAD_STAT "$plusHighNum{$ecolname}\n" : print KAD_STAT "0\n";
}

close KAD_STAT;
close LOG;

