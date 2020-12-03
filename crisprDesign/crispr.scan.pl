#!/usr/bin/perl -w
# ==========================================================================
# genecrispr.pl
# Sanzhen Liu
# update 2/27/2020
# ==========================================================================
use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /;
use POSIX qw(strftime);
use FindBin;
use Getopt::Long;
use Term::ANSIColor;

my $script_version = "v0.21";
my $scriptPath = $FindBin::Bin; # script path

# initiate parameters:
my $nthreads = 1;
my $ext = 0;
my $feature = "exon";
my $pattern = "G[ATGC]{20}GG";
my $strand = "both";
my $motif = "GCGATG,GGTCTC";
my $max_perfect = 1;
my $offtarget_cutoff_score = 10;

sub errINF {
	print <<EOF;
	Usage: perl genecrispr.pl --fasta <fasta> --gff <gff> [options]
	genecrispr.pl: to design crispr oligos based on the specified gene name
	
	[Options]
	--fasta <fasta>      :FASTA file on which GFF3 was built; required
	--gtf <GTF>          :GTF file; required
	--gene <file>        :gene list; one gene per line; required
	--prefix <str>       :prefix for outputs (by default, --gene input used as the prefix)
	--pattern <str>	     :regular expression for search pattern ("$pattern")
	--strand <str>       :strand to design; options are plus, minus, and both ($strand)
	--feature <str>	     :gene, exon, CDS, or others ($feature)
	--ext <int>          :bp length for the interval extension ($ext)
	--max_perfect <int>  :maximum number of perfectly matched targets ($max_perfect)
	--offcutoff <int>    :cutoff score of off-target crRNA ($offtarget_cutoff_score);
	--crispr_out <str>   :CRISPR oligo output filename; default = <gene>_design.txt
	--fasta_out <str>    :CRISPR oligo output in a FASTA format default = <gene>_design.fas
	--bowtie_db <indexDB>:bowtie database; default is to build the index files from scratch
	--nthreads <num>     :number of threads to run bowtie ($nthreads)
	--motif <str>        :DNA sequence separated by a comma ($motif)
	                      Motif(s) and their reversed complementary sequences are NOT allowed on gRNA sequences.
                          Default sequences match restriction sites of BtgZI and BsaI, respectively.
                          If no motifs are used, --motif can be set to be "NNN".
	--force:             :force to overwrite outputs if specified
	--version:           :version information
	--help:              :help information
EOF
	exit;
}

my (@line, $chromosome, $chr, %chr_seq);
my ($gtf, $fasta, $glist, $seq);
my ($crispr_out, $fasta_out);
my ($prefix, $bowtie_db);
my ($force, $version, $help);

GetOptions("gtf=s"        => \$gtf,
           "fasta=s"      => \$fasta,
           "glist=s"      => \$glist,
           "pattern=s"    => \$pattern,
           "prefix=s"     => \$prefix,
           "strand=s"     => \$strand,
           "feature=s"    => \$feature,
           "ext=i"        => \$ext, 
           "max_perfect=i"=> \$max_perfect,
		   "offcutoff=i"  => \$offtarget_cutoff_score,
		   "crispr_out=s" => \$crispr_out,
           "fasta_out=s"  => \$fasta_out,
           "bowtie_db=s"  => \$bowtie_db,
		   "nthreads=i"   => \$nthreads,
		   "motif=s"      => \$motif,
		   "force"        => \$force,
		   "version"      => \$version,
		   "help"         => \$help);

&errINF if $help;

if ($version) {
	print "$script_version\n";
	exit;
}

if (!defined $gtf and !defined $fasta and !defined $glist) {
	print STDERR color('red');
	print STDERR "--gtf, --fasta, and --glist are required\n"; 
	print STDERR color('reset');
	&errINF;
}

$prefix = "crispr" if (!defined $prefix);
if (-d $prefix) {
	if (!defined $force) {
		print "- $prefix exist!\n";
		print "Add --force to overwrite outputs.\n";
		exit;
	}
} else {
	`mkdir $prefix`;
}

$crispr_out = $prefix."/2.gRNA.txt" if (!defined $crispr_out);
$fasta_out = $prefix."/3.gRNA.fas" if (!defined $fasta_out);
my $design_out = $prefix."/4.gRNA.design";

#####################
### 0. check dependency
#####################
my $bowtiepath=`which bowtie`;
chomp $bowtiepath;
if ($bowtiepath eq "") {
	print STDERR "bowtie commands need to be in an executable path\n";
	exit;
}

my $rscript=`which Rscript`;
chomp $rscript;
if ($rscript eq "") {
	print STDERR "R needs to be installed or loaded\n";
	exit;
}

&logmessage("Dependency was detected");

#####################
### 1. read reference
#####################

&logmessage("Read reference FASTA sequences");

open(IN, "<", $fasta) || die;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		if (defined $seq) { 
			$chr_seq{$chr} = $seq;
		}       
		$chr = $1;
		$seq = "";
	} else {
		$seq .= $_; 
	}
}
$chr_seq{$chr} = $seq;
close IN;


#####################
# genelist
#####################
&logmessage("Read gene list");
my $gcount=0;
my %goi;
open(IN, "<", $glist) || die;
while (<IN>) {
	chomp;
	if ($_ !~ /^$/) {
		$goi{$_}++;  # genes of interest
		$gcount++;
	}
}
close IN;

# exit if no genes found
if ($gcount < 1) {
	print STDERR "no genes in the $glist\n";
	exit;
}


#####################
# 2. extract gtf data
#####################

&logmessage("Extract GTF data for genes");

my $genegtf = $prefix."/1.gtf";
`grep -f $glist $gtf > $genegtf`; # extract gtf information for the gene

#####################
# guideRNA design
#####################

&logmessage("Design crRNA on $feature");

my (%oligo_info);
# 1   gramene exon    46229   46342   .   +   .   gene_id "Zm00001d027230"; transcript_id "Zm00001d027230_T001";
open (IN, $genegtf) || die;
while (<IN>) {
	if (!/^$/ and !/^#/) { # not empty line or comment line
		chomp;
		my $gene = "";
		my $subtarget = "";
		@line = split(/\t/,$_);
		my $info = $line[8];
		if ($line[2] eq $feature) {
			if ($info =~ /^gene_id \"(.+?)\"; transcript_id \"(.+?)\"/) {
				$gene = $1;
				$subtarget = $2;
			} elsif ($info =~ /^gene_id \"(.+?)\";/) {
				$gene = $1;
				$subtarget = $gene;
			}
			my $chr = $line[0];
			my $start;
			
			# go through every gene
			if (exists $goi{$gene}) {
				# assign start point:
				if (($line[3] - $ext) <= 0) {
					$start = 1;
				} else {
					$start = $line[3] - $ext;
				}
				# assign end point:
				my $end = $line[4] + $ext;
				my $completed_forward = 0;
				my $completed_reverse = 0;
				my $orientation = "plus";
				if (exists $chr_seq{$chr}) {
					my $current_seq = substr($chr_seq{$chr}, ($start-1), ($end-$start+1));
					while (!$completed_forward | !$completed_reverse) {
						my $design;
						if ($strand eq "plus") {
							$completed_reverse = 1;
						} elsif ($strand eq "both") {
							if ($completed_forward) {
								$orientation = "minus";
								$completed_reverse = 1;
							}
						} elsif ($strand eq "minus") {
							$orientation = "minus";
							$completed_reverse = 1;
						} else {
							die "Only plus, minus, or both can be specified for --strand\n";
						}
						$completed_forward = 1;
						$design = &pattern_search($current_seq, $start, $orientation, $pattern);
						my @design = @{$design};
						foreach (@design) {
							my ($design_oligo, $design_ori_pos, $design_ref_pos) = split(/\t/, $_);
							my $search_info = $gene."-".$subtarget."-".$orientation."-".$design_ori_pos."-".$chr."_".$design_ref_pos;
							
							# chech motif if it is required to check:
							my $is_has_motif = 0;
							if (defined $motif) {
								$is_has_motif = &motif_check($design_oligo, $motif);
							}
							if (!$is_has_motif) { # motif was not identified
								push(@{$oligo_info{$design_oligo}}, $search_info);
							}
						}
					}
				} else {
					print STDERR "$chr WAS NOT FOUND!\n";
				}
			}
		}
	}
}
close IN;

open (TXTOUT, ">$crispr_out") || die;
open (FASOUT, ">$fasta_out") || die;
my @alloligos = keys %oligo_info;
my $count = 0;
my %gRNAs;
my %genecheck;
if ($#alloligos > 0) {
	foreach (sort {$oligo_info{$a} cmp $oligo_info{$b}} keys  %oligo_info) {
		$count++;
		my @oligo_info = @{$oligo_info{$_}};
		my $gRNAseq = $_;
		my $gRNAinfo = join(";", @oligo_info);
		my $genename = $gRNAinfo;
		$genename =~ s/\-.*//g;
		if (!exists $genecheck{$genename}) {
			$genecheck{$genename}++;
			$count=0;
		}
		my $gRNAname = $genename.":".$count;
		print TXTOUT "$gRNAname\t$gRNAseq\t$gRNAinfo\n";
		$gRNAs{$gRNAname} = $gRNAseq."\t".$gRNAinfo;
		print FASOUT ">$genename:$count\n$_\n";
	}
}
close TXTOUT;
close FASOUT;


#####################
# align gRNAs to ref
#####################

&logmessage("Align crRNAs on the reference");

my $aln = File::Temp->new(TEMPLATE => 'tempXXXXX', SUFFIX => '.aln.tmp');
# if bowtie DB is not provided, build it.
if (! defined $bowtie_db) {
	&logmessage("no Bowtie DB found; Build it now");
	my $bowtie_db_dir = $prefix."/".".4.refdb";
	`mkdir $bowtie_db_dir`;
	$bowtie_db="./".$bowtie_db_dir."/refdb";
	`ln $fasta $bowtie_db`;
	`bowtie-build $bowtie_db $bowtie_db`; 
}

# run bowtie alignment
`bowtie --quiet -v 3 -m 100 -n 2 -l 18 -a -f --sam-nohead -B 1 -n 2 -p $nthreads $bowtie_db $fasta_out > $aln`;

#####################
# evaluate gRNAs
#####################

&logmessage("Evaluate crRNAs");

# score rule
# 1: -5
# 2-5: -2
# 6-8: -3
# 9-20:-5
# 22-23: -10
my %entry_record;
open(IN, $aln) || die;
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	my $entry = $line[0];
	my $ori = $line[1];
	my $chr = $line[2];
	my $pos = $line[3];
	my $seq = $line[4];
	if ($ori eq "-") {
		$seq = &revcom($seq);
	}
	my $score = length($seq) - 1;
	if ($score != 22) {
		print "WARNINGS: oligo $seq is not 23 nt long\n";
	}

	if (exists $line[7]) {
		my $var = $line[7];
		my @var = split(",", $var);
		for (@var) {
			my ($var_pos, $var_base) = split(/:/, $_);
			my $var_penalty = &score_cal($var_pos, $seq, $ori);
			$score -= $var_penalty;
		}
	}

	my $aln_record = $chr."(".$ori.")".$pos;
	push(@{$entry_record{$entry}{$score}}, $aln_record);
}


######################
# outpu crRNA
#####################
#GRMZM2G101523:44	22	1	chr5(-)178831410	20	1	chr10(-)39135768	-	0	-	GRMZM2G101523_T02-CDS.394837-minus-945_923_178831432_178831410
&logmessage("Output crRNA design");

open(OUT, ">$design_out") || die;
print OUT "crOligo_name\tBestScore\tBestNum\tBestHits\tSecondBestScore\tSecondBestNum\tSecondBestHits\tThirdBestScore\tThirdBestNum\tThirdBestHits\tcrOligo\tSource\n";
foreach my $each_entry (keys %entry_record) {
	print OUT "$each_entry";
	my %allscores = %{$entry_record{$each_entry}};
	my @score_list = sort {$b <=> $a} keys %allscores;
	for (my $i = 0; $i <= 2; $i++) {
		my ($hit_num, $hits);
		my $cur_score = 0;
		if (exists $score_list[$i]) {
			my @hits = @{$allscores{$score_list[$i]}};
			$hit_num = $#hits + 1;
			$hits = join(";", @hits);
			$cur_score = $score_list[$i];
		} else {
			$hit_num = 0;
			$hits = 0;
		}

		print OUT "\t$cur_score\t$hit_num\t$hits";
	}
	print OUT "\t$gRNAs{$each_entry}\n";
}
close IN;

###############################################
# Rmd for report
###############################################

&logmessage("Generate HTML report");

# Rmd
my $rmd = $scriptPath."\/utils\/crRNAonGene.Rmd";
my $wdir = `pwd`;
chomp($wdir);


#html report:
my $htmlout = $prefix.".gRNAdesign.report.html";
my $tmpRscript = $prefix."/.".$prefix.".Rmd.render.R";
open(TMPR, ">$tmpRscript") || die;
print TMPR "library\(rmarkdown\)\n";
print TMPR "library\(knitr\)\n";
print TMPR "\n";
print TMPR "render\(\'$rmd\',\n";
print TMPR "  params = list(\n";
print TMPR "    wd=\"$wdir\",\n";
print TMPR "    crdesign=\"$design_out\",\n";
print TMPR "    gtf=\"$genegtf\",\n";
print TMPR "    prefix=\"$prefix\",\n";
print TMPR "    nperfect=\"$max_perfect\",\n";
print TMPR "    offcutoff=\"$offtarget_cutoff_score\"\),\n";
print TMPR "  knit_root_dir=getwd\(\),\n";
print TMPR "  output_dir=getwd\(\),\n";
print TMPR "  output_format=\"html_document\",\n";
print TMPR "  output_file=\"$htmlout\"\)\n";
close TMPR;

# run
`Rscript $tmpRscript`;
#`rm $tmpRscript`;


#############
# modules
#############
### function for formatted output:
sub format_print {
	my ($inseq, $formatlen) = @_; 
	while (my $chunk = substr($inseq, 0, $formatlen, "")) {
		print "$chunk\n";
	}   
}

### script for pattern searching
sub pattern_search {
	my ($search_fasta, $ref_start, $ori, $search_pattern) = @_; 
	my $nextpos = 0;
	my @match_out;
	if ($ori eq "minus") {
		$search_fasta = &revcom($search_fasta);
	}
	while ($search_fasta =~ /($search_pattern)/gi) {
		my $match = $1;
		my $fas_len = length($search_fasta);
		my $matchstart = index($search_fasta, $match, $nextpos);
		my $matchend = $matchstart + length($match);
		my ($original_matchstart, $original_matchend);
		if ($ori eq "minus") {
			$original_matchstart = $fas_len - $matchstart;
			$original_matchend = $original_matchstart - length($match) + 1;
		} else {
			$original_matchstart = $matchstart + 1;
			$original_matchend = $matchend;
		}
		my $ref_matchstart = $original_matchstart + $ref_start - 1;
		my $ref_matchend = $original_matchend + $ref_start - 1;
		my @match_ori_pos = ($original_matchstart, $original_matchend);
		my @match_ref_pos = ($ref_matchstart, $ref_matchend);
		my $match_ori_info = join("_", @match_ori_pos);
		my $match_ref_info = join("_", @match_ref_pos);
		my $match_info = $match."\t".$match_ori_info."\t".$match_ref_info;
		push(@match_out, $match_info);
		$nextpos = $matchstart;
	}
	return \@match_out;
}

### script for reverse and complementary
sub revcom {
	my $inseq = shift @_; 
	my $revcom = reverse($inseq);
	$revcom =~ tr/AGCTagct/TCGAtcga/;
	return $revcom;
}

### check particular sequence on gRNA sequences
sub motif_check {
	my $has_motif = 0;
	my ($inseq, $inmotif) = @_;
	my $rc_inseq = &revcom($inseq);
	my @inmotif = split(/,/, $inmotif);
	foreach my $emotif (@inmotif) {
		if ($inseq =~ /$emotif/ or $rc_inseq =~ /$emotif/) {
			$has_motif = 1;	
		}
	}
	return($has_motif);
}

### determine penalty for mismatch for guide RNA design
sub score_cal {
	my $penalty = 0;
	my ($vp, $in_seq, $in_ori) = @_;
	
	if ($in_ori eq "-") {
		$vp = length($in_seq) - $vp + 1;
	}

	if ($vp == 1 or ($vp >= 9 and $vp <= 20)) {
		$penalty = 5;
	} elsif ($vp >= 2 and $vp <= 5) {
		$penalty = 2;
	} elsif ($vp >= 6 and $vp <= 8) {
		$penalty = 3;
	} elsif ($vp >= 22 and $vp <= 23) {
		$penalty = 10;
	}	
	return $penalty;
}

### date_time
sub logmessage {
	my $inmessage = shift;
	my $datestring = strftime "[%b-%d-%Y %H:%M:%S] ", localtime;
	print STDERR color('red');
	print STDERR "$datestring";
	print STDERR color('reset');
	print STDERR "$inmessage\n";
}

