#!/use/bin/perl -w
#
#======================================================================
# kad.pm
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 8/10/2019
#======================================================================

package kad;
use strict;
use warnings;

sub kad_cal {
# calculate KAD
	my ($cok_data, $cok_asm, $cmode_val) = @_;
	# cok_data: counts of k-mers in the data
	# cok_asm: counts of k-mers in the assembly
	# cmode_val: mode of cok_data
	my $kad_value = &log2(($cok_data + $cmode_val) / ($cok_asm + 1) / $cmode_val);
	my $kad_value_d3 = &d3($kad_value); # 3 digits
	return $kad_value_d3;
}

sub log2 {
# log2 function
	my $n = shift;
	return log($n)/log(2);
}

sub d3 {
# keep 3 digits
	my $infloat = shift;
	my $outfloat = sprintf("%.3f", $infloat);
	$outfloat =~ s/.000//g;
	return $outfloat;
}

sub num2bin {
# round a number to binsize;
# 1.54 and binsize 0.5; return 1.5
# 1.04 and binsize 0.5; return 1.0
	my ($innum, $inbinsize) = @_;
	my $binned_num = int($innum / $inbinsize + 0.5) * $inbinsize;
	return $binned_num;
}

sub cmode {
# mode value from a set of numbers
	my ($infile, $col) = @_;
	my $outmode;
	my %intmpcounts;
	open(INTMP, $infile) || die;
	$_ = <INTMP>;
	while (<INTMP>) {
		chomp;
		my @intmpline = split;
		my $intmpcount = $intmpline[$col - 1];
		if ($intmpcount > 0) {
			$intmpcounts{$intmpcount}++;
		}
	}
	close INTMP;
	
	foreach (sort {$intmpcounts{$b} <=> $intmpcounts{$a}} keys %intmpcounts) {
		$outmode = $_;
		last;
	}
	return $outmode;
}

sub hash_min_val {
# the min value of a hash
	my @inhash = @_;
	my $outval;
	foreach my $inhash (@inhash) {
		my %inhash = %{$inhash};
		foreach (sort {$inhash{$a} <=> $inhash{$b}} keys %inhash) {
			if (!defined $outval) {
				$outval = $inhash{$_};
			} elsif ($inhash{$_} < $outval) {
				$outval = $inhash{$_};
			}
			last;
		}
	}
	return $outval;
}

sub hash_max_val {
# the max value of a hash
	my @inhash = @_;
	my $outval;
	foreach my $inhash (@inhash) {
		my %inhash = %{$inhash};
		foreach (sort {$inhash{$b} <=> $inhash{$a}} keys %inhash) {
			if (!defined $outval) {
				$outval = $inhash{$_};
			} elsif ($inhash{$_} > $outval) {
				$outval = $inhash{$_};
			}
			last;
		}
	}
	return $outval;
}

sub hash_max_val_key {
# the key of the max value of a hash
	my @inhash = @_;
	my $outkey=undef;
	foreach my $inhash (@inhash) {
		my %inhash = %{$inhash};
		foreach (sort {$inhash{$b} <=> $inhash{$a}} keys %inhash) {
			if (!defined $outkey) {
				$outkey = $_;
			} elsif ($_ > $outkey) {
				$outkey = $_;
			}
			last;
		}
	}
	return $outkey;
}

sub dist_print {
# print counts distribution
	my ($inname, $inhash, $inmaxcount, $inmax_starnum) = @_;
	my %inhash = %{$inhash};

	print STDERR "\n===\n";
	print STDERR "Summary of $inname\n";
	foreach (sort {$a <=> $b} keys %inhash) {
		my $start_num = int($inmax_starnum * ( $inhash{$_} / $inmaxcount));
		print STDERR "$_\t";
	
		for (my $i=0; $i<$start_num; $i++) {
			print STDERR "\*";
		}

		for (my $i=$inmax_starnum; $i>$start_num; $i--) {
			print STDERR " ";
		}
		print STDERR " $inhash{$_}\n";
	}
}

1;

