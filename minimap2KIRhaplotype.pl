#!/usr/bin/perl -w
use strict;
use warnings;

my $usage = "minimap2KIRHaplotype.pl myfile.paf\n";

# file input
my $paf_file = shift(@ARGV) or die $usage;
#KIR2DL1 15119   12      15106   +       KP420442.1      141674  53205   68297   14607   15095   60      tp:A:P  cm:i:1432       s1:i:14607      s2:i:7610      dv:f:0.0023
#KIR2DL2 14540   15      14489   +       KP420442.1      141674  22508   36988   11152   14509   60      tp:A:P  cm:i:1051       s1:i:11137      s2:i:8334      dv:f:0.0187
#KIR2DL3 14544   5       14528   +       KP420442.1      141674  22490   37021   13920   14535   60      tp:A:P  cm:i:1365       s1:i:13918      s2:i:7573      dv:f:0.0034


open (PAF, $paf_file);
my %PAF_hash;
# key = geneID
# value = 
# 	PID
#	alnProp
#	refStart
#	refEnd
#	print
while (my $line = <PAF>) {
	chomp $line;
	my @line_array = split("\t", $line);
	my $PID = $line_array[9] / $line_array[10];
	my $alnProp = $line_array[10] / $line_array[1];
	if (($PID > 0.80 ) && ($alnProp > 0.40)) {
		my $geneID = $line_array[0];
		$geneID =~ s/KIR//g;
		my $refStart = $line_array[7];
		my $refEnd = $line_array[8];
		if (exists $PAF_hash{$geneID} ) {
		
			# add duplicate copy if PID and alnProp high for both copies		
			if 	(($PID > 0.85) && ($alnProp > 0.9) &&
				($PAF_hash{$geneID}{'copy1'}{'PID'} > 0.85) &&
				($PAF_hash{$geneID}{'copy1'}{'alnProp'} > 0.90)) {
				
				$PAF_hash{$geneID}{'copy2'}{'alnProp'} = $alnProp;
				$PAF_hash{$geneID}{'copy2'}{'PID'} = $PID;
				$PAF_hash{$geneID}{'copy2'}{'refStart'} = $refStart;
				$PAF_hash{$geneID}{'copy2'}{'refEnd'} = $refEnd;
				$PAF_hash{$geneID}{'copy2'}{'print'} = 1;
			}
			
			# update better aln
			elsif ($alnProp > $PAF_hash{$geneID}{'copy1'}{'alnProp'}) {
				$PAF_hash{$geneID}{'copy1'}{'alnProp'} = $alnProp;
				$PAF_hash{$geneID}{'copy1'}{'PID'} = $PID;
				$PAF_hash{$geneID}{'copy1'}{'refStart'} = $refStart;
				$PAF_hash{$geneID}{'copy1'}{'refEnd'} = $refEnd;
				$PAF_hash{$geneID}{'copy1'}{'print'} = 1;
			}
		}
		else {
		# add first copy
			$PAF_hash{$geneID}{'copy1'}{'alnProp'} = $alnProp;
			$PAF_hash{$geneID}{'copy1'}{'PID'} = $PID;
			$PAF_hash{$geneID}{'copy1'}{'refStart'} = $refStart;
			$PAF_hash{$geneID}{'copy1'}{'refEnd'} = $refEnd;
			$PAF_hash{$geneID}{'copy1'}{'print'} = 1;
		}
	}
}


foreach my $key1 (keys %PAF_hash) {
	foreach my $copy1 (keys %{$PAF_hash{$key1}}) {
		foreach my $key2 (keys %PAF_hash) {
			foreach my $copy2 (keys %{$PAF_hash{$key2}}) {
				unless ($key1 eq $key2) {
					my @coords = ();
					push (@coords, $PAF_hash{$key1}{$copy1}{'refStart'});
					push (@coords, $PAF_hash{$key1}{$copy1}{'refEnd'});
					push (@coords, $PAF_hash{$key2}{$copy2}{'refStart'});
					push (@coords, $PAF_hash{$key2}{$copy2}{'refEnd'});
					my $overlap = overlap(@coords);
					if ($overlap > 0) {
						if ($PAF_hash{$key1}{$copy1}{'PID'} > $PAF_hash{$key2}{$copy2}{'PID'}) {
							$PAF_hash{$key2}{$copy2}{'print'} = 0;
						}
						elsif ($PAF_hash{$key2}{$copy2}{'PID'} > $PAF_hash{$key1}{$copy1}{'PID'}) {
							$PAF_hash{$key1}{$copy1}{'print'} = 0;
						}
					}
				}
			}
		}
	}
}


foreach my $key1 (keys %PAF_hash) {
	foreach my $copy1 (keys %{$PAF_hash{$key1}}) {
		if ($PAF_hash{$key1}{$copy1}{'print'} == 1){
			print $key1, "\n";
		}
	}
}

sub overlap {
	my @input = @_;
	my $a_start = shift @input;
	my $a_end = shift @input;
	my $b_start = shift @input;
	my $b_end = shift @input;
	# no overlap
	if ((($a_start < $b_start) && ($a_end < $b_start)) ||
		(($b_start < $a_start) && ($b_end < $a_start))) {
		return 0;
	}
	# nested
	elsif (($a_start < $b_end) && ($a_start > $b_start) &&
			($a_end > $b_start) && ($a_end < $b_end)) {
		return 1;
	} 
	# nested
	elsif (($b_start < $a_end) && ($b_start > $a_start) &&
			($b_end > $a_start) && ($b_end < $a_end)) {
		return 1;
	}
	# overlap
	elsif (($b_start > $a_start) && ($b_start < $a_end) &&
			($a_end > $b_start) && ($a_end < $b_end)) {
		return 2;
	}
	# overlap
	elsif (($a_start > $b_start) && ($a_start < $b_end) &&
			($b_end > $a_start) && ($b_end < $a_end)) {
		return 2;
	}
	else {
		return 3;
	}
}
	
