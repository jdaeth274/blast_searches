#!/usr/bin/perl

use strict;
use warnings;
my $usage = "\n Takes input from velvet and creates a single long dna file, shamelessy copied from Nick's compare genome script";

if (@ARGV <= 0){
  print STDERR "$usage";
  exit(0);
}


my $filein = shift(@ARGV);

chomp $filein;

#open(fazza, $filein) or die "\nSorry could not open that file\n";

#$filename = <fazza>;

#close fazza;

## based on the first line in fasta, lets alter this to use the file name
#$filename =~ s/>\.//g;
#$filename =~ s/\.1\s//g;
#$filename =~ s/>.*\.//g;
#$filename =~ s/>//g;
#$filename =~ s/\s+$//;

my $filename = $filein;
# Get the basename of the file
my @files_B = split('\/([^\/]+)$', $filename);

my $size = @files_B;
my $split_file;
# If just given basename use this
if($size == 2) {
	$split_file = $files_B[1];
}else{
	$split_file = $files_B[0];
}



my @files = ($filein, $split_file);

my $dna_out = make_dna(@files);

sub make_dna {
	my $in = shift;
	my $isolate = shift;
	my $out = $isolate;

	if ($isolate =~ /\.fa$|\.fasta$|\.fa$|\.cons$|\.seq$|\.fna$/) {
		$out =~ s/\.fa$|\.fasta$|\.fa$|\.cons$|\.seq$|\.fna$/\.dna/g;
	} else {
		$out = $out;
	}
	print "\n The file's name is $out \n";
	my $genome;
	open IN, $in or die print STDERR "Unable to open input file $in\n";
	while (my $line = <IN>) {
		unless ($line =~ /^>/) {
			chomp $line;
			$genome.=$line;
		}
	}
	close IN;
	open OUT, "> $out" or die print STDERR "Unable to open output file $out\n";
	print OUT ">$isolate\n";
	my @lines = unpack("(A60)*",$genome);
	foreach my $line (@lines) {
		print OUT "$line\n";
	}
	close OUT;
	return $out;
}
