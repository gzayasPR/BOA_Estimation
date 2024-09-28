#!/usr/bin/perl -w

#use strict;

if ((scalar @ARGV) != 2){
    print "USAGE::./thisprog.pl <LAMPLD-outfile> <postprocessed-file>\n";
    die;
}
my $infile = shift @ARGV;
my $outfile = shift @ARGV;
my @hap1;
my @hap2;

open IN, $infile or die "print cannot open $infile\n";
open OUT, ">$outfile" or die "print cannot open $outfile\n"; 

$count = 0; 
while (my $line = <IN>) {
   $line =~ s/^\s+//; 
   chomp $line; 
  # print "Ind: $count\n";
    my @temp =  split /\s+/, $line;
    my $start = 0;
    undef @hap1;
    undef @hap2;
    foreach my $bkpt (@temp){
#	print $bkpt."\n";
        my ($pop, $end ) = split /:/, $bkpt;
	for (my $i =$start; $i < $end; $i++){
	    my @ancs = split //, $pop ;
	    push(@hap1, $ancs[0]);
	    push(@hap2, $ancs[1]);
	}
      $start = $end;
    }
    print OUT join("",@hap1)."\n";
    print OUT join("",@hap2)."\n";
$count++;
}
close(IN);
close(OUT);
