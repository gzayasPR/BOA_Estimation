#!/usr/local/bin/perl -w
use strict;

if(@ARGV != 6){
    print "USAGE::perl run_LAMPPED.pl <posfile> <EURHaps> <NAHaps> <YRIHaps> <genfile> <outfile>\n";
    die;
}


my $numStatesHMM = 15;
my $winSize = 300;
my $posfile = shift @ARGV;
my $EURfile = shift @ARGV;
my $NAfile = shift @ARGV;
my $YRIfile = shift @ARGV;
my $genfile = shift @ARGV;
my $outfile = shift @ARGV;



my $cmd ="./bin/gettriophase $posfile $genfile $outfile";
print "Running haplotype inference based on Mendelian Inheritance rules::$cmd\n";
`$cmd`;

$cmd ="./bin/haplanc $winSize $numStatesHMM $posfile $EURfile $NAfile $YRIfile $outfile.phase $outfile.lanc 1";
print "Running LAMPED::$cmd\n";

`$cmd`;


