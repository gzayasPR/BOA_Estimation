#!/usr/local/bin/perl -w
use strict;

if(@ARGV != 7){
    print "USAGE::perl run_LAMPLD.pl <numStatesHMM> <windowsize> <posfile> <EURHaps> <YRIHaps> <genfile> <outfile>\n";
    die;
}

my $numStatesHMM = shift @ARGV;
my $winSize = shift @ARGV;
my $posfile = shift @ARGV;
my $EURfile = shift @ARGV;
my $YRIfile = shift @ARGV;
my $genfile = shift @ARGV;
my $outfile = shift @ARGV;

my $cmd ="./bin/unolanc2way $winSize $numStatesHMM $posfile $EURfile $YRIfile $genfile $outfile";
print "Running LAMPLD::$cmd\n";

`$cmd`;


