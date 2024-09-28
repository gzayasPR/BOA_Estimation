#!/usr/local/bin/perl -w
use strict;

if(@ARGV != 6){
    print "USAGE::perl run_LAMPLD.pl <posfile> <EURHaps> <YRIHaps> <CHBHaps> <genfile> <outfile>\n";
    die;
}

my $numStatesHMM = 15;
my $winSize = 15;
my $posfile = shift @ARGV;
my $EURfile = shift @ARGV;
my $CHBfile = shift @ARGV;
my $YRIfile = shift @ARGV;
my $genfile = shift @ARGV;
my $outfile = shift @ARGV;

my $cmd ="./bin/unolanc $winSize $numStatesHMM $posfile $EURfile $YRIfile $CHBfile $genfile $outfile";
print "Running LAMPLD::$cmd\n";

`$cmd`;
