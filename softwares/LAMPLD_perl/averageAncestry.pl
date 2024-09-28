#!/usr/bin/perl

###########################
#Daniel Hui################
#University of Pittsburgh##
#9/8/15####################
###########################

use strict;
use warnings;

if ( @ARGV != 4){
        print "USAGE:: perl averageAncestry.pl <phased or unphased> <ways_admixed> <standardized input> <avg ancestry>\n";
        print "There are " . scalar @ARGV . " arguments instead of 4.\n";
        die;
}

open my $in, "$ARGV[2]" or die $!;
open my $out, ">$ARGV[3]" or die $!;

if (lc($ARGV[0]) eq "phased"){
        if($ARGV[1] eq "2"){
                &phased2way($in,$out);
        } elsif ($ARGV[1] eq "3"){
                &phased3way($in,$out);
        } else{
                die $!;
        }
} elsif (lc($ARGV[0]) eq "unphased"){
         if($ARGV[1] eq "2"){
                &unphased2way($in,$out);
        } elsif ($ARGV[1] eq "3"){
                &unphased3way($in,$out);
        } else{
                die $!;
        }
} else{
        print "input format is:: perl averageAncestry.pl <phasing> <ways_admixed> <infile> <outfile>\n";
	print "The ways admixed are <2> or <3>.\n";
	print "The phasing is <phased> or <unphased>, phased has haplotype information unphased does not.\n";
        die;
}


######SUBS#####
sub phased2way{
	my ($in, $out) = @_;

	print $out "SAMPLE POP1 POP2\n";

	my $count = 0;
	my $line1 = 0;
	my @lLine;

	while(<$in>){
       		#take first line
		if($line1==0){
                	@lLine = split /\s+/;
                	$line1++;
        	}

        	#on the second line
        	else{
                	my @line = split /\s+/;
                	my $sum1 =0;
		
			#one column of header
                	for(my $i = 1 ; $i < @line; $i+=2){
                        	$sum1 = $sum1 + $line[$i] + $lLine[$i];
                	}

               		$sum1 = $sum1 / (scalar (@line-2));
                	my $sum2 = 1 - $sum1;

                	$count++;
                	$line1 = 0;
                	@lLine = "";

                	print $out "$count: $sum1 $sum2\n";
        	}
	}
}


sub phased3way{
	my ($in, $out) = @_;

	print $out "SAMPLE POP1 POP2 POP3\n";

	my $count;
	my $line1 = 0;
	my @lLine;

	while(<$in>){
        	#take first line
        	if($line1==0){
               		@lLine = split /\s+/;
                	$line1++;
        	}

        	#on the second line
        	else{
                	my @line = split /\s+/;
                	my $sum1 =0;
                	my $sum2 =0;
			my $sum3 =0;
		
			#one column of header
                	for(my $i = 1 ; $i < @line; $i+=3){
                        	$sum1 = $sum1 + $line[$i] + $lLine[$i];
                	}

                	for(my $i = 2; $i < @line; $i+=3){
                        	$sum2 = $sum2 + $line[$i] + $lLine[$i];
                	}
			
			for(my $i = 3; $i < @line; $i+=3){
				$sum3 = $sum3 + $line[$i] + $lLine[$i];
			}

                	$sum1 = $sum1 / (scalar (2*(@line-2)/3));
                	$sum2 = $sum2 / (scalar (2*(@line-2)/3));
			$sum3 = $sum3 / (scalar (2*(@line-2)/3));

                	$count++;
                	$line1 = 0;
                	@lLine = "";

                	print $out "$count: $sum1 $sum2 $sum3\n";
        	}
	}
} 


sub unphased2way{
	my ($in,$out) = @_;

	print $out "SAMPLE POP1 POP2\n";

	while(<$in>){
        	my @line = split /\s+/;
       		my $sum1=0;

		#one line of header
        	for(my $i = 1 ; $i < @line; $i+=2){
                	$sum1 = $sum1 + $line[$i];
        	}

        	$sum1 = $sum1 / (scalar (@line-1)  );
        	my $sum2 = 1 - $sum1;

        	print $out "Sample$.: $sum1 $sum2\n";
	}
}


sub unphased3way{
	my ($in,$out) = @_;

	print $out "SAMPLE CEU YRI CHB\n";

	while(<$in>){
        	my @line = split /\s+/;

        	my $sum1=0;
        	my $sum2=0;

		#one column of header
        	for(my $i = 1 ; $i < @line; $i+=3){
                	$sum1 = $sum1 + $line[$i];
        	}

        	for(my $j = 2 ; $j < @line; $j+=3){
                	$sum2 = $sum2 + $line[$j];
        	}

        	$sum1 = $sum1 / (scalar (@line-1) - (@line/3) );
        	$sum2 = $sum2 / (scalar (@line-1) - (@line/3) );
        	my $sum3 = 1 - ($sum1 + $sum2);

        	print $out "$.: $sum1 $sum2 $sum3\n";
	}
}
