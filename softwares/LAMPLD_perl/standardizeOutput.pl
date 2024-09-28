#!/usr/bin/perl

###########################
#Daniel Hui################
#University of Pittsburgh##
#9/4/15####################
###########################

#use strict;
use warnings;

if ( @ARGV != 4){
	print "USAGE:: perl lait.pl <software_name> <ways_admixed (phasing for ELAI)> <software_output> <standard_output>\n";
       	print "There are " . scalar @ARGV . " arguments instead of 4.\n";
       	die;
}

open my $in, "$ARGV[2]" or die $!;
open my $out, ">$ARGV[3]" or die $!;

if (lc($ARGV[0]) eq "lamp"){
	if($ARGV[1] eq "2"){ 
		&lamp2way($in,$out);
        } elsif ($ARGV[1] eq "3"){
		&lamp3way($in,$out); 
        } else{
                print "Sorry, only 2-way or 3-way for $ARGV[1]-way admixture for LAMP.\n";
                die;
        }
} elsif (lc($ARGV[0]) eq "lampld" || lc($ARGV[0]) eq "lamp-ld"){
	 if($ARGV[1] eq "2"){
                &lampld2way($in,$out);
        } elsif ($ARGV[1] eq "3"){
                &lampld3way($in,$out);
        } else{
                print "Sorry, only 2-way or 3-way for $ARGV[1]-way admixture for LAMP-LD.\n";
                die;
        }
} elsif (lc($ARGV[0]) eq "hapmix"){
	 if($ARGV[1] eq "2"){
                &hapmix($in,$out);
        } else{
                print "Sorry, only 2-way for $ARGV[1]-way admixture for HAPMIX.\n";
                die;
        }
} elsif (lc($ARGV[0]) eq "elai"){
	 if(lc($ARGV[1]) eq "unphased"){
                &elaiUnphased($in,$out);
        } elsif (lc($ARGV[1]) eq "phased"){
                &elaiPhased($in,$out);
        } else{
                print "Phased or unphased  $ARGV[1]-way admixture for ELAI.\n";
                die;
        }
} else{
        print "input format is:: perl lait.pl <software> <ways_admixed (phasing for ELAI)> <infile> <outfile>\n";
        print "The software choices are lamp, lampld, elai, or hapmix.\n";
        print "The parameter for elai is the output format <phased> or <unphased>. For the other software it is the ways-admixed, either <2> or <3>.\n";
        die;
}


#####################SUBS######################

########
##LAMP##
########
sub lamp2way{
my ($in, $out) = @_;
my $secondLine = 0;
my @firstLine;
my $sampleCount = 0;
while (<$in>){
        chomp;
        #put first line of two into array, then exit this iteration of while loop
        if ($secondLine == 0){
                @firstLine = split("\t", $_);

                #multiple each value by 2
                for (my $i = 1; $i < @firstLine; $i++){
                       	$firstLine[$i] = ($firstLine[$i] * 2);
                }

                $secondLine++;
        }else{
                #print out sample no.
                $sampleCount++;
                print $out "Sample_$sampleCount\t";

                #multiply each value by 2
                my @secondLine = split("\t", $_);

                for (my $i = 1; $i < @secondLine; $i++){
                        $secondLine[$i] = ($secondLine[$i] * 2);
                }

                #print values of first and second line
                for (my $i = 1; $i < @secondLine ; $i++){
                        print $out "$firstLine[$i] $secondLine[$i]\t";
                }

                print $out "\n";

                #reset values for the next pair of lines
                @firstLine = ();
                $secondLine = 0;
        }
}

close $in;
close $out;

}


sub lamp3way{
	my ($in, $out) = @_;
	my $secondLine = 0;
	my $thirdLine = 0;
	my $firstLine = 1;

	my @firstLine;
	my @secondLine;

	my $sampleCount = 0;

	while (<$in>){
       		chomp;

        	if ($firstLine == 1){
                	@firstLine = split /\s+/;

                	for(my $i=1; $i<@firstLine; $i++){
                        	$firstLine[$i] = ($firstLine[$i] * 2);
                	}

                	$firstLine = 0;
                	$secondLine = 1;
        	}

        	#put first line of two into array, then exit this iteration of while loop
        	elsif  ($secondLine == 1){
                	@secondLine = split /\s+/;

                	#multiple each value by 2
                	for (my $i =1; $i<@secondLine;$i++){
                                $secondLine[$i] = ($secondLine[$i] * 2);
                	}

                	$secondLine = 0;
                	$thirdLine = 1;

        	}elsif ($thirdLine == 1){

                	#print out sample no.
                	$sampleCount++;
                	print $out "Sample_$sampleCount\t";

                	#multiply each value by 2
                	my @thirdLine = split /\s+/;

                	for (my $i=1; $i<@thirdLine;$i++){
                        	$thirdLine[$i] = ($thirdLine[$i] * 2);
                	}


                	#print values of first and second line
                	for (my $i = 1; $i < @thirdLine ; $i++){
                        	print $out "$firstLine[$i] $secondLine[$i] $thirdLine[$i]\t";     
                	}
                	print $out "\n";

                	#reset values for the next pair of lines
                	@firstLine = ();
                	@secondLine = ();
                	@thirdLine = ();
                	$secondLine = 0;
                	$thirdLine = 0;
                	$firstLine = 1;
        	}
	}
	close $in;
	close $out;
}

###########
##LAMP-LD##
###########

sub lampld2way{
	my ($in, $out) = @_;

	my $sampleNo = 1;
	my $hapNo = 0; 
	while (<$in>){
		chomp;
        	#sample and hap no.
        	$hapNo++;
        	print $out "Sample_$sampleNo-Hap_$hapNo\t";
                
        	if ($hapNo == 2){
                	$hapNo = 0;
                	$sampleNo++;
        	}
                
        	my @line = split(//, $_);
                
        	foreach (@line){
                	if ($_ == 0){
                      		print $out "1 0 \t";
                	} elsif ($_ == 1){
                     		print $out "0 1 \t";
                	} else{
                     		print $out "? ?  ";
                	}
        	}
   	     	print $out "\n";
	}
	close $in;
	close $out;
}


sub lampld3way{
	my ($in, $out) = @_;

	my $sampleNo = 1;
	my $hapNo = 0;
	while (<$in>){
		chomp;

       		#sample and hap no.
       		$hapNo++;
        	print $out "Sample_$sampleNo-Hap_$hapNo\t";

        	if ($hapNo == 2){
        		$hapNo = 0;
                	$sampleNo ++;
        	}

        	#processing
        	my @line = split(//, $_);

        	foreach (@line){
        		if ($_ == 0){
              			print $out "1 0 0\t";
                	} elsif ($_ == 1){
                		print $out "0 1 0\t";
                	} elsif ($_ == 2){
               			print $out "0 0 1\t";
                	} else{
                		print $out  "9 9 9 ";
                	}
        	}
		print $out  "\n";
	}
	close $in;
	close $out;
}

##########
##HAPMIX##
##########
sub hapmix{
	my ($in, $out) = @_;

	#columns to rows
	while (<$in>){
        	chomp;
        	@line = split ("", $_);

        	$lastcol = 0;
        	$lastcol = $#line if $#line > $lastcol;
        	$oldlastcol = $lastcol;

        	for (my $i = $oldlastcol; $i < $lastcol; $i++) {
                	$outline[$i] = "" x $oldlastcol;
        	}

        	for (my $i=0; $i <=$lastcol; $i++) {
                	$outline[$i] .= "$line[$i]"
                	}
	}

	#conversion and print
	my $count = 1;

	for (my $i=0; $i <= $lastcol; $i++) {

		print $out "Sample_".$count."\t";
       		#sample no
        	$count++;
		

        	#conversion
        	my @line =  split(//, $outline[$i]);

        	foreach (@line){
        		if ($_ eq 0){
                        	print $out ("0 2\t");
                	} elsif ($_ eq 1){
                        	print $out ("1 1\t");
                	} elsif ($_ eq 2){
                		print $out ("2 0\t");
                	} else{
                		print $out ("? ?\t");
                	}
        	}
        	print $out "\n";
	}

	close $in;
	close $out;
}

########
##ELAI##
########
sub elaiPhased{
	my ($in,$out) = @_;
	my $sampleNo = 1;
	my $hapNo = 0;
	while (<$in>){
        	chomp;
        	#sample and hap no.
        	$hapNo++;
        	print $out "Sample$sampleNo-Hap_$hapNo\t";

        	if ($hapNo == 2){
                	$hapNo = 0;
                	$sampleNo ++;
        	}
        	print $out "$_\n";
	}       
	close $in;
	close $out;
}

	sub elaiUnphased{
	my ($in,$out) = @_;
	my $sampleNo = 1;
	while (<$in>){
        	chomp;
       		#sample and hap no.
        	print $out "Sample_$sampleNo $_\n";
        	$sampleNo ++;
	} 
	close $in;
	close $out;
}

