#!/usr/bin/perl

###########################
#Daniel Hui
#dah124@pitt.edu
#University of Pittsburgh
#April 2017
###########################

use strict;
use warnings;

my $software = lc($ARGV[0]);

print "Please be sure .ped has 6 columns of header and .hap has two columns of header!\n";

###################LAMP##################
if($software eq "lamp"){
	if($ARGV[1] eq "2"){ 
		if ( @ARGV != 5){
        	print "USAGE:: perl lait.pl <lamp> <2> <map> <ped> <output_path>\n";
		print "There are " . scalar @ARGV . " arguments instead of 5.\n";
        	die;
		}	
	} elsif ($ARGV[1] eq "3"){ 
		if ( @ARGV != 8){
		print "USAGE:: perl lait.pl <lamp> <3> <map> <ped> <freqs_pop1> <freqs_pop2> <freqs_pop3> <output_path>\n";
	   	print "There are " . scalar @ARGV . " arguments instead of 8.\n";
		die;
		}
	} else{
        	print "Sorry, LAIT does not currently support $ARGV[1]-way admixture for LAMP.\n";
        	die;
        }

	print "Creating files...\n";

	open my $MAP, "$ARGV[2]", or die "mapfile open error";
        open my $PED, "$ARGV[3]", or die "pedfile open error";
	my $path = "$ARGV[-1]";
	
	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";
	
	open my $GENO, "+>$path/geno.txt";
	open my $CONFIG, ">$path/config.txt";
	open my $SNPS, "+>$path/chr_.pos";

	my %snps;
	my %alts;
	
	if ($ARGV[1] eq "3"){
		open my $FREQ1, "$ARGV[4]" or die "freq1 error";
		open my $FREQ2, "$ARGV[5]" or die "freq2 error";
		open my $FREQ3, "$ARGV[6]" or die "freq3 error";
	
		open my $OUTFREQ1, ">$path/freqs_pop1.txt" or die "outfreq1 error";
		open my $OUTFREQ2, ">$path/freqs_pop2.txt" or die "outfreq2 error";
		open my $OUTFREQ3, ">$path/freqs_pop3.txt" or die "outfreq3 error";

		%snps = &lampfreq($FREQ1, $FREQ2, $FREQ3, $MAP, $OUTFREQ1, $OUTFREQ2, $OUTFREQ3, %snps);

		#get alternate alleles from any frequency file (used pop1)
		#can choose any since has to be subset anyway
		#ref and alt alleles are same for each ref pop
		seek($FREQ1, 0, 0);
		while(<$FREQ1>){
			chomp;
			my @line = split /\s+/;
	  		if (exists $snps{$line[0]}){
	    			$alts{$snps{$line[0]}} = $line[4];
			}
		}
	}

	&lampconfig($CONFIG, $ARGV[1]);
	seek ($MAP, 0, 0) or die "can't go back: $!";

	#lampsnps, didn't work when passed hash into sub
	#"Odd number of elements in hash assignment"
	#edit* you have to use hash references
	#fix it later i suppose, it all works now...
	if ($ARGV[1] eq "3"){

	#SNP file
	#lampsnps()
	my %exists;
	while (<$MAP>){
		chomp;
     		my @line = split /\s+/;
        	if (exists $snps{$line[1]} && $line[3] ne "0"){
			if (!exists $exists{$line[1]}){
				$exists{$line[1]} = $line[1];
       		 		print $SNPS "$line[3]\n";
			}
		}
	}
	print "SNP file done.\n";

	#Genotype file
	#lampgeno()
	while (<$PED>){
        	my @line = split /\s+/;

		for(my $i = 0; $i < 6; $i++){
			shift @line;
		}

        	#alt/variant is 0, reference is 1
        	for (my $i = 0; $i < scalar(@line) - 1; $i+=2){
			my $count = $i/2;
			if (exists $alts{$count}){
				if ( ($line[$i] eq $alts{$count}) && ($line[$i+1] eq $alts{$count}) ){
                        		print $GENO "0\t";
		 		}elsif ( ($line[$i] ne $alts{$count}) && ($line[$i+1] ne $alts{$count}) ) {
                        		print $GENO "2\t";
                  		}elsif ( ($line[$i] eq $alts{$count}) && ($line[$i+1] ne $alts{$count} )  || ($line[$i] ne $alts{$count}) && ($line[$i+1] eq $alts{$count} )  ) {
                        		print $GENO "1\t";
                  		}else{
                        		print $GENO "-1\t";
                        	}
			}
		}
		print $GENO "\n";

	}
	print "Genotype file done.\n";

	#2way
	}else{
		open my $INTER, "+>$path/interfile.txt" or die "interfile: $!";
		open my $INTER0, "+>$path/inter0.txt" or die "inter0: $!";
		my %lines = &lampsnps($MAP, $SNPS);
		&lampgeno($PED, $GENO, $INTER, $INTER0, $path, $SNPS, \%lines);
	}

	print "All done! Check out $path\n";
	
}


####################LAMPLD###################
elsif($software eq "lampld" || $software eq "lamp-ld"){
	my $path = "$ARGV[-1]";
	if($ARGV[1] eq "2"){
	        if(@ARGV != 8){
       			print "USAGE:: perl lait.pl <lampld> <2> <map> <ped> <ref_snps>".
				" <haps_pop1> <haps_pop2> <output_path>\n";
		   	print "There are " . scalar @ARGV . " arguments instead of 8.\n";
       	     		die;
        	}	

	#3way
       	} elsif($ARGV[1] eq "3"){
		if(@ARGV != 9){
      			print "USAGE:: perl lait.pl <lamp-ld> <3> <map> <ped> <snps>". 
				" <haps_pop1> <haps_pop2> <haps_pop3> <output_path>\n";
    			print "There are " . scalar @ARGV . " arguments instead of 9.\n";
          		die;
        	}
  		print "Creating files...\n";
        } else {
        	print "Sorry, $ARGV[-1]-way admixture is not supported";
		die;
	}

        #END 3WAY
 
	if($ARGV[1] eq "2"){
        	print "Creating files...\n";
        }

	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";

	open my $MAP, "$ARGV[2]", or die "mapfile open error";
        open my $PED, "$ARGV[3]", or die "pedfile open error";
        open my $SNPS, "$ARGV[4]", or die "snps open error";
        open my $HAPS1, "$ARGV[5]", or die "pop1 haps open error";
        open my $HAPS2, "$ARGV[6]", or die "pop2 haps open error";
	
	open my $POS, ">$path/chr.pos";
	open my $GENO, ">$path/genofile.gen";
	open my $POP1, ">$path/pop1.hap", or die "pop1 geno open error";
        open my $POP2, ">$path/pop2.hap", or die "pop2 geno open error";
	
	my ($pedCount, $hapCount) = &ldpos($MAP, $SNPS, $POS);

	my $coding = &ldgeno($PED, $GENO, $pedCount);
	&ldhap($HAPS1, $POP1, $coding);
	print "Reference haplotypes pop1 done.\n";
	&ldhap($HAPS2, $POP2, $coding);
	print "Reference haplotypes pop2 done.\n";
	
	if($ARGV[1] eq "3"){
        	open my $HAPS3, "$ARGV[7]", or die "pop3 haps open error";
            	open my $POP3, ">$path/pop3.hap", or die "pop3 geno error: $!";
            	&ldhap($HAPS3, $POP3, $coding);
	    	print "Reference haplotypes pop3 done.\n";
	}
	
	print "All done! Check out $path\n";
	
}

######################ELAI#########################3
elsif($software eq "elai"){
	
	my $path = "$ARGV[-1]";
	
	#ELAI 2WAY
	if($ARGV[1] eq "2"){
		print "You have chosen ELAI for 2way admixed samples.\n";
	
		if(@ARGV != 9){
            		print "USAGE:: perl lait.pl <elai> <2> <map> <ped> ". 
			"<haps_pop1> <snps_pop1> <haps_pop2> <snps_pop2> <output_path>\n";
	    		print "There are " . scalar @ARGV . " arguments instead of 9.\n";
            		die;
            	}	
            
        #END 2WAY
        
        #ELAI 3WAY
        } elsif($ARGV[1] eq "3"){	
		print "You have chosen ELAI for 3way admixed samples.\n";	
		if(@ARGV != 11){
            		print "USAGE:: perl lait.pl <elai> <3> <map> <ped>".
			" <haps_pop1> <snps_pop1> <haps_pop2> <snps_pop2>".
			" <haps_pop3> <snps_pop3> <output_path>\n";
	    		print "There are " . scalar @ARGV . " arguments instead of 11.\n";
            		die;
            	}
        }
        #END 3WAY
        
	else{
        	print "Sorry, LAIT does not currently support $ARGV[1]-way admixture.";		
        	die;
	}
	
	#ELAI ANY WAY

	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";
	print "Creating files...\n";
	
	open my $MAP, "$ARGV[2]", or die "mapfile open error";
        open my $PED, "$ARGV[3]", or die "pedfile open error";
        open my $HAPS1, "$ARGV[4]", or die "pop1 haps open error";
        open my $SNPS1, "$ARGV[5]", or die "pop1 snps open error";
        open my $HAPS2, "$ARGV[6]", or die "pop2 haps open error";
        open my $SNPS2, "$ARGV[7]", or die "pop2 snps open error";
	
	open my $POS, ">$path/SNP.pos", or die "posfile open error";
        open my $GENO, ">$path/admix.geno", or die "genofile open error";
        open my $POP1, ">$path/pop1.geno", or die "pop1 geno open error";
        open my $POP2, ">$path/pop2.geno", or die "pop2 geno open error";
	
	&admixandpos($MAP, $PED, $POS, $GENO);
	&sourcegeno($HAPS1, $SNPS1, $POP1);
	print "Reference genotype pop1 file done.\n";
	&sourcegeno($HAPS2, $SNPS2, $POP2);
	print "Reference genotype pop2 file done.\n";

	if($ARGV[1] eq "3"){
        	open my $HAPS3, "$ARGV[8]", or die "pop3 haps open error";
        	open my $SNPS3, "$ARGV[9]", or die "pop3 snps open error";
        	open my $POP3, ">$path/pop3.geno", or die "pop3 geno open error";
        	&sourcegeno($HAPS3, $SNPS3, $POP3);
		print "Reference genotype pop3 file done.\n";
	}

	print "All done! Check out $path\n";

}


###############HAPMIX#####################
elsif($software eq "hapmix"){

	if($ARGV[1] eq "2"){ 
		if (@ARGV != 10){
      			print "USAGE:: perl lait.pl <hapmix> <2> <chr#>".
			" <map> <ped> <ref_snps> <recombmap> <haps_pop1>".
			" <haps_pop2> <output_path>\n";
			print "There are " . scalar @ARGV . " arguments instead of 10.\n";
        		die;
		}
       	}else{
        	print "HAPMIX can only do 2-way admixture.\n";       
		die;
	}
        
	print "Creating files...\n";
        
	my $chrNo = $ARGV[2];
	open my $MAP, "$ARGV[3]", or die "mapfile open error";
        open my $PED, "$ARGV[4]", or die "pedfile open error";
        open my $SNPS, "$ARGV[5]", or die "snps open error";
        open my $RECMAP, "$ARGV[6]" or die "recmap open error";
        open my $HAPS1, "$ARGV[7]", or die "pop1 haps open error";
        open my $HAPS2, "$ARGV[8]", or die "pop2 haps open error";
	my $path = "$ARGV[-1]";
	
	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";

	my $dir2 = "$path/RUN";
  	unless(mkdir $dir2) {
     		print "$dir2 already exists\n";
    	}	
	
	open my $GENO, ">$path/AAgenofile.$chrNo" or die "genofile open error";
	open my $AASNPS, ">$path/AAsnpfile.$chrNo" or die "AAsnpfile open error";
	open my $POP1, ">$path/pop1genofile.$chrNo", or die "pop1 geno open error";
        open my $POP2, ">$path/pop2genofile.$chrNo", or die "pop2 geno open error";
        open my $SNPS1, ">$path/pop1snpfile.$chrNo", or die;
        open my $SNPS2, ">$path/pop2snpfile.$chrNo" or die;
        open my $INTER, ">$path/AAinterfile.txt" or die "interfile open error";
        open my $INTER1, ">$path/interfile1.txt" or die "interfile open error";
        open my $INTER2, ">$path/interfile2.txt" or die "interfile open error";
        open my $RATES, ">$path/rates.$chrNo" or die "rates open error";
	open my $PAR, ">$path/parameters.$chrNo.par" or die ".par error";

	#get smallest and largest positions
	my $last;
	my $minPos;
	my $lc;
	while(<$RECMAP>){ 
		$lc++;
		$last = $_; 
		if($lc == 2){
			my @fline = split /\s+/, $_;
			$minPos = $fline[0];	
		}
	}
	my @line =  split /\s+/, $last;
	my $maxPos = $line[0];

	#make files
	&par($PAR, $chrNo, $path);
        &hapmixadmixsnps($MAP, $SNPS, $AASNPS, $SNPS1, $SNPS2, $maxPos, $minPos);
	print "SNP files done.\n";	
	
	seek ($RECMAP, 0, 0) or die "can't go back";
	&hapmixrec($RECMAP, $RATES);
	print "Rates file done.\n";
	
	seek ($SNPS, 0, 0) or die "can't go back";
	seek ($MAP, 0, 0) or die "can't go back";
	
	
	my ($coding, $SNPlist) = &hapmixadmixgeno($SNPS, $MAP, $PED, $INTER, $GENO, $path, $maxPos, $minPos);
	print "AAgenofile done.\n";
	seek ($SNPS, 0, 0);
	seek ($MAP, 0, 0);
	&hapmixrefhaps($SNPS, $MAP, $HAPS1, $INTER1, $POP1, $path, "interfile1.txt", $maxPos, $minPos, $coding, $SNPlist);
	print "Reference haps pop1 done.\n";

        seek ($SNPS, 0, 0);
        seek ($MAP, 0, 0);

	&hapmixrefhaps($SNPS, $MAP, $HAPS2, $INTER2, $POP2, $path, "interfile2.txt", $maxPos, $minPos, $coding, $SNPlist);
	print "Reference haps pop2 done.\n";	

	print "All done! Check out $path\n";
	

} else{
	print "See inputs by entering:: perl lait.pl <software> <ways_admixed>\n";
	print "The software choices are lamp, lampld, elai, or hapmix.\n";
	print "The ways_admixed are 2 or 3.\n";
	die;
}




#<--------------------------#######SUBS########---------------------------->


#****************************LAMP****************************

################lamp freq
sub lampfreq{
	my ($FREQ1, $FREQ2, $FREQ3, $MAP, $OUT1, $OUT2, $OUT3, %all) = @_;
	my %freq1;
	my %freq2;
	my %freq3;

	while(<$FREQ1>){
		chomp;
		my @line = split /\s+/;
		$freq1{$line[0]}=$line[0];
	}

	while(<$FREQ2>){
		chomp;
		my @line = split /\s+/;
		$freq2{$line[0]}=$line[0];
	}

	while(<$FREQ3>){
		chomp;
		my @line = split /\s+/;
		$freq3{$line[0]}=$line[0];
	}

	#put each position from POSFILE into a hash as both the key and value
	while (<$MAP>){
		chomp;
		my @line = split /\s+/;
		if( exists $freq1{$line[1]} && exists $freq2{$line[1]} && exists $freq3{$line[1]} && $line[3] ne "0") {
        		$all{$line[1]} = $.;
	        }
	}

	#get the position value for each line in freqfile, and 
	#if it exists in positions hash then print
	#note that this is MINOR allele frequency

	seek($FREQ1,0,0);
	seek($FREQ2,0,0);
	seek($FREQ3,0,0);

	#rm dups
	my $ind = 0;
	while (<$FREQ1>){
		chomp;
	        my @line = split /\s+/;

		my $count = 0;
		foreach(@line){
			if ($_ eq "otherallele_freq"){
				$ind = $count;
				#last;		
			}	
			$count++;		
		}

       		if (exists $all{$line[0]} )  {
        		print $OUT1 "$line[$ind]\n";
        	}
	}
	print "Frequency file pop1 done.\n";

	while (<$FREQ2>){
		chomp;
		my @line = split /\s+/;
       		if (exists $all{$line[0]} )  {
        		print $OUT2 "$line[$ind]\n";
        	}
	}
	print "Frequency file pop2 done.\n";

	while (<$FREQ3>){
		chomp;
		my @line = split /\s+/;
		if (exists $all{$line[0]} )  {
        		print $OUT3 "$line[$ind]\n";
	        }
	}
	print "Frequency file pop3 done.\n";

	return %all;
}




################lamp snp list
sub lampsnps{
	my ($MAP, $SNPS) = @_;
	my $count = 0;
	my %lines;
	my %snps;
	
	while (<$MAP>){
		chomp;
		my @line = split /\s+/;

		if($line[3] ne "0"){
			if(!exists $snps{$line[3]}){
       		 		$snps{$line[3]} = $line[3];
				$lines{$count} = $count;  
        			print $SNPS "$line[3]\n"; 
        		}
		}
		$count++;
	}

	print "SNP file done.\n";
	return %lines;
}



################lamp ref geno
sub lampgeno{
	my ($PED, $GENO, $INTER, $INTER0, $path, $SNPS, $lines_ref) = @_;
	my %lines = %$lines_ref;

	my $firstline = 1;
	my @refarr;
	while(<$PED>){
		my @line = split(/\s+/, $_);
		
		#remove 6 columns of header info
		for(my $i = 0; $i < 6; $i++){
			shift @line;
		}
	
		if($firstline == 1){
			for(my $i = 0; $i < @line ; $i+=2){
				$refarr[$i/2] = $line[$i];
			}
		}	

		for (my $i = 0; $i < @line - 1; $i+=2){
	  		if(exists $lines{$i/2}){
				if (($line[$i] eq $refarr[$i/2]) && ($line[$i+1] eq $refarr[$i/2]) ){
					print $INTER0 "0";
		 		} elsif(($line[$i] eq $refarr[$i/2]) && ($line[$i+1] ne $refarr[$i/2])  || ($line[$i] ne $refarr[$i/2] && $line[$i+1] eq $refarr[$i/2]) ){
					print $INTER0 "1";
		 		} elsif (($line[$i] ne $refarr[$i/2]) && ($line[$i+1] ne $refarr[$i/2])) {
					print $INTER0 "2";
		 		}else{
					print $INTER0 "9";
				}	
	  		}
		}	
		$firstline = 0;
		print $INTER0 "\n";
	}

	#rm all mono het sites
	seek ($INTER0, 0, 0);
	my @outline;
	my $lastcol;

	while (<$INTER0>){
		chomp;
		my @line = split "";

		$lastcol = 0;
		$lastcol = $#line if $#line > $lastcol;
		my $oldlastcol = $lastcol;

		for (my $i = $oldlastcol; $i < $lastcol; $i++) {
			if (defined($oldlastcol)){
				$outline[$i] = "" x $oldlastcol;
        		}
		}

		for (my $i=0; $i <=$lastcol; $i++) {
			$outline[$i] .= "$line[$i]"
		}	

	}

	#get mono het sites and rm, add line no. to hash
	my %rmmono;
	
	for (my $i=0; $i <= $lastcol; $i++) {
		my @line = split "", $outline[$i];
		foreach(my $j = 0; $j < @line; $j++){
			if ($line[$j] ne "1"){
				print $INTER $outline[$i]."\n";
				last;
			} 

			if ($j == scalar (@line - 1)){
				$rmmono{$i} = $i;	
			}
		}
	}

	#col back to row again
	seek ($INTER, 0, 0);
	my @outline2;
        my $lastcol2;

        while (<$INTER>){
                chomp;
                my @line = split "";

                $lastcol2 = 0;
                $lastcol2 = $#line if $#line > $lastcol2;
                my $oldlastcol = $lastcol2;

                for (my $i = $oldlastcol; $i < $lastcol2; $i++) {
                        if (defined($oldlastcol)){
                                $outline2[$i] = "\t" x $oldlastcol;
                        }
                }

                for (my $i=0; $i <=$lastcol2; $i++) {
			if($line[$i] eq "9"){
				$line[$i] = "-1";
			}
			
                        $outline2[$i] .= "$line[$i]\t"
                }

        }

	for (my $i=0; $i <= $lastcol2; $i++) {
             print $GENO $outline2[$i]."\n";                           
        }

        unlink "$path/interfile.txt" or print "oops: $!\n";
	unlink "$path/inter0.txt" or print "goshdarnit: $!\n";
	
	seek ($SNPS, 0, 0);
	my $cnt = 0;

	open my $POS, ">$path/chr.pos" or die "sorry\n"; 

	while (<$SNPS>){
		chomp;
		if (!exists $rmmono{$.}){
			print $POS "$_\n";
		}		
		$cnt++;
	}
		
	unlink "$path/chr_.pos" or print "couldnt del temp.pos b/c: $!\n";
	print "Genotype file done.\n";
}


################lamp config file
sub lampconfig{
	my ($CONFIG, $POPS) = @_;
					#no. of populations
	print $CONFIG "# Number of populations\npopulations=$POPS\n\n"
."# To use LAMP with ancestral allele frequencies, provide files for allele frequencies for the pure populations\n";

	#allele frequencies for the pure populations
	if ($POPS == 3){ print $CONFIG "pfile=freqs_pop1.txt,freqs_pop2.txt,freqs_pop3.txt\n\n#######################################################################\n";}
	else{ print $CONFIG "#pfile = freqs_pop1.txt,freqs_pop2.txt\n\n#######################################################################\n";}
	print $CONFIG "# Required files\n#######################################################################\n"
		
		#genotype file			    #SNP pos file
	."# Genotypes\ngenofile=geno.txt\n# SNP positions.\n";
	if ($POPS == 2){print $CONFIG "posfile=chr.pos\n"
	}else { print $CONFIG "posfile=chr_.pos\n"}
	print $CONFIG "# Output file of inferred ancestries.\n"

						#output file
	."# Defaults to 'ancestry.txt'\noutputancestryfile=ancestry.txt\n\n"
	."#######################################################################\n"
	."# Parameters for the windowing scheme\n#######################################################################\n\n"

	#offset for adjacent windows				      #recombination rate
	."# The offset of adjacent windows\noffset=0.2\n# Recombination rate\nrecombrate=1e-6\n";

	#no. of generations					#alpha
	if ($POPS == 3) { print $CONFIG "# Number of generations\ngenerations=10\n# Alpha (must sum to one)\nalpha=0.6,0.2,0.2\n\n\n"}
	else {
	print $CONFIG "# Number of generations\ngenerations=7\n# Alpha (must sum to one)\nalpha=0.2,0.8\n\n\n"
	}
	print $CONFIG "#######################################################################\n"
	."#######################################################################\n"

	#R^2 cutoff
	."# R^2 Cutoff to  use for pruning SNPs\nldcutoff=0.1";
	print "Configuration file done.\n";	
}


#****************************LAMP-LD****************************

################lamp-ld ref haps
sub ldhap{
	my ($INHAP, $OUTHAP, $count_ref) = @_;
	my %count = %$count_ref;

	while (<$INHAP>){
		chomp;

		my @orig = split /\s+/;
		shift(@orig);
		shift(@orig);	

		my @line = split("", $orig[0]);

		for (my $i = 0; $i < @line ; $i++){
			if ( exists $count{$i}){
				if($line[$i] eq "A" || $line[$i] eq "T" || $line[$i] eq "C" || $line[$i] eq "G"){
					if($line[$i] eq $count{$i}){
						print $OUTHAP "0";
					}else{ 
						print $OUTHAP "1";
					}
				}else{
					print $OUTHAP "?";
				}
			}
		}
		print $OUTHAP "\n";
	}
	close $INHAP;
	close $OUTHAP;
}


################lamp-ld snp pos
sub ldpos{
	my ($MAP, $SNPS, $POS) = @_;
	my %snps;
	my %subSnps;
	my %pc;
	my %hc;

	while(<$SNPS>){
		chomp;
		$snps{$_} = $_;
	}

	#subset positions and
	#line counts for .ped file
	my $mc = 0;
	while (<$MAP>){
		chomp;
		my @line = split /\s+/;
	
		if( exists $snps{$line[1]} && !exists $subSnps{$line[1]}){
			print $POS "$line[3]\n";
			$pc{$mc} = $mc;
			$subSnps{$line[1]} = $mc;
		}
		$mc++;
	}

	#line counts for hap file
	seek($SNPS, 0, 0) or die "can't go back: $!";

	my %pc_to_hc;
	my $sc=0;
	while (<$SNPS>){
		chomp;
		if(exists $subSnps{$_}){
			$hc{$sc} = $sc;
			$pc_to_hc{$subSnps{$_}} = $sc;
		}	
		$sc++;
	}

	print "Position file done.\n";		
	return (\%pc_to_hc, \%hc);
}


################lamp-ld ref geno
sub ldgeno{
	my ($PED, $GENO, $count_ref) = @_;
	my %count = %$count_ref;	
	
	my %hapCoding; #hash for coding in reference haplotype files

	while (<$PED>){
		chomp;
		my @line = split /\s+/;
		
		#rm headers
		for(my $i = 0; $i < 6; $i++){
			shift @line;
		}
	
		#do coding
		for (my $i = 0 ; $i < @line ; $i+=2){
			if (exists $count{$i/2}){
				if($. == 1){
					#hash for hap file
					$hapCoding{$count{$i/2}} = $line[$i];

					#store ref and alt allele
					$count{$i/2} = $line[$i];

					#letter to number
					if ($line[$i] eq $line[$i+1]){
						print $GENO "0";
					}elsif ($line[$i] ne $line[$i+1] && ($line[$i+1] eq "T" || $line[$i+1] eq "C" || $line[$i+1] eq "G" || $line[$i+1] eq "A" )){
						print $GENO "1";
					}else{
						print $GENO "?";
					}
				}else{
					if($line[$i] ne "A" || $line[$i] ne "T" || $line[$i] ne "C" || $line[$i] ne "G" || $line[$i+1] ne "A" || $line[$i+1] ne "T" || $line[$i+1] ne "C" || $line[$i+1] ne "G"){
						if($line[$i] eq $count{$i/2} && $line[$i] eq $count{$i/2}){
							print $GENO "0";
						}elsif($line[$i] ne $count{$i/2} && $line[$i+1] ne $count{$i/2}){
							print $GENO "2";
						}else{
							print $GENO "1";
						}
					}else{
						print $GENO "?";
					}
				}		
			}
		}
		print $GENO "\n"; 
	}
	print "Admixed genotype file done.\n";
	return (\%hapCoding);
}


#****************************ELAI****************************

################elai, admix geno and snp pos
sub admixandpos{
	my ($map, $ped, $pos, $geno) = @_;

	my $individuals;
	my $SNPs;
	my @genoArray;

	#lists no. of SNPs and rs no.
	while (<$map>){
		chomp;
		my @line = split /\s+/;
		print $pos "$line[1] $line[3] $line[0]\n";
		push (@genoArray, "$line[1]");
	
		$SNPs++;
	}
	print "Position file done.\n";

	#lists no. of individuals and geno
	while (<$ped>){
		chomp;
		$individuals += 1;
		my @line = split /\s+/;

		#rm headers if they're there
		for (my $i = 0; $i < 6; $i++){
			shift(@line);
		}

		#pair each 2 indices after that and delete old line
		my @newLine;
		my $lineCount=0;
		for (my $i = 0; $i < scalar @line; $i+=2){
			$newLine[$lineCount] = "$line[$i]$line[$i+1]";
			$lineCount++;
		}
	
		@line = ();
	
		#add new pairs to each SNP
		for (my $i = 0; $i < scalar @newLine; $i++){
			$genoArray[$i] = "$genoArray[$i] $newLine[$i]";
		}
	}
	
	#print statements
	say $geno "$individuals ";
	say $geno $SNPs;

	foreach (@genoArray){
		print $geno "$_\n";
	}

	print "Admixed genotype file done.\n";

	close $map;
	close $ped;
	close $geno;
}


################elai ref geno
sub sourcegeno{
	my ($hap, $snp, $pop) = @_;

	my $individuals;
	my $SNPs;
	my @genoArray;

	#lists no. of SNPs and rs no.
	while (<$snp>){
		chomp;
		push(@genoArray, $_);
		$SNPs++;
	}
	
	#lists no. of individuals and genotype
	my $count = 0;
	my @oldLine;

	while (<$hap>){
		chomp;
		$individuals += 1;
	
		#do ops depending on every 2 lines
		if ($count == 0){
			my @orig = split /\s+/;
			shift(@orig);
			shift(@orig);

			@oldLine = split("", $orig[0]);
		}
	
		if ($count == 1){
			my @orig = split /\s+/;
			shift(@orig);
			shift(@orig);

			my @line = split("", $orig[0]);


			#add new pairs to each SNP
			for (my $i = 0; $i < @line; $i++){
				$genoArray[$i] = "$genoArray[$i],$oldLine[$i]$line[$i]";
			}
		}
	
		$count++;
		if ($count == 2){
				$count = 0;
		}
	}

	say $pop $individuals/2 ." =" ;
	say $pop $SNPs;

	foreach (@genoArray){
		print $pop "$_\n";
	}
}



#****************************HAPMIX****************************
###############hapmix parameter
sub par{
	my ($PAR, $chrNo, $path) = @_;
	print $PAR
"GENOTYPE:1
OUTPUT_SITES:0
SITE_POSITIONS: 1 1000000000
THETA:0.2
LAMBDA:6.0
RECOMBINATION_VALS:600 900
MUTATION_VALS:0.2 0.2 0.01
MISCOPYING_VALS:0.05 0.05
REFPOP1GENOFILE:$path/pop1genofile.$chrNo
REFPOP2GENOFILE:$path/pop2genofile.$chrNo
REFPOP1SNPFILE:$path/pop1snpfile.$chrNo
REFPOP2SNPFILE:$path/pop2snpfile.$chrNo
ADMIXSNPFILE:$path/AAsnpfile.$chrNo
ADMIXGENOFILE:$path/AAgenofile.$chrNo
REF1LABEL:CEU
REF2LABEL:YRI
RATESFILE: $path/rates.$chrNo
ADMIXPOP: AA
CHR:$chrNo
OUTDIR:$path/RUN
HAPMIX_MODE:LOCAL_ANC
OUTPUT_DETAILS:ANC_INT_THRESH
THRESHOLD:0.0
KEEPINTFILES:0";

	print "Parameter file done.\n";
}


################hapmix ref haps
sub hapmixrefhaps{
        my ($SNPS, $MAP, $HAP, $INTER, $OUTHAP, $path, $str, $maxPos, $minPos, $coding_ref, $crossSNPsref) = @_;
	my %coding = %$coding_ref;
	my %crossSNPs = %$crossSNPsref;

	#now that we have hash %crossSNPs w/ keys
	# and values as no. of column of SNP that we want in .ped, 
	#we skip each outline that is not the key/value that we want

	while (<$HAP>) {
		chomp;

		my @orig = split /\s+/;
		shift(@orig);
		shift(@orig);

		my @line = split("", $orig[0]);

		for (my $i = 0; $i < scalar(@line) ; $i++){

			if (exists $crossSNPs{$i}){
				if($line[$i] eq "A" || $line[$i] eq "T" || $line[$i] eq "C" || $line[$i] eq "G"){
					if($line[$i] eq $coding{$i}){
						print $INTER "0";
					}else{ 
						print $INTER "1";
					}
				}else{
					print $INTER "?";
				}
			}
        	}
		print $INTER "\n";
	}
		

	close $INTER;
	open $INTER, "<$path/$str" or die "error";

	#seek($INTER,0,0); instead^

	my @outline;
	my $lastcol;

	while (<$INTER>){
		chomp;
		my @line = split "";

		$lastcol = 0;
		$lastcol = $#line if $#line > $lastcol;
		my $oldlastcol = $lastcol;

		for (my $i = $oldlastcol; $i < $lastcol; $i++) {
			if (defined($oldlastcol)){
				$outline[$i] = "\t" x $oldlastcol;
        		}
		}

		for (my $i=0; $i <=$lastcol; $i++) {
			$outline[$i] .= "$line[$i]\t"
		}	

	}

	#rm whitespace and print
	for (my $i=0; $i <= $lastcol; $i++) {
		$outline[$i] =~ s/\s+//g;
		print $OUTHAP $outline[$i]."\n";                                                    
	}

	close $INTER;
	unlink "$path/$str" or print "can't delete ref interfile at $INTER: $!\n";
}


################hapmix admixed geno
sub hapmixadmixgeno{
	my ($SNPS, $MAP, $PED, $INTER, $GENO, $path, $maxPos, $minPos) = @_; 
	
	#read in and put admixed SNPs in hash
	my %admixSNP;
	my $lc = 0;
	while (<$SNPS>){
		chomp;
        	$admixSNP{$_} = $lc;
		$lc++;
	}

	#cross admixed SNPs with those in simu .map 
	my %crossSNPs;
	my $pedCol = 0;
	my %snpPosInHapfile;
	my %AAtoHap;

	while(<$MAP>){
        	chomp;
        	my @line = split /\s+/;
        	if ($line[3] < $maxPos && $line[3] > $minPos){ #check with rate file
        		if (exists $admixSNP{$line[1]}){
				$crossSNPs{$pedCol} = $pedCol;
				$snpPosInHapfile{$admixSNP{$line[1]}} = $admixSNP{$line[1]};
				$AAtoHap{$pedCol} = $admixSNP{$line[1]};
               		}
		}
        	$pedCol++;
	}

	#now that we have hash %crossSNPs w/ keys and 
	# values as no. of column of SNP that we want in .ped, 
	#we skip each outline that is not the key/value that we want

	my %coding; #ref/alt allele coding
	my %refCoding; #coding for reference haps

	while (<$PED>) {        
		chomp;
		my @line = split /\s+/;

		for(my $i = 0; $i < 6; $i++){
			shift @line;
		}

		#change letters to numbers
		for (my $i = 0 ; $i < @line ; $i+=2){
			if (exists $crossSNPs{$i/2}){
				if($. == 1){
					$coding{$i/2} = $line[$i]; #this is the line in the admixed data
					#need another hash for the reference data
					$refCoding{$AAtoHap{$i/2}} = $line[$i];

					if($line[$i] eq $line[$i+1]){
						print $INTER "0";
					}elsif ($line[$i] ne $line[$i+1] && ($line[$i+1] eq "T" || $line[$i+1] eq "C" || $line[$i+1] eq "G" || $line[$i+1] eq "A" )){
						print $INTER "1";
					}else{
						print $INTER "?";
					}
				}else{
					if($line[$i] ne "A" || $line[$i] ne "T" || $line[$i] ne "C" || $line[$i] ne "G" || $line[$i+1] ne "A" || $line[$i+1] ne "T" || $line[$i+1] ne "C" || $line[$i+1] ne "G"){
						if($line[$i] eq $coding{$i/2} && $line[$i+1] eq $coding{$i/2}){
							print $INTER "0";
						}elsif($line[$i] ne $coding{$i/2} && $line[$i+1] ne $coding{$i/2}){
							print $INTER "2";
						}else{
							print $INTER "1";
						}
					}else{
						print $INTER "?";
					}
				}		
			}
		}
		print $INTER "\n";
	}

	close $INTER;
	open $INTER, "$path/AAinterfile.txt" or die "erorr";

	my @outline;
	my $lastcol;

	while (<$INTER>){
        	chomp;
	        my @line = split "";
       
        	$lastcol = 0;
	        $lastcol = $#line if $#line > $lastcol;
	        my $oldlastcol = $lastcol;

	        for (my $i = $oldlastcol; $i < $lastcol; $i++) {
        	        $outline[$i] = "\t" x $oldlastcol;
	        }

        	for (my $i=0; $i <=$lastcol; $i++) {
	                $outline[$i] .= "$line[$i]"
        	        }
       	 }


	for (my $i=0; $i <= $lastcol; $i++) {
        	print $GENO $outline[$i]."\n";
	}

	close $INTER;
	unlink "$path/AAinterfile.txt" or print "can't delete AA interfile\n";
	return (\%refCoding, \%snpPosInHapfile);
}





################hapmix recombination rates
sub hapmixrec{
	my ($RECMAP, $RATES) = @_;
	
	#read in RECMAPFILE, ignore first line
	my @outline;
	my $lineCount = 1;
	while (<$RECMAP>) {
		if ($lineCount > 1){
			chomp;
			my @line = split " ";	

			#columns to rows
			my $lastcol = 0;
			$lastcol = $#line if $#line > $lastcol;
			my $oldlastcol = $lastcol;

			for (my $i = $oldlastcol; $i < $lastcol; $i++) {
				$outline[$i] = "\t" x $oldlastcol;
			}

			for (my $i=0; $i <=$lastcol; $i++) {
				$outline[$i] .= "$line[$i] " #<-change to space or \t or w/e
			}
		}
		$lineCount++;
	}	

	print $RATES ":sites:";
	print $RATES ($lineCount - 2) . "\n";
	print $RATES $outline[0]."\n";
	print $RATES $outline[2]."\n";
}


################hapmix admixed snps
sub hapmixadmixsnps{
	my ($MAP, $SNPS, $OUT, $OUT1, $OUT2, $maxPos, $minPos) = @_;
	#read in SNP file and put each value into %SNPhash, with key and value the same
	my %SNPhash;
	while (<$SNPS>){
		chomp;
		$SNPhash{$_} = $_;
	}

	#test if SNP is in mapfile and if so, print
	while (<$MAP>){
		chomp;
		
		my @line = split /\s+/;
		
		if ($line[3] < $maxPos && $line[3] > $minPos){
			if (exists($SNPhash{$line[1]})) {
			 	print $OUT "\t$line[1]\t$line[0]\t$line[2]\t$line[3]\n";
				print $OUT1 "\t$line[1]\t$line[0]\t$line[2]\t$line[3]\n";
				print $OUT2 "\t$line[1]\t$line[0]\t$line[2]\t$line[3]\n";
			}
		}
	}	
}
