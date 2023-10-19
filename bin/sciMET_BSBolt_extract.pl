#!/usr/bin/perl

BEGIN {
	use FindBin '$RealBin';
	push @INC, $RealBin;
}

use Getopt::Std; %opt = ();

getopts("O:t:sm:G:C:", \%opt);

# dfaults
$minSize = 10000000;

$die = "

sciMET_BSBolt_extract.pl (options) [rmdup & filtered bam file, name-sorted]

Extracts methylation from a BSBolt bam and outputs bismark annotated
calls in the chrom folder format.

Options:
   -O   [STR]   Out prefix
   -m   [INT]   Minimum chromosome size to retain (def = $minSize)
                  Used to exclude random and other small contigs
   -G   [BED]   Bed file of GpCs - will create additional folders
                  for GpC or HpC contexts.
   -t   [INT]   Max number of concurrent threads (def = 1)
   -C   [STR]   List of cellIDs to include (use instead of pre-
                  filtering bam file)
   
";

if (!defined $ARGV[0]) {die $die};

if (defined $opt{'C'}) {
	open IN, "$opt{'C'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$INCLUDE{$cellID} = 1;
	} close IN;
}

if (!defined $opt{'s'}) { # main thread
	#defaults
	if (!defined $opt{'t'}) {$opt{'t'} = 1};
	%CALL_CONV = ("x" => "z", "X" => "Z", "y" => "x", "Y" => "X", "z" => "h", "Z" => "H");
	
	# open log file
	#open LOG, ">$opt{'O'}.log";
	$ts = localtime(time);
	print "$ts\tProgram called.\n\n============== Phase 1: Parsing bam and processing per-cell calls ==============\n\n";

	# setup output directory
	system("mkdir $opt{'O'}");
	$ts = localtime(time);
	print "$ts\tOutput directory created: $opt{'O'}\n";
	
	# setup threading stuff
	%QUEUE = ();
	
	# start parsing bam
	$ts = localtime(time);
	print "$ts\tParsing bam file: $ARGV[0]\n";
	
	open IN, "samtools view $ARGV[0] |";
	$currentBarc = "null";
	$methCol = "null";
	
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$barc = $P[0]; $barc =~ s/:.+$//;
		
		# check if barc is included, ortherwise skip line
		if (!defined $opt{'C'} || defined $INCLUDE{$barc}) {
		
			# check to close out barcode
			if ($currentBarc ne $barc) {
				# close out prev barc
				if ($currentBarc ne "null") {
					close OUT;
					$QUEUE{$currentBarc} = 0;
				}
				# check thread stats
				check_threads();
				
				# set next barcode & open outfile
				$currentBarc = $barc;
				open OUT, ">$opt{'O'}/$currentBarc.meth";
			}
			# parse read
			$chr = $P[2];
			$pos = $P[3];
			# ID column in bam that has the meth col if not yet found
			if ($methCol eq "null") {
				for ($i = 11; $i < @P; $i++) {
					if ($P[$i] =~ /^XB/) {
						$methCol = $i;
						$i+=100;
					}
				}
			}
			$meth = $P[$methCol];
			$meth =~ s/^XB:Z://;
			print OUT "$chr\t$pos\t$meth\n";
			
		}
	} close IN;
	
	# finish last one
	close OUT;
	$QUEUE{$currentBarc} = 0;
	
	# process rest of queue
	$incomplete = 1;
	while ($incomplete > 0) {
		$incomplete = 0;
		sleep(10); # wait
		check_threads();
		foreach $queued_barc (keys %QUEUE) {
			if ($QUEUE{$queued_barc} < 2) { # active or waiting
				$incomplete++;
			}
		}
	}
	$ts = localtime(time);
	print "$ts\tAll threads complete!\n\n============== Phase 2: Generating chrom folders and summary ==============\n\n";
	#done with threads!
	
	# cleanup and finish processing......
	
	# first info file
	open OUT, ">$opt{'O'}.cellInfo.txt";
	print OUT "#CellID\tCoverage\tCG_Cov\tCG_mC_Pct\tCH_Cov\tCH_mC_Pct\n";
	close OUT;
	foreach $barc (keys %QUEUE) {
		system("cat $opt{'O'}/$barc.info >> $opt{'O'}.cellInfo.txt");
	}
	
	# now chrom methylaiton files...
	# get chromosome names to prep for making files
	$chromCT = 0;
	open HEADER, "samtools view -H $ARGV[0] |";
	while ($l = <HEADER>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[0] =~ /^\@SQ/) {
			$P[2] =~ s/^LN://;
			if ($P[2] >= $minSize) {
				$P[1] =~ s/^SN://;
				$CHROMS{$P[1]} = "$P[1]";
				$chromCT++;
			}
		}
	} close HEADER;
	$ts = localtime(time);
	print "$ts\tFound $chromCT chromosomes that met min size requirement.\n";
	
	# read in GC sites
	if (defined $opt{'G'}) {
		$ts = localtime(time);
		print "$ts\tReading in GC sites from $opt{'G'}.\n";
		if ($opt{'G'} =~ /gz$/) {
			open GC, "zcat $opt{'G'} |";
		} else {
			open GC, "$opt{'G'}";
		}
		while ($l = <GC>) {
			chomp $l;
			($chrom,$start,$end,$strand) = split(/\t/, $l);
			$GPC_CHR_POS{$chrom}{$start} = $strand;
			$inGpC++;
		} close GC;
		$ts = localtime(time);
		print "\t$ts\tDone, $inGpC found.\n";
	}
	
	# open out folders
	if (defined $opt{'G'}) {
		system("mkdir $opt{'O'}.HCH.chroms");
		system("mkdir $opt{'O'}.GCH.chroms");
		system("mkdir $opt{'O'}.HCG.chroms");
		system("mkdir $opt{'O'}.GCG.chroms");
		foreach $chrom (keys %CHROMS) {
			$handle = $CHROMS{$chrom}."HCH";
			$HANDLES{$handle} = 1;
			open $handle, ">$opt{'O'}.HCH.chroms/$chrom.bed";
			$handle = $CHROMS{$chrom}."GCH";
			$HANDLES{$handle} = 1;
			open $handle, ">$opt{'O'}.GCH.chroms/$chrom.bed";
			$handle = $CHROMS{$chrom}."HCG";
			$HANDLES{$handle} = 1;
			open $handle, ">$opt{'O'}.HCG.chroms/$chrom.bed";
			$handle = $CHROMS{$chrom}."GCG";
			$HANDLES{$handle} = 1;
			open $handle, ">$opt{'O'}.GCG.chroms/$chrom.bed";
		}
	} else {
		system("mkdir $opt{'O'}.CH.chroms");
		system("mkdir $opt{'O'}.CG.chroms");
		foreach $chrom (keys %CHROMS) {
			$handle = $CHROMS{$chrom}."CH";
			$HANDLES{$handle} = 1;
			open $handle, ">$opt{'O'}.CH.chroms/$chrom.bed";
			$handle = $CHROMS{$chrom}."CG";
			$HANDLES{$handle} = 1;
			open $handle, ">$opt{'O'}.CG.chroms/$chrom.bed";
		}
	}
	
	# print to files
	foreach $barc (keys %QUEUE) {
		open IN, "$opt{'O'}/$barc.calls";
		while ($l = <IN>) {
			chomp $l;
			($chrom,$start,$call) = split(/\t/, $l);
			if (defined $opt{'G'}) {
				if (defined $GPC_CHR_POS{$chrom}{$start}) {
					if ($call =~ /x/i) { # CG
						$handle = $CHROMS{$chrom}."GCG";
					} else { # CH
						$handle = $CHROMS{$chrom}."GCH";
					}
				} else {
					if ($call =~ /x/i) { # CG
						$handle = $CHROMS{$chrom}."HCG";
					} else { # CH
						$handle = $CHROMS{$chrom}."HCH";
					}
				}
			} else {
				if ($call =~ /x/i) { # CG
					$handle = $CHROMS{$chrom}."CG";
				} else { # CH
					$handle = $CHROMS{$chrom}."CH";
				}
			}
			
			if (defined $CALL_CONV{$call}) {
				print $handle "$chrom\t$start\t$start\t$barc\t$CALL_CONV{$call}\n";
			}
		}
	}
	
	# close folder files
	foreach $handle (keys %HANDLES) {
		close $handle;
	}
	
	$ts = localtime(time);
	print "$ts\tDONE!\n"
	
} else { # subthread
	open IN, "$ARGV[0]/$ARGV[1].meth";
	while ($l = <IN>) {
		chomp $l;
		($chr,$pos,$meth) = split(/\t/, $l);
		@M = split(//, $meth);
		for ($i = 0; $i < @M; $i++) {
			if ($M[$i] !~ /[0-9]/) {
				$pos++;
				$coord = "$chr\t$pos";
				if (!defined $COORD_meth{$coord}) { # only new cov
					$cov++;
					$COORD_meth{$coord} = $M[$i];
					$METH_ct{$M[$i]}++;
				}
			} elsif ($M[$i] =~ /[0-9]/) {
				$add = $M[$i];
				while ($M[$i+1] !~ /xyz/i && $M[$i+1] =~ /[0-9]/) {
					$i++;
					$add .= $M[$i];
				}
				$pos += $add;
			}
		}
	} close IN;
	
	open OUT, ">$ARGV[0]/$ARGV[1].calls";
	foreach $coord (keys %COORD_meth) {
		print OUT "$coord\t$COORD_meth{$coord}\n";
	} close OUT;
	
	$CG_cov = $METH_ct{"x"} + $METH_ct{"X"};
	$CH_cov = $METH_ct{"y"} + $METH_ct{"Y"} + $METH_ct{"z"} + $METH_ct{"Z"};
	if ($CG_cov > 0) {
		$CG_pct = sprintf("%.2f", ($METH_ct{"X"}/$CG_cov)*100);
	} else {$CG_pct = "0.00"};
	if ($CH_cov > 0) {
		$CH_pct = sprintf("%.2f", (($METH_ct{"Y"} + $METH_ct{"Z"})/$CH_cov)*100);
	} else {$CH_pct  = "0.00"};
	
	open OUT, ">$ARGV[0]/$ARGV[1].info";
	print OUT "$ARGV[1]\t$cov\t$CG_cov\t$CG_pct\t$CH_cov\t$CH_pct\n";
	close OUT;
	
	exit;
}

sub check_threads {
	$running = 0; @WAITING = ();
	foreach $queued_barc (keys %QUEUE) {
		if ($QUEUE{$queued_barc} == 1) { # active
			$running++;
			if (-e "$opt{'O'}/$queued_barc.info") { # finished
				$running--;
				$QUEUE{$queued_barc} = 2;
				print "\t\t\t$queued_barc COMPLETED!\n";
			}
		} elsif ($QUEUE{$queued_barc} == 0) { # store to kick off
			push @WAITING, $queued_barc;
		}
	}
	if ($running < $opt{'t'}) { # if not at full threading
		$ts = localtime(time);
#		print LOG "\t$ts\tThreads active: $running, adding...\n";
		for ($waitingID = 0; $waitingID < ($opt{'t'} - $running); $waitingID++) { # see how many short
			if (defined $QUEUE{$WAITING[$waitingID]}) { # if wiating is occupied
				system("$RealBin/sciMET_BSBolt_extract.pl -s $opt{'O'} $WAITING[$waitingID] &"); # start the thread
				$QUEUE{$WAITING[$waitingID]} = 1; # log it as active
				print "\t\t$WAITING[$waitingID] now running.\n";
			}
		}
	}
}