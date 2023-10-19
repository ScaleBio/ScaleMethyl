#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("F:C:O:B:", \%opt);

#defaults


$die = "
Usage:

sciMET_meth2mtx.pl [options] -B [window bed file] -O [output prefix] -F [methylation chrom folder] ...

Options:

-F   [DIR]   Directory with methylation by cell for each chromosome (req)
             Can be a comma separated list to process multiple experiments.
             Or can be a file ending in 'txt' with aline for each folder
-B   [STR]   Bed file of windows for computing methylation (req)
-O   [STR]   Output prefix (req)

-C   [STR]   File of list of passing cellIDs (will filter to them)

";

# PARSE OPTS

if (!defined $opt{'F'} ||
	!defined $opt{'O'} ||
	!defined $opt{'B'}) {die $die};


if (defined $opt{'C'}) {
	open IN, "$opt{'C'}";
	while ($l = <IN>) {
		chomp $l;
		$l =~ s/\s.+$//;
		$PASSING_CELLS{$l} = 1;
	} close IN;
}

$winID = 0;
open BED, "$opt{'B'}";
while ($l = <BED>) {
	chomp $l;
	@P = split(/\t/, $l);
	$win_name = "$P[0]_$P[1]_$P[2]";
	$winID++;
	$WINID_name{$winID} = $win_name;
	$WINNAME_winID{$win_name} = $winID;
	$CHR_START_win{$P[0]}{$P[1]} = $winID;
	$WINID_end{$winID} = $P[2];

	$CHR_HANDLES{$P[0]} = 1;
	
}
$lastWin = $winID;

$CONTEXT_total{'H'}=0;
$CONTEXT_total{'h'}=0;
$CONTEXT_total{'X'}=0;
$CONTEXT_total{'x'}=0;
$CONTEXT_total{'Z'}=0;
$CONTEXT_total{'z'}=0;
$not_in_window=0;
$win_with_cov=0;
@CELLIDs = ();

if ($opt{'F'} =~ /txt$/) {
	@FOLDERS = ();
	open F, "$opt{'F'}";
	while ($l = <F>) {chomp $l; push @FOLDERS, $l};
	close F;
} else {
	@FOLDERS = split(/,/, $opt{'F'});
}
foreach $chr_dir (@FOLDERS) {
	foreach $chr (keys %CHR_HANDLES) {
		
		if (-e "$chr_dir/$chr.bed.gz") {
			$ts = localtime(time);
			print STDERR "$ts\tFound $chr_dir/$chr.bed.gz, processing...\n";
			
			open INT, "bedtools intersect -a $chr_dir/$chr.bed.gz -b $opt{'B'} -wa -wb |";
			# 0   1    2     3      4   5    6    7
			#chr pos1 pos2 cellID meth chr start end ...
			
			while ($l = <INT>) {
				chomp $l;
				@P = split(/\t/, $l);
				$cellID = $P[3];
				if (!defined $opt{'C'} || defined $PASSING_CELLS{$cellID}) {
					$winName = "$P[5]_$P[6]_$P[7]";
					$winID = $WINNAME_winID{$winName};
					if (!defined $WINID_totalCov{$winID}) {
						$win_with_cov++;
					}
					$WINID_totalCov{$winID}++;
					
					if (!defined $WINID_CELLID_covstr{$winID}{$cellID}) {
						$WINID_cellCT{$winID}++;
					}
					$WINID_totalCov{$winID}++;
					$WINID_CELLID_covstr{$winID}{$cellID} .= $P[4];
					

					if (!defined $CELLID_totalCov{$cellID}) {
						$CELLID_totalCov{$cellID} = 1;
						$cellCT++;
						push @CELLIDs, $cellID;
						$CELLID_CONTEXT_cov{$cellID}{'H'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'h'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'X'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'x'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'Z'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'z'} = 0;
					} else {
						$CELLID_totalCov{$cellID}++;
					}
					$CELLID_CONTEXT_cov{$cellID}{$P[4]}++;
					$CONTEXT_total{$P[4]}++;
					
				}
				
			} close INT;
			
		} else {
			print STDERR "WARNING: Cannot find $chr_dir/$chr.bed.gz!\n";
		}
	}
}

$CH_status = 0;
$CG_status = 0;

if ($CONTEXT_total{'h'} > 0 ||
	$CONTEXT_total{'H'} > 0 ||
	$CONTEXT_total{'X'} > 0 ||
	$CONTEXT_total{'x'} > 0) {
	$CH_status = 1;
}

if ($CONTEXT_total{'z'} > 0 ||
	$CONTEXT_total{'Z'} > 0) {
	
	$CG_status = 1;
}

$header = "";
for ($i = 0; $i < @CELLIDs; $i++) {
	$cellID = $CELLIDs[$i];
	$header .= "$cellID\t";
	
	# global meth calc
	if ($CH_status>0) {
		if (($CELLID_CONTEXT_cov{$cellID}{'h'}+$CELLID_CONTEXT_cov{$cellID}{'x'}+$CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})<1) {$CELLID_CONTEXT_cov{$cellID}{'h'}=100};
		if (($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})<1) {$CELLID_CONTEXT_cov{$cellID}{'H'}=1}; # pseudocount if 0
		$CELLID_mCH{$cellID} = sprintf("%.3f", ($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})/($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'}+$CELLID_CONTEXT_cov{$cellID}{'h'}+$CELLID_CONTEXT_cov{$cellID}{'x'}));
	}
	
	if ($CG_status>0) {
		if (($CELLID_CONTEXT_cov{$cellID}{'z'}+$CELLID_CONTEXT_cov{$cellID}{'Z'})<1) {$CELLID_CONTEXT_cov{$cellID}{'z'}=10};
		if (($CELLID_CONTEXT_cov{$cellID}{'Z'})<1) {$CELLID_CONTEXT_cov{$cellID}{'Z'}=1}; # pseudocount if 0
		$CELLID_mCG{$cellID} = sprintf("%.3f", ($CELLID_CONTEXT_cov{$cellID}{'Z'})/($CELLID_CONTEXT_cov{$cellID}{'Z'}+$CELLID_CONTEXT_cov{$cellID}{'z'}));
	}

}
$header =~ s/\t$//;

if ($CH_status>0) {
	open CHMeth, ">$opt{'O'}.CH.mtx";
	open CHRat, ">$opt{'O'}.CH.ratio.mtx";
	open CHCov, ">$opt{'O'}.CH.cov.mtx";
	open CHDiff, ">$opt{'O'}.CH.score.mtx";
	print CHMeth "$header\n"; print CHRat "$header\n"; print CHCov "$header\n"; print CHDiff "$header\n";
	for ($winID = 1; $winID <= $lastWin; $winID++) {
		print CHMeth "$WINID_name{$winID}";
		print CHRat "$WINID_name{$winID}";
		print CHCov "$WINID_name{$winID}";
		print CHDiff "$WINID_name{$winID}";
		foreach $cellID (@CELLIDs) {
			$X = $WINID_CELLID_covstr{$winID}{$cellID} =~ tr/X//;
			$x = $WINID_CELLID_covstr{$winID}{$cellID} =~ tr/x//;
			$H = $WINID_CELLID_covstr{$winID}{$cellID} =~ tr/H//;
			$h = $WINID_CELLID_covstr{$winID}{$cellID} =~ tr/h//;
			$cov = $X+$x+$H+$h;
			if ($cov > 1) {
				$meth = sprintf("%.3f", ($X+$H)/$cov);
				$rat = sprintf("%.3f", $meth/$CELLID_mCH{$cellID});
				$diff = $meth - $CELLID_mCH{$cellID};
				if ($diff>0) {
					$score = sprintf("%.3f", $diff/(1-$CELLID_mCH{$cellID}));
				} else {
					$score = sprintf("%.3f", $diff/$CELLID_mCH{$cellID});
				}
			} else {
				$meth = 0;
				$rat = 1;
				$score = 0;
			}
			print CHMeth "\t$meth";
			print CHRat "\t$rat";
			print CHCov "\t$cov";
			print CHDiff "\t$score";
		}
		print CHMeth "\n"; print CHRat "\n"; print CHCov "\n"; print CHDiff "\n";
	}
	close CHMeth; close CHRat; close CHCov; close CHDiff;
}

if ($CG_status>0) {
	open CGMeth, ">$opt{'O'}.CG.mtx";
	open CGRat, ">$opt{'O'}.CG.ratio.mtx";
	open CGCov, ">$opt{'O'}.CG.cov.mtx";
	open CGDiff, ">$opt{'O'}.CG.score.mtx";
	print CGMeth "$header\n"; print CGRat "$header\n"; print CGCov "$header\n"; print CGDiff "$header\n";
	for ($winID = 1; $winID <= $lastWin; $winID++) {
		print CGMeth "$WINID_name{$winID}";
		print CGRat "$WINID_name{$winID}";
		print CGCov "$WINID_name{$winID}";
		print CGDiff "$WINID_name{$winID}";
		foreach $cellID (@CELLIDs) {
			$Z = $WINID_CELLID_covstr{$winID}{$cellID} =~ tr/Z//;
			$z = $WINID_CELLID_covstr{$winID}{$cellID} =~ tr/z//;
			$cov = $Z+$z;
			if ($cov > 1) {
				$meth = sprintf("%.3f", $Z/$cov);
				$rat = sprintf("%.3f", $meth/$CELLID_mCG{$cellID});
				$diff = $meth - $CELLID_mCG{$cellID};
				if ($diff>0) {
					$score = sprintf("%.3f", $diff/(1-$CELLID_mCG{$cellID}));
				} else {
					$score = sprintf("%.3f", $diff/$CELLID_mCG{$cellID});
				}
			} else {
				$meth = 0;
				$rat = 1;
				$score = 0;
			}
			print CGMeth "\t$meth";
			print CGRat "\t$rat";
			print CGCov "\t$cov";
			print CGDiff "\t$score";
		}
		print CGMeth "\n"; print CGRat "\n"; print CGCov "\n"; print CGDiff "\n";
	}
	close CGMeth; close CGRat; close CGCov; close CGDiff;
}

exit;
