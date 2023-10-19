#!/usr/bin/perl

$die = "

sciMET_mergeMtx.pl [output name] [matrix 1] [matrix 2] (matrix 3) ...

Must have no overlapping cells unless adding annot. Must be same rows.

Can specify [output annot file],[output matrix] [matrix1 annot],[matrix 1] etc...
This will add the annot ID to the cell name and also
output an annot file with new names and annots.

Designed for merging the cells split into chroms
files during the extract process or across experiments
if the annot oiptions are used.

";

if (!defined $ARGV[2]) {die $die};

if ($ARGV[0] =~ /,/) {
	($out_annot,$out_mtx) = split(/,/, $ARGV[0]);
	open ANNOT, ">$out_annot";
} else {
	$out_mtx = "$ARGV[0]";
}

open OUT, ">$out_mtx";
$header_out = "";
for ($i = 1; $i < @ARGV; $i++) {
	if (defined $out_annot) {
		if ($ARGV[$i] !~ /,/) {die "ERROR: When specifying an annot output, each matrix must have an accompanying annot name!\n"};
		($annot,$in_mtx) = split(/,/, $ARGV[$i]);
		$ID_annot{$i} = $annot;
	} else {
		if ($ARGV[$i] =~ /,/) {
			($annot,$in_mtx) = split(/,/, $ARGV[$i]);
			$ID_annot{$i} = $annot;
		} else {
			$in_mtx = $ARGV[$i];
		}
	}
	open $i, "$in_mtx";
	$l = <$i>;
	chomp $l;
	if (defined $out_annot) {
		@CELLS = split(/\t/, $l);
		for ($j = 0; $j < @CELLS; $j++) {
			$header_out .= "$ID_annot{$i}_$CELLS[$j]\t";
			if (defined $out_annot) {
				print ANNOT "$ID_annot{$i}_$CELLS[$j]\t$ID_annot{$i}\n";
			}
		}
	} else {
		$header_out .= $l."\t";
	}
}
if (defined $out_annot) {close ANNOT};
$header_out =~ s/\t$//;
print OUT "$header_out\n";

$first = 1;
while ($l = <$first>) {
	chomp $l;
	$out_line = $l."\t";
	for ($i = 2; $i < @ARGV; $i++) {
		$l = <$i>;
		chomp $l;
		@P = split(/\t/, $l);
		shift @P;
		$l = join("\t", @P);
		$out_line .= "$l\t";
	}
	$out_line =~ s/\t$//;
	print OUT "$out_line\n";
}

close OUT;
