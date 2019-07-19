#!/usr/bin/perl 

#6/29/12

$path = $ARGV[0];
open (BP, "<$path/Human_BP_pairwise_spec_sens_score.txt") || die "cannot open BP";
open (CC, "<$path/Human_CC_pairwise_spec_sens_score.txt")|| die "cannot open CC";
open (MF, "<$path/Human_MF_pairwise_spec_sens_score.txt")|| die "cannot open MF";


open (OUT, ">Human_funsim_scores.txt");
open (OUT_BP_MF, ">Human_funsim_scores_BP+MF.txt");

my %bp = ();
$bpmax;
while ($line = <BP>){
	chomp $line;
	#print "get_funsim2.pl : BP $line\n";
	@line = ();
	@line = split /\t/, $line;

	$prot1 = $line[0];
	$prot2 = $line[1];

	if ($line[2] > $line[3]){
		$bpmax = $line[2];
	}
	if ($line[3] > $line[2]){
		$bpmax = $line[3];
	}
	if ($line[2] == $line[3]){
		$bpmax = $line[2];
	}
	
	$line = <MF>;
	chomp $line;
	#print "get_funsim2.pl : MF $line\n";
	@line = ();
	@line = split /\t/, $line;

	if ($line[2] > $line[3]){
		$mfmax = $line[2];
	}
	if ($line[3] > $line[2]){
		$mfmax = $line[3];
	}
	if ($line[2] == $line[3]){
		$mfmax = $line[2];
	}

	$line = <CC>;
	chomp $line;
	#print "get_funsim2.pl : CC $line\n";
	@line = ();
	@line = split /\t/, $line;

	if ($line[2] > $line[3]){
		$ccmax = $line[2];
	}
	if ($line[3] > $line[2]){
		$ccmax = $line[3];
	}
	if ($line[2] == $line[3]){
		$ccmax = $line[2];
	}

	#compute funsim

	$bpsq = $bpmax * $bpmax;
	$mfsq = $mfmax * $mfmax;
	$ccsq = $ccmax * $ccmax;

	$gettingThere = $bpsq + $ccsq + $mfsq;
	$funsim = $gettingThere * (1/3);
	$bpmf = ($bpsq + $mfsq) * (1/2);

	print OUT "$prot1\t$prot2\t$funsim\n";
	print OUT_BP_MF "$prot1\t$prot2\t$bpmf\n";
		
}
print "Done Funsim\n";

close(BP);
close(MF);
close(CC);

close(OUT);

close(OUT_BP_MF);



