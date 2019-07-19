#!/usr/bin/perl

my %avgprec = ();
my %avgrecall = (); 
my %avgfscore = (); 
my @unknownprot = ();
my %prec = ();
my %recall = ();
my %fscore = ();
my %tp = ();

$iter = $ARGV[0];

open (BP, "<data/BPparents.txt") || die "Cannot open BP\n";
open (MF, "<data/MFparents.txt") || die "Cannot open MF\n";
open (CC, "<data/CCparents.txt") || die "Cannot open CC\n";

%parents = ();

while (<BP>)
{

	chomp $_;
	@line = split(" ", $_);		#GO:0000015	GO:0000015,GO:0005575,GO:0005622,GO:0005623,GO:0005737,GO:0005829,GO:0032991,GO:0043234,GO:0044424,GO:0044444,GO:0044445,GO:0044464,
	$parents{$line[0]} = $line[1];
}
close(BP);
print "accuracy_relaxed.pl : Done loading BP parents\n";

$count=1;
while (<MF>)
{

	chomp $_;
	@line = split(" ", $_);		#GO:0000015	GO:0000015,GO:0005575,GO:0005622,GO:0005623,GO:0005737,GO:0005829,GO:0032991,GO:0043234,GO:0044424,GO:0044444,GO:0044445,GO:0044464,
	$parents{$line[0]} = $line[1];
}
close(MF);
print "accuracy_relaxed.pl : Done loading MF parents\n";

$count=1;
while (<CC>)
{

	chomp $_;
	@line = split(" ", $_);		#GO:0000015	GO:0000015,GO:0005575,GO:0005622,GO:0005623,GO:0005737,GO:0005829,GO:0032991,GO:0043234,GO:0044424,GO:0044444,GO:0044445,GO:0044464,
	$parents{$line[0]} = $line[1];
}
close(CC);
print "accuracy_relaxed.pl : Done loading CC parents\n";

my %true = ();
	
open(TRUE, "</net/kihara/khan27/Desktop/Ishita_Projects/Group_Function_Prediction/GO_enrichment/CustomSlim/RA-MAPK_full_enriched_annot_CustomSlim0.3.txt") || die "Cannot open TRUE\n";
my $np = 0;
while(<TRUE>)
{
	chomp $_;
	@line = split(" ", $_);
	$prot = $line[0];
	$np++;
	$true{"trueGO"} = $line[1];
}
close(TRUE);
my @trueGO = split(",", $true{"trueGO"});


open(DIST, "<data/GO_Distance_max2014.txt") || die "Cannot open DIST\n";
my %godist = ();
my %gocat = ();
while(<DIST>)
{
	chomp $_;
	@line = split("\t", $_);	#GO:0000001	mitochondrion inheritance	biological_process	5
	$go = $line[0];
	$dist = $line[3];
	if($go =~ /GO:/)
	{
		$godist{$go} = $dist;
		$gocat{$go} = $line[2];
	}
}
close(DIST);

open(GO, "<global_GO.txt") || die "cannot open GO\n";
my %func = ();
while(<GO>)
{
	chomp $_;
	@line = split(" ", $_);		#0 GO:0006811
	$func{$line[0]} = $line[1];
}
close(GO);

open(TL, "<cluster_training_labels_withGeneFlag_$iter.txt") || die "cannot open TL\n";
while(<TL>)
{
	chomp $_;
	@line = split(" ", $_);		#cluster_1	UNKNOWN	2
	if($line[2] > 0)	#cluster containes disease genes
	{
		push(@unknownprot, $line[0]);
	}
}
close(TL);

sub relaxedTP
{
	$tpFlag = 0;
	my ($GOlist, $goterm) = @_;
#	print "relaxedTP Function pred $goterm\n";
	foreach $go_(@$GOlist)
	{

#		print "relaxedTP Function golist $go_ go $goterm\n";
		if($gocat{$goterm} ne $gocat{$go_})
		{
			next;
		}

		my @trueAnc = split(",", $parents{$go_});
		my @predAnc = split(",", $parents{$goterm});
#		my @isect = intersect(@trueAnc, @predAnc);
		my %trueAnc_hash = ();	
		%trueAnc_hash = map{$_ =>1} @trueAnc;
		my @common = ();				
		foreach $p_anc(@predAnc)	#compute common parents		
		{
			if(exists($trueAnc_hash{$p_anc}))
			{
				push(@common, $p_anc);
			}
		}
		if(scalar(@common) == 0)	#no common parent
		{
			next;
		}
		foreach $anc(@common)
		{
#			print "common anc $anc\n";
			if((exists($godist{$anc})) && (exists($godist{$go_})))
			{
				$gap = abs($godist{$anc} - $godist{$go_});
				if($gap <= 2)
				{
#					print "relaxed golist $go_ PredAnc $anc gap $gap\n";
					$tpFlag = 1;
					last;
				}
			} 
		}
		if($tpFlag == 1)
		{
			last;
		}
	}
	return $tpFlag;
}

sub accuracy
{
	$class = $_[0];
	$iter = $_[1];
	$mode = $_[2];

	my %pred = ();
	if($mode == 1)
	{
		$predfile = "clusterlabel_crf_$iter.txt";	
	}
	if($mode == 2){		$predfile = "clusterlabel_$iter.txt";	}
	
	open(PRED, "<$predfile") || die "Cannot open PRED\n";
	while(<PRED>)
	{
		chomp $_;
		@line = split(" ", $_);
		$prot = $line[0];
		$pred{$prot} = $_;

		@predLab = split(" ", $pred{$prot});
		$nGO = scalar(@predLab) - 1;
	}
	close(PRED);

	foreach $prot(sort keys %pred)
	{
		$tp = 0;
		$tpr = 0;
		$fp = 0;
		$fn = 0;
		$prec = 0;
		$recall = 0;
		@predLab = split(" ", $pred{$prot});
		$nGO = scalar(@predLab) - 1;
		@predGO = ();
		for($i = 1; $i < $nGO; $i++)
		{
			$pl = $predLab[$i];
			$go = $func{$i - 1};
			if($pl == 1){	push(@predGO, $go);	}
		}
#		print "T = " . scalar(@trueGO) . "\n";
#		foreach $true(@trueGO)
#		{
#			print "true $true\n"
#		}
#		foreach $pred(@predGO)
#		{
#			print "pred $pred\n"
#		}

		#check false positives with relaxed TP constraint
		#(PFP_proteins - 2008) A term is considered correct here if it shares a common ancestor at a GO depth >= 1 (GO category root depth = 0) and is within two edges of a known annotation.
		
		#compute tp and fp
		foreach $pred(@predGO)
		{
#			print "Compute for $prot pred $pred\n";
			if(grep(/$pred/, @trueGO))	#exact
			{
#				print "$prot $pred exact TP\n";
				$tp++;
			}
			
			else				#relaxed
			{
				$tpFlag = relaxedTP(\@trueGO, $pred);
				if($tpFlag == 1)
				{
#					print "pred $pred relaxed TP\n";
					$tp++;
				}
				else
				{
					$fp++;
#					print "pred $pred relaxed FP\n";
				}
			}
		}

		#compute tpr and fn
		foreach $true(@trueGO)
		{
#			print "Compute for $prot pred $pred\n";
			if(grep(/$true/, @predGO))	#exact
			{
#				print "$prot $true exact TP\n";
				$tpr++;
			}
			
			else				#relaxed
			{
				$tpFlag = relaxedTP(\@predGO, $true);
#				$tpFlag = 0;
				if($tpFlag == 1)
				{
#					print "true $true relaxed TP\n";
					$tpr++;
				}
				else
				{
					$fn++;
#					print "true $true relaxed FP\n";
				}
			}
		}

		#compute prec recall
		if(($tp + $fp) > 0)
		{
			$prec = $tp / ($tp + $fp);
		}
	
		if(($tp + $fn) > 0)
		{
			$recall = $tpr / ($tpr + $fn);
		}
		
		if(($prec + $recall) > 0)
		{
			$fscore{$mode}{$prot} = (2*($prec*$recall)/($prec+$recall));
		}

		$prec{$prot} = $prec;
		$recall{$prot} = $recall;
		$tp{$prot} = $tp;
		$T = scalar(@trueGO);
		$P = scalar(@predGO);
		print "$prot prec $prec recall $recall fscore " . $fscore{$mode}{$prot} . " tp $tp fp $fp fn $fn T $T P $P\n";	
#		last;
	} 
	$avgPrec = 0;
	$avgRecall = 0;
	$avgFscore = 0;
	$avgTp = 0;

	foreach $prot(keys %prec)
	{
		$avgPrec += $prec{$prot};
		$avgRecall += $recall{$prot};
		$avgFscore += $fscore{$prot};
		$avgTp += $tp{$prot};
	}

	$np = keys %pred;
	$avgPrec = $avgPrec / $np;
	$avgRecall = $avgRecall / $np;
	$avgFscore = $avgFscore / $np;
	$avgTp = $avgTp / $np;

	$avgprec{$iter}{$class} = $avgPrec;
	$avgrecall{$iter}{$class} = $avgRecall;
	$avgfscore{$iter}{$class} = $avgFscore;
	$avgtp{$iter}{$class} = $avgTp;

	return $np;
};

#----------- MAIN CV -----------#

print "************ Accuracy Relaxed **************\n";
$mode = 1;	#1 = CRF 2 = MajorityVote
$np = accuracy(1, $iter, $mode);
$mode = 2;	#1 = CRF 2 = MajorityVote
$np = accuracy(1, $iter, $mode);

open(out, ">CRF_fscore_relaxed_$iter.txt");
print out "Cluster CRF-Pred MajorityVote\n";
foreach $prot(sort @unknownprot)
{
	$prot_ = $prot;
	$prot_ =~ s/_/-/g;
	print out "$prot_ ". $fscore{1}{$prot} . " " . $fscore{2}{$prot} . "\n";
}
close(out);

