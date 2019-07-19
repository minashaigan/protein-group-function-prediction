#!/usr/bin/perl

sub compute_agreement
{
	$iter1 = $_[0];
	$iter2 = $_[1];

	my %enriched_annot = ();
	open(in, "<enrichment_disease_genes/DiseaseGene_CRF_pred_Enriched_$iter1.txt") || die "Cant open annot DiseaseGene_CRF_pred_Enriched_$iter1.txt\n";
	while(<in>)
	{
		chomp $_;	#iter0 GO:0019901,GO:0009897,GO:0050794,GO:0050896,
		@line = split(" ", $_);
		$ac = $line[0];
		@goline = split(",", $line[1]);
	
		my %seen = ();
		foreach $go(@goline)
		{
			unless($seen{$go})
			{
				$seen{$go} = 1;
				push(@{$enriched_annot{$iter1}}, $go);
			}
		}
	}
	close(in);
	open(in, "<enrichment_disease_genes/DiseaseGene_CRF_pred_Enriched_$iter2.txt") || die "Cant open annot DiseaseGene_CRF_pred_Enriched_$iter1.txt\n";
	while(<in>)
	{
		chomp $_;	#iter0 GO:0019901,GO:0009897,GO:0050794,GO:0050896,
		@line = split(" ", $_);
		$ac = $line[0];
		@goline = split(",", $line[1]);
	
		my %seen = ();
		foreach $go(@goline)
		{
			unless($seen{$go})
			{
				$seen{$go} = 1;
				push(@{$enriched_annot{$iter2}}, $go);
			}
		}
	}
	close(in);

	$numerator = 0;
	$denominator = scalar(@{$enriched_annot{$iter1}}) + scalar(@{$enriched_annot{$iter2}});
	#find agreement = (A intersection B) / (A union B) between annot{$iter1} and $annot{$iter2}
	my %enriched_annot_hash = ();
	foreach $go(@{$enriched_annot{$iter1}})
	{
		$enriched_annot_hash{$go} = 1;
	}
	foreach $go(@{$enriched_annot{$iter2}})
	{
		if(exists($enriched_annot_hash{$go}))
		{
			$numerator++;
		}
	}
	$score = 0;
	if($denominator > 0)
	{
		$score = $numerator / $denominator;
	}
	return $score;
}


#---------------------------------Pipeline Call Starts Here-----------------------------------------#
 


my $count_args = $#ARGV + 1;
my $networkName = "";
my $geneGroupFileName = "";
my $pipeline_iter = 0;

if ($count_args == 2) {
	$geneGroupFileName = $ARGV[0];
	$networkName = $ARGV[1];
} else{
	print "usage: perl $0 <geneGroupFileName> <networkName>\n";
	print "example: perl $0 ../Data/Nets/RA/MAPK/RA_Panoga_MAPK.txt ../Data/Nets/RA/MAPK/RA_SNF_W.txt\n";
	exit;
}

print "********************Post Processing Input Network $networkName********************\n";
system("perl postprocessSNF.pl $geneGroupFileName $networkName\n");

$pipeline_iter = 0;

while(1)
{
	print "*********************perl clustering.pl $pipeline_iter $geneGroupFileName $networkName********************\n";
	system("perl clustering.pl $pipeline_iter $geneGroupFileName $networkName\n");

	print "*********************perl crf-funsim.pl $pipeline_iter********************\n";
	system("perl crf-funsim.pl $pipeline_iter\n");

	print "*********************perl postCRF2.pl $pipeline_iter********************\n";
	system("perl postCRF2.pl $pipeline_iter\n");

	print "**********************perl enrichment_disease_genes/get_GO.pl $pipeline_iter***************************\n";
	system("perl enrichment_disease_genes/get_GO.pl $pipeline_iter $geneGroupFileName\n");


	if($pipeline_iter >= 5)
	{
		$aScore1 = compute_agreement($pipeline_iter, $pipeline_iter-1);
		$aScore2 = compute_agreement($pipeline_iter-1, $pipeline_iter-2);
		$aScore3 = compute_agreement($pipeline_iter-2, $pipeline_iter-3);
		$diff1 = $aScore1 - $aScore2;
		$diff2 = $aScore1 - $aScore3;  
		
		print "iter = $pipeline_iter Scores $aScore1 $aScore2 $aScore3 diff1 $diff1 diff2 $diff2\n";
		if($diff1 <= 0.05 && $diff2 <= 0.08)
		{
			last;
		}
	}

	if($pipeline_iter >= 10)
	{
		last;
	}
	$pipeline_iter++;
}

$cmd = "cat ";
for($i = 0; $i <= $pipeline_iter; $i++)
{
	$cmd = $cmd . "enrichment_disease_genes/DiseaseGene_CRF_pred_Enriched_$i.txt ";
}
$cmd = $cmd . "> enrichment_disease_genes/DiseaseGene_CRF_pred_Enriched_comb.txt";
system("$cmd\n");

print("\n\nDone!\nOutput Group Function Prediction for the input gene group at end of each iteration is creaded in enrichment_disease_genes/DiseaseGene_CRF_pred_Enriched_comb.txt\n")
