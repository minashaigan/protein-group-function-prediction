#!/usr/bin/perl

open(GO, "<global_GO.txt") || die "cannot open GO\n";
my %func = ();
while(<GO>)
{
	chomp $_;
	@line = split(" ", $_);		#0 GO:0006811
	$func{$line[0]} = $line[1];
}
close(GO);

my %genes = ();
$geneGroupFileName = $ARGV[1];
open(FILE, "<$geneGroupFileName") || die "cannot open FILE $geneGroupFileName\n";
while(<FILE>)
{
	chomp $_;
	if(/#/)
	{
		next;
	}
	@line = split(" ", $_);		#CBL	P22681	1
	$genes{$line[1]} = 1;
}


$iter = $ARGV[0];

open(PRED, "<nodelabel_$iter.txt") || die "cannot open PRED\n";

my %geneannot = ();
while(<PRED>)
{
	chomp $_;
	@predLab = split(" ", $_);
	$prot = $predLab[0];
	if(exists($genes{$prot}))
	{
		$nGO = scalar(@predLab) - 1;
	
		for($i = 1; $i < $nGO; $i++)
		{
			$pl = $predLab[$i];
			$go = $func{$i - 1};
			if($pl == 1)
			{
				push(@{$geneannot{$prot}}, $go);
			}
		}
	}
}
close(PRED);

open(out, ">enrichment_disease_genes/DiseaseGene_CRF_pred_$iter.txt");
foreach $prot (sort keys %geneannot)
{
	print out "$prot ";
	foreach $go(sort @{$geneannot{$prot}})
	{
		print out "$go,";
	}
	print out "\n";
}
close(out);


system("rm GO-above-pval0.01.txt\n");
system("python enrichment_disease_genes/enrich_offline.py -o 9606 -f enrichment_disease_genes/DiseaseGene_CRF_pred_$iter.txt\n");

open(ENR, ">enrichment_disease_genes/DiseaseGene_CRF_pred_Enriched_$iter.txt");

open(enrich, "<GO-above-pval0.01.txt") || die "cant open enrich\n";
my $enrichterms = <enrich>;
close(enrich);
print ENR "iter$iter $enrichterms\n";
print "iter$iter $enrichterms\n";

close(ENR);
