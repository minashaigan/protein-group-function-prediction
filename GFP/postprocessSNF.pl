#!/usr/bin/perl

#takes input ---> 
#	Nets/RA/MAPK/RA_SNF_W.txt
#	data/SlimGOterms.txt
#	data/Human_fullannot_GO_to_Slim_MAP.txt
#	data/Human_uniprot_annot_full.txt
#generates output --->
#	global_GO.txt			#0 GO:0006811
#	nodelabel.txt			#p1 1 0 0 1 0 0 1 0 1
#	node_training_labels.txt	#p1 1
#	node_graph.txt			#node_1 node_2 0.567

my $groupGenefile = $ARGV[0];
my $SNFfile = $ARGV[1];
print "postprocessSNF: $SNFfile\n";

#----------------- Load Integrated Network --------------------#
open(SNF, "<$SNFfile") || die "Cannot open SNF\n";
$l = 0;
my @nodes = (); 
my %network = ();
while(<SNF>)
{
	chomp $_;
	s/"//g;
	@line = split(" ", $_);
	if($l == 0)		#a b c d e f g
	{
		@nodes = @line;
		$l++;	
		next;
	}
	
	$prot1 = $line[0];	#a 1.0 0.9 0.8 0.0 0.0 0.0
	for($i = 1; $i < scalar(@line); $i++)
	{
		$prot2 = $nodes[$i - 1];
		#print "$prot1 $prot2 $line[$i]\n";
		$network{$prot1}{$prot2} = $line[$i];

	}
}
close(SNF);
print "Done Loading Network\n";

#----------------------Load Slimmed GO set------------------------------#

$slimfile = "data/Term_list_for_HumanSlim0.3_GOSlim";
open(SLIM, "<$slimfile") || die "Cannot open SLIM\n";
my @slimGO = ();

while(<SLIM>)
{
	chomp $_;
	push(@slimGO, $_);
}
close(SLIM);

print "Loaded Slim SET; Total Slim GO ". scalar(@slimGO) . "\n";


my %slimID = ();
$id = 0;
foreach $go(@slimGO)
{
	$slimID{$go} = $id;
	$id++;	
}

my @sortedSlimGO = sort { $slimID{$a} <=> $slimID{$b} } keys %slimID;

#------------------ Load GO Slim mapping for entire Human annotation ----------------#

$slimMAPfile = "data/CustomSlim_Human_UP_IC0.3.txt";
open(SLIM, "<$slimMAPfile") || die "Cannot open SLIM\n";

my %slimMAP = ();
my %seen = ();
while(<SLIM>)
{
	chomp $_;
	@line = split(" ", $_);		#GO:0000012 GO:0006259,GO:0051716,GO:0006950,GO:0044763,GO:0044710,
	@slimset = split(",", $line[1]);
	
	$go = $line[0];
	foreach $goS(@slimset)
	{
		if(grep(/$goS/, @slimGO))
		{
			push(@{$slimMAP{$go}}, $goS);
		}
	}	
}
close(SLIM);

#----------------------Load GO annotation for nodes--------------------#
$annotfile = "data/Human_uniprot_annot_full.txt";
open(ANNOT, "<$annotfile") || die "Cannot open full human ANNOT\n";
my %annot = ();

while(<ANNOT>)
{
	chomp $_;
	@line = split(" ", $_);			#A0M8Q6 GO:0006956,GO:0006958,GO:0038095,GO:0038096,GO:0045087,GO:0050776,GO:0003823,GO:0005576,GO:0005886,
	@goset = split(",", $line[1]);

	$prot = $line[0];
	foreach $go(@goset)
	{
		push(@{$annot{$prot}}, $go);
	}
}
close(ANNOT);
print "Done loading Annot\n";

#----------------------ReLoad GO annotation for disease related nodes--------------------#
system("perl get_annot.pl $groupGenefile\n");
$annotfile = "data/group_annot_full.txt";
open(ANNOT, "<$annotfile") || die "Cannot open group gene ANNOT\n";
while(<ANNOT>)
{
	chomp $_;
	@line = split(" ", $_);			#A0M8Q6 GO:0006956,GO:0006958,GO:0038095,GO:0038096,GO:0045087,GO:0050776,GO:0003823,GO:0005576,GO:0005886,
	@goset = split(",", $line[1]);

	$prot = $line[0];
	print "removed annot for diseasegene $prot before length = " . scalar(@{$annot{$prot}}) . "\n";
	$annot{$prot} = ();
	foreach $go(@goset)
	{
		push(@{$annot{$prot}}, $go);
	}
	print "removed annot for diseasegene $prot after length = " . scalar(@{$annot{$prot}}) . "\n";
}
close(ANNOT);
print "Done loading Annot\n";

#---------------MAP annpt to slim and create binary GO vector for nodes in the graph-------------------
my %nodeLabels = ();
$annot_go = 0;
my %missed = ();
my $nKnown = 0;
foreach $prot(@nodes)
{
	$known = 0;
	if(exists($annot{$prot}))
	{
		foreach $go(@{$annot{$prot}})
		{
			if(exists($slimMAP{$go}))
			{
#				print "SLIM-MAP $go->";
				foreach $goS(@{$slimMAP{$go}})
				{
					$nodeLabels{$prot}{$slimID{$goS}} = 1; 
#					print "$goS:$slimID{$goS},";
				
				}
#				print "\n";
				$missed{$go} = 2;
				$known = 1;
			}
			else
			{
#				print "SLIM-MAP EMPTY for $go\n";
				$missed{$go} = 1;
			}
		}
	}
	else
	{
		print "prot $prot NO annotations in Uniprot\n";
	}
	if($known == 1)
	{
		$nKnown++;
	}
}

#unique GO terms in the network
my %seen = ();
my @uniqGO = ();
foreach $prot(@nodes)
{
	if(exists($annot{$prot}))
	{
		foreach $go(@{$annot{$prot}})
		{
			unless($seen{$go})
			{
				$seen{$go} = 1;
				push(@uniqGO, $go);
			}
		}
	}
}

#how many GO terms missed slim mapping

my $missed_map = 0;
foreach $go(@uniqGO)
{
	if($missed{$go} == 1)
	{
		$missed_map++;
	}
}

open(OUT, ">nodelabel_0.txt");
open(OUT1, ">nodeGO.txt");
open(OUT3, ">originalAnnot.txt");
open(OUT2, ">node_training_labels.txt");

print "Done creating Binary GO vectors; #GO terms missing Slim Map = $missed_map out of " . scalar(@uniqGO) . "\n";

my %emptyGO = ();

#count occurances of GO terms in proteins

foreach $go(@sortedSlimGO)
{
	$id = $slimID{$go};
	$emptyGO{$id} = 1;		#GO column is empty
	foreach $prot(@nodes)
	{
		if(exists($nodeLabels{$prot}{$id}))
		{
			$emptyGO{$id}++;	#GO column is non empty
		}	
	}
	#print "$go " . $emptyGO{$id} . "\n";
}

my $nProt = $nKnown;
my $nempty = 0;
my $nempty1 = 0;
$occCutoff = 0;

open(out, ">go_perc.txt");
foreach $go(@sortedSlimGO)
{
	$goId = $slimID{$go};
	if($emptyGO{$goId} == 1)
	{
		$nempty++;
		next;
	}
	
	$occ = $emptyGO{$goId} - 1;
	$perc = $occ / $nProt;

	print out "$go occ $occ out of $nProt perc $perc\n";
	if($perc < $occCutoff)
	{
		$emptyGO{$goId} = 1;
		$nempty++;
		next;
	}

	if($perc > (1-$occCutoff))
	{
#		print "too many 1s in $go\n";
		$emptyGO{$goId} = 1;
		$nempty1++;
		next;
	}
	
}
close(out);

foreach $prot(@nodes)
{
	print OUT "$prot";
	print OUT1 "$prot";
	print OUT2 "$prot";
	print OUT3 "$prot";
	$train = 0;
	foreach $go(@sortedSlimGO)
	{
		$id = $slimID{$go};
		
		if($emptyGO{$id} >= 1)		#non empty GO column
		{
			if(exists($nodeLabels{$prot}{$id}))
			{
				print OUT " 1";
				print OUT1 " $go";
				$train = 1;
			}
			else
			{
				print OUT " 0";
			}
		}
		
	}

	foreach $go(@{$annot{$prot}})
	{
		print OUT3 " $go";
	}

	print OUT "\n";
	print OUT1 "\n";
	print OUT2 " $train\n";
	print OUT3 "\n";
}

close(OUT);
close(OUT1);
close(OUT2);
close(OUT3);

$perc = ($nKnown / scalar(@nodes)) * 100;
$nUnknown = scalar(@nodes) - $nKnown;
print "Percentage of Known protiens in Network $perc % out of ". scalar(@nodes) . " known " . "$nKnown unknown $nUnknown\n". "#GO with under occ = $nempty #GO with over occ = $nempty1 occ cutoff $occCutoff\n";

open(out, ">node_graph_ori.txt") || die "Cant open node_graph.txt\n";

for($i = 0; $i < scalar(@nodes); $i++)
{
	for($j = $i + 1; $j < scalar(@nodes); $j++)
	{
		$node1 = $nodes[$i];
		$node2 = $nodes[$j];
		$score = $network{$node1}{$node2};
		print out "$node1 $node2 $score\n";
	}
}
close(out);


open(OUT, ">global_GO.txt") || die "Cant open OUT\n";
$index = 0;
foreach $go(@sortedSlimGO)
{
	$id = $slimID{$go};
	if($emptyGO{$id} >= 1)		#non empty GO column
	{
		print OUT $index . " ". $go . "\n";
		$index++;
	}
}
close(OUT);

print "Total Slim GO = ". scalar(@slimGO) . " Slim GO after filter wtih occ cutoff $occCutoff = $index\n";
