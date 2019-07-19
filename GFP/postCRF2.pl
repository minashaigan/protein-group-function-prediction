#!/usr/bin/perl

#takes input ---> 
#	nodelabel.txt			#p1 1 0 0 1 0 0 1 0 1
#	clusterlabel_crf.txt		#p1 1 0 0 1 0 0 1 0 1
#	current_cluster.txt		#cluster_2 Q15084, P13667, Q6UXH1, Q13795,
#	node_training_labels.txt	#p1 1
#	global_GO.txt			#0 GO:0006811

#generates output --->
#	nodelabel.txt			#p1 1 0 0 1 0 0 1 0 1

my %labels = ();
my %ori_label = ();
my $cls_label = ();
my @nodes = ();
my %cluster = ();
my %training = ();
my %func = ();
$funsim_code_filepath = "slim-funsim";
my %allsim = ();
my $ss_cutoff = 0.3;
my $tot_training_prot = 0;
$iter = $ARGV[0];

sub load_data
{
	#load node labels
	open(LAB, "<nodelabel_$iter.txt") || die "Cannot open LAB\n";
	while(<LAB>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1 0 0 1 0 0 1 0 1
		$p = $line[0];
		for($i = 1; $i < scalar(@line); $i++)
		{
			$l = $line[$i];
			$label{$p}{$i-1} =  $l;
		}
		push(@nodes, $p);
	}
	close(LAB);
	
	print "postCRF.pl : Done loading node lables\n";

	#load original node labels
	open(LAB, "<nodelabel_0.txt") || die "Cannot open LAB0\n";
	while(<LAB>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1 0 0 1 0 0 1 0 1
		$p = $line[0];
		for($i = 1; $i < scalar(@line); $i++)
		{
			$l = $line[$i];
			$ori_label{$p}{$i-1} =  $l;
		}
	}
	close(LAB);
	
	print "postCRF.pl : Done loading original node lables\n";

	#load cluster labels
	open(LAB, "<clusterlabel_crf_$iter.txt") || die "Cannot open cluster LAB\n";
	while(<LAB>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1 0 0 1 0 0 1 0 1
		$p = $line[0];
		for($i = 1; $i < scalar(@line); $i++)
		{
			$l = $line[$i];
			$cls_label{$p}{$i-1} =  $l;
		}
	}
	close(LAB);
	
	print "postCRF.pl : Done loading cluster labels iter = $iter\n";

	#load cluster
	open(CLS, "<current_cluster_$iter.txt") || die "Cannot open current cluster\n";
	while(<CLS>)
	{
		chomp $_;
		@line = split(" ", $_);		#cluster_2 Q15084, P13667, Q6UXH1, Q13795,
		@members = split(",", $line[1]);
		$cls = $line[0];
		foreach $member(@members)
		{
			push(@{$cluster{$cls}}, $member);
		}
	}
	close(CLS);
	print "postCRF.pl : Done loading cluster memberships\n";

	open(TL, "<node_training_labels.txt") || die "cannot open TL\n";
	
	while(<TL>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1
		$training{$line[0]} = $line[1];
		if($line[1] == 0)
		{
			$training{$line[0]} = 2;
		}
		else
		{
			$tot_training_prot++;
		}
	}
	close(TL);
	print "postCRF.pl : Done Loading node training labels\n";

	open(GO, "global_GO.txt") || die "cannot open GO\n";

	while(<GO>)
	{
		chomp $_;
		@line = split(" ", $_);		#0 GO:0006811
		$func{$line[0]} = $line[1];
	}
	close(GO);
	print "postCRF.pl : Done loading GO indexes\n";

	open (BPSIM, "<$funsim_code_filepath/BP_cursim_table_Rel.txt") || die "Cannot open BPSIM\n";
	while (<BPSIM>)
	{
		chomp $_;
		@line = split(" ", $_);
		$allsim{$line[0]}{$line[1]} = $line[2];
		$allsim{$line[1]}{$line[0]} = $line[2];
	}
	close(BPSIM);

	open (MFSIM, "<$funsim_code_filepath/MF_cursim_table_Rel.txt") || die "Cannot open MFSIM\n";
	while (<MFSIM>)
	{
		chomp $_;
		@line = split(" ", $_);
		$allsim{$line[0]}{$line[1]} = $line[2];
		$allsim{$line[1]}{$line[0]} = $line[2];
	}
	close(MFSIM);

	open (CCSIM, "<$funsim_code_filepath/CC_cursim_table_Rel.txt") || die "Cannot open CCSIM\n";
	while (<CCSIM>)
	{
		chomp $_;
		@line = split(" ", $_);
		$allsim{$line[0]}{$line[1]} = $line[2];
		$allsim{$line[1]}{$line[0]} = $line[2];
	}
	close(CCSIM);

	open(FAM, "<data/FAM_customslimGO_IC0.3.txt") || die "cannot open FAM\n";
	while(<FAM>)
	{
		chomp $_;
		@line = split(" ", $_);		#GO:0000228 GO:0000003 0.0942163552970131
		$score = sprintf("%.3f", $line[2]);
		if(exists($allsim{$line[0]}{$line[1]}))
		{
			next;
		}
		$allsim{$line[0]}{$line[1]} = $score;
		$allsim{$line[1]}{$line[0]} = $score;
		
	}
	close(FAM);
};

#for a given GO term and an array of GO terms, find max SS score between the given term to any term in the array (except with itself)
sub get_max_SS
{
	my ($GOlist, $goterm) = @_;
	$max = -1;
	foreach $go_(@$GOlist)
	{
		if($goterm eq $go_)	
		{
			next;
		}
		if(exists($allsim{$goterm}{$go_}))
		{
			$ss = $allsim{$goterm}{$go_};
#			print "SS of pair $goterm $go_ = $ss\n";
			if($ss >= $max)
			{
				$max = $ss;
			}
		}
		else
		{
#			print "pair $goterm $go_ dont exist\n";
		}
	}
	return $max;
};
sub reassign_graph_nodes
{
	#output nodelabel.txt	p1 1 0 0 1 0 0 1 0 1

	$n_changed = 0;
	$nGO = keys %func;
	foreach $cls(sort keys %cluster)
	{
		foreach $prot(@{$cluster{$cls}})
		{
			if($training{$prot} == 2)		#unknown protein
			{
				for($i = 0; $i < $nGO; $i++)		#reassign the label of the prot using old_label = new_group_label
				{
					$ori_label{$prot}{$i} = $cls_label{$cls}{$i};
				}
			}
	
			else					#known protein
			{
#				print "postCRF.pl : Reassigning lables of $prot in $cls\n";
				
				#create array of existing GO terms
				my @oldgo = ();
				for($i = 0; $i < $nGO; $i++)
				{
					if($ori_label{$prot}{$i} == 1)
					{
						$go = $func{$i};
						push(@oldgo, $go);
					}
				}
				for($i = 0; $i < $nGO; $i++)
				{
					if($cls_label{$cls}{$i} == 1 && $ori_label{$prot}{$i} == 0)		#check whether should add new GO term to the original training label
					{
						$newgo = $func{$i};
						$maxSS = get_max_SS(\@oldgo, $newgo);
						if($maxSS >= $ss_cutoff)
						{
							$n_changed++;
							$ori_label{$prot}{$i} = $cls_label{$cls}{$i}; 	#update group function
						}
					}
				}
			}
			
		}
	}
	
	$niter = $iter + 1;
	open(out, ">nodelabel_$niter.txt") || die "cannot open nodeLabels\n";
	foreach $prot(@nodes)
	{
		print out "$prot";
		for($i = 0; $i < $nGO; $i++)		
		{
			print out " " . $ori_label{$prot}{$i};
		}
		print out "\n";
	}
	close(out);

	print "postCRF.pl : Done reassigning labels of Graph nodes based on Cluster Lables\n";
	
	$tot_prot = 0;
	$n_changed = $n_changed / ($nGO * $tot_training_prot);

	print "postCRF.pl : label added for fraction of known proteins labels = $n_changed\n";
};

sub clear
{
	%cluster = ();
	%training = ();
	%func = ();
	%allsim = ();
};

load_data();
reassign_graph_nodes();		#change node labels of unknown proteins
clear();
%allsim = ();
