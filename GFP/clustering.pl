#!/usr/bin/perl

#takes input ---> 
#	nodelabel.txt			#p1 1 0 0 1 0 0 1 0 1
#	global_GO.txt			#0 GO:0006811
#	SNFfile / toy_SNF_graph.txt
#	node_training_labels.txt	#p1 1
#generates output --->
#	cluster_graph.txt		#cluster_1 cluster_2 0.567
#	node_graph.txt			#node_1 node_2 0.567
#	current_cluster.txt		#cluster_2 Q15084, P13667, Q6UXH1, Q13795,
#	clusterlabel.txt		#p1 1 0 0 1 0 0 1 0 1
#	cluster_training_labels.txt	#cluster_1 0

my %label = ();
my %cls_label = ();
my %func = ();
my @nodes = ();
my %funsim = ();
my %cluster = ();
my @clusternodes = ();
my %network = ();
my %training = ();

$iter = $ARGV[0];
$geneGroupFileName = $ARGV[1]; 
$networkName = $ARGV[2];

$funsim_code_filepath = "slim-funsim";
sub load_data
{	
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
		#push(@nodes, $p);
	}
	close(LAB);

	open(GO, "<global_GO.txt") || die "cannot open GO\n";

	while(<GO>)
	{
		chomp $_;
		@line = split(" ", $_);		#0 GO:0006811
		$func{$line[0]} = $line[1];
	}
	close(GO);
	print "Clustering.pl : Done loading node labels iter = $iter\n";

	#load the network from SNF
	$SNFfile = $networkName;
	open(SNF, "<$SNFfile") || die "Cannot open SNF at load_data $SNFfile\n";
	$l = 0;
	
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
	print "Clustering.pl : Done Loading SNF Network\n";

	open(TL, "<node_training_labels.txt") || die "cannot open TL\n";
	
	while(<TL>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1
		$training{$line[0]} = $line[1];
	}
	close(TL);
	print "Clustering.pl : Done Loading node training labels\n";
}


sub get_funsim_matrix
{	
	my $nGO = keys %func;
	$one = 1.0;
	$score1 = sprintf("%.3f", $one);
	
	open(annot, ">$funsim_code_filepath/annotation.txt") || die "cannot open annot\n";
	for($i = 0; $i < scalar(@nodes); $i++)
	{
		$prot = $nodes[$i];
		$funsim{$prot}{$prot} = $score1;
				
		print annot "$prot ";
		for($j = 0; $j < $nGO; $j++)
		{
			if($label{$prot}{$j} == 1)
			{
				$go = $func{$j};
				print annot "$go,";
			}
		}
		print annot "\n";
	}
	close(annot);

	print "Clustering.pl : Calling Funsim Computation\nperl $funsim_code_filepath/run.pl\n";
	system("perl $funsim_code_filepath/run.pl $funsim_code_filepath");
	print "Clustering.pl : DONE computing Funsim\n";
	
	open(funsim, "<Human_funsim_scores.txt") || die "cannot open funsim\n";
	while(<funsim>)
	{
		chomp $_;
		@line = split(" ", $_);		#P55145	Q15084	0.041329807605584
		$score = sprintf("%.3f", $line[2]);
		$funsim{$line[0]}{$line[1]} = $score;
#		print "Clustering.pl : Reading funsim $line[0] $line[1] $funsim{$line[0]}{$line[1]}\n";
	}
	close(funsim);		
};

sub print_funsim_matrix
{
	open(out, ">funsim_matrix_$iter");
	
	print out "\"" . $nodes[0] . "\"";
	for($i = 1; $i < scalar(@nodes); $i++)
	{
		$prot1 = $nodes[$i];
		print out " \"$prot1\"";
	}
	print out "\n";

	for($i = 0; $i < scalar(@nodes); $i++)
	{
		$prot1 = $nodes[$i];
		print out "\"$prot1\"";
		for($j = 0; $j < scalar(@nodes); $j++)
		{
			$prot2 = $nodes[$j];
			
#			print "Clustering.pl : Printing funsim $prot1 $prot2\n";
			if(exists($funsim{$prot1}{$prot2}))	
			{
				print out " " . $funsim{$prot1}{$prot2};
				next;
			}
			if(exists($funsim{$prot2}{$prot1}))	
			{
				print out " " . $funsim{$prot2}{$prot1};
				next;
			}
		}
		print out "\n";
	}

	close(out);
	print "Clustering.pl : DONE printing Funsim matrix iter = $iter\n";
};

sub parsecluster
{
	open(AP, "<apcluster_result_$iter");
	$nline = 0;
	$n = 1;
	$prefix = "cluster";
	while(<AP>)
	{
		chomp $_;
		if($nline == 0)			#save number of clusters
		{
			$ncluster = $_;
			$nline++;
			next;
		}
		if(/"Var1" "Freq"/)		#new cluster starts
		{
			$clustername = $prefix . "_" . "$n";
			push(@clusternodes, $clustername);
			$n++;	
			next;
		}

		#else parse this to get cluster members 		"1" "P55145" 1
		s/"//g;
		@line = split(" ", $_);
		$member = $line[1];
#		print "Clustering.pl : Found $member in $clustername\n";
		push(@{$cluster{$clustername}}, $member);
		
	}
	
	print "Clustering.pl : Done parsing cluster results from R iter = $iter; Tolal $ncluster clusters\n";
	close(AP);

	open(out, ">current_cluster_$iter.txt");
	foreach $cls(sort keys %cluster)
	{
		print out "$cls ";
		foreach $prot(@{$cluster{$cls}})
		{
			print out "$prot,";
		}
		print out "\n";
	}
	close(out);
};


sub construct_clustergraph
{
	#output cluster_graph.txt	p1 p2 0.2
	
	open(out, ">cluster_graph_$iter.txt") || die "Cant open cluster_graph.txt\n";
	$nClus = keys %cluster;
	for($i = 0; $i < $nClus; $i++)	
	{
		for($j = $i+1; $j < $nClus; $j++)	
		{
			$cluster1 = $clusternodes[$i];
			$cluster2 = $clusternodes[$j];
			
			$edge_weight = 0;
			$edge_count = 0;
			foreach $prot1(@{$cluster{$cluster1}})
			{
				foreach $prot2(@{$cluster{$cluster2}})
				{
					if(exists($network{$prot1}{$prot2}))
					{
						if($network{$prot1}{$prot2} > 0)
						{
#							print "Clustering.pl : within cluster edge $prot1 $prot2 ". $network{$prot1}{$prot2} . "\n";
							$edge_weight += $network{$prot1}{$prot2};
							$edge_count++;
						}
					}
					elsif(exists($network{$prot2}{$prot1}))
					{
						if($network{$prot1}{$prot2} > 0)
						{
#							print "Clustering.pl : within cluster edge $prot1 $prot2 ". $network{$prot1}{$prot2} . "\n";
							$edge_weight += $network{$prot2}{$prot1};
							$edge_count++;
						}
					}
				}
			}
			if($edge_count > 0)
			{
				$edge_weight = $edge_weight / $edge_count;
			}
			$score = sprintf("%.3f", $edge_weight);
			print out "$cluster1 $cluster2 $score\n";
		}
	}
	close(out);
	
};

sub construct_nodegraph
{
	open(out, ">node_graph_$iter.txt") || die "Cant open node_graph.txt\n";

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
}

sub construct_clusterlabels
{
	#output clusterlabel.txt	P55145 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0

	$nGO = keys %func;
	
	for($i = 0; $i < $nGO; $i++)
	{
		foreach $cls(sort keys %cluster)
		{
			$onevote = $zerovote = 0;
			foreach $prot(@{$cluster{$cls}})
			{
				#if($i == 0){	print "$cls $prot "."label ". $label{$prot}{$i} . "\n" ;	}
				if($training{$prot} == 1)		#put cluster labels based on known member proteins
				{
					if($label{$prot}{$i} == 1)
					{
						$onevote++;
					}
					if($label{$prot}{$i} == 0)
					{
						$zerovote++;
					}	
				}
			}
			if($onevote >= $zerovote){ $cls_label{$cls}{$i} = 1;	}
			else{	$cls_label{$cls}{$i} = 0;	}
			
		}
	}
	open(out, ">clusterlabel_$iter.txt") || die "Cannot open clusterlabel.txt\n";
	foreach $cls(sort keys %cluster)
	{
		print out $cls;
		for($i = 0; $i < $nGO; $i++)
		{
			print out " ". $cls_label{$cls}{$i}; 
		}
		print out "\n";
	}
	close(out);
	print "Clustering.pl : Done printing cluster labels iter = $iter\n";
};

sub mark_cluster_training
{
	#output cluster_training_labels.txt	P55145 1

	open(out, ">cluster_training_labels_$iter.txt") || die "Cannot open cluster_training_labels\n";
	open(out1, ">clustering_log_$iter.txt") || die "Cannot open cluster_training_labels\n";

	$frac_cutoff = 0.1;
	my $unknown_cls = 0;
	foreach $cls(sort keys %cluster)
	{
		$nMember = 0;
		$known = 0;
		foreach $prot(@{$cluster{$cls}})
		{
			$nMember++;
			if($training{$prot} == 1)
			{
				$known++;
			}			
		}
		
		$frac = $known / $nMember; 
		$unknown_frac = 1 - $frac;
		print out1 "Clustering.pl : $cls known $known total $nMember fraction unknown $unknown_frac\n";
		if($unknown_frac >= $frac_cutoff)
		{
			print out "$cls 0\n";
			$unknown_cls++;
		}
		else
		{
			print out "$cls 1\n";
		}
	}
	print out1 "Clustering.pl : Done printing cluster Training labels; iter = $iter; Total $ncluster Cluster #unknown cluster = $unknown_cls member cutoff $frac_cutoff\n";
	close(out);
	close(out1);
};

sub mark_cluster_genes
{
	#output cluster_training_labels_withGeneFlag.txt	clustername traininglabel geneflag		(geneflag >= 1 is the cluster has any member from the input gene group)


	my %genes = ();
	open(FILE, "<$geneGroupFileName") || die "cannot open FILE\n";
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

	open(out, ">cluster_training_labels_withGeneFlag_$iter.txt") || die "Cannot open cluster_training_labels\n";

	$frac_cutoff = 0.1;
	foreach $cls(sort keys %cluster)
	{
		$nMember = 0;
		$known = 0;
		$geneflag = 0;
		foreach $prot(@{$cluster{$cls}})
		{
			$nMember++;
			if($training{$prot} == 1)
			{
				$known++;
			}	
			if(exists($genes{$prot}))
			{
				$geneflag++;
			}
		}
		$frac = $known / $nMember; 
		$unknown_frac = 1 - $frac;
		if($unknown_frac >= $frac_cutoff)
		{
			print out "$cls\tUNKNOWN\t$geneflag\n";
		}
		else
		{
			print out "$cls\tKNOWN\t$geneflag\n";
		}
	}
	
	close(out);
	print "Clustering.pl : Done printing cluster Training labels with geneflag iter = $iter\n";

};

sub clear
{
	%label = ();
	%cls_label = ();
	%func = ();
	@nodes = ();
	%funsim = ();
	%cluster = ();
	@clusternodes = ();
	%network = ();
	%training = ();
};


load_data();
print "Clustering.pl : nGO : $nGO iter = $iter #nodes ". scalar(@nodes) . "\n";
get_funsim_matrix();
print_funsim_matrix();

#call clustering
system("Rscript cluster.R $iter $networkName\n");
print "Clustering.pl : Done with clustering iteration $iter\n";

#make cluster graph
parsecluster();
construct_clustergraph();
construct_nodegraph();
construct_clusterlabels();
mark_cluster_training();
mark_cluster_genes();	#print cluster Training labels with geneflag
clear();





