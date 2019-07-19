#!/usr/bin/perl

#takes input ---> 
#	cluster_graph.txt		#cluster_1 cluster_2 0.567
#	clusterlabel.txt		#p1 1 0 0 1 0 0 1 0 1
#	global_GO.txt			#0 GO:0006811
#	cluster_training_labels.txt	#cluster_1 0

#generates output --->
#	clusterlabel.txt		#p1 1 0 0 1 0 0 1 0 1


use Graph::Undirected;
use Math::Random qw(random_normal);
use Math::Random qw(random_uniform);
use Statistics::R;

my $net = Graph::Undirected->new; # An undirected graph.
my %label = ();
my %func = ();
my %gomap = ();
my %training = ();
my %features = ();
my %fam = ();
my %mean_prob = ();
my $t = 1;		#iterations
$funsim_code_filepath = "slim-funsim";
my %edges = ();
my %cat = ();
my %bpgo = ();
my %mfgo = ();
my %ccgo = ();
my %bpsim = ();
my %mfsim = ();
my %ccsim = ();

$pIter = $ARGV[0];	#pipeline iteration number

sub get_graph()
{
	

	open(GRAPH, "<cluster_graph_$pIter.txt") || die "Cannot open GRAPH\n";
#	open(GRAPH, "<node_graph_known.txt") || die "Cannot open GRAPH\n";

	while(<GRAPH>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 p2 0.2
		if($line[2] > 0)
		{
			$net -> add_weighted_edge($line[0], $line[1], $line[2]);
			$edges{$line[0]}{$line[1]} = $line[2];
			$edges{$line[1]}{$line[0]} = $line[2];
		}
	}
	close(GRAPH);

	print "-----------Network Loaded-------------\n";
	print $net."\n";

	my @node = $net->vertices;
	for($i = 0; $i < scalar(@node) - 1; $i++)
	{
		for($j = $i + 1; $j < scalar(@node); $j++)
		{
			if($net -> has_edge($node[$i], $node[$j]))
			{
				print $node[$i] . "-" . $node[$j] . " " . $net -> get_edge_weight($node[$i], $node[$j]) . "\n";
			}
		}
	}

	open(LAB, "<clusterlabel_$pIter.txt") || die "Cannot open LAB\n";
#	open(LAB, "<nodelabel_rand.txt") || die "Cannot open LAB\n";

	while(<LAB>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1 0 0 1 0 0 1 0 1
		$p = $line[0];
		for($i = 1; $i <= scalar(@line); $i++)
		{
			$l = $line[$i];
			$label{$p}{$i-1} =  $l;
		}
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

	open(TL, "<cluster_training_labels_$pIter.txt") || die "cannot open TL\n";
#	open(TL, "<node_training_labels_rand.txt") || die "cannot open TL\n";

	while(<TL>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1
		$training{$line[0]} = $line[1];
		if($line[1] == 0)
		{
			$training{$line[0]} = 2;
		}
	}
	close(TL);

	open(FAM, "<data/FAM_customslimGO_IC0.3.txt") || die "cannot open FAM\n";
	while(<FAM>)
	{
		chomp $_;
		@line = split(" ", $_);		#GO:0000228 GO:0000003 0.0942163552970131
		$score = sprintf("%.3f", $line[2]);
		$fam{$line[0]}{$line[1]} = $score;
		
	}
	close(FAM);

	open(CAT, "<$funsim_code_filepath/GO_201401_category.txt") || die "split_annot.pl : Cant open CAT\n";
	while(<CAT>)
	{
		chomp $_;
		@line = split(" ", $_);		#8049	p
		$len = length($line[0]);
		$prefix = "GO:";
		$n_zeros = 7 - $len;
	
		for($i = 0; $i < $n_zeros; $i++)
		{
			$prefix = $prefix . "0";
		}
		$go = $prefix.$line[0];
		$cat{$go} = $line[1];
		#print "get_rel_ancestors.pl : CAT $go $line[1]\n";
	}
	close(CAT);

	open (BPSIM, "<$funsim_code_filepath/BP_cursim_table_Rel.txt") || die "Cannot open BPSIM\n";
	while (<BPSIM>)
	{
		chomp $_;
		@line = split(" ", $_);
		$bpsim{$line[0]}{$line[1]} = $line[2];
		$bpsim{$line[1]}{$line[0]} = $line[2];
	}
	close(BPSIM);

	open (MFSIM, "<$funsim_code_filepath/MF_cursim_table_Rel.txt") || die "Cannot open MFSIM\n";
	while (<MFSIM>)
	{
		chomp $_;
		@line = split(" ", $_);
		$mfsim{$line[0]}{$line[1]} = $line[2];
		$mfsim{$line[1]}{$line[0]} = $line[2];
	}
	close(MFSIM);

	open (CCSIM, "<$funsim_code_filepath/CC_cursim_table_Rel.txt") || die "Cannot open CCSIM\n";
	while (<CCSIM>)
	{
		chomp $_;
		@line = split(" ", $_);
		$ccsim{$line[0]}{$line[1]} = $line[2];
		$ccsim{$line[1]}{$line[0]} = $line[2];
	}
	close(CCSIM);
	
};

sub get_fam					#fam($cur_go|$y_go)
{
	$y_go = $_[0];			
	$cur_go = $_[1];

	if(exists($fam{$y_go}{$cur_go}))
	{
		$score = $fam{$y_go}{$cur_go};
		return $score;
	}
	if(exists($fam{$cur_go}{$y_go}))	#reminder!! FAM is not symmetric
	{
		$score = $fam{$cur_go}{$y_go};
		return $score;
	}
	else
	{
		$zero = 0;
		$score = sprintf("%.3f", $zero);
		return $score;
	}
	
};

sub compute_node_feature
{
	$go_i = $_[0];
	$prot = $_[1];
	$nGO = keys %func;
	
	$y_go = $func{$go_i};
	$FAM_pos = 0;
	$FAM_neg = 0;
	$fam_cutoff = 0.25;		#PFP_proteins paper uses 5 FAM cutoffs - 0.1, 0.25, 0.5, 0.75, 0.9
#	$fam_cutoff = 0;		#PFP_proteins paper uses 5 FAM cutoffs - 0.1, 0.25, 0.5, 0.75, 0.9

	for($j = 0; $j < $nGO; $j++)
	{
		$cur_go = $func{$j};
		$golabel = $label{$prot}{$j};
		$fam = get_fam($y_go, $cur_go);

		if($j == $go_i)
		{
			$fam = 1.0;
			#next;
		}

		if($fam >= $fam_cutoff)
		{
			if($golabel == 1)	
			{
				$FAM_pos += $fam;
				#print "crf.pl : $go_i prot $prot FAM POS $prot $y_go $cur_go-$j $fam\n";
			} 	
			elsif($golabel == 0)
			{
				$FAM_neg += $fam;
				#print "crf.pl : $go_i prot $prot FAM NEG $prot $y_go $cur_go-$j $fam\n";
			}
		}
	}
	
	#print "crf.pl : $go_i prot $prot TOTAL FAM POS  $FAM_pos TOTAL FAM NEG  $FAM_neg\n";

	$score = sprintf("%.4f", $FAM_pos);	# + random_normal(1, 0, 0.0001));
	$features{$prot}{"FAM_pos"} = $score;
	$score = sprintf("%.4f", $FAM_neg);	# + random_normal(1, 0, 0.0001));
	$features{$prot}{"FAM_neg"} = $score; 
};


sub compute_spec_sens
{
	$prot1 = $_[0];
	$prot2 = $_[1];
	$mode = $_[2];
	
	@go1 = ();
	@go2 = ();
	@rowMax = ();
	@colMax = ();
	@rMax = ();
	@cMax = ();
	$finR = 0;
	$finC = 0;
	%matrix = ();
	$normColMax = $normMax = 0;	

	if($mode == 1){		@go1 = @{$bpgo{$prot1}};	@go2 = @{$bpgo{$prot2}};		}
	if($mode == 2){		@go1 = @{$mfgo{$prot1}};	@go2 = @{$mfgo{$prot2}};		}
	if($mode == 3){		@go1 = @{$ccgo{$prot1}};	@go2 = @{$ccgo{$prot2}};		}
	#print "mode $mode\n";
	for ($q=1; $q <= scalar(@go1);$q++){
		$term1 = $go1[$q-1];
		#print "TERM 1:$term1\tTERM2: ";	#check matrix row
		for($z=1; $z <= scalar(@go2); $z++){
			$term2 = $go2[$z-1];
			#print "$term2 ";		#checks matrix column
			if($mode == 1){		$score = $bpsim{$term1}{$term2};	}
			if($mode == 2){		$score = $mfsim{$term1}{$term2};	}
			if($mode == 3){		$score = $ccsim{$term1}{$term2};	}
			
			$matrix{$q}{$z} = $score;
			#print "score $matrix{$q}{$z} \n";
		}
		#print TEST "\n\n";#debugging	
	}
	for($w=1;$w<=scalar(@go1);$w++){
		for($x=1;$x<=scalar(@go2);$x++){
			push (@rowMax, $matrix{$w}{$x});
		}
		@rowMax = sort {$a <=> $b} @rowMax;
		$Rmax = $rowMax[$#rowMax];
		#print "rmax $Rmax\n";
		push (@rMax, $Rmax);
		@rowMax = ();
	}
	foreach (@rMax) {
		if (scalar(@rMax) == 1 ){
			$finR = $rMax[0];
		} else {
			$finR += $_;
		}
	}
	#print "finR $finR\n";
	if (scalar(@go1) != 0) {		$normMax = $finR/scalar(@go1);		}
	#print "normMax $normMax\n";

	for ($b=1;$b<=scalar(@go2);$b++){
		for($f=1;$f<=scalar(@go1);$f++){
			push (@colMax, $matrix{$f}{$b});
		}
		@colMax = sort {$a <=> $b} @colMax;
		$Cmax = $colMax[$#colMax];
		#print "cMax $Cmax\n";
		push (@cMax, $Cmax);
		@colMax = ();
	}
	foreach (@cMax) {
		if (scalar(@cMax) == 1 ) {
		$finC = $cMax[0];
		} else {
			$finC += $_;
		}
	}
	#print "finC $finC\n";
	if (scalar(@go2) != 0 ) {		$normColMax = $finC/scalar(@go2);	}
	#print "normColMax $normColMax\n";

	if($normMax >= $normColMax){		return $normMax;	}
	else{		return $normColMax;	}
	
};

sub get_funsim_score
{
	$prot1 = $_[0];
	$prot2 = $_[1];
	my $nGO = keys %func;	
	%bpgo = ();
	%mfgo = ();
	%ccgo = ();
	$funsim_cutoff = 0.4;		#top 5% human funsim cutoff
#	$funsim_cutoff = 0;		#top 5% human funsim cutoff
	for($i = 0; $i < $nGO; $i++)
	{
		if($label{$prot1}{$i} == 1)
		{
			$go = $func{$i};
			#print "$prot1 $i $go\n";	
			if($cat{$go} eq 'p'){	push(@{$bpgo{$prot1}}, $go); 		}
			if($cat{$go} eq 'f'){	push(@{$mfgo{$prot1}}, $go);		}
			if($cat{$go} eq 'c'){	push(@{$ccgo{$prot1}}, $go);		}
		}
	}
	
	for($i = 0; $i < $nGO; $i++)
	{
		if($label{$prot2}{$i} == 1)
		{
			$go = $func{$i};
			#print "$prot2 $i $go\n";	
			if($cat{$go} eq 'p'){	push(@{$bpgo{$prot2}}, $go); 		}
			if($cat{$go} eq 'f'){	push(@{$mfgo{$prot2}}, $go);		}
			if($cat{$go} eq 'c'){	push(@{$ccgo{$prot2}}, $go);		}
		}
	}
	if((scalar(@bpgo{$prot1}) == 0) || (scalar(@bpgo{$prot2}) == 0)){	$bpmax = 0;	}
	else{	$bpmax = compute_spec_sens($prot1, $prot2, 1);		}
	if((scalar(@mfgo{$prot1}) == 0) || (scalar(@mfgo{$prot2}) == 0)){	$mfmax = 0;	}
	else{	$mfmax = compute_spec_sens($prot1, $prot2, 2);		}		
	if((scalar(@ccgo{$prot1}) == 0) || (scalar(@ccgo{$prot2}) == 0)){	$ccmax = 0;	}
	else{	$ccmax = compute_spec_sens($prot1, $prot2, 3);		}		
	
	$bpsq = $bpmax * $bpmax;
	$mfsq = $mfmax * $mfmax;
	$ccsq = $ccmax * $ccmax;

	$gettingThere = $bpsq + $ccsq + $mfsq;
	$funsim = $gettingThere * (1/3);
	$bpmf = ($bpsq + $mfsq) * (1/2);

	#print "funsim $prot1 $prot2 $funsim\n";
	if($funsim < $funsim_cutoff)		#discard small associations
	{
		$funsim = 0;
	}
	return $funsim;
};

sub get_vectorSim_score
{
	$go_i = $_[0];
	$prot1 = $_[1];
	$prot2 = $_[2];

	$nGO = keys %func;
	$match = 0;
	for($i = 0; $i < $nGO; $i++)
	{
		if($label{$prot1}{$i} == $label{$prot2}{$i})
		{
			$match++;
		}
	}
	$sim = $match / $nGO;
	$score = sprintf("%.3f", $sim);
	return $score;
};

sub compute_edge_feature_funsim
{
	$go_i = $_[0];
	$prot = $_[1];
	$nGO = keys %func;
	
	my $neiE_pos = 0;
	my $neiE_neg = 0;
	my $funsim_pos = 0;
	my $funsim_neg = 0;

#	my @neighbors = $net->neighbours($prot);
	foreach $nei(keys %{$edges{$prot}})
	{
		#print "crf.pl : $go_i edge feature prot $prot nei $nei\n";
		if($label{$nei}{$go_i} == 1)		#positive neighbors
		{
			$fs = get_funsim_score($prot, $nei);
#			$fs = get_vectorSim_score($go_i, $prot, $nei);
			$funsim_pos += $fs;
			#print "crf.pl : $go_i edge feature prot $prot nei $nei Positive edge $edge FS $fs\n";
		}
		
		elsif($label{$nei}{$go_i} == 0)		#negative neighbors
		{
			$fs = get_funsim_score($prot, $nei);
#			$fs = get_vectorSim_score($go_i, $prot, $nei);
			$funsim_neg += $fs;
			#print "crf.pl : $go_i edge feature prot $prot nei $nei Negative edge $edge FS $fs\n";
		}
	}

	#print "crf.pl : $go_i feature neiE_pos -> prot $prot GO $go_i score $neiE_pos FS $funsim_pos\n";
	#print "crf.pl : $go_i feature neiE_neg -> prot $prot GO $go_i score $neiE_neg FS $funsim_neg\n";

	$score = sprintf("%.4f", $funsim_pos); 	#+ random_normal(1, 0, 0.01));
	$features{$prot}{"funsim_pos"} = $score;
	$score = sprintf("%.4f", $funsim_neg); 	#+ random_normal(1, 0, 0.01));
	$features{$prot}{"funsim_neg"} = $score;
};

sub compute_edge_feature
{
	$go_i = $_[0];
	$prot = $_[1];
	$nGO = keys %func;
	
	my $neiE_pos = 0;
	my $neiE_neg = 0;
	my $funsim_pos = 0;
	my $funsim_neg = 0;

#	my @neighbors = $net->neighbours($prot);
	foreach $nei(keys %{$edges{$prot}})
	{
		#print "crf.pl : $go_i edge feature prot $prot nei $nei\n";
		$edge = $edges{$prot}{$nei};		#$net->get_edge_weight($prot, $nei);
		if($label{$nei}{$go_i} == 1)		#positive neighbors
		{
			$neiE_pos += $edge;

			$fs = get_funsim_score($prot, $nei);
#			$fs = get_vectorSim_score($go_i, $prot, $nei);
			$funsim_pos += $fs;
			#print "crf.pl : $go_i edge feature prot $prot nei $nei Positive edge $edge FS $fs\n";
		}
		
		elsif($label{$nei}{$go_i} == 0)		#negative neighbors
		{
			$neiE_neg += $edge;

			$fs = get_funsim_score($prot, $nei);
#			$fs = get_vectorSim_score($go_i, $prot, $nei);
			$funsim_neg += $fs;
			#print "crf.pl : $go_i edge feature prot $prot nei $nei Negative edge $edge FS $fs\n";
		}
	}

	#print "crf.pl : $go_i feature neiE_pos -> prot $prot GO $go_i score $neiE_pos FS $funsim_pos\n";
	#print "crf.pl : $go_i feature neiE_neg -> prot $prot GO $go_i score $neiE_neg FS $funsim_neg\n";

	$score = sprintf("%.4f", $neiE_pos);	# + random_normal(1, 0, 0.01));
	$features{$prot}{"neiE_pos"} = $score;
	$score = sprintf("%.4f", $neiE_neg);	# + random_normal(1, 0, 0.01));
	$features{$prot}{"neiE_neg"} = $score;
	$score = sprintf("%.4f", $funsim_pos); 	#+ random_normal(1, 0, 0.01));
	$features{$prot}{"funsim_pos"} = $score;
	$score = sprintf("%.4f", $funsim_neg); 	#+ random_normal(1, 0, 0.01));
	$features{$prot}{"funsim_neg"} = $score;
};

sub get_features
{
	%features = ();
	$go_i = $_[0];		#mode = 0:training  mode = 1:unknown  mode = 2:all
	$md = $_[1];
#	print "********mode $md go $go_i*********\n";
	foreach $prot(sort keys %training)
	{
		if($md == 0)
		{	
			if($training{$prot} == 1)	#training data
			{
				compute_node_feature($go_i, $prot);
				compute_edge_feature($go_i, $prot);
#				compute_edge_feature_funsim($go_i, $prot);
			}
		}
		if($md == 1)
		{
			if($training{$prot} == 2)	#unknown data
			{
				#$label{$prot}{$go_i} = 1;
				compute_node_feature($go_i, $prot);
				compute_edge_feature($go_i, $prot);	
#				compute_edge_feature_funsim($go_i, $prot);
			}
		}
		if($md == 2)
		{
			#print "crf.pl : get-features for $prot GO $go_i\n";
			compute_node_feature($go_i, $prot);
			compute_edge_feature($go_i, $prot);	
#			compute_edge_feature_funsim($go_i, $prot);
		}
	}

};

sub print_features
{
	$go_i = $_[0];
	open(f, ">features.txt");
	print f "prot neiE_pos neiE_neg funsim_pos funsim_neg FAM_pos FAM_neg label\n";
#	print f "prot neiE_pos neiE_neg funsim_pos funsim_neg label\n";
#	print f "prot neiE_pos neiE_neg label\n";
	foreach $prot(sort keys %features)
	{
		print f "$prot ". " " . $features{$prot}{"neiE_pos"} . " " . $features{$prot}{"neiE_neg"} . " " . $features{$prot}{"funsim_pos"} . " " . $features{$prot}{"funsim_neg"} . " " . $features{$prot}{"FAM_pos"} . " " .  $features{$prot}{"FAM_neg"} . " " ;
#		print f "$prot ". $features{$prot}{"neiE_pos"} . " " . $features{$prot}{"neiE_neg"} . " " . $features{$prot}{"funsim_pos"} . " " . $features{$prot}{"funsim_neg"} . " " ;
#		print f "$prot ". $features{$prot}{"neiE_pos"} . " " . $features{$prot}{"neiE_neg"} . " " ;
		print f $label{$prot}{$go_i} . "\n";
	}
	close(f);
};

sub read_updated_labels
{
	$go_i = $_[0];
	open(LAB, "<upd_label.txt") || die "Cannot open LAB\n";
	while(<LAB>)
	{
		chomp $_;
		@line = split(" ", $_);		#p1 1
		$p = $line[0];
		$label{$p}{$go_i} =  $line[1];
		#print "crf.pl : read updated label prot $p label $line[1]\n";
	}
	close(LAB);
};

sub initialize_unknown_labels
{	
	my $nGO = keys %func;
	$seed = 42 + $pIter;
	srand($seed);		#set seed for random generation
	for($i = 0; $i < $nGO; $i++)
	{
		$c1 = $c2 = $np = 0;
		foreach $prot(keys %training)
		{
			if($training{$prot} == 1)
			{
				$np++;
				if($label{$prot}{$i} == 1)	#protein has go_i
				{
					$c1 = $c1 + 1;
				}
			}
		}
		$c2 = $np - $c1;
		foreach $prot(sort keys %training)
		{
			if($training{$prot} == 2)	#unknown data
			{
				$prob = 0;
				if(($c1+$c2) > 0)
				{
					$prob = ($c1/($c1+$c2));
				}
				$random_label = 0;
#				$randn = random_uniform(1, 0, 1);
				$randn = rand();
				if($prob >= $randn)
				{
					$random_label = 1;	#random_binomial(1,1,$prob);
				}
				$label{$prot}{$i} = $random_label;
#				print "crf.pl : initialize unknown labels for $prot GO $i prob $prob label randn $randn label $random_label\n";
			}
		}
	}#end for loop i

};

sub store_mean_prob
{
	$go_i = $_[0];
	
	open(mp, "<mean_prob") || die "crf.pl: cannot open mean prob\n";
	while(<mp>)
	{
		chomp $_;
		@line = split(" ", $_);		#O14933 0.7525079
		$mean_prob{$line[0]}{$go_i} = $line[1];
	}
	close(mp);
};
sub crf_main
{
	$go_i = $_[0];

	my $R = Statistics::R->new();
	#$R->startR; 

	#calculate features for training
	get_features($go_i, 0);			#get features for training data
	print_features($go_i);
	system("Rscript myglm.R\n");		#call R module for glmfit
#	$R->run_from_file( "myglm.R" );
	print "crf.pl : Done initial training for $go_i\n";

#	#set random values for unknown labels
#	initialize_unknown_labels($go_i);
#	print "crf.pl : Done initializing random labels of unknown for $go_i\n";

	#CRF loop
	$t = 25;
	system("rm metropolis.txt\n");
	system("rm prob1.txt\n");
	for($iter = 0; $iter < $t; $iter++)
	{
		get_features($go_i, 1);			#get features for unknown data
		print_features($go_i);
		system("Rscript inference.R\n");	#call R module for inference
#		$R->run_from_file( "inference.R" );
#		print "crf.pl : iteration $iter for $go_i - DONE inference\n";

		read_updated_labels($go_i);		#unknown data will have updated labels. 
		get_features($go_i, 2);			#get features for all data
		print_features($go_i);
		system("Rscript training.R\n");		#call R module for training
#		$R->run_from_file( "training.R" );
#		print "iter $iter\n";
	}	
	
	#call R module to discard burn in iterations, calculate mean of prob1 of lag period iterations, and make final inference 
	
	print "\nGO $go_i piter $pIter\n";
	get_features($go_i, 1);			#get features for unknown data
	print_features($go_i);
	system("Rscript gibbs.R $t\n");
	read_updated_labels($go_i);	
	store_mean_prob($go_i);
	#$R->stopR();
};

sub clear
{
	%label = ();
	%func = ();
	%gomap = ();
	%training = ();
	%features = ();
	%fam = ();
	%mean_prob = ();
	%edges = ();
	%cat = ();
	%bpgo = ();
	%mfgo = ();
	%ccgo = ();
	%bpsim = ();
	%mfsim = ();
	%ccsim = ();
};


#-----------------Main driver----------------------# 

get_graph();
my $nGO = keys %func;
initialize_unknown_labels();

open(RANDOUT, ">Initial_FlatRandom_slimMatrix.txt") || die "Cant open PRED\n";
foreach $prot(sort keys %training)
{
	if($training{$prot} == 2)
	{
		print RANDOUT $prot;
		print $prot;
		for($i = 0; $i < $nGO; $i++)
		{
			print RANDOUT " ". $label{$prot}{$i};
			print " ". $label{$prot}{$i};
		}
		print RANDOUT "\n";
		print "\n";
	}	
}
close(RANDOUT);

#$nGO = 1;

print "Total $nGO GO iterations of CRF starting....\n";
for($c = 0; $c < $nGO; $c++)
{
	print "----------main : Calling CRF for $c ". $func{$c} . "-----------\n";
	crf_main($c);
}

#output file cluster labels 	 clusterlabel.txt	p1 1 0 0 1 0 0 1 0 1

open(out, ">clusterlabel_crf_$pIter.txt") || die "Cannot open out\n";
foreach $cls(sort keys %training)
{
	print out $cls;
	for($i = 0; $i < $nGO; $i++)
	{
		print out " ".$label{$cls}{$i};
	}
	print out "\n";
	
}
close(out);
clear();
