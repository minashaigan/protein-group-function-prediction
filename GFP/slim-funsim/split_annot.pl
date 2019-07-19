#!/usr/bin/perl

$file = $ARGV[0];
$path = $ARGV[1];

open(CAT, "<$path/GO_201401_category.txt") || die "split_annot.pl : Cant open CAT\n";
my %cat = ();
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

#print "get_rel_ancestors.pl : Done loading Category\n";

#separate annotations into BP,MF,CC 
open(ANNOT, "<$path/$file.txt") || die "split_annot.pl : Cant open ANNOT $path/$file.txt\n";

open(bp, ">$path/annotation_bp.txt") || die "split_annot.pl : Cant open annot_bp\n";
open(mf, ">$path/annotation_mf.txt") || die "split_annot.pl : Cant open annot_mf\n";
open(cc, ">$path/annotation_cc.txt") || die "split_annot.pl : Cant open annot_cc\n";

while(<ANNOT>)
{
	chomp $_;
	my @line = split(" ", $_);		#A0AVK6 GO:0033301,GO:0008283,GO:0060718,GO:0070365,GO:0032466,GO:0000122,GO:0001890,GO:0032877,GO:0045944
	my @goterms = split(/,/, $line[1]);
	$prot = $line[0];
	print bp "$prot ";
	print mf "$prot ";
	print cc "$prot ";

	foreach $go(@goterms)
	{
		if($cat{$go} eq 'p')
		{
			print bp "$go,";
			next;
		}
		if($cat{$go} eq 'f')
		{
			print mf "$go,";
			next;
		}
		if($cat{$go} eq 'c')
		{
			print cc "$go,";
			next;
		}
		print "split_annot.pl : $go not found in CAT !!!\n" 
		
	}
	print bp "\n";
	print mf "\n";
	print cc "\n";

}
close(ANNOT);
close(bp);
close(mf);
close(cc);

#print "split_annot.pl : DONE\n";

