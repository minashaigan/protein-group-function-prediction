$data = $ARGV[0];
open(FILE, "<$data") || die "get_annot.pl: Cant open $data\n";
my %genes = ();
while(<FILE>)
{
	chomp $_;
	if(/#/)
	{
		next;
	}
	@line = split(" ", $_);		#CACNA1A	O00555	1
	$genes{$line[1]} = 1;
}
close(FILE);

$n = keys %genes;
print "From $data $n genes loaded\n";

#----------------------Load GO annotation for human--------------------#

$annotfile = "data/Human_uniprot_annot_full.txt";
open(ANNOT, "<$annotfile") || die "get_annot.pl: Cannot open ANNOT\n";
my %annot = ();

while(<ANNOT>)
{
	chomp $_;
	@line = split(" ", $_);			#A0M8Q6 GO:0006956,GO:0006958,GO:0038095,GO:0038096,GO:0045087,GO:0050776,GO:0003823,GO:0005576,GO:0005886,
	@goset = split(",", $line[1]);

	$prot = $line[0];
	if(exists($genes{$prot}))
	{
		foreach $go(@goset)
		{
			push(@{$annot{$prot}}, $go);
		}
	}
}
close(ANNOT);

$n = keys %annot;
print "From $data $n genes have annotation in uniprot\n";


#print

open(out, ">data/group_annot_full.txt");
foreach $prot(sort keys %genes)
{
	print out "$prot ";
	if(exists($annot{$prot}))
	{
		foreach $go(@{$annot{$prot}})
		{
			print out "$go,"
		}
	}
	print out "\n";
}
close(out);
