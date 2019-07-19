#$file = "/home/khan27/Desktop/Ishita_Projects_May201/Group_Function_Prediction/Networks/GO_association/Data_Files/Human_uniprot_annot_full";				#filename = annotation.txt
$file = "annotation";
$path = $ARGV[0];

print "slim-funsim/run.pl: perl $path/split_annot.pl $file $path\n";
system("perl $path/split_annot.pl $file $path\n");

print "slim-funsim/run.pl: perl $path/get_go_bp_score.pl $path\n";
system("perl $path/get_go_bp_score.pl $path\n");		#pairwise spec sens BP
print "slim-funsim/run.pl: perl $path/get_go_mf_score.pl $path\n";	
system("perl $path/get_go_mf_score.pl $path\n");		#pairwise spec sens MF
print "slim-funsim/run.pl: perl $path/get_go_cc_score.pl $path\n";		
system("perl $path/get_go_cc_score.pl $path\n");		#pairwise spec sens CC

print "slim-funsim/run.pl: perl $path/get_fun_sim2.pl $path\n";
system("perl $path/get_fun_sim2.pl $path\n");			#pairwise funsim score of input proteins

