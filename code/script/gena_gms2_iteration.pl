#!/usr/bin/perl -w
use strict;


USAGE() if scalar(@ARGV) != 8 && scalar(@ARGV) != 9;
my ($seq_f, $anno_f,  $MGM_mod_f, $min_gene_len, $ORDER, $genetic_code, $start_class, $motif_order, $thr ) = @ARGV;
unless($thr){
	$thr = 0;
}

if($genetic_code != 4 && $genetic_code != 11){
	USAGE();
}
my @temp_file;

#my $GENA_DIR = "/storage2/gms2/bin/gena/ncbi";
my $GENA_DIR = "/storage2/starter/experiments/gms2_old_vs_new/gms2_old";
my $HMM_DIR = "/storage2/gms2/bin/gms"; 


my $min_seq_len = 200000;
my $GC = `$GENA_DIR/gc_of_seqs.pl $seq_f`;
$GC =~/^(\d*)/;
$GC = $1;


my $seq_size = -s $seq_f;
if ($seq_size < $min_seq_len){
	print STDERR "Error: sequence too short: $seq_size!\n";
	exit;
}
#--------------------
# run initial MGM prediction using whole-genome GC
# currently does not support genetic code 4
#--------------------
`rm -f performance`;
`rm -f log`;
my $iter = 0;

print "start iteration $iter\n";
run("$HMM_DIR/gmhmmp -g 11  $seq_f  -m $MGM_mod_f  -o itr_$iter.lst");
run("echo \"iteration $iter\" >> performance");
run("$GENA_DIR/compredict.pl $anno_f itr_$iter.lst  3 4 2 3 4 2 0 | tail -1 >> performance");
$iter ++;

#---------------------------
#run iterations until convergence
#---------------------------
## prepare ORFs
run("$GENA_DIR/parse_orfs_for_gms2.pl --in $seq_f --out all.orfs  --min $min_gene_len  --upstr 3 --gcode $genetic_code");

while(1){
	my $pre_iter = $iter - 1;
	if($iter >= 10){
		finish($pre_iter);
	}
	print "Enter interation $iter\n";
	## -----------trainining--------------

	## train coding model
	print "training ...\n";
	
	## find genes predicted by native
	my $PRE_START_LENGTH = 40;
	
	if($pre_iter == 0){
		run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default  --mkmod itr_$iter.mod --NDEC 150  --MINC $min_gene_len  --seq $seq_f --geneset itr_$pre_iter.lst --ORDM 2 --order_non 2 --revcomp_non 1 --pre_start startseq.$iter --PRE_START_WIDTH $PRE_START_LENGTH ");
	}
    elsif($pre_iter >= 4){
		run("$GENA_DIR/classify.pl _all.orfs.MGM.local.$pre_iter _all.orfs.native.$pre_iter   itr_$pre_iter.lst > _itr_$pre_iter.lst");
		run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default  --mkmod itr_$iter.mod --CDEC 800 --NDEC 100 --MINC $min_gene_len --SEPARATE_COD 1  --seq $seq_f --geneset _itr_$pre_iter.lst --ORDM $ORDER --order_non 2 --revcomp_non 1  --pre_start startseq.$iter.40 --PRE_START_WIDTH 40 ");
		run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default  --mkmod itr_$iter.mod --CDEC 800 --NDEC 100 --MINC $min_gene_len --SEPARATE_COD 1  --seq $seq_f --geneset _itr_$pre_iter.lst --ORDM $ORDER --order_non 2 --revcomp_non 1  --pre_start startseq.$iter.20 --PRE_START_WIDTH 22 ");
	}
	else{
		run("$GENA_DIR/classify.pl _all.orfs.MGM.local.$pre_iter _all.orfs.native.$pre_iter   itr_$pre_iter.lst > _itr_$pre_iter.lst");
		run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default  --mkmod itr_$iter.mod  --NDEC 300  --CDEC 150  --MINC $min_gene_len --SEPARATE_COD 1  --seq $seq_f --geneset _itr_$pre_iter.lst --ORDM 2 --order_non 2 --revcomp_non 1 --pre_start startseq.$iter.20 --PRE_START_WIDTH 22");
		run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default  --mkmod itr_$iter.mod  --NDEC 300  --CDEC 150  --MINC $min_gene_len --SEPARATE_COD 1  --seq $seq_f --geneset _itr_$pre_iter.lst --ORDM 2 --order_non 2 --revcomp_non 1 --pre_start startseq.$iter.40 --PRE_START_WIDTH 40");
	}
	
	run("$GENA_DIR/mod_prok_to_euk.pl  itr_$iter.mod > itr_$iter.mod.log");
	push @temp_file, "itr_$iter.mod.log";
	
	## train start context model
	run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default --mkmod prestart_$iter.mod --seq $seq_f --geneset itr_$pre_iter.lst --ORDM 0 --order_non 0 --revcomp_non 1 --fixmotif --PRE_START_WIDTH 3 --width 18");
	
	#######################################
	#    train rbs model
	#######################################
	# input: prestart sequences are in file "startseq.$iter"
	#
	# steps: train rbs model, calculate rbs score
	#
	# output: output score to a file named 'rbs.$iter'
	my $MOTIF_LENGTH = -20;
	if(  $iter > 1 ){	
#		if( 0 ){
			if($start_class == 2){
				#use mixed motif for Synechosystis
				
				run("LD_LIBRARY_PATH=/usr/local/gcc/lib64/ $GENA_DIR/biogem  mfinder_sc  -g $seq_f  -G fasta -l itr_$pre_iter.lst  -L lst -m AGGAGG -M 4 -E 0 -W 20 -v 20 -u 22 -w 6 -o $motif_order -c 2 -p 1 -a nucl -x -n itr_$iter.mod.log -O $iter.karl");	
				run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default --gibbs $iter.karl.motif  --mod rbs_$iter.mod --seq $iter.karl.seqs");
				run("cat itr_$iter.mod rbs_$iter.mod >itr_$iter.rbs.mod");
				#run("$GENA_DIR/hmm -K rbs.$iter.rbs -m itr_$iter.rbs.mod -w -r -z score.$pre_iter  $seq_f");
				run("$GENA_DIR/hmm -K rbs.$iter.rbs -m itr_$iter.rbs.mod  -r -z score.$pre_iter  $seq_f");
				
				run("merge_start.pl   $iter.karl.sc rbs.$iter.rbs  label $iter.karl.class  > rbs.$iter");
				
				
			}
            elsif($start_class == 3){
				#use mixed motif for Halobacterium
				run("LD_LIBRARY_PATH=/usr/local/gcc/lib64/ $GENA_DIR/biogem mfinder_mixed -g $seq_f  -G fasta -w 6 -o $motif_order -r $motif_order -p 1 -a nucl -b itr_$pre_iter.lst -B lst -O $iter.karl -T 22 -u 40 -v 22  -U 22 ");
				
				run("LD_LIBRARY_PATH=/usr/local/gcc/lib64/ $GENA_DIR/biogem  score_rbs -g $seq_f -G fasta -u 40 -w 6 -o $motif_order -p 1 -m itr_$iter.mod.log -r $iter.karl.first >  rbs.$iter.first");
				run("LD_LIBRARY_PATH=/usr/local/gcc/lib64/ $GENA_DIR/biogem  score_rbs -g $seq_f -G fasta -u 22 -w 6 -o $motif_order -p 1 -m itr_$iter.mod.log -r $iter.karl.notfirst > rbs.$iter.notfirst");

				run("classify_first_nofirst.pl itr_$pre_iter.lst  3 4 2 rbs.$iter.first rbs.$iter.notfirst   > rbs.$iter");
			}
			else{
				#motif exists
				run("LD_LIBRARY_PATH=/usr/local/gcc/lib64/ $GENA_DIR/biogem mfinder -s startseq.$iter.20  -S fasta -w 6 -o $motif_order -p 1 -a nucl > gibbs_out.$iter "); #Gibbs4
				
				run("$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default --gibbs gibbs_out.$iter  --mod rbs_$iter.mod --seq startseq.$iter.20");
				run("cat itr_$iter.mod rbs_$iter.mod >itr_$iter.rbs.mod");
				run("$GENA_DIR/hmm -K rbs.$iter -m itr_$iter.rbs.mod  -r -z score.$pre_iter  $seq_f");
				
			}

		push @temp_file, "startseq.$iter";
	}
	#######################################

	## ---------------predicting--------------------
	## prepare model
	run("cat  itr_$iter.mod.log > native.mod.log");
	run("tail $MOTIF_LENGTH  prestart_$iter.mod >> native.mod.log"); 
	run("$GENA_DIR/replace_start_codon.pl  $MGM_mod_f.log itr_$iter.mod  > MGM.mod.log");
	#run("cp  $MGM_mod_f.log  MGM.mod.log");
	run("tail $MOTIF_LENGTH  prestart_$iter.mod >> MGM.mod.log");
	## calculate logodd score

	print "calculate logodd score ...\n";
	run("$GENA_DIR/cal_logodd.pl  all.orfs  native.mod.log  0  0 $GC > all.orfs.native.$iter"); 
	run("$GENA_DIR/merge_score_and_rbs.pl  all.orfs.native.$iter  rbs.$iter > _all.orfs.native.$iter");
	run("$GENA_DIR/cal_logodd.pl  all.orfs  MGM.mod.log  0 0 $GC > all.orfs.MGM.local.$iter");
	run("$GENA_DIR/merge_score_and_rbs.pl  all.orfs.MGM.local.$iter  rbs.$iter | awk '\$5>$thr'  > _all.orfs.MGM.local.$iter");
	run("$GENA_DIR/max_score.pl  1 _all.orfs.native.$iter  _all.orfs.MGM.local.$iter | awk '\$5>0'  > score.$iter"); ## need to change this if $thr<0
	
	## predict genes
	print "predict genes ...\n";
	run("$GENA_DIR/hmm -m itr_$iter.mod  -z score.$iter -x 0.5 -o itr_$iter.lst  $seq_f");
	## record performance
	run("echo \"iteration $iter\" >> performance");
	run("$GENA_DIR/compredict.pl $anno_f itr_$iter.lst  3 4 2 3 4 2 0 | tail -1 >> performance");
	
	push @temp_file, ("all.orfs.native.$iter","_all.orfs.native.$iter","all.orfs.MGM.local.$iter","_all.orfs.MGM.local.$iter");
	
	##--------------------------------compare---------------------------------
	my $sim = `$HMM_DIR/probuild --par $HMM_DIR/par_$genetic_code.default  --compare --source itr_$iter.lst --target itr_$pre_iter.lst`;
	if( ($sim > 0.99 && $iter > 4) || $iter >=9 ){
		print "Finish at iteration $iter!\n";
		finish($iter);
	}

	$iter ++;
}

sub finish{
	my $iter = shift;
	run("cp itr_$iter.mod  GMS2.mod");
	run("$GENA_DIR/classify.pl _all.orfs.MGM.local.$iter _all.orfs.native.$iter  itr_$iter.lst > GMS2.lst");
	if($start_class == 1 || $start_class == 2){
		run("cp itr_$iter.rbs.mod rbs.mod");
	}else{
		run("cp itr_$iter.rbs.first.mod rbs.mod");
	}
	run("rm -f _all");
	exit;
}


sub run{
	my $command = shift;
	open(FILE, ">>log");
	print FILE $command."\n";
	system($command);
	close(FILE);
}
sub USAGE{
	print "Usage: gms2_iteration.pl <seq.fna> <ptt.lst> <MGM.mod> <min_gene_len> <order> <genetic_code(4|11)> <start_class> <motif_order> (optional: threshold_for_score)\n";
	print "start class: 1) RBS 2)RBS+upstream_signal 3)RBS+promoter\n";
	exit;
}
