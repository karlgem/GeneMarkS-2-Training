#!/usr/bin/perl -w


use strict;
use warnings;

USAGE() if scalar(@ARGV) != 8 && scalar(@ARGV) != 9;
my ($seq_f, $anno_f, $MGM_mod_f, $min_gene_len, $ORDER, $genetic_code, $start_class, $motif_order, $thr) = @ARGV;
unless($thr) {
    $thr = 0;
}

# genetic code check
if ($genetic_code != 4 && $genetic_code != 11) {
    USAGE();
}

my $HMM_DIR = "";       # path to HMM
my $GMS2_DIR = "";      # path to gms2 training procedures

# minimum sequence length
my $min_seq_len = 200000;
my $seq_size = -s $seq_f;
if ($seq_size < $min_seq_len) {
    print STDERR "Error: sequence too short: $seq_size!\n";
    exit;
}


#--------------------------------------------------
# Parameters
#--------------------------------------------------
my $ITER_SWITCH = 4;        # iteration after which set 2 parameters are used

#--------------------------------------------------
# Run initial MGM prediction using whole-genome GC
# currently does not support genetic 4
#--------------------------------------------------

`rm -f performance`;
`rm -f log`;

my $iter = 0;

print "Start: iteration 0\n";
run("$HMM_DIR/gmhmmp -g 11 $seq_f -m $MGM_mod_f -o itr_$iter.lst");
run("echo \"iteration $iter\" >> performance");
run("$GMS2_DIR/compredict.pl $anno_f itr_$iter.lst 3 4 2 3 4 2 0 | tail -1 >> performance");

# get genome class
my $genomeClass = `$GMS2_DIR/gmsuite genome-class itr_$iter.lst`;

$iter++;


#--------------------------------------------------
# Run iterations until convergence
#--------------------------------------------------
my $MAX_ITER = 10;

for ($iter = 1; $iter < $MAX_ITER; $iter++) {

    print "Start: iteration $iter\n";

    ## ---------- training ---------- ##
    print "Training ... \n"


    # first 4 iterations are run
    if ($iter <= $ITER_SWITCH) {
        # set parameters
    }
    # run remaining iterations by increasing threshold
    elsif ($iter > $ITER_SWITCH) {
        # set parameters
    }

    # train models based on set parameters and class genome

    ## ---------- prediction ---------- ##
    print "Predicting ... \n"

    # run gmhmm


    ## ---------- check for convergence ---------- ##

}








