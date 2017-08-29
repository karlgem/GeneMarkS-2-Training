#!/usr/bin/perl

# Author: Karl Gemayel
# Created: November 30, 2016
#
# Run the GeneMarkS-2 gene-finder.

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;

# get script name
my $scriptName = basename($0);

# get path of current script
my $scriptPath = abs_path($0);
$scriptPath =~ /^(.*)\//;
$scriptPath = $1;

# my $trainer = "$scriptPath/biogem";        # training of parameters 
my $trainer = "/home/karl/repos/biogem-cpp/code/bin/biogem";
my $predictor = "$scriptPath/gmhmmp2";      # predicting genes

my $comparePrediction = "$scriptPath/compp";    # compare prediction files to check for convergence

# ------------------------------ #
#      Modes for iterations      #
# ------------------------------ # 
my $modeNoMotif             = "no-motif";
my $modePureRBS             = "group-d";
my $modeArchaeaLeaderless   = "group-a";
my $modeBacteriaLeaderless  = "group-b";
my $modeUpstreamSignature   = "group-e";
my $modeClassC              = "group-c";
my @validIterationModes = ($modeNoMotif, $modePureRBS, $modeArchaeaLeaderless, $modeBacteriaLeaderless, $modeUpstreamSignature);

# Class Names
my $classA_mode = "group-a";
my $classB_mode = "group-b";
my $classC_mode = "group-c";
my $classD_mode = "group-d";
my $classE_mode = "group-e";
my $classA2_mode = "group-a2";


# ------------------------------ #
#    Default Variable Values     #
# ------------------------------ # 
my $D_GENETIC_CODE                      = 11                                ;
my $D_FNOUTPUT                          = "gms2.lst"                        ;
my $D_FORMAT_OUTPUT                     = "lst"                             ;
my $D_MGMTYPE                           = "auto"                            ;
my $D_PROM_WIDTH_A                      = 12                                ;
my $D_PROM_WIDTH_B                      = 6                                 ;
my $D_RBS_WIDTH                         = 6                                 ;
my $D_PROM_UPSTR_LEN_A                  = 40                                ;
my $D_PROM_UPSTR_LEN_B                  = 20                                ;
my $D_RBS_UPSTR_LEN                     = 20                                ;
my $D_SPACER_SCORE_THRESH_A             = 0.1                               ;
my $D_SPACER_SCORE_THRESH_B             = 0.25                              ;
my $D_SPACER_DIST_THRESH                = 14                                ;
my $D_SPACER_WINDOW_SIZE                = 1                                 ;
my $D_16S                               = "TAAGGAGGTGA"                     ;
my $D_MIN_MATCH_16S                     = 4                                 ;
my $D_MIN_MATCH_RBS_PROM                = 3                                 ;
my $D_MIN_FRAC_RBS_16S_MATCH            = 0.5                               ;
my $D_UPSTR_SIG_LENGTH                  = 35                                ;
my $D_UPSTR_SIG_ORDER                   = 2                                 ;
my $D_MAX_ITER                          = 10                                ;
my $D_CONV_THRESH                       = 0.99                              ;
my $D_COD_ORDER                         = 5                                 ;
my $D_NONCOD_ORDER                      = 2                                 ;
my $D_START_CONTEXT_ORDER               = 2                                 ;
my $D_FGIO_DIST_THRESH                  = 25                                ;

# ------------------------------ #
#    Command-line variables      #
# ------------------------------ # 
my $fn_genome                                                               ;       # Name of file containing genome sequence
my $genomeType                                                              ;       # Type of genome: Options: archaea, bacteria, auto
my $geneticCode                         = $D_GENETIC_CODE                   ;       # Genetic code
my $fnoutput                            = $D_FNOUTPUT                       ;       # Name of final output file
my $formatOutput                        = $D_FORMAT_OUTPUT                  ;       # Format for output file
my $fnAA                                                                    ;       # amino acid sequences
my $fnNN                                                                    ;       # nucleotide sequences

# Class-A
my $classA_widthPromoter                = $D_PROM_WIDTH_A                   ;
my $classA_widthRBS                     = $D_RBS_WIDTH                      ;
my $classA_promoterUpstreamLength       = $D_PROM_UPSTR_LEN_A               ;
my $classA_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $classA_spacerScoreThresh            = $D_SPACER_SCORE_THRESH_A          ;
my $classA_spacerDistThresh             = $D_SPACER_DIST_THRESH             ;
my $classA_spacerWindowSize             = $D_SPACER_WINDOW_SIZE             ;

# Class-B
my $classB_widthPromoter                = $D_PROM_WIDTH_B                   ;
my $classB_widthRBS                     = $D_RBS_WIDTH                      ;
my $classB_promoterUpstreamLength       = $D_PROM_UPSTR_LEN_B               ;
my $classB_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $classB_spacerScoreThresh            = $D_SPACER_SCORE_THRESH_B          ;
my $classB_spacerWindowSize             = $D_SPACER_WINDOW_SIZE             ;
my $classB_spacerDistThresh             = $D_SPACER_DIST_THRESH             ;
my $classB_tail16S                      = $D_16S                            ;
my $classB_minMatchToTail               = $D_MIN_MATCH_16S                  ;

# Class-C
my $classC_widthRBS                     = $D_RBS_WIDTH                      ;
my $classC_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $classC_minMatchPromoterRBS          = $D_MIN_MATCH_RBS_PROM             ;
my $classC_minMatchRBS16S               = $D_MIN_MATCH_16S                  ;

# Class-D
my $classD_widthRBS                     = $D_RBS_WIDTH                      ;
my $classD_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $classD_percentMatchRBS              = $D_MIN_FRAC_RBS_16S_MATCH         ;
my $classD_minMatchRBS16S               = $D_MIN_MATCH_16S                  ;
my $classD_tail16S                      = $D_16S                            ;

# Class-E
my $classE_widthRBS                     = $D_RBS_WIDTH                      ;
my $classE_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $classE_upstreamSignatureLength      = $D_UPSTR_SIG_LENGTH               ;
my $classE_upstreamSignatureOrder       = $D_UPSTR_SIG_ORDER                ;
my $classE_tail16S                      = $D_16S                            ;

# Iteration control
my $MAX_ITER                            = $D_MAX_ITER                       ;       # number of max iterations in main cycle
my $CONV_THRESH                         = $D_CONV_THRESH                    ;       # convergence threshold
my $numIterWithoutRBS                   = 1                                 ;

# Model Hyperparameters
my $orderCod                            = $D_COD_ORDER                      ;       # order for coding model
my $orderNon                            = $D_NONCOD_ORDER                   ;       # order for noncoding model
my $scOrder                             = $D_START_CONTEXT_ORDER            ;       # start context order
my $fgioDistThresh                      = $D_FGIO_DIST_THRESH               ;
#
# Misc Variables
my $toMgmProb                           = 0.15                              ;
my $toNativeProb                        = 0.85                              ;
my $fixedNativeAtypicalProb;
my $trainNonCodingOnFullGenome;
my $minAtypicalProb                     = 0.02                              ;
my $runMFinderWithoutSpacer;
my $showAdvancedOptions;            
my $mgmType = $D_MGMTYPE                                                    ;       # Type of MGM model: options: "bac, arc, auto"
my $verbose                                                                 ;       # verbose mode
my $keepAllFiles                                                            ;
my $twoStepClassA                                                           ;
my $cutPromTrainSeqs                                                        ;       # train promoter on fragment of upstream region

# Parse command-line options
GetOptions (
    'seq=s'                                 =>  \$fn_genome,
    'genome-type=s'                         =>  \$genomeType,
    'gcode=i'                               =>  \$geneticCode,
    'output=s'                              =>  \$fnoutput,
    'format=s'                              =>  \$formatOutput,
    'faa=s'                                 =>  \$fnAA,
    'fnn=s'                                 =>  \$fnNN,
    # Class-A
    'class-a-width-promoter=i'              =>  \$classA_widthPromoter,
    'class-a-width-rbs=i'                   =>  \$classA_widthRBS,
    'class-a-promoter-upstream-length=i'    =>  \$classA_promoterUpstreamLength,
    'class-a-rbs-upstream-length=i'         =>  \$classA_rbsUpstreamLength,
    'class-a-spacer-score-thresh=f'         =>  \$classA_spacerScoreThresh,
    'class-a-spacer-dist-thresh=i'          =>  \$classA_spacerDistThresh,
    'class-a-spacer-window-size=i'          =>  \$classA_spacerWindowSize,
    # Class-B
    'class-b-width-promoter=i'              =>  \$classB_widthPromoter,
    'class-b-width-rbs=i'                   =>  \$classB_widthRBS,
    'class-b-promoter-upstream-length=i'    =>  \$classB_promoterUpstreamLength,
    'class-b-rbs-upstream-length=i'         =>  \$classB_rbsUpstreamLength,
    'class-b-spacer-score-thresh=f'         =>  \$classB_spacerScoreThresh,
    'class-b-spacer-window-size=i'          =>  \$classB_spacerWindowSize,
    'class-b-tail-16s=s'                    =>  \$classB_tail16S,
    'class-b-min-match-to-tail=i'           =>  \$classB_minMatchToTail,
    # Class-C
    'class-c-width-rbs=i'                   =>  \$classC_widthRBS,
    'class-c-rbs-upstream-length=i'         =>  \$classC_rbsUpstreamLength,
    'class-c-min-match-promoter-rbs=i'      =>  \$classC_minMatchPromoterRBS,
    # Class-D
    'class-d-width-rbs=i'                   =>  \$classD_widthRBS,
    'class-d-rbs-upstream-length=i'         =>  \$classD_rbsUpstreamLength,
    'class-d-percent-match-rbs=f'           =>  \$classD_percentMatchRBS,
    # Class-E
    'class-e-width-rbs=i'                   =>  \$classE_widthRBS,
    'class-e-rbs-upstream-length=i'         =>  \$classE_rbsUpstreamLength,
    'class-e-upstream-signature-length=i'   =>  \$classE_upstreamSignatureLength,
    'class-e-upstream-signature-order=i'    =>  \$classE_upstreamSignatureOrder,
    'class-e-tail-16s=s'                    =>  \$classE_tail16S,
    # Iteration control
    'max-iter=i'                            =>  \$MAX_ITER,
    'conv-thresh=f'                         =>  \$CONV_THRESH,
    # Model Hyperparameters: Orders
    'order-cod=i'                           =>  \$orderCod,
    'order-non=i'                           =>  \$orderNon,
    'order-sc=i'                            =>  \$scOrder,
    # Model Hyperparameters: lengths
    'fgio-dist-thresh=i'                    =>  \$fgioDistThresh,
    # Misc
    'fixed-native-atypical-prob'            =>  \$fixedNativeAtypicalProb,
    'train-noncoding-on-full-genome'        =>  \$trainNonCodingOnFullGenome,
    'min-atypical-prob=f'                   =>  \$minAtypicalProb,
    'run-mfinder-without-spacer'            =>  \$runMFinderWithoutSpacer,
    'v'                                     =>  \$verbose,
    'advanced-options'                      =>  \$showAdvancedOptions,
    'mgm-type=s'                            =>  \$mgmType,
    'keep-all-files'                        =>  \$keepAllFiles,
    'two-step-class-a'                      =>  \$twoStepClassA,
    'cut-prom-train-seqs'                   =>  \$cutPromTrainSeqs,
);

Usage($scriptName) if (!defined $fn_genome or !defined $genomeType or !isValidGenomeType($genomeType));


# setup temporary file collection
my @tempFiles;

my $extraGroupA = "";
if (defined $cutPromTrainSeqs) {
    $extraGroupA = "--cut-prom-train-seqs";
}

# create "single fasta format" from multifasta file
my $fnseq = "tmpseq.fna";
MultiToSingleFASTA($fn_genome, $fnseq);

# add temporary files
push @tempFiles, ($fnseq) unless $keepAllFiles;

my $mgmMod = "$scriptPath/mgm_$geneticCode.mod";        # name of MGM mod file (based on genetic code)
my $modForFinalPred = "tmp.mod";                        # used to keep a version of the model at every iteration 

my $alignmentInMFinder = "right";
if (defined $runMFinderWithoutSpacer) {
    $alignmentInMFinder = "none";
}


my $runClassA = ($genomeType eq "archaea" or $genomeType eq "auto");
my $runClassB = ($genomeType eq "bacteria" or $genomeType eq "auto");
my $runClassATwoStep = ($runClassA and defined $twoStepClassA);

#----------------------------------------
# Run initial MGM prediction
#----------------------------------------
my $mgmPred = CreatePredFileName("0");                  # create a prediction filename for iteration 0
#run("$scriptPath/gmhmmp2 -M $mgmMod -s $fnseq -o $mgmPred --mgm_type $mgmType ");       # Run MGM
run("$predictor -M $mgmMod -s $fnseq -o $mgmPred --mgm_type $mgmType ");       # Run MGM

# add temporary files
push @tempFiles, ($mgmPred) unless $keepAllFiles;


# Compute probability of bacteria #bac/total; add probability to native mod file
my ($bacProb, $arcProb) = EstimateBacArc($mgmPred);

#----------------------------------------
# Main Cycle: Run X iterations 
#----------------------------------------
my $prevPred = $mgmPred;        # Previous iteration prediction: start with MGM predictions
my $prevMod = $mgmMod;          # Previous iteration model:      start with MGM model

# Run iterations without RBS
my $iterBegin = 1;
my $iterEnd = $numIterWithoutRBS;
my $prevIter = RunIterations( { "mode" => $modeNoMotif, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );

$prevPred = "itr_$prevIter.lst";
$prevMod = "itr_$prevIter.mod";


# run single class A iteration
if ($runClassA and !$runClassA) {           # HERE
    $iterBegin = $prevIter + 1;
    $iterEnd   = $iterBegin;            # run single iteration
    $prevIter = RunIterations( { "mode" => $classA_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd} );

    $prevPred = CreatePredFileName($prevIter);
    $prevMod = CreateModFileName($prevIter);

    print "Ran class A iteration: $prevIter\n" if defined $verbose;
}


# If: class A testing enabled, and FGIO have a singnal localized > 14nt         HERE
if ($runClassA and !$runClassA and FGIOHaveSignalAfterThresh($prevIter, $classA_spacerDistThresh, $classA_spacerScoreThresh, $classA_spacerWindowSize)) {
    my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);
    $iterBegin = $prevIter + 1;
    $iterEnd = $iterBegin + $numIterRemain - 1;
    $prevIter = RunIterations( { "mode" => $classA_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd, "extra" => $extraGroupA });
}
# Otherwise, either class A disabled, or FGIO don't have localized signal > 14nt
else {

    # go back one iteration (to cancel archaea)
    my $prevModCancel = CreateModFileName($prevIter);
    my $prevPredCancel = CreatePredFileName($prevIter);

    if ($runClassA and !$runClassA) {
        `mv $prevModCancel classA.mod`;
        `mv $prevPredCancel classA.lst`;

        print "Entering class B: Moved model $prevModCancel\n" if defined $verbose;

        $prevIter -= 1;
        $prevPred = CreatePredFileName($prevIter);
        $prevMod = CreateModFileName($prevIter);

        print "Moved back one iteration to: $prevIter\n" if defined $verbose;
    }

    $prevPred = CreatePredFileName($prevIter);
    $prevMod = CreateModFileName($prevIter);

    if ($runClassATwoStep) {
        # single iteration
        $iterBegin = $prevIter + 1;
        $iterEnd   = $iterBegin;            # run single iteration
        $prevIter = RunIterations( { "mode" => $classA2_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );

        $prevPred = CreatePredFileName($prevIter);
        $prevMod = CreateModFileName($prevIter);

        print "Ran class A2 iteration: $prevIter\n" if defined $verbose;
    }


    if ($runClassATwoStep and FGIOHaveSignalAfterThresh($prevIter, $classA_spacerDistThresh, $classA_spacerScoreThresh, $classA_spacerWindowSize)) {
        my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);
        $iterBegin = $prevIter + 1;
        $iterEnd = $iterBegin + $numIterRemain - 1;
        $prevIter = RunIterations( { "mode" => $classA2_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd, "extra" => $extraGroupA } );

        $prevPred = CreatePredFileName($prevIter);
        $prevMod = CreateModFileName($prevIter);
    }
    else {


        if ($runClassATwoStep) {
            # go back one iteration (to cancel archaea)
            my $prevModCancel = CreateModFileName($prevIter);
            my $prevPredCancel = CreatePredFileName($prevIter);


            `mv $prevModCancel classA2.mod`;
            `mv $prevPredCancel classA2.lst`;

            print "Entering class B: Moved model $prevModCancel\n" if defined $verbose;

            $prevIter -= 1;
            $prevPred = CreatePredFileName($prevIter);
            $prevMod = CreateModFileName($prevIter);

            print "Moved back one iteration to: $prevIter\n" if defined $verbose;
        }

        # run single iteration of class B
        # if ($runClassB) {
            $iterBegin = $prevIter + 1;
            $iterEnd = $iterBegin;          # Run single iteration
            $prevIter = RunIterations( { "mode" => $classB_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );

            $prevPred = CreatePredFileName($prevIter);
            $prevMod = CreateModFileName($prevIter);
        # }
        
        # If: class B testing enabled, and FGIO-16S have signal localized <= 14nt, and promoter/RBS don't match
        if($runClassB and FGIONotMatching16SHaveSignalBeforeThresh($prevIter, $classB_spacerDistThresh, $classB_spacerScoreThresh, $classB_spacerWindowSize) and not PromoterAndRBSConsensusMatch($prevIter, $classC_minMatchPromoterRBS) ) {
            my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);
            $iterBegin = $prevIter + 1;
            $iterEnd = $iterBegin + $numIterRemain - 1;
            $prevIter = RunIterations( { "mode" => $classB_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );

            $prevPred = CreatePredFileName($prevIter);
            $prevMod = CreateModFileName($prevIter);
        }
        # Otherwise, either class B disabled, or FGIO-16S don't have localized signal <= 14nt, or promoter/RBS match
        else {

            my $promoterAndRBSMatched = PromoterAndRBSConsensusMatch($prevIter, $classC_minMatchPromoterRBS);;
            # go back one iteration (to cancel class C)
            $prevModCancel = CreateModFileName($prevIter);
            $prevPredCancel = CreatePredFileName($prevIter);

            # if ($runClassB) {
                `mv $prevModCancel classB.mod`;
                `mv $prevPredCancel classB.lst`;

                # add temporary files
                push @tempFiles, ("classB.mod", "classB.lst") unless $keepAllFiles;

                $prevIter -= 1;
            # }
            # elsif (!$runClassB and $runClassA) {
            #     `mv $prevModCancel classA.mod`;
            #     `mv $prevPredCancel classA.lst`;   

            #     push @tempFiles, ("classA.mod", "classA.lst") unless $keepAllFiles;
            # }

            

            $prevPred = CreatePredFileName($prevIter);
            $prevMod = CreateModFileName($prevIter);

            # run single iteration
            $iterBegin = $prevIter + 1;
            $iterEnd = $iterBegin;
            $prevIter = RunIterations( { "mode" => $classC_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );

            $prevPred = CreatePredFileName($prevIter);
            $prevMod = CreateModFileName($prevIter);

            # if RBS consensus doesn't match 16S tail, then this is non-canonical rbs
            # if ($promoterAndRBSMatched and !RBSConsensusAnd16SMatch($prevIter, $classC_minMatchRBS16S)) {     # matchpr
            if (!RBSConsensusAnd16SMatch($prevIter, $classC_minMatchRBS16S)) {
                my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);
                $iterBegin = $prevIter + 1;
                $iterEnd = $iterBegin + $numIterRemain - 1;
                $prevIter = RunIterations( { "mode" => $classC_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
            }
            # class D (canonical rbs)
            else {

                # go back one iteration (to cancel class C)
                $prevModCancel = CreateModFileName($prevIter);
                $prevPredCancel = CreatePredFileName($prevIter);
                `mv $prevModCancel classC.mod`;
                `mv $prevPredCancel classC.lst`;

                push @tempFiles, ("classC.mod", "classC.lst") unless $keepAllFiles;

                $prevIter -= 1;

                $prevPred = CreatePredFileName($prevIter);
                $prevMod = CreateModFileName($prevIter);
                
                # run single iteration
                $iterBegin = $prevIter + 1;
                $iterEnd = $iterBegin;
                $prevIter = RunIterations( { "mode" => $classD_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd  } );

                $prevPred = CreatePredFileName($prevIter);
                $prevMod = CreateModFileName($prevIter);

                if (PredictedRBSMatch16S($prevPred, $classD_tail16S, $classD_minMatchRBS16S)) {
                    my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);
                    $iterBegin = $prevIter+1;
                    $iterEnd = $iterBegin + $numIterRemain - 1;
                    $prevIter = RunIterations( { "mode" => $classD_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
                }
                # class E
                else {
                    # go back one iteration (to cancel class D)
                    $prevModCancel = CreateModFileName($prevIter);
                    $prevPredCancel = CreatePredFileName($prevIter);
                    `mv $prevModCancel classD.mod`;
                    `mv $prevPredCancel classD.lst`;

                    push @tempFiles, ("classD.mod", "classD.lst") unless $keepAllFiles;

                    $prevIter -= 1;
                    $prevPred = CreatePredFileName($prevIter);
                    $prevMod = CreateModFileName($prevIter);

                    my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);
                    $iterBegin = $prevIter + 1;
                    $iterEnd =  $iterBegin + $numIterRemain - 1;
                    $prevIter = RunIterations( { "mode" => $classE_mode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
                } 
            }
        }
    }
}


$prevPred = CreatePredFileName($prevIter);      # Prediction file: get name of previous iteration
$prevMod = CreateModFileName($prevIter);        # Model file: get name of previous iteration;



#----------------------------------------
# Clean up and get scores
#----------------------------------------


my $finalPred = $fnoutput;
my $finalMod = "GMS2.mod";
my $finalMGM = $mgmMod;  
run("mv $modForFinalPred $finalMod");



# add bacteria and archaea probability to modfile
AddToModel($finalMod, "TO_ATYPICAL_FIRST_BACTERIA", $bacProb);
AddToModel($finalMod, "TO_ATYPICAL_SECOND_ARCHAEA", $arcProb);

if (not $fixedNativeAtypicalProb) {
    ($toNativeProb, $toMgmProb) = EstimateNativeAtypical($prevPred);
}
# add mgm and native probabilities to modfile
AddToModel($finalMod, "TO_MGM", $toMgmProb);
AddToModel($finalMod, "TO_NATIVE", $toNativeProb);
# AddToModel($finalMod, "TO_NATIVE", 15);

# check for fnn and faa output options
my $extraOutput = "";
$extraOutput .= "--AA $fnAA " if defined $fnAA;
$extraOutput .= "--NT $fnNN " if defined $fnNN;


run("$predictor -m $finalMod -M $finalMGM -s $fn_genome -o $finalPred --format $formatOutput $extraOutput");
# run("$predictor -m $finalMod -M $finalMGM -s $fn_genome -o $finalPred.mgm.algo_test   --format $formatOutput --native_prob 0.0000000001 --mgm_prob 1");
# run("$predictor -m $finalMod -M $finalMGM -s $fn_genome -o $finalPred.norm.algo_test  --format $formatOutput --native_prob $toNativeProb ");


run ("rm -f @tempFiles");


# Calculate number of iterations remaining (until we reach max number of allowed iterations)
sub NumOfIterRemaining {
    my ($prevIter, $maxIter) = @_;

    return $maxIter - $prevIter;
}

# Return true if the FGIO that don't match 16S have a localized signal located before the distance threshold
sub FGIONotMatching16SHaveSignalBeforeThresh {
    my ($prevIter, $distThresh, $scoreThresh, $windowSize) = @_;

    $prevPred = "itr_$prevIter.lst";
    $prevMod = "itr_$prevIter.mod";
    #my $isBacteriaProm = run("$trainer experiment promoter-is-valid-for-bacteria --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh");
    my $isBacteriaProm = run("$trainer experiment promoter-is-valid-for-bacteria --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh --window-size $windowSize --min-leaderless-percent 11 --min-leaderless-count 100 --fnlabels $prevPred --fnseq $fnseq");

    chomp $isBacteriaProm;
    return $isBacteriaProm eq "yes";
}

sub RBSSignalLocalized {
    my ($prevIter, $distThresh, $scoreThresh, $windowSize) = @_;

    my $fnmod = "itr_$prevIter.mod";

    my $rbsIsLocalized = run("$trainer experiment rbs-is-localized --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh --window-size $windowSize");
    chomp $rbsIsLocalized;

    return $rbsIsLocalized eq "yes";
}

# Return true if the FGIO have a localized motif signal located further than a distance threshold
sub FGIOHaveSignalAfterThresh {
    my ($prevIter, $distThresh, $scoreThresh, $windowSize) = @_;

    $prevMod = "itr_$prevIter.mod";
    #my $isArchaea = run("$trainer experiment promoter-is-valid-for-archaea --fnmod $prevMod --dist-thresh $distThresh");
    my $isArchaea = run("$trainer experiment promoter-is-valid-for-archaea --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh --window-size $windowSize");

    chomp $isArchaea;

    return $isArchaea eq "yes";
}

# Return true if the fraction of predicted RBS is greater than the threshold
sub PredictedRBSMatch16S {
    my ($fnpred, $seq16S, $minMatch) = @_;

    my $rbsMatchedOutput = run("$trainer experiment match-rbs-to-16s --fnlabels $fnpred --match-to $seq16S --min-match $minMatch --allow-ag-sub");

    my @matchInfo = split(' ', $rbsMatchedOutput);
    my $percentMatched = $matchInfo[1] / $matchInfo[0];

    # my $denom=system('cat $fnpred | grep -E "(native|atypical)[[:space:]]+[ACGT]+[[:space:]]+[[:digit:]]+[[:space:]]+1[[:space:]]*" | wc -l');
    # $percentMatched = $matchInfo[1] / $denom;

    print "Percent of matched RBS: $percentMatched\n" if defined $verbose;

    return $percentMatched >= $classD_percentMatchRBS;
}

# Return true if the Promoter and RBS model consensus sequences match each other
sub PromoterAndRBSConsensusMatch {
    my ($prevIter, $minMatch) = @_;

    my $fnmod = "itr_$prevIter.mod";

    my $isClassC = run("$trainer experiment promoter-and-rbs-match --fnmod $fnmod --match-thresh $minMatch");
    chomp $isClassC;

    return $isClassC eq "yes";
}

# Return true if the RBS model consensus matches the 16S tail
sub RBSConsensusAnd16SMatch {
    my ($prevIter, $minMatch) = @_;

    my $fnmod = "itr_$prevIter.mod";

    my $isMatched = run("$trainer experiment rbs-consensus-and-16s-match --fnmod $fnmod --allow-ag-sub");
    chomp $isMatched;

    return $isMatched eq "yes";
}

# Run GMS2 iterations in a particular mode
sub RunIterations {
    my %params = %{ $_[0] };

    my $mode        = $params{"mode"};
    my $iterBegin   = $params{"iteration-begin"};
    my $iterEnd     = $params{"iteration-end"};

    my $extraOptions  = "";
    if (exists $params{"extra"}) {
        $extraOptions = $params{"extra"};
    }

    my $iter = $iterBegin;
    while ($iter <= $iterEnd) {

        print "Mode $mode: Entering iteration $iter...\n" if defined $verbose;

        # train on native genes
        my $nativeOnly = 0;
        if ($iter > 1) {
            $nativeOnly = 1;
        }

        my $currMod  = CreateModFileName($iter);        # model file for current iteration
        my $currPred = CreatePredFileName($iter);       # prediction file for current iteration

        my $currMFinderResult = CreateMFinderResultFileName($iter);
        my $prevMFinderResult = CreateMFinderResultFileName($iter-1);


        # Training step: use prediction of previous iteration
        # my $trainingCommand = "$trainer gms2-training -s $fnseq -l $prevPred -m $currMod --coding-order $orderCod --noncoding-order $orderNon --train-on-native-only $nativeOnly --genetic-code $geneticCode --sc-order $scOrder $noncodingOnFullGenome --fgio-dist-thresh $fgioDistThresh";

       my $trainingCommand = "$trainer gms2-training -s $fnseq -l $prevPred -m $currMod --order-coding $orderCod --order-noncoding $orderNon --only-train-on-native $nativeOnly --genetic-code $geneticCode --order-start-context $scOrder --fgio-dist-thr $fgioDistThresh";

        $trainingCommand .= " $extraOptions ";

        # if version of mfinder supports reusing previous iteration
        if ($trainer eq "$scriptPath/kgsuite_mfinder_reuse") {
            if ($iter == $iterBegin) {
                $trainingCommand .= " --fn-mfinder-output $currMFinderResult";
            }
            elsif ($iter > $iterBegin) {
                $trainingCommand .= " --fn-mfinder-output $currMFinderResult  --fn-mfinder-input $prevMFinderResult";
            }
        }



        if ($mode eq $modeNoMotif) {
            $trainingCommand .= " --run-motif-search false";
            $trainingCommand .= " --genome-group D";        # FIXME: remove requirement from training
        }
        elsif ($mode eq $classA_mode) {
            $trainingCommand .= " --genome-group A --ga-upstr-len-rbs $classA_rbsUpstreamLength --align $alignmentInMFinder --ga-width-rbs $classA_widthRBS --ga-upstr-len-prom $classA_promoterUpstreamLength --ga-width-prom $classA_widthPromoter";
        }
        elsif ($mode eq $classA2_mode) {
            $trainingCommand .= " --genome-group A2 --ga-upstr-len-rbs $classB_rbsUpstreamLength --align $alignmentInMFinder --ga-width-rbs $classB_widthRBS --ga-upstr-len-prom $classA_promoterUpstreamLength --ga-width-prom $classA_widthPromoter --ga-extended-sd $classB_tail16S";            
        }
        elsif ($mode eq $classB_mode) {
            $trainingCommand .= " --genome-group B --gb-upstr-len-rbs $classB_rbsUpstreamLength --align $alignmentInMFinder --gb-width-rbs $classB_widthRBS --gb-upstr-len-prom $classB_promoterUpstreamLength --gb-width-prom $classB_widthPromoter --gb-extended-sd $classB_tail16S";
        }
        elsif ($mode eq $classC_mode) {
            $trainingCommand .= " --genome-group C --gc-upstr-len-rbs $classC_rbsUpstreamLength --align $alignmentInMFinder --gc-width-rbs $classC_widthRBS";
        }
        elsif ($mode eq $classD_mode) {
            $trainingCommand .= " --genome-group D --gd-upstr-len-rbs $classD_rbsUpstreamLength --align $alignmentInMFinder --gd-width-rbs $classD_widthRBS";
        }
        elsif ($mode eq $classE_mode) {
            $trainingCommand .= " --genome-group E --ge-upstr-len-rbs $classE_rbsUpstreamLength --align $alignmentInMFinder --ge-width-rbs $classE_widthRBS --ge-len-upstr-sig $classE_upstreamSignatureLength --ge-order-upstr-sig $classE_upstreamSignatureOrder --ge-extended-sd $classE_tail16S";
        }
        else {
            die "Mode invalid: should not reach this point";
        }

        run("$trainingCommand");

        # if ($mode eq $classA2_mode) {
        #     run("sed -i 's/GENOME_TYPE bacteria-promoter/GENOME_TYPE archaea-promoter/g' $currMod");
        # }

        `cp $currMod $modForFinalPred`;

        # add bacteria and archaea probability to model file
        AddToModel($currMod, "TO_ATYPICAL_FIRST_BACTERIA", $bacProb);
        AddToModel($currMod, "TO_ATYPICAL_SECOND_ARCHAEA", $arcProb);

        if (not $fixedNativeAtypicalProb and $iter > 1) {
            ($toNativeProb, $toMgmProb) = EstimateNativeAtypical($prevPred);
        }
        # add mgm and native probabilities to modfile
        AddToModel($currMod, "TO_MGM", $toMgmProb);
        AddToModel($currMod, "TO_NATIVE", $toNativeProb);

        # Prediction step: using current model file
        my $errCode = run("$predictor -m $currMod -M $mgmMod -s $fnseq -o $currPred --format train");

        # Check for convergence
        my $similarity = run("$comparePrediction -n -a $prevPred -b $currPred -G");

        print "Iteration : $similarity\n" if defined $verbose;


        # add temporary files
        push @tempFiles, ($currMod, $currPred) unless $keepAllFiles;
        

        if ( ($similarity > 99 && $iter > 2) ) {
            print "Converged at iteration $iter\n" if defined $verbose;
            return $iter;
        }

        # set previous prediction (before exiting loop)
        $prevPred = $currPred;
        $prevMod = $currMod;


        $iter++;
    }

    return $iter-1;
}

# Run a system command and log it
sub run {
    my $command = shift;
    open(FILE, ">>log");
    print FILE $command . "\n";
    my $value = `$command`;
    chomp($value);
    return $value;
}

# Estimate bacteria and archaea probabilities based on the counts in the prediction file
sub EstimateBacArc {
    my $fname = shift;
    my $counts_all = 0;
    my $counts_bac = 0;

    my $min_gene_length_bac_arc = 600;

    open( my $IN , $fname ) or die "error on open file $fname: $!\n";
    while( my $line = <$IN>)
    {
        next if ( $line =~ /^\s*$/ );
        next if ( $line =~ /^#/ );
        next if ( $line =~ /SequenceID:/ );

        if ( $line =~ /^\s*\d+\s+[+-]\s+\S+\s+\S+\s+(\d+)\s+bac\s*/ )
        {
            if ( $1 >= $min_gene_length_bac_arc )
            {
                ++$counts_bac;
                ++$counts_all;
            }
        }
        elsif ( $line =~ /^\s*\d+\s+[+-]\s+\S+\s+\S+\s+(\d+)\s+arc\s*/ )
        {
            if ( $1 >= $min_gene_length_bac_arc )
            {
                ++$counts_all;
            }
        }
        else {die;}
    }
    close $IN;

    if (!$counts_all) {print "error, unexpected format foundin file: $fname"; exit 1; }

    if (defined $verbose) {
        my $counts_arc = $counts_all - $counts_bac;
        my $bacProb = $counts_bac / $counts_all;
        my $arcProb = $counts_arc / $counts_all;
        print "NumBac = $counts_bac\n";
        print "NumArc = $counts_arc\n";
        print "Bacteria Probability: $bacProb\n";
        print "Archaea Probability: $arcProb\n";
    }

    return ( sprintf( "%.5f", $counts_bac/$counts_all ), sprintf( "%.5f", ($counts_all - $counts_bac)/$counts_all ) );
}

# Esitmate native and atypical probabilities based on the counts in the prediction file
sub EstimateNativeAtypical {
    my $fname = shift;
    my $counts_all = 0;
    my $counts_native = 0;

    my $min_gene_length_native_atypical = 600;

    open( my $IN , $fname ) or die "error on open file $fname: $!\n";
    while( my $line = <$IN>)
    {
        next if ( $line =~ /^\s*$/ );
        next if ( $line =~ /^#/ );
        next if ( $line =~ /SequenceID:/ );

        my $currLength = -1;
        if ( $line =~ /^\s*\d+\s+[+-]\s+\d+\s+\d+\s+(\d+)\s*/ ) {
            $currLength = $1;
        }

        # if ( $line =~ /^\s*\d+\s+[+-]\s+\d+\s+\d+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+native\s*/ )
        if ( $line =~ /\s+native\s*$/ ) 
        {
            if ( $currLength >= $min_gene_length_native_atypical )
            {
                ++$counts_native;
                ++$counts_all;
            }
        }
        #elsif ( $line =~ /^\s*\d+\s+[+-]\s+\d+\s+\d+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+atypical\s*/ )
        elsif ( $line =~ /\s+atypical\s*$/ ) 
        {
            if ( $currLength >= $min_gene_length_native_atypical )
            {
                ++$counts_all;
            }
        }
        #else {die;}
    }
    close $IN;

    if (!$counts_all) {print "error, unexpected format foundin file: $fname"; exit 1; }

    if (defined $verbose) {
        my $counts_atypical = $counts_all - $counts_native;
        my $nativeProb = $counts_native / $counts_all;
        my $atypicalProb = $counts_atypical / $counts_all;
        print "NumNative = $counts_native\n";
        print "NumAtypical = $counts_atypical\n";
        print "Native Probability: $nativeProb\n";
        print "Atypical Probability: $atypicalProb\n";
    }

    my $probNative = $counts_native/$counts_all;
    my $probAtypical = ($counts_all - $counts_native)/$counts_all;
    if ($probAtypical < $minAtypicalProb) {
        $probAtypical = $minAtypicalProb;
        $probNative = 1 - $probAtypical;
    }

    return ( sprintf( "%.5f", $probNative ), sprintf( "%.5f", $probAtypical) );
}

# Converts a multifasta file to single fasta by concatenating all sequences into one
sub MultiToSingleFASTA {
    # input and output filenames
    my ($fnin, $fnout) = @_;
    
    run ("echo '>anydef' > $fnout");
    run ("grep -v '>' $fnin >> $fnout");
    return;
}

# Add label/value pair to a model file
sub AddToModel {
    my ( $fname, $label, $value ) = @_;
    open (my $fout, ">>", $fname) or die "Error: Could not open file $fname\n";

    print $fout "\$" . $label . " $value\n";
    close $fout;
}

# Create name for model file based on iteration number
sub CreateModFileName {
    my $iter = $_[0];
    return "itr_$iter.mod";
}

# Create name for prediction file based on iteration number
sub CreatePredFileName {
    my $iter = $_[0];
    return "itr_$iter.lst";
}

# Create name for Motif Finder Output file based on iteration number
sub CreateMFinderResultFileName {
    my $iter = $_[0];
    return "itr_$iter.mfinderresult";
}

# Returns true if the genome type is valid: archaea, bacteria, auto
sub isValidGenomeType {
    my $gt = $_[0];
    return ($gt eq "archaea" or $gt eq "bacteria" or $gt eq "auto");
}



# FIXME: fnn, faa, check gcode 11,4, keep-all-files




# Usage function: print usage message and exit script
sub Usage {
    my $name = $_[0];
    print "Usage: $name --seq SEQ --genome-type TYPE";
    print
"
Basic Options: 
--seq                                   File containing genome sequence in FASTA format
--genome-type                           Type of genome: archaea, bacteria, auto 
--gcode                                 The genetic code number (default: $D_GENETIC_CODE. Choices: 11 and 4)
--output                                Name of output file (default: $D_FNOUTPUT)
--format                                Format of output file (default: $D_FORMAT_OUTPUT)
--fnn                                   Name of output file that will hold nucleotide sequences of predicted genes
--faa                                   Name of output file that will hold protein sequences of predicted genes
--advanced-options                      Show the advanced options
";

    if (defined $showAdvancedOptions) {
        print 
"
Advanced Options:
# Iteration control
--max-iter                              Number of max iterations in the main cycle (default: $D_MAX_ITER)
--conv-thresh                           The convergence threshold (in range [0,1]) (default: $D_CONV_THRESH)

# Misc
fixed-native-atypical-prob              Fix the native and atypical prior probabilities
train-noncoding-on-full-genome          Train the non-coding model on the full genome
min-atypical-prob                       Set minimum prior probability for atypical genes
run-mfinder-without-spacer              Disable the \"location distribution\" in the motif finders
mgm-type                                Type of genome model to use for MGM predictions
                                        Option: bac, arc, auto. Default: (default: $D_MGMTYPE)
keep-all-files                          Keep all intermediary files 
two-step-class-a                        If set, the two-step classification is applied to class A
fgio-dist-thresh                        Distance threshold for FGIO identification

# Class-A
class-a-width-promoter                  Width of the promoter motif model (default: $D_PROM_WIDTH_A)
class-a-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
class-a-promoter-upstream-length        Upstream length for promoter training (default: $D_PROM_UPSTR_LEN_A)
class-a-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
class-a-spacer-score-thresh             Minimum peak threshold for the spacer distribution (default: $D_SPACER_SCORE_THRESH_A)
class-a-spacer-dist-thresh              Minimum distance threshold for the spacer distribution (default: $D_SPACER_DIST_THRESH)
class-a-spacer-window-size              Window size for calculating the \"peak value\" to compare
                                        against the score threshold. (default: $D_SPACER_WINDOW_SIZE)
# Class-B
class-b-width-promoter                  Width of the promoter motif model (default: $D_PROM_WIDTH_B)
class-b-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
class-b-promoter-upstream-length        Upstream length for promoter training (default: $D_PROM_UPSTR_LEN_B)
class-b-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
class-b-spacer-score-thresh             Minimum peak threshold for the spacer distribution (default: $D_SPACER_SCORE_THRESH_B)
class-b-spacer-window-size              Window size for calculating the \"peak value\" to compare
                                        against the score threshold (default: $D_SPACER_WINDOW_SIZE)
class-b-tail-16s                        The 16S rRNA tail used for selecting training sequences for
                                        the promoter model (default: $D_16S)
class-b-min-match-to-tail               Minimum number of consecutive nucleotide matches to the 16S (default: $D_MIN_MATCH_16S)

# Class-C
class-c-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
class-c-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
class-c-min-match-promoter-rbs          Minimum number of consecutive nucleotide matches between the 
                                        promoter and RBS (default: $D_MIN_MATCH_RBS_PROM)

# Class-D
class-d-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
class-d-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
class-d-percent-match-rbs               Minimum percentage of predicted RBS sequences that match to 16S (default: $D_MIN_FRAC_RBS_16S_MATCH)

# Class-E
class-e-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
class-e-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
class-e-upstream-signature-length       Length of the upstream-signature Nonuniform Markov model (default: $D_UPSTR_SIG_LENGTH)
class-e-upstream-signature-order        Order of the upstream-signature Nonuniform Markov model (default: $D_UPSTR_SIG_ORDER)
class-e-tail-16s                        The 16S rRNA tail used for selecting training sequences for
                                        the RBS model (default: $D_16S)
";
    }

    exit;
}
