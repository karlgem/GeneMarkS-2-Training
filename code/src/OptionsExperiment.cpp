//
//  OptionsExperiment.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/12/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsExperiment.hpp"
#include "NumSequence.hpp"
#include <iostream>

using namespace std;
using namespace gmsuite;

#define STR_MATCH_SEQ_TO_UPSTREAM   "match-seq-to-upstream"
#define STR_MATCH_SEQ_TO_NONCODING  "match-seq-to-noncoding"
#define STR_BUILD_START_MODELS      "build-start-models"
#define STR_BUILD_START_MODELS2     "build-start-models2"
#define STR_BUILD_START_MODELS3     "build-start-models3"
#define STR_SCORE_STARTS            "score-starts"
#define STR_MATCH_RBS_TO_16S        "match-rbs-to-16s"
#define STR_SCORE_LABELED_STARTS    "score-labeled-starts"
#define STR_PROMOTER_IS_VALID_FOR_ARCHAEA    "promoter-is-valid-for-archaea"
#define STR_PROMOTER_IS_VALID_FOR_BACTERIA    "promoter-is-valid-for-bacteria"
#define STR_START_MODEL_STRATEGY_2  "start-model-strategy-2"
#define STR_PROMOTER_AND_RBS_MATCH "promoter-and-rbs-match"
#define STR_RBS_CONSENSUS_AND_16S_MATCH "rbs-consensus-and-16s-match"
#define STR_RBS_IS_LOCALIZED        "rbs-is-localized"
#define STR_COMPUTE_MOTIF_SCORE_FOR_STARTS "compute-motif-score-for-starts"

namespace gmsuite {
    // convert string to experiment_t
    istream& operator>>(istream& in, OptionsExperiment::experiment_t& unit) {
        string token;
        in >> token;
        
        if      (token == STR_MATCH_SEQ_TO_UPSTREAM)    unit = OptionsExperiment::MATCH_SEQ_TO_UPSTREAM;
        else if (token == STR_MATCH_SEQ_TO_NONCODING)   unit = OptionsExperiment::MATCH_SEQ_TO_NONCODING;
        else if (token == STR_BUILD_START_MODELS)       unit = OptionsExperiment::BUILD_START_MODELS;
        else if (token == STR_BUILD_START_MODELS2)      unit = OptionsExperiment::BUILD_START_MODELS2;
        else if (token == STR_BUILD_START_MODELS3)      unit = OptionsExperiment::BUILD_START_MODELS3;
        else if (token == STR_SCORE_STARTS)             unit = OptionsExperiment::SCORE_STARTS;
        else if (token == STR_MATCH_RBS_TO_16S)         unit = OptionsExperiment::MATCH_RBS_TO_16S;
        else if (token == STR_SCORE_LABELED_STARTS)     unit = OptionsExperiment::SCORE_LABELED_STARTS;
        else if (token == STR_PROMOTER_IS_VALID_FOR_ARCHAEA)     unit = OptionsExperiment::PROMOTER_IS_VALID_FOR_ARCHAEA;
        else if (token == STR_PROMOTER_IS_VALID_FOR_BACTERIA)     unit = OptionsExperiment::PROMOTER_IS_VALID_FOR_BACTERIA;
        else if (token == STR_START_MODEL_STRATEGY_2)   unit = OptionsExperiment::START_MODEL_STRATEGY_2;
        else if (token == STR_PROMOTER_AND_RBS_MATCH)   unit = OptionsExperiment::PROMOTER_AND_RBS_MATCH;
        else if (token == STR_RBS_CONSENSUS_AND_16S_MATCH) unit = OptionsExperiment::RBS_CONSENSUS_AND_16S_MATCH;
        else if (token == STR_RBS_IS_LOCALIZED)         unit = OptionsExperiment::RBS_IS_LOCALIZED;
        else if (token == STR_COMPUTE_MOTIF_SCORE_FOR_STARTS)   unit = OptionsExperiment::COMPUTE_MOTIF_SCORE_FOR_STARTS;
        else
            throw boost::program_options::invalid_option_value(token);
        
        return in;
    }
}

// constructor
OptionsExperiment::OptionsExperiment(string mode) : Options(mode) {
    
}


bool OptionsExperiment::parse(int argc, const char **argv) {
    try {
        
//        GenericOptions genericOptions;
        
        // Group of options allowed only on command line
        po::options_description generic("General Options");
        Options::addProcessOptions_GenericOptions(genericOptions, generic);
        
        // Group of options allowed on both CML and in config files
        po::options_description config("Configuration");
        config.add_options()
        ;
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
            ("mode", po::value<string>(&mode)->required(), "Program Mode")
            ("experiment", po::value<experiment_t>(&experiment)->required(), "Experiment")
            ("subargs", po::value<std::vector<std::string> >(), "Arguments for experiment")
        ;
        
        
        // Congregate options into further groups
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);
        
        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);
        
        po::options_description visible("Allowed options");
        visible.add(generic).add(config);
        
        // Specify positional arguments
        po::positional_options_description pos;
        pos.add("mode",1);
        pos.add("experiment",1);
        pos.add("subargs",-1);
        
        
        
        // create storage component for storing names and values of arguments
        po::variables_map vm;
        
        // store command-line options
        po::parsed_options parsed = po::command_line_parser(argc, argv).          // pass in input
        options(cmdline_options).                     // specify options list
        positional(pos).                              // specify which are positional
        allow_unregistered().                         // allow unregistered args since we haven't yet set-up the utility-based ones
        run();
        
        po::store(parsed,vm);

        if (!vm.count("experiment")){
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // try parsing arguments.
        po::notify(vm);
      
        po::options_description expDesc;
        
        // Experiment: match sequence to upstream regions
        if (experiment == MATCH_SEQ_TO_UPSTREAM) {
            addProcessOptions_MatchSeqToUpstreamOptions(matchSeqToUpstream, expDesc);
        }
        // Experiment: match sequence to simulated noncoding regions
        if (experiment == MATCH_SEQ_TO_NONCODING) {
            addProcessOptions_MatchSeqToNoncodingOptions(matchSeqToNoncoding, expDesc);
        }
        // Experiment: build start models
        if (experiment == BUILD_START_MODELS) {
            addProcessOptions_BuildStartModelsOptions(buildStartModels, expDesc);
        }
        // Experiment: build start models2
        if (experiment == BUILD_START_MODELS2) {
            addProcessOptions_BuildStartModels2Options(buildStartModels2, expDesc);
        }
        // Experiment: build start models3
        if (experiment == BUILD_START_MODELS3) {
            addProcessOptions_BuildStartModels3Options(buildStartModels3, expDesc);
        }
        // Experiment: score starts
        if (experiment == SCORE_STARTS){
            addProcessOptions_ScoreStarts(scoreStarts, expDesc);
        }
        // Experience: match rbs to 16s
        if (experiment == MATCH_RBS_TO_16S) {
            addProcessOptions_MatchRBSTo16SOptions(matchRBSTo16S, expDesc);
        }
        // Experiment: get start-model type
        if (experiment == PROMOTER_IS_VALID_FOR_ARCHAEA)
            addProcessOptions_PromoterIsValidForArchaea(promoterIsValidForArchaea, expDesc);
        // Experiment: get start-model bacteria
        if (experiment == PROMOTER_IS_VALID_FOR_BACTERIA)
            addProcessOptions_PromoterIsValidForBacteria(promoterIsValidForBacteria, expDesc);
        // Experiment: start model strategy 2
        if (experiment == START_MODEL_STRATEGY_2)
            addProcessOptions_StartModelStrategy2Options(startModelStrategy2, expDesc);
        if (experiment == PROMOTER_AND_RBS_MATCH)
            addProcessOptions_PromoterAndRBSMatchOptions(promoterAndRBSMatch, expDesc);
        if (experiment == RBS_CONSENSUS_AND_16S_MATCH)
            addProcessOptions_RBSConsensusAnd16SMatch(rbsConsensusAnd16SMatch, expDesc);
        if (experiment == RBS_IS_LOCALIZED)
            addProcessOptions_RBSIsLocalized(rbsIsLocalized, expDesc);
        if (experiment == COMPUTE_MOTIF_SCORE_FOR_STARTS)
            addProcessOptions_ComputeMotifScoreForStarts(computeMotifScoreForStarts, expDesc);
        
        cmdline_options.add(expDesc);
        
        // Collect all the unrecognized options from the first pass. This will include the
        // (positional) mode and command name, so we need to erase them
        vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());       // erase mode
        opts.erase(opts.begin());       // erase command name
        
        // Parse again...
        po::store(po::command_line_parser(opts).options(expDesc).run(), vm);
        
        
        // if help specified, print usage message and quit
        if (vm["help"].as<bool>()) {
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // update all values and make sure required are provided
        po::notify(vm);
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
}



void OptionsExperiment::addProcessOptions_MatchSeqToUpstreamOptions(MatchSeqToUpstreamOptions &options, po::options_description &processOptions) {
    
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
        ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ;
    
}


void OptionsExperiment::addProcessOptions_MatchSeqToNoncodingOptions(MatchSeqToNoncodingOptions &options, po::options_description &processOptions) {
    
    Options::addProcessOptions_GenReadSeqAndLabelsOptions(options, processOptions);
    processOptions.add_options()
        ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
        ("length", po::value<size_t>(&options.length)->default_value(40), "Length of noncoding sequences")
        ("num-noncoding", po::value<size_t>(&options.numNoncoding)->default_value(1000), "Number of non-coding sequences.")
        ("order-nonc", po::value<unsigned>(&options.order)->default_value(0), "Order of non-coding model")
        ("pcounts", po::value<double>(&options.pcounts)->default_value(1), "Pseudocounts for noncoding model")
    ;
}

void OptionsExperiment::addProcessOptions_BuildStartModelsOptions(BuildStartModelsOptions &options, po::options_description &processOptions) {
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
    ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ("min-match", po::value<unsigned>(&options.min16SMatch)->default_value(4), "Minimum accepted match length from upstream to 16S tail")
    ("allow-ag-sub", po::bool_switch(&options.allowAGSubstitution)->default_value(false), "Allow G to be substituted for A when matching to 16S tail")
    ;
    
    // mfinder options
    po::options_description mfinder ("Motif Finder");
    OptionsMFinder::addProcessOptions(options.mfinderOptions, mfinder);
    processOptions.add(mfinder);
}

void OptionsExperiment::addProcessOptions_BuildStartModels2Options(BuildStartModels2Options &options, po::options_description &processOptions) {
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
    ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ("min-match", po::value<unsigned>(&options.min16SMatch)->default_value(4), "Minimum accepted match length from upstream to 16S tail")
    ("allow-ag-sub", po::bool_switch(&options.allowAGSubstitution)->default_value(false), "Allow G to be substituted for A when matching to 16S tail")
    ;
    
    // mfinder options
    po::options_description mfinder ("Motif Finder");
    OptionsMFinder::addProcessOptions(options.mfinderOptions, mfinder);
    processOptions.add(mfinder);
}



void OptionsExperiment::addProcessOptions_BuildStartModels3Options(BuildStartModels3Options &options, po::options_description &processOptions) {
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
    ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ("min-match", po::value<unsigned>(&options.min16SMatch)->default_value(4), "Minimum accepted match length from upstream to 16S tail")
    ("allow-ag-sub", po::bool_switch(&options.allowAGSubstitution)->default_value(false), "Allow G to be substituted for A when matching to 16S tail")
    ("nfgio-thresh", po::value<size_t>(&options.nfgioThresh)->default_value(20), "Used to extract second genes in operon")
    ;
    
    // mfinder options
    po::options_description mfinder ("Motif Finder");
    OptionsMFinder::addProcessOptions(options.mfinderOptions, mfinder);
    processOptions.add(mfinder);
}


void OptionsExperiment::addProcessOptions_ScoreStarts(ScoreStarts &options, po::options_description &processOptions) {
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
    ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ("min-match", po::value<unsigned>(&options.min16SMatch)->default_value(4), "Minimum accepted match length from upstream to 16S tail")
    ("allow-ag-sub", po::bool_switch(&options.allowAGSubstitution)->default_value(false), "Allow G to be substituted for A when matching to 16S tail")
    ("nfgio-thresh", po::value<size_t>(&options.nfgioThresh)->default_value(20), "Used to extract second genes in operon")
    ("fgio-thresh", po::value<size_t>(&options.fgioThresh)->default_value(40), "Used to extract first genes in operon")
    ("search-len", po::value<size_t>(&options.searchUpstrLen)->default_value(20), "Determines search window for 16S match")
    ("upstr-len-rbs", po::value<size_t>(&options.upstreamLenRBS)->default_value(20), "Length of upstream for RBS motif search")
    ("upstr-len-promoter", po::value<size_t>(&options.upstreamLenPromoter)->default_value(40), "Length of upstream for Promoter motif search")
    ("gms2-mod", po::value<string>(&options.gms2mod)->required(), "Name of GMS2 mod file containing noncoding model")
    ;
    
    
    // mfinder options
    po::options_description mfinder ("Motif Finder ");
    OptionsMFinder::addProcessOptions(options.mfinderOptions, mfinder);
    processOptions.add(mfinder);
    
    
}

void OptionsExperiment::addProcessOptions_MatchRBSTo16SOptions(MatchRBSTo16S &options, po::options_description &processOptions) {
//    Options::addProcessOptions_GenericOptions(options, processOptions);
    processOptions.add_options()
    ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ("fnlabels", po::value<string>(&options.fnlabels)->required(), "File containing gene labels with predicted RBS.")
    ("min-match", po::value<unsigned>(&options.min16SMatch)->default_value(4), "Minimum number of consecutively matched nucleotides for a match to be considered as a match.")
    ("allow-ag-sub", po::bool_switch(&options.allowAGSubstitution)->default_value(false), "Allow G to be substituted for A when matching to 16S tail")
    ;
    
}


void OptionsExperiment::addProcessOptions_ScoreLabeledStarts(ScoreLabeledStarts &options, po::options_description &processOptions) {
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of GMS2 mod file containing RBS model.")
    ;
}




void OptionsExperiment::addProcessOptions_PromoterIsValidForArchaea(PromoterIsValidForArchaea &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of mod file containing RBS spacer.")
    ("dist-thresh", po::value<size_t>(&options.distanceThresh)->default_value(22), "Distance threshold after which spacer indicates promoter.")
    ("score-thresh", po::value<double>(&options.scoreThresh)->default_value(0.1), "Minimum score above which spacer is considered localized.")
    ("window-size", po::value<size_t>(&options.windowSize)->default_value(1), "Size of window in which to determine localization")
    ;
}

void OptionsExperiment::addProcessOptions_PromoterIsValidForBacteria(PromoterIsValidForBacteria &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of mod file containing RBS spacer.")
    ("dist-thresh", po::value<size_t>(&options.distanceThresh)->default_value(15), "Distance threshold before which spacer indicates promoter.")
    ("score-thresh", po::value<double>(&options.scoreThresh)->default_value(0.1), "Minimum score above which spacer is considered localized.")
    ("min-leaderless-percent", po::value<double>(&options.minLeaderlessPercent)->default_value(0.0), "Minimum percentage of leaderless transcripts")
    ("min-leaderless-count", po::value<size_t>(&options.minLeaderlessCount)->default_value(0), "Minimum number of leaderless transcripts")
    ("window-size", po::value<size_t>(&options.windowSize)->default_value(1), "Size of window in which to determine localization")
    ("fnlabels", po::value<string>(&options.fnlabels)->default_value(""), "Labels file, if provided, is used to calculate number of leaderless and first-genes-in-operon")
    ("fnseq", po::value<string> (&options.fnseq)->default_value(""), "Sequence file")
    ("min-gene-length", po::value<NumSequence::size_type>(&options.minGeneLength)->default_value(300), "Minimum gene length allowed in training")
    ("match-to", po::value<string>(&options.matchTo)->default_value("TAAGGAGGTGA"), "16S tail")
    ("allow-ag-substitution", po::bool_switch(&options.allowAGSubstitution)->default_value(true), "Allow AG substitution.")
    ("match-thresh", po::value<unsigned>(&options.matchThresh)->default_value(4), "Match threshold for 16S tail.")
    ("fgio-distance-thresh", po::value<size_t>(&options.fgioDistThresh)->default_value(25), "Minimum distance between genes classified as first-genes-in-operon")
    ;
}


void OptionsExperiment::addProcessOptions_PromoterAndRBSMatchOptions(PromoterAndRBSMatch &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of mod file containing RBS spacer.")
    ("match-thresh", po::value<size_t>(&options.numberOfMatches)->default_value(4), "Match threshold for 16S tail.")
    ;
}


void OptionsExperiment::addProcessOptions_StartModelStrategy2Options(StartModelStrategy2Options &options, po::options_description &processOptions) {
    
    Options::addProcessOptions_GenReadSeqAndLabelsOptions(options, processOptions);

    processOptions.add_options()
    ("seq-16s",         po::value<string>(&options.seq16S)->default_value("TAAGGAGGTGA"), "Sequence to match to.")
    ("min-match",       po::value<size_t>(&options.min16SMatch)->default_value(4), "Minimum accepted match length from upstream to 16S tail")
    ("allow-ag-sub",    po::bool_switch(&options.allowAGSubstitution)->default_value(false), "Allow G to be substituted for A when matching to 16S tail")
    ("fgio-thresh",     po::value<size_t>(&options.fgioDistanceThresh)->default_value(25), "Used to extract second genes in operon")
    ("ig-thresh",       po::value<size_t>(&options.igDistanceThresh)->default_value(20), "Used to extract first genes in operon")
    ("fgio-upstrLen",   po::value<size_t>(&options.fgioUpstreamLength)->default_value(40), "Upstream length for first genes in operon")
    ("ig-upstrLen",     po::value<size_t>(&options.igUpstreamLength)->default_value(20), "Upstream length for interior genes")
    ("min-gene-length", po::value<size_t>(&options.minGeneLength)->default_value(0), "Minimum gene length")
    ("match-to-upstream-of-length", po::value<size_t> (&options.matchToUpstreamOfLength)->default_value(20), "Length of upstream region we're matching against")
    ("fn_out",         po::value<string>(&options.fn_out)->required(), "Name of output file")
    ("fgio-matched-upstream-length",    po::value<size_t> (&options.upstreamLengthFGIOMatched)->default_value(40), "FGIO Match Upstream Length")
    ("fgio-unmatched-upstream-length",  po::value<size_t> (&options.upstreamLengthFGIOUnmatched)->default_value(40), "FGIO Unmatch Upstream Length")
    ("ig-matched-upstream-length",      po::value<size_t> (&options.upstreamLengthIGMatched)->default_value(20), "IG Match Upstream Length")
    ("ig-unmatched-upstream-length",    po::value<size_t> (&options.upstreamLengthIGUnmatched)->default_value(20), "IG Unmatch Upstream Length")
    ;
    
    
    // RBS mfinder options
    po::options_description mfinderFGIOMatched ("FGIO Matched - Motif Finder ");
    OptionsMFinder::addProcessOptions(options.mfinderFGIOMatchedOptions, mfinderFGIOMatched, false, "fgio-matched");
    processOptions.add(mfinderFGIOMatched);
    
    po::options_description mfinderFGIOUnmatched ("FGIO Unmatched - Motif Finder ");
    OptionsMFinder::addProcessOptions(options.mfinderFGIOUnmatchedOptions, mfinderFGIOUnmatched, false, "fgio-unmatched");
    processOptions.add(mfinderFGIOUnmatched);
    
    po::options_description mfinderIGMatched ("IG Matched - Motif Finder ");
    OptionsMFinder::addProcessOptions(options.mfinderIGMatchedOptions, mfinderIGMatched, false, "ig-matched");
    processOptions.add(mfinderIGMatched);
    
    po::options_description mfinderIGUnmatched ("IG Unmatched - Motif Finder ");
    OptionsMFinder::addProcessOptions(options.mfinderIGUnmatchedOptions, mfinderIGUnmatched, false, "ig-unmatched");
    processOptions.add(mfinderIGUnmatched);
}


void OptionsExperiment::addProcessOptions_RBSConsensusAnd16SMatch(RBSConsensusAnd16SMatch &options, po::options_description &processOptions){
    processOptions.add_options()
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of mod file containing RBS.")
    ("match-thresh", po::value<unsigned>(&options.matchThresh)->default_value(4), "Match threshold for 16S tail.")
    ("match-to", po::value<string>(&options.matchTo)->default_value("TAAGGAGGTGA"), "16S tail")
    ("allow-ag-sub",    po::bool_switch(&options.allowAGSubstitution)->default_value(true), "Allow G to be substituted for A when matching to 16S tail")
    ;
}




void OptionsExperiment::addProcessOptions_RBSIsLocalized(RBSIsLocalized &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of mod file containing RBS spacer.")
    ("dist-thresh", po::value<size_t>(&options.distanceThresh)->default_value(15), "Distance threshold before which spacer indicates promoter.")
    ("score-thresh", po::value<double>(&options.scoreThresh)->default_value(0.1), "Minimum score above which spacer is considered localized.")
    ("window-size", po::value<size_t>(&options.windowSize)->default_value(1), "Size of window in which to determine localization")
    ;
}

void OptionsExperiment::addProcessOptions_ComputeMotifScoreForStarts(ComputeMotifScoreForStarts &options, po::options_description &processOptions) {
    
    processOptions.add_options()
    ("fnsequence", po::value<string>(&options.fnsequence)->required(), "Name of sequences file")
    ("fnmod", po::value<string>(&options.fnmod)->required(), "Name of mod file containing motif and spacer.")
    ("fnlabels", po::value<string>(&options.fnlabels)->required(), "File containing gene labels.")
    ("upstream-len", po::value<size_t>(&options.upstreamLength)->default_value(20), "Length of upstream region to search in")
    ;
}
