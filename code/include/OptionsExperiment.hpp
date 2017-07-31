//
//  OptionsExperiment.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/12/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsExperiment_hpp
#define OptionsExperiment_hpp

#include <stdio.h>

#include "Options.hpp"
#include "OptionsMFinder.hpp"

namespace gmsuite {
    
    /**
     * @class OptionsExperiment
     * @brief A class that deals with parsing command-line options for Experiment module
     */
    class OptionsExperiment : public Options {
        
    public:
        
        // experiment types
        typedef enum {
            MATCH_SEQ_TO_UPSTREAM,
            MATCH_SEQ_TO_NONCODING,
            BUILD_START_MODELS,
            BUILD_START_MODELS2,
            BUILD_START_MODELS3,
            SCORE_STARTS,
            MATCH_RBS_TO_16S,
            SCORE_LABELED_STARTS,
            PROMOTER_IS_VALID_FOR_ARCHAEA,
            PROMOTER_IS_VALID_FOR_BACTERIA,
            START_MODEL_STRATEGY_2,
            PROMOTER_AND_RBS_MATCH,
            RBS_CONSENSUS_AND_16S_MATCH
        }
        experiment_t;
        
        /**
         * Constructor: initialize an experiment
         *
         * @param mode the options mode (placeholder name)
         */
        OptionsExperiment(string mode = "experiment");
        
        /**
         * Parse the command line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        bool parse(int argc, const char *argv[]);
        
        experiment_t experiment;        // the chosen experiment
        
        
        
        
        
        /**********************************************\
         *          Experiment-based options          *
        \**********************************************/
        
        GenericOptions genericOptions;
        
        // match-seq-to-upstream options
        struct MatchSeqToUpstreamOptions : public GenExtractUpstreamsOptions {
            string matchTo;                 // the sequence to be matched
            
        }
        matchSeqToUpstream;
        
        
        // match-seq-to-noncoding options
        struct MatchSeqToNoncodingOptions : public GenReadSeqAndLabelsOptions {
            string matchTo;                 // the sequence to be matched
            size_t length;                  // length of noncoding sequences
            size_t numNoncoding;            // number of simulated non-coding
            unsigned order;                 // order of noncoding model
            double pcounts;                 // pseudocounts
        }
        matchSeqToNoncoding;
        
        // build-start-models options
        struct BuildStartModelsOptions : public GenExtractUpstreamsOptions {
            string matchTo;                 // the sequence to be matched
            unsigned min16SMatch;           // the minimum accepted match length with 16S tail
            bool allowAGSubstitution;       // whether A and G can be substituted while matching
            OptionsMFinder mfinderOptions;  // options for motif finder
        }
        buildStartModels;
        
        
        // build-start-models options
        struct BuildStartModels2Options : public GenExtractUpstreamsOptions {
            string matchTo;                 // the sequence to be matched
            unsigned min16SMatch;           // the minimum accepted match length with 16S tail
            bool allowAGSubstitution;       // whether A and G can be substituted while matching
            OptionsMFinder mfinderOptions;  // options for motif finder
        }
        buildStartModels2;
        
        // build-start-models3 options
        struct BuildStartModels3Options : public GenExtractUpstreamsOptions {
            string matchTo;                 // the sequence to be matched
            unsigned min16SMatch;           // the minimum accepted match length with 16S tail
            size_t nfgioThresh;             // non-first-gene-in-operon threshold
            bool allowAGSubstitution;       // whether A and G can be substituted while matching
            OptionsMFinder mfinderOptions;  // options for motif finder
        }
        buildStartModels3;
        
        // score-starts options
        struct ScoreStarts : public GenExtractUpstreamsOptions {
            string matchTo;                     // sequence to be matched
            unsigned min16SMatch;               // the minimum accepted match length with 16S tail
            size_t nfgioThresh;                 // non-first-gene-in-operon threshold
            size_t fgioThresh;                  // first-gene-in-operon threshold
            size_t searchUpstrLen;              // the upstream length in which to search for 16S match
            bool allowAGSubstitution;           // whether A and G can be subsituted while matching
            size_t upstreamLenRBS;              // upstream length for RBS search
            size_t upstreamLenPromoter;         // upstream length for promoter search
            string gms2mod;                     // gms2 mod file
            OptionsMFinder mfinderOptions;      // RBS
        }
        scoreStarts;
        
        
        // match RBS to 16S
        struct MatchRBSTo16S : public GenericOptions {
            string fnlabels;                    // file containing labels and RBS sequence
            string matchTo;                     // 16S sequence to be matched to
            unsigned min16SMatch;               // the minimum accepted match length with 16S tail
            bool allowAGSubstitution;       // whether A and G can be substituted while matching
        }
        matchRBSTo16S;
        
        struct ScoreLabeledStarts : public GenExtractUpstreamsOptions {
            string fnmod;
        }
        scoreLabeledStarts;
        
        // get-start-model-type
        struct PromoterIsValidForArchaea : public GenericOptions {
            string fnmod;
            size_t distanceThresh;
            double scoreThresh;
            size_t windowSize;
        }
        promoterIsValidForArchaea;
        
        // get-start-model-type
        struct PromoterIsValidForBacteria : public GenericOptions {
            string fnmod;
            size_t distanceThresh;
            double scoreThresh;
            double minLeaderlessPercent;
            size_t minLeaderlessCount;
            size_t windowSize;
            
            string fnlabels;
            string fnseq;
            size_t fgioDistThresh;
            size_t minGeneLength;
            string matchTo;
            bool allowAGSubstitution;
            unsigned matchThresh;
        }
        promoterIsValidForBacteria;
        
        // promoter/rbs match
        struct PromoterAndRBSMatch : public GenericOptions {
            string fnmod;
            size_t numberOfMatches;
        }
        promoterAndRBSMatch;
        
        struct RBSConsensusAnd16SMatch : public GenericOptions {
            string fnmod;
            string matchTo;             // 16S
            bool allowAGSubstitution;
            unsigned matchThresh;
        }
        rbsConsensusAnd16SMatch;
        
        
        struct StartModelStrategy2Options : public GenReadSeqAndLabelsOptions {
            string seq16S;                      // 16S rRNA tail match
            size_t min16SMatch;               // the minimum accepted match length with 16S tail
            bool allowAGSubstitution;           // whether A and G can be subsituted while matching
            size_t fgioDistanceThresh;          // FGIO distance threshold
            size_t igDistanceThresh;            // interior genes distance treshold
            size_t fgioUpstreamLength;          // upstream length for fgio
            size_t igUpstreamLength;            // upstream length for ig
            size_t minGeneLength;               // minimu gene length
            size_t matchToUpstreamOfLength;     // length of the upstream region we're matching against
            string fn_out;                      // name of output file
            
            size_t upstreamLengthFGIOMatched;
            size_t upstreamLengthFGIOUnmatched;
            size_t upstreamLengthIGMatched;
            size_t upstreamLengthIGUnmatched;
            
            
            OptionsMFinder mfinderFGIOMatchedOptions;
            OptionsMFinder mfinderFGIOUnmatchedOptions;
            OptionsMFinder mfinderIGMatchedOptions;
            OptionsMFinder mfinderIGUnmatchedOptions;
        }
        startModelStrategy2;
        
        
        /**********************************************\
         *              Option Processing             *
        \**********************************************/
        
        static void addProcessOptions_MatchSeqToUpstreamOptions(MatchSeqToUpstreamOptions &options, po::options_description &processOptions);
        static void addProcessOptions_MatchSeqToNoncodingOptions(MatchSeqToNoncodingOptions &options, po::options_description &processOptions);
        static void addProcessOptions_BuildStartModelsOptions(BuildStartModelsOptions &options, po::options_description &processOptions);
        static void addProcessOptions_BuildStartModels2Options(BuildStartModels2Options &options, po::options_description &processOptions);
        static void addProcessOptions_BuildStartModels3Options(BuildStartModels3Options &options, po::options_description &processOptions);
        static void addProcessOptions_MatchRBSTo16SOptions(MatchRBSTo16S &options, po::options_description &processOptions);
        
        static void addProcessOptions_ScoreStarts(ScoreStarts &options, po::options_description &processOptions);
        static void addProcessOptions_ScoreLabeledStarts(ScoreLabeledStarts &options, po::options_description &processOptions);
        static void addProcessOptions_PromoterIsValidForArchaea(PromoterIsValidForArchaea &options, po::options_description &processOptions);
        static void addProcessOptions_PromoterIsValidForBacteria(PromoterIsValidForBacteria &options, po::options_description &processOptions);
        
        static void addProcessOptions_StartModelStrategy2Options(StartModelStrategy2Options &options, po::options_description &processOptions);
        static void addProcessOptions_PromoterAndRBSMatchOptions(PromoterAndRBSMatch &options, po::options_description &processOptions);
        static void addProcessOptions_RBSConsensusAnd16SMatch(RBSConsensusAnd16SMatch &options, po::options_description &processOptions);
        
        
    };
}

#endif /* OptionsExperiment_hpp */
