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
            BUILD_START_MODELS3
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
            bool allowAGSubstitution;           // whether A and G can be subsituted while matching
            OptionsMFinder mfinderRBSOptions;   // RBS
            OptionsMFinder mfinderPromoterOptions;  // Promoter
        }
        scoreStarts;
        
        
        
        /**********************************************\
         *              Option Processing             *
        \**********************************************/
        
        static void addProcessOptions_MatchSeqToUpstreamOptions(MatchSeqToUpstreamOptions &options, po::options_description &processOptions);
        static void addProcessOptions_MatchSeqToNoncodingOptions(MatchSeqToNoncodingOptions &options, po::options_description &processOptions);
        static void addProcessOptions_BuildStartModelsOptions(BuildStartModelsOptions &options, po::options_description &processOptions);
        static void addProcessOptions_BuildStartModels2Options(BuildStartModels2Options &options, po::options_description &processOptions);
        static void addProcessOptions_BuildStartModels3Options(BuildStartModels3Options &options, po::options_description &processOptions);
        
        static void addProcessOptions_ScoreStarts(ScoreStarts &options, po::options_description &processOptions);
        
        
        
        
        
    };
}

#endif /* OptionsExperiment_hpp */
