//
//  OptionsUtilities.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsUtilities_hpp
#define OptionsUtilities_hpp

#include <stdio.h>
#include <string>

#include "Options.hpp"
#include "OptionsGMS2Training.hpp"

using std::string;

namespace gmsuite {
    
    /**
     * @class OptionsUtilities
     * @brief A class that deals with parsing command-line options for Utilities module
     */
    class OptionsUtilities : public Options {
        
        
    public:
        
        typedef enum {
            EXTRACT_UPSTR,
            START_MODEL_INFO,
            MATCH_SEQ_TO_UPSTREAM,
            MATCH_SEQ_TO_NONCODING,
            LABELS_SIMILARITY_CHECK,
            EMIT_NON_CODING,
            COUNT_NUM_ORF,
            EXTRACT_SC_PER_OPERON_STATUS,
            EXTRACT_SC_PER_MOTIF_STATUS
        }
        utility_t;
        
        OptionsUtilities(string mode);
        
        /**
         * Parse the command-line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        bool parse(int argc, const char *argv[]);
        
        
        // Below, create a variable for each parameter, to make for easy access
    public:
        
        utility_t utility;             // the chosen utility
        
        // below you'll find utility-based options
        
        // generic options
        struct GenericOptions {
            bool debug;
        };
        
        // extract-upstream options
        struct ExtractUpstreamUtility : public GenericOptions {
            string fn_sequence;             // sequence filename
            string fn_label;                // label filename
            string fn_output;               // output filename
            size_t length;                  // length of upstream regions
            bool allowOverlaps;             // allow upstream region to overlap coding region
            size_t minimumGeneLength;       // minimum gene length associated with upstream
        }
        extractUpstreamUtility;
        
        
        // start-model-info options
        struct StartModelInfoUtility : public GenericOptions {
            string fn_sequence;             // sequence filename
            string fn_label;                // label filename
            bool allowOverlaps;             // allow upstream region to overlap coding region
            size_t numOfSimNonCoding;       // number of simulated non-coding sequences
            OptionsGMS2Training optionsGMS2Training;            // options for running gms2 training
        }
        startModelInfoUtility;
        
        
        // match-seq-with-upstream
        struct MatchSeqWithUpstream : public ExtractUpstreamUtility {
            string matchTo;                 // the sequence to be matched
        }
        matchSeqWithUpstream;
        
        // match-seq-to-noncoding
        struct MatchSeqWithNoncoding : public StartModelInfoUtility {
            string matchTo;                 // the sequence to be matched
        }
        matchSeqWithNoncoding;
        
        struct LabelsSimilarityCheck : public GenericOptions {
            string fn_labelsA;              // label file A
            string fn_labelsB;              // label file B
        } labelsSimilarityCheck;
        
        struct EmitNonCoding : public GenericOptions {
            string fn_mod;                  // input model file containing noncoding model
            string fn_out;                  // sequence output file
            NumSequence::size_type length;  // length of non-coding sequence
        } emitNonCoding;
        
        struct CountNumORF : public GenericOptions {
            string fn_sequence;             // sequence file
            string fn_mod;                  // input model file containing genetic code and other params
            bool printSeq;                  // print sequences instead of number
        } countNumORF;
        
        struct ExtractStartContextPerOperonStatus : public GenericOptions {
            string fn_sequence;             // sequence filename
            string fn_label;                // label filename
            string fn_output;               // output filename
            size_t length;                  // length of upstream regions
            bool allowOverlaps;             // allow upstream region to overlap coding region
            size_t minimumGeneLength;       // minimum gene length associated with upstream
            size_t distThreshFGIO;
            size_t distThreshIG;
        } extractStartContextPerOperonStatus;
        
        struct ExtractStartContextPerMotifStatus : public GenericOptions {
            string fn_sequence;             // sequence filename
            string fn_label;                // label filename
            string fn_output;               // output filename
            size_t length;                  // length of upstream regions
            bool allowOverlaps;             // allow upstream region to overlap coding region
            size_t minimumGeneLength;       // minimum gene length associated with upstream
            size_t distThreshFGIO;
            size_t distThreshIG;
            string matchTo;
            size_t matchThresh;
            bool allowAGSubstitution;
        } extractStartContextPerMotifStatus;
        
        static void addProcessOptions_ExtractUpstream(ExtractUpstreamUtility &options, po::options_description &processOptions);
        static void addProcessOptions_StartModelInfo(StartModelInfoUtility &options, po::options_description &processOptions);
        static void addProcessOptions_LabelsSimilarityCheck(LabelsSimilarityCheck &options, po::options_description &processOptions);
        static void addProcessOptions_EmitNonCoding(EmitNonCoding &options, po::options_description &processOptions);
        static void addProcessOptions_CountNumORF(CountNumORF &options, po::options_description &processOptions);
        static void addProcessOptions_ExtractStartContextPerOperonStatus(ExtractStartContextPerOperonStatus &options, po::options_description &processOptions);
        static void addProcessOptions_ExtractStartContextPerMotifStatus(ExtractStartContextPerMotifStatus &options, po::options_description &processOptions);
        
    };
}


#endif /* OptionsUtilities_hpp */
