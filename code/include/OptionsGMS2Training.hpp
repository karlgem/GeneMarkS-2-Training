//
//  OptionsGMS2Training.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsGMS2Training_hpp
#define OptionsGMS2Training_hpp

#include <stdio.h>
#include <string>

#include "Options.hpp"
#include "NumSequence.hpp"
#include "GeneticCode.hpp"
#include "OptionsMFinder.hpp"
#include "ProkGeneStartModel.hpp"

using std::string;

namespace gmsuite {
    
    /**
     * @class OptionsGMS2Training
     * @brief A class that deals with parsing command-line options for GMS2Training module
     */
    class OptionsGMS2Training : public Options {
        
        
    public:
        
        OptionsGMS2Training(string mode="gms2-training");
        
        /**
         * Parse the command-line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        bool parse(int argc, const char *argv[]);
        
        
        static void addProcessOptions(OptionsGMS2Training &options, po::options_description &processOptions);
        
        // Below, create a variable for each parameter, to make for easy access
    public:
        
        typedef ProkGeneStartModel::genome_class_t genome_class_t;
        typedef GeneticCode::gcode_t gcode_t;
        
        string fn_sequence;             /**< Input filename containing DNA sequence */
        string fn_labels;               /**< Input filename containing labels */
        string fn_outmod;               /**< Output model file */
        size_t upstrLength;             /**< Length of upstream sequences for motif search */
        genome_class_t genomeClass;     /**< The genome's class */
        gcode_t geneticCode;            /**< Genetic code */
        string fn_settings;             /**< Settings to place in output mod file */
        bool runMotifSearch;            /**< Run motif search */
        size_t upstrLenFGIO;            /**< Upstream length for first-genes-in-operon */ 
        unsigned widthArchaeaPromoter;
        // prediction parameters
        double nonProbN;
        double codProbN;
        double nonDurationDecay;
        double codDurationDecay;
        NumSequence::size_type geneMinLengthPrediction;
        string matchTo;
        bool allowAGSubstitution;
        unsigned matchThresh;
        NumSequence::size_type upstreamSignatureLength;
        unsigned upstreamSignatureOrder;
        bool trainNonCodingOnFullGenome;            /**< If set, non-coding model is trained on full genome (instead of non-coding regions) */
        unsigned fgioDistThresh;
        bool cutPromTrainSeqs;
        
        // GMS2 model parameters
        double pcounts;                             /**< Pseudocounts */
        unsigned codingOrder;                       /**< Coding model's Markov order */
        unsigned noncodingOrder;                    /**< Noncoding model's Markov order */
        unsigned startContextOrder;                 /**< Start-context Markov order */
        int startContextMargin;                     /**< Start-context margin (distance from start) */
        NumSequence::size_type upstreamLength;      /**< Upstream length used for motif search */
        NumSequence::size_type startContextLength;  /**< length of start context model */
        NumSequence::size_type MIN_GENE_LEN;        /**< minimum gene length required for using gene in training */ 
        bool trainOnNativeOnly;                     /**< train on genes predicted by native model */
        OptionsMFinder optionsMFinder;              /**< Options for Motif Finder */
        
    };
}


#endif /* OptionsGMS2Training_hpp */
