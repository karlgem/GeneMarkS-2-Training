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
        
        OptionsGMS2Training(string mode);
        
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
        
        typedef ProkGeneStartModel::genome_class_t genome_class_t;
        
        string fn_sequence;             /**< Input filename containing DNA sequence */
        string fn_labels;               /**< Input filename containing labels */
        string fn_outmod;               /**< Output model file */
        size_t upstrLength;             /**< Length of upstream sequences for motif search */
        genome_class_t genomeClass;     /**< The genome's class */
        
        
        // GMS2 model parameters
        double pcounts;                             /**< Pseudocounts */
        unsigned codingOrder;                       /**< Coding model's Markov order */
        unsigned noncodingOrder;                    /**< Noncoding model's Markov order */
        unsigned startContextOrder;                 /**< Start-context Markov order */
        NumSequence::size_type upstreamLength;      /**< Upstream length used for motif search */
        NumSequence::size_type startContextLength;  /**< length of start context model */
        
        OptionsMFinder optionsMFinder;              /**< Options for Motif Finder */
        
    };
}


#endif /* OptionsGMS2Training_hpp */
