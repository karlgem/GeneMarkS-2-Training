//
//  OptionsMFinder.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsMFinder_hpp
#define OptionsMFinder_hpp

#include "Options.hpp"
#include "MFinderModelParams.hpp"

#include <stdio.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace gmsuite {
    
    /**
     * @class OptionsMFinder
     * @brief A class that deals with parsing command-line options for MFinder module
     */
    class OptionsMFinder : public Options {
        
        
    public:
        
        OptionsMFinder(string mode = "mfinder");
        
        /**
         * Parse the command-line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        bool parse(int argc, const char *argv[]);
        
        
        static void addProcessOptions(OptionsMFinder &optionsMFinder, po::options_description &processOptions, bool allowSingleLetter = true, string prefix="");
        
        
        // Below, create a variable for each parameter, to make for easy access
    public:
        
        typedef MFinderModelParams::align_t align_t;
        
        string fname_in;                /**< Input filename containing DNA sequence */
        double pcounts;                 /**< The pseudocount value */
        unsigned width;                 /**< The width of the motif */
        unsigned motifOrder;            /**< Order for motif Markov model */
        unsigned bkgdOrder;             /**< Order of background Markov model */
        align_t align;                  /**< Set to left, right, or none, to indicate whether positional distribution should be used, and from which direction */
        unsigned tries;                 /**< Number of restarts; i.e. number of times the algorithm is executed with new random initialization */
        unsigned maxIter;               /**< Number of Gibbs iterations per single try */
        unsigned maxEMIter;             /**< Number of EM iterations per single try */
        unsigned shiftEvery;            /**< Number of iterations before attempting to shift the motif left and right */
        double filterThresh;            /**< Threshold used to filter unwanted */
        
        
    };
    
}

#endif /* OptionsMFinder_hpp */
