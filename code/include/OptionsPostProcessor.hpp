//
//  OptionsPostProcessor.hpp
//  Biogem CPP
//
//  Created by Karl Gemayel on 9/18/18.
//  Copyright © 2018 Georgia Institute of Technology. All rights reserved.
//

#ifndef OptionsPostProcessor_hpp
#define OptionsPostProcessor_hpp

#include <stdio.h>
#include <string>

#include "Options.hpp"
#include "GeneticCode.hpp"

using std::string;

namespace gmsuite {
    
    /**
     * @class OptionsPostProcessor
     * @brief A class that deals with parsing command-line options for the GMS2 post processor module.
     */
    class OptionsPostProcessor : public Options {
        
    public:
        
        OptionsPostProcessor(string mode="post-processor");
        
        /**
         * Parse the command-line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        bool parse(int argc, const char *argv[]);
        
        static void addProcessOptions(OptionsPostProcessor &options, po::options_description&processOptions);
        
    
        string fnsequence;
        string fnlabels;
        string fnmod;
        GeneticCode::gcode_t gcode;
        
        size_t windowUpstream;
        size_t windowDownstream;
        
        size_t neighborhoodUpstream;
        size_t neighborhoodDownstream;
        
        bool printWindow;
        
    };
}

#endif /* OptionsPostProcessor_hpp */
