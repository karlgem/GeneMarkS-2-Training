//
//  OptionsMFinder.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsMFinder_hpp
#define OptionsMFinder_hpp

#include <stdio.h>
#include "Options.hpp"

namespace gmsuite {
    
    /**
     * @class OptionsMFinder
     * @brief A class that deals with parsing command-line options for MFinder module
     */
    class OptionsMFinder : public Options {
        
        
    public:
        
        OptionsMFinder(string mode);
        
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
        
        string fname_in;                /**< Input filename containing DNA sequence */
        double pcounts;                 /**< The pseudocount value */
    };
    
}

#endif /* OptionsMFinder_hpp */
