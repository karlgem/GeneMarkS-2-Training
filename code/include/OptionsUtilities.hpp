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

using std::string;

namespace gmsuite {
    
    /**
     * @class OptionsUtilities
     * @brief A class that deals with parsing command-line options for Utilities module
     */
    class OptionsUtilities : public Options {
        
        
    public:
        
        #define EXTRACT_UPSTR "extract_upstream"
        
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
        
        string utility;             // the chosen utility
        
        // below you'll fine utility-based options
        
        // generic options
        struct GenericOptions {
            bool debug;
        };
        
        // extract-upstream options
        struct ExtractUpstreamUtility : public GenericOptions {
            size_t length;                  // length of upstream regions
            bool allowOverlaps;             // allow upstream region to overlap coding region
        }
        extractUpstreamUtility;
        
    };
}


#endif /* OptionsUtilities_hpp */
