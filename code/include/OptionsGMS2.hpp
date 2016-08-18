//
//  OptionsGMS2.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsGMS2_hpp
#define OptionsGMS2_hpp

#include <stdio.h>
#include <string>
#include <boost/program_options.hpp>

#include "Options.hpp"

using std::string;

namespace po = boost::program_options;

namespace gmsuite {
    
    /**
     * @class OptionsGMS2
     * @brief A class that deals with parsing command-line options for GMS2 module
     */
    class OptionsGMS2 : public Options {
        
        
    public:
        
        OptionsGMS2(string mode) { this->mode = mode; }
        
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
    };
}

#endif /* OptionsGMS2_hpp */
