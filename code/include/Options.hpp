//
//  Options.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/11/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Options_hpp
#define Options_hpp

#include <string>
#include <vector>

#include <boost/program_options.hpp>

using std::string;
using std::vector;

namespace po = boost::program_options;

namespace gmsuite {
    
    /**
     * @class Options
     * @brief A class that deals with parsing command-line options
     */
    class Options {
        
    public:
        
        /**
         * Parse the command-line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        virtual bool parse(int argc, const char *argv[]) = 0;
        
        
        
        
        
    protected:
        
        
        /**
         * Create three groups of parameters: (1) all, (2) hidden, (3) positional
         */
        
        /**
         * Retrieve the program name from the path (i.e. the word after the 
         * final '/' character). E.g.
         *
         * If the input is "/home/prog", the output will be "prog"
         *
         * @param path a path ending with the program name
         * @return the program name
         */
        string basename(const string &path) const;
        
        /**
         * Get string representations of the positional arguments, to be used
         * in the usage message.
         * 
         * @param pod contains information on positional arguments (e.g. names)
         * @return a vector of strings containing positional argument names.
         */
        vector<string> get_positional_args(const po::positional_options_description &pod) const;
        
        /**
         * Get a string representation of the usage message.
         *
         * @param program_name the program's (base)name
         * @param desc contains information on arguments and their values
         * @param pod contains information on positional arguments
         *
         * @return a usage message.
         */
        string make_usage_string(const string &program_name,
                                 const po::options_description& desc,
                                 const po::positional_options_description& pod) const;
        
    // Below, create a variable for each parameter, to make for easy access
    public:
        
        // basic config
        int verbose;
        
        // program parameters
        string mode;                    // parameters for which mode? e.g. GMS2
        
        
        // developer options
        
        
    };
}

#endif /* Options_hpp */
