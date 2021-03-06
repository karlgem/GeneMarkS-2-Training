//
//  OptionsMFinder.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsMFinder.hpp"

#include <vector>
#include <limits.h>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

using namespace std;
using namespace gmsuite;
namespace po = boost::program_options;


namespace gmsuite {
    // convert string to genome_group_t
    std::istream& operator>>(std::istream& in, OptionsMFinder::align_t& unit) {
        string token;
        in >> token;
        
        if (token == "none")
            unit = MFinderModelParams::NONE;
        else if (token == "left")
            unit = MFinderModelParams::LEFT;
        else if (token == "right")
            unit = MFinderModelParams::RIGHT;
        else
            throw po::invalid_option_value(token);
        
        return in;
        
    }
}


OptionsMFinder::OptionsMFinder(string mode) : Options(mode) {
    
}


// parse CMD options
bool OptionsMFinder::parse(int argc, const char *argv[]) {
    
    try {
        vector<string> config_fnames;                           // holds names of all configuration files (specified by user)
        
        // Declare a group of options that will be allowed
        // only on the command line (CML)
        po::options_description generic("General Options");
        generic.add_options()
        ("version", "Print version string")
        ("help,h", "Display help message")
        ("config", po::value(&config_fnames), "Config file where options may be specified (can be specified more than once)")
        ;
        
        // Declare a group of options that will be allowed
        // on both CML and in the config files
        po::options_description config("Configuration");
        config.add_options()
            ("verbose,v", po::value<int>(&verbose)->default_value(0), "Verbose level")
        ;
        
        addProcessOptions(*this, config);
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
        ("mode", po::value<string>(&mode)->required(), "Program Mode")
        ("fname", po::value<string>(&fname_in)->required(), "Name of sequence file");
        ;
        
        // Congregate options into further groups
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);
        
        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);
        
        po::options_description visible("Allowed options");
        visible.add(generic).add(config);
        
        // Specify positional arguments
        po::positional_options_description pos;
        pos.add("mode",1);
        pos.add("fname",-1);
        
        
        
        // create storage component for storing names and values of arguments
        po::variables_map vm;
        
        // store command-line options
        po::store(po::command_line_parser(argc, argv).          // pass in input
                  options(cmdline_options).                     // specify options list
                  positional(pos).                              // specify which are positional
                  run(),                                        // parse options
                  vm);                                          // specify storage container
        
        // if help specified, print usage message and quit
        if (vm.count("help")) {
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // try parsing arguments.
        po::notify(vm);
        
        // if config file specified, read it and get values
        if (vm.count("config") > 0) {
            config_fnames = vm["config"].as<vector<string> >();
            
            // for each config filename, make sure file exists. Then try parsing it
            for (size_t i = 0; i < config_fnames.size(); i++) {
                
                // open it
                ifstream ifs(config_fnames[i].c_str());
                
                // if can't open file, print error and return false
                if (ifs.fail()) {
                    cerr << "Error opening config file: " << config_fnames[i] << endl;
                    return false;
                }
                
                // Otherwise, parse config file and store the parameter values.
                // NOTE: store doesn't change the value of an option if it is already assigned.
                //       Therefore, values specified in the command line take precedence over those
                //       in the config file, as expected.
                po::store(po::parse_config_file(ifs, config_file_options), vm);
                
                // try parsing arguments. Exception is returned if require option not provided.
                po::notify(vm);
            }
        }
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    
    return true;
    
    
}



void OptionsMFinder::addProcessOptions(OptionsMFinder &optionsMFinder, po::options_description &processOptions, bool allowSingleLetter, string prefix) {
    
    if (!prefix.empty())
        prefix += "-";
    
    string opt_width        = prefix + "width"          + (allowSingleLetter ? ",w" : "");
    string opt_motifOrder   = prefix + "motif-order"    + (allowSingleLetter ? ",o" : "");
    string opt_bkgdOrder    = prefix + "bkgd-order"     + (allowSingleLetter ? ",b" : "");
    string opt_align        = prefix + "align"          ;
    string opt_tries        = prefix + "tries"          ;
    string opt_maxIter      = prefix + "max-iter"       ;
    string opt_maxEMIter    = prefix + "max-em-iter"    ;
    string opt_shiftEvery   = prefix + "shift-every"    ;
    string opt_pcount       = prefix + "pcount"         ;
    string opt_filterThresh = prefix + "filter-thresh"  ;
    
    
    processOptions.add_options()
    (opt_width.c_str(),         po::value<unsigned> (&optionsMFinder.width)->default_value(6), "Width of motif")
    (opt_motifOrder.c_str(),    po::value<unsigned> (&optionsMFinder.motifOrder)->default_value(0), "Order of the motif's Markov model")
    (opt_bkgdOrder.c_str(),     po::value<unsigned> (&optionsMFinder.bkgdOrder)->default_value(0), "Order of the background's Markov model")
    (opt_align.c_str(),         po::value<align_t>  (&optionsMFinder.align)->default_value(MFinderModelParams::NONE), "If set, positional information is considered by the model")
    (opt_tries.c_str(),         po::value<unsigned> (&optionsMFinder.tries)->default_value(10), "Number of restarts")
    (opt_maxIter.c_str(),       po::value<unsigned> (&optionsMFinder.maxIter)->default_value(60), "Number of Gibbs iterations per single try")
    (opt_maxEMIter.c_str(),     po::value<unsigned> (&optionsMFinder.maxEMIter)->default_value(10), "Number of EM iterations per single try")
    (opt_shiftEvery.c_str(),    po::value<unsigned> (&optionsMFinder.shiftEvery)->default_value(10), "Number of iterations before shifting motif")
    (opt_pcount.c_str(),        po::value<double>   (&optionsMFinder.pcounts)->default_value(1), "Pseudocounts")
    (opt_filterThresh.c_str(),  po::value<double>   (&optionsMFinder.filterThresh)->default_value(-std::numeric_limits<double>::infinity()), "Value for filtering out motifs with low score")
    ;
}

























