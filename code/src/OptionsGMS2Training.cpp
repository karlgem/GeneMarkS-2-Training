//
//  OptionsGMS2Training.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsGMS2Training.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

using namespace std;
using namespace gmsuite;
namespace po = boost::program_options;

OptionsGMS2Training::OptionsGMS2Training(string mode) : Options(mode), optionsMFinder(mode) {
    
}

// parse CMD options
bool OptionsGMS2Training::parse(int argc, const char *argv[]) {
    
    
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
        ("fn-sequence,s", po::value<string>(&fn_sequence)->required(), "Name of sequence file")
        ("fn-labels,l", po::value<string>(&fn_labels)->required(), "Name of labels file")
        ("fn-mod,m", po::value<string>(&fn_outmod)->required(), "Name of output model file")
        ("genome-class", po::value<int>()->required(), "The genome's class: 1,2,3")
        
        // GMS2 model parameters
        ("pcounts", po::value<double>(&pcounts)->default_value(1), "Pseudocounts for gms2 models")
        ("coding-order", po::value<unsigned>(&codingOrder)->default_value(4), "Order for coding Markov model")
        ("noncoding-order", po::value<unsigned>(&noncodingOrder)->default_value(2), "Order for noncoding Markov model")
        ("sc-order", po::value<unsigned>(&startContextOrder)->default_value(0), "Order for start-context model")
        ("sc-length", po::value<NumSequence::size_type>(&startContextLength)->default_value(12), "Length of start-context model")
        ("upstream-length", po::value<NumSequence::size_type>(&upstreamLength)->default_value(40), "Length of upstream region for motif search")
        
//        // MFinder options
//        ("pcounts-mfinder", po::value<double>(&optionsMFinder.pcounts)->default_value(1), "Pseudocounts for mfinder models")
//        ("width", po::value<unsigned>(&optionsMFinder.width)->default_value(6), "Width of motif in MFinder")
//        ("motif-order", po::value<unsigned>(&optionsMFinder.motifOrder)->default_value(0), "Order for motif Markov model")
//        ("bkgd-order", po::value<unsigned>(&optionsMFinder.bkgdOrder)->default_value(0), "Order for background model in MFinder")
//        ("align", po::value<string>(&optionsMFinder.align)->default_value("none"), "Set to left or right to allow use of spacer distribution in mfinder")
//        ("tries", po::value<unsigned>(&optionsMFinder.tries)->default_value(10), "Number of tries in mfinder")
        ;
        
        // mfinder options
        po::options_description mfinder("Motif Finder");
        OptionsMFinder::addProcessOptions(optionsMFinder, mfinder);
        config.add(mfinder);
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
        ("mode", po::value<string>(&mode)->required(), "Program Mode")
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
        pos.add("mode",-1);
        
        
        
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
        
        // get genome sequence
        switch (vm["genome-class"].as<int>()) {
            case 1:
                genomeClass = ProkGeneStartModel::C1;
                break;
            case 2:
                genomeClass = ProkGeneStartModel::C2;
            case 3:
                genomeClass = ProkGeneStartModel::C3;
                
            default:
                throw invalid_argument("Unknown genome class");
        }
        
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
    
    
}
