//
//  OptionsPostProcessor.cpp
//  Biogem CPP
//
//  Created by Karl Gemayel on 9/18/18.
//  Copyright © 2018 Georgia Institute of Technology. All rights reserved.
//

#include "OptionsPostProcessor.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>


using namespace std;
using namespace gmsuite;

namespace po = boost::program_options;

//namespace gmsuite {
//    std::istream& operator>>(std::istream& in, GeneticCode::gcode_t& unit)
//    {
//        std::string token;
//        in >> token;
//        if (token == "11")
//            unit = GeneticCode::ELEVEN;
//        else if (token == "4")
//            unit = GeneticCode::FOUR;
//        //    else
//        //        throw boost::program_options::validation_error("Invalid genome class");
//        
//        return in;
//    }
//}

OptionsPostProcessor::OptionsPostProcessor(string mode) : Options(mode) {
    
}

bool OptionsPostProcessor::parse(int argc, const char **argv) {
    
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
        
        if (vm.count("config") > 0) {
            config_fnames = vm["config"].as<std::vector<std::string> >();
            
            
            for (size_t n = 0; n < config_fnames.size(); n++) {
                std::ifstream file(config_fnames[n].c_str());
                
                if(file.fail())
                {
                    std::cerr << "Error opening config file: " << config_fnames[n] << std::endl;
                    return false;
                }
                
                po::store(po::parse_config_file(file, config_file_options), vm);
                file.close();
            }
        }
        
        
        // if help specified, print usage message and quit
        if (vm.count("help")) {
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // try parsing arguments.
        po::notify(vm);
        
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
}


void OptionsPostProcessor::addProcessOptions(OptionsPostProcessor &options, po::options_description&processOptions) {
    
    processOptions.add_options()
    ("sequences", po::value<string> (&options.fnsequence)->required(), "File containing sequences")
    ("labels", po::value<string> (&options.fnlabels)->required(), "File containing labels")
    ("mod", po::value<string> (&options.fnmod)->required(), "file containing models")
    ("window-upstream", po::value<size_t>(&options.windowUpstream)->default_value(60), "Upstream length of window")
    ("genetic-code",        po::value<GeneticCode::gcode_t>      (&options.gcode                              )->default_value(GeneticCode::ELEVEN), "Genetic code")
    ("window-downstream", po::value<size_t> (&options.windowDownstream)->default_value(60), "Downstream length of window")
    ("neigh-downstream", po::value<size_t> (&options.neighborhoodDownstream)->default_value(150), "Downstream neighborhood in which to search for candidate starts")
    ("neigh-upstream", po::value<size_t> (&options.neighborhoodUpstream)->default_value(150), "Upstream neighborhood in which to search for candidate starts")
    ;
}
