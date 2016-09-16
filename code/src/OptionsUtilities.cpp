//
//  OptionsUtilities.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsUtilities.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

using namespace std;
using namespace gmsuite;
namespace po = boost::program_options;

OptionsUtilities::OptionsUtilities(string mode) : Options(mode) {
    
}

// parse CMD options
bool OptionsUtilities::parse(int argc, const char *argv[]) {
    
    
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
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
        ("mode", po::value<string>(&mode)->required(), "Program Mode")
        ("utility", po::value<string>(&utility)->required(), "Utility function")
        ("subargs", po::value<std::vector<std::string> >(), "Arguments for command")
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
        pos.add("utility",1);
        pos.add("subargs",-1);
        
        
        
        // create storage component for storing names and values of arguments
        po::variables_map vm;
        
        // store command-line options
        po::parsed_options parsed = po::command_line_parser(argc, argv).          // pass in input
            options(cmdline_options).                     // specify options list
            positional(pos).                              // specify which are positional
            allow_unregistered().                         // allow unregistered args since we haven't yet set-up the utility-based ones
            run();
        
        po::store(parsed,vm);
        
        // try parsing arguments.
        po::notify(vm);
        
        // Utility: Extract upstream
        if (utility == EXTRACT_UPSTR) {
            
            po::options_description utilDesc(string(EXTRACT_UPSTR) + " options");
            utilDesc.add_options()
                ("sequence,s", po::value<string>(&extractUpstreamUtility.fn_sequence)->required(), "Sequence filename")
                ("label,l", po::value<string>(&extractUpstreamUtility.fn_label)->required(), "Label filename")
                ("output,o", po::value<string>(&extractUpstreamUtility.fn_output)->required(), "Output filename")
                ("length", po::value<size_t>(&extractUpstreamUtility.length)->required(), "Upstream length")
                ("allow-overlap-with-cds", "If set, then upstream (non-coding) regions are allowed to overlap with coding regions. If not set, these sequences are ignored.")
                ("min-gene-length", po::value<size_t>(&extractUpstreamUtility.minimumGeneLength)->default_value(0), "Minimum gene length")
            ;
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            extractUpstreamUtility.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
        }
        else                                                                        // unrecognized utility
            throw po::invalid_option_value(utility);
        
        // if help specified, print usage message and quit
        if (vm.count("help")) {
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // update all values and make sure required are provided
        po::notify(vm);
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
    
    
}
