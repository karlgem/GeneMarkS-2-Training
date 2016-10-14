//
//  OptionsExperiment.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/12/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsExperiment.hpp"
#include <iostream>

using namespace std;
using namespace gmsuite;

#define STR_MATCH_SEQ_TO_UPSTREAM  "match-seq-to-upstream"
#define STR_MATCH_SEQ_TO_NONCODING "match-seq-to-noncoding"

namespace gmsuite {
    // convert string to experiment_t
    istream& operator>>(istream& in, OptionsExperiment::experiment_t& unit) {
        string token;
        in >> token;
        
        if      (token == STR_MATCH_SEQ_TO_UPSTREAM)    unit = OptionsExperiment::MATCH_SEQ_TO_UPSTREAM;
        else if (token == STR_MATCH_SEQ_TO_NONCODING)   unit = OptionsExperiment::MATCH_SEQ_TO_NONCODING;
        else
            throw boost::program_options::invalid_option_value(token);
        
        return in;
    }
}

// constructor
OptionsExperiment::OptionsExperiment(string mode) : Options(mode) {
    
}


bool OptionsExperiment::parse(int argc, const char **argv) {
    try {
        
        GenericOptions genericOptions;
        
        // Group of options allowed only on command line
        po::options_description generic("General Options");
        Options::addProcessOptions_GenericOptions(genericOptions, generic);
        
        // Group of options allowed on both CML and in config files
        po::options_description config("Configuration");
        config.add_options()
        ;
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
            ("mode", po::value<string>(&mode)->required(), "Program Mode")
            ("experiment", po::value<experiment_t>(&experiment)->required(), "Experiment")
            ("subargs", po::value<std::vector<std::string> >(), "Arguments for experiment")
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
        pos.add("experiment",1);
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

        if (!vm.count("experiment")){
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // try parsing arguments.
        po::notify(vm);
      
        po::options_description expDesc;
        
        // Experiment: match sequence to upstream regions
        if (experiment == MATCH_SEQ_TO_UPSTREAM) {
            addProcessOptions_MatchSeqToUpstreamOptions(matchSeqToUpstream, expDesc);
        }
        // Experiment: match sequence to simulated noncoding regions
        if (experiment == MATCH_SEQ_TO_NONCODING) {
            addProcessOptions_MatchSeqToNoncodingOptions(matchSeqToNoncoding, expDesc);
        }
        
        cmdline_options.add(expDesc);
        
        // Collect all the unrecognized options from the first pass. This will include the
        // (positional) mode and command name, so we need to erase them
        vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());       // erase mode
        opts.erase(opts.begin());       // erase command name
        
        // Parse again...
        po::store(po::command_line_parser(opts).options(expDesc).run(), vm);
        
        
        // if help specified, print usage message and quit
        if (vm["help"].as<bool>()) {
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



void OptionsExperiment::addProcessOptions_MatchSeqToUpstreamOptions(MatchSeqToUpstreamOptions &options, po::options_description &processOptions) {
    
    Options::addProcessOptions_GenExtractUpstreamsOptions(options, processOptions);
    processOptions.add_options()
        ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
    ;
    
}


void OptionsExperiment::addProcessOptions_MatchSeqToNoncodingOptions(MatchSeqToNoncodingOptions &options, po::options_description &processOptions) {
    
    Options::addProcessOptions_GenReadSeqAndLabelsOptions(options, processOptions);
    processOptions.add_options()
        ("match-to", po::value<string>(&options.matchTo)->required(), "Sequence to match to.")
        ("length", po::value<size_t>(&options.length)->default_value(40), "Length of noncoding sequences")
        ("num-noncoding", po::value<size_t>(&options.numNoncoding)->default_value(1000), "Number of non-coding sequences.")
        ("order-nonc", po::value<unsigned>(&options.order)->default_value(0), "Order of non-coding model")
    ;
}
