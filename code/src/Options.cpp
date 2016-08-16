//
//  Options.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/11/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Options.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;
using namespace gmsuite;

// parse cmd arguments
bool Options::parse(int argc, const char **argv) {
    
    vector<string> config_fnames;                           // holds names of all configuration files (specified by user)
    
    po::options_description desc("General Options");        // holds general arguments and values
    desc.add_options()
        ("help,h", "Display help message")
        ("config", po::value(&config_fnames), "Config file where options may be specified (can be specified more than once)")
        ("verbose", po::value<int>(&verbose)->default_value(0), "Verbose level")
    ;
    
    // create set of hidden arguments which correspond to positional arguments. This is used
    // to add positional arguments, while not putting their description in the "options" section
    po::options_description positional_hidden;
    positional_hidden.add_options()
        ("mode", po::value<string>(&mode)->required(), "Program Mode")
    ;
    
    // congregate all options
    po::options_description all_options;
    all_options.add(desc);
    all_options.add(positional_hidden);
    
    // specify positional arguments
    po::positional_options_description pos;
    pos.add("mode", -1);
    
    
    // create storage component for storing names and values of arguments
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).          // pass in input
              options(all_options).                         // specify options list
              positional(pos).                              // specify which are positional
              run(),                                        // parse options
              vm);                                          // specify storage container
    
    // if help specified, print usage message and quit
    if (vm.count("help")) {
        cout << make_usage_string(basename(argv[0]), desc, pos) << endl;
        return false;
    }
    
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
            po::store(po::parse_config_file(ifs, all_options), vm);
        }
    }
    
    
    // try parsing arguments. Exception is returned if require option not provided.
    try {
        po::notify(vm);
    }
    catch (po::required_option &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
}


// extract program name from path
string Options::basename(const string &path) const {
    
#ifdef HAVE_BOOST_FILESYSTEM
    return boost::filesystem::path(p).stem().string();
#else
    size_t start = path.find_last_of("/");
    if(start == std::string::npos)
        start = 0;
    else
        ++start;
    return path.substr(start);
#endif
    
}


// get positional argument names for usage message
vector<string> Options::get_positional_args(const po::positional_options_description &pod) const {
    
    // reasonable upper limit for number of positional options:
    const int MAX = 1000;
    string last = pod.name_for_position(MAX);
    
    vector<string> pos_args;     // will contain positional argument names
    
    unsigned int i = 0;
    
    // keep looping until all arguments are read
    while (true) {
        string curr = pod.name_for_position(i);
    
        pos_args.push_back(curr);           // add argument
        
        // if reached last, break
        if (curr == last) {
            break;
        }
    }
    
    return pos_args;
}


string Options::make_usage_string(const string &program_name,
                                  const po::options_description &desc,
                                  const po::positional_options_description &pod) const {
    
    vector<string> components;
    
    components.push_back("Usage: ");
    components.push_back(program_name);
    
    // add positional arguments
    vector<string> pos_args = get_positional_args(pod);
    components.insert(components.end(), pos_args.begin(), pos_args.end());
    
    // add option descriptions
    if (desc.options().size() > 0) {
        components.push_back("[options]");
    }
    
    // create final usage string
    ostringstream oss;
    copy(components.begin(),
         components.end(),
         ostream_iterator<string>(oss, " "));       // add space between each member of "components"
    
    // add new line and print option descriptions
    oss << endl << desc;
    
    // return final string
    return oss.str();
        
    
    
}

