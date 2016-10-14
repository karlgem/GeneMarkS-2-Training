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

Options::Options(string name) {
    this->mode = name;
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

void Options::addProcessOptions_GenericOptions(GenericOptions &options, po::options_description &processOptions) {
    processOptions.add_options()
        ("help,h", po::bool_switch(&options.help)->default_value(false), "Display help message")
        ("debug", po::bool_switch(&options.debug)->default_value(false), "Debug")
        ("verbose", po::bool_switch(&options.verbose)->default_value(false), "Verbose")
    ;
}

void Options::addProcessOptions_GenReadSeqAndLabelsOptions(GenReadSeqAndLabelsOptions &options, po::options_description &processOptions) {
    
    processOptions.add_options()
        ("sequence,s", po::value<string>(&options.fn_seqeuence)->required(), "Sequence filename")
        ("label,l", po::value<string>(&options.fn_labels)->required(), "Labels filename")
    ;
}

void Options::addProcessOptions_GenExtractUpstreamsOptions(GenExtractUpstreamsOptions &options, po::options_description &processOptions) {
    
    addProcessOptions_GenReadSeqAndLabelsOptions(options, processOptions);
    processOptions.add_options()
        ("length", po::value<size_t>(&options.length)->required(), "Upstream length")
        ("allow-overlap-with-cds", po::bool_switch()->default_value(false), "Allow upstream sequences to overlap with previous CDS")
        ("min-gene-length", po::value<size_t>(&options.minGeneLength)->default_value(0), "Minimum gene length")
    ;
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
        
        ++i;
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

