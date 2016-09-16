//
//  ModuleUtilities.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleUtilities.hpp"

#include "SequenceParser.hpp"

using namespace std;
using namespace gmsuite;


// constructor
ModuleUtilities::ModuleUtilities(const OptionsUtilities& opt) : options(opt) {
    
}

// Run module according to the provided options.
void ModuleUtilities::run() {
    
    if (options.utility == EXTRACT_UPSTR) {
        runExtractUpstream();
    }
    
    // unrecognized utility to run
    throw invalid_argument("Unknown utility function " + options.utility);
    
}



void ModuleUtilities::runExtractUpstream() {
    
    
    
}
