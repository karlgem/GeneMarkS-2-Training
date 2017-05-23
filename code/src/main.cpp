//
//  main.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/21/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <map>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>

#include "Options.hpp"
#include "OptionsGMS2.hpp"
#include "OptionsMFinder.hpp"
#include "OptionsUtilities.hpp"
#include "OptionsExperiment.hpp"
#include "OptionsGMS2Training.hpp"

#include "Module.hpp"
#include "ModuleGMS2.hpp"
#include "ModuleMFinder.hpp"
#include "ModuleUtilities.hpp"
#include "ModuleExperiment.hpp"
#include "ModuleGMS2Training.hpp"

using namespace std;
using namespace gmsuite;

#define MOD_GMS2 "gms2"
#define MOD_MFINDER "mfinder"
#define MOD_UTILITIES "utilities"
#define MOD_EXPERIMENT "experiment"
#define MOD_GMS2_TRAINING "gms2-training"

#define GMSUITE_VERSION 0.01

string usage_message(string progName) {
    stringstream ssm;
    
    ssm << "KGSuite, version " << GMSUITE_VERSION << endl;
    ssm << "Usage: " << progName << " mode" << endl;
    ssm << "The valid modes are:" << endl;
    ssm << "\t" << MOD_GMS2 << "\t" << "GeneMarkS2" << endl;
    ssm << "\t" << MOD_MFINDER << "\t" << "MFinder" << endl;
    ssm << "\t" << MOD_UTILITIES << "\t" << "Utilities" << endl;
    ssm << "\t" << MOD_EXPERIMENT << "\t" << "Experiment" << endl;
    ssm << "\t" << MOD_GMS2_TRAINING << "\t" << "GMS2 Training step" << endl;
    
    return ssm.str();
}

int main(int argc, const char * argv[]) {

    // first argument should be a valid "mode" (i.e. which module to run)
    if (argc < 2) {
        cerr << usage_message(argv[0]) << endl;
        return 1;
    }
    
    // get mode
    string aMode = argv[1];
    
    // create option and module objects based on mode
    if (aMode == MOD_GMS2) {                            // GeneMarkS2
        OptionsGMS2 options(aMode);                     // create option object
        
        if (!options.parse(argc, argv))                 // parse input arguments
            return 1;
        
        ModuleGMS2 module (options);                    // create module with options
        module.run();                                   // run module
    }
    else if (aMode == MOD_MFINDER) {
        OptionsMFinder options(aMode);                  // Motif Finder
        
        if (!options.parse(argc, argv))                 // parse input arguments
            return 1;
        
        ModuleMFinder module (options);                 // create module with options
        module.run();                                   // run module
    }
    else if (aMode == MOD_UTILITIES) {
        OptionsUtilities options(aMode);                // Utilities
        
        if (!options.parse(argc, argv))                 // parse input arguments
            return 1;
        
        ModuleUtilities module (options);               // create module with options
        module.run();                                   // run module
    }
    else if (aMode == MOD_EXPERIMENT) {
        OptionsExperiment options(aMode);               // Experiment
        
        if (!options.parse(argc, argv))                 // parse input arguments
            return 1;
        
        ModuleExperiment module (options);               // create module with options
        module.run();                                   // run module
    }
    else if (aMode == MOD_GMS2_TRAINING) {
        OptionsGMS2Training options(aMode);             // GMS2 Training
        
        if (!options.parse(argc, argv))                 // parse input arguments
            return 1;
        
        ModuleGMS2Training module (options);            // create module with options
        module.run();                                   // run module
    }
    else if (aMode == "--version") {
        cout << "Version: " << GMSUITE_VERSION << endl;
    }
    
    return 0;
}
