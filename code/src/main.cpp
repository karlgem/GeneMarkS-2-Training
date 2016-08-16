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

using namespace std;
using namespace gmsuite;

// Define helper methods to create OptionsXXX objects, depending on the mode set
template<typename T> Options* createOptionsFromMode() { return new T; };

typedef map<string, Options*(*)()> mode_options_t;      // map a key (mode) to corresponding OptionsXXX object

mode_options_t init_mode_options_map() {
    mode_options_t map;
    
    map["gms2"] = &createOptionsFromMode<OptionsGMS2>;          // add gms2 options pair
    
    return map;
}


int main(int argc, const char * argv[]) {
    

    Options options;
    
    // get "mode" from options
    if (!options.parse(argc, argv))
        return 1;
    
    /***** Now that mode is parsed, we create the corresponding OptionsXXX object to handle its parameters *****/
    
    // start by creating a map that maps modes to OptionsXXX objects
    mode_options_t mode_map = init_mode_options_map();
    
    // use that map to create the corresponding options object
    Options* optionsMode;
    
    // depending on mode, compile new options object
    try {
        optionsMode = mode_map.at(options.mode)();
    }
    catch (out_of_range &ex) {
        cerr << "Error" << "Invalid mode: " << options.mode << endl;
        cerr << "Use the -help command for more info" << endl;
        return 1;
    }
    
    // parse CML with new mode-specific options object
    optionsMode->parse(argc, argv);
    
    // do something with options object
    
    
    // free options object
    delete optionsMode;
    
    
    return 0;
}