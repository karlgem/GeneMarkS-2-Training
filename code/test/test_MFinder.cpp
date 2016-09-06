//
//  test_MFinder.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

#include "catch.hpp"
#include "ModuleMFinder.hpp"
#include "OptionsMFinder.hpp"


using namespace gmsuite;

TEST_CASE("Testing Motif Finder") {
 
    srand(1);
    
    OptionsMFinder options("mfinder");
    options.fname_in = "/Users/Karl/repos/GeneMarkS-2/code/tmp/ecoli-motifs.fa";
//    options.fname_in = "/Users/Karl/repos/GeneMarkS-2/code/tmp/sample-motifs.fa";
//    options.fname_in = "/Users/Karl/repos/GeneMarkS-2/code/tmp/sample-rand-motifs.fa";
    options.pcounts = 0;
    
    ModuleMFinder mfinder (options);
    
    mfinder.run();
    
}