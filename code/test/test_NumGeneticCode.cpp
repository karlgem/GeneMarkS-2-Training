//
//  test_NumGeneticCode.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/22/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>

#include "catch.hpp"
#include "GeneticCode.hpp"
#include "NumGeneticCode.hpp"
#include "CharNumConverter.hpp"

#include <set>
#include <string>
#include <iostream>

using namespace std;
using namespace gmsuite;

TEST_CASE("Testing NumGeneticCode - ") {
    
    
    
    SECTION("Numeric Genetic Code 11") {
        
        // create genetic code and numeric converter
        AlphabetDNA alphabet;
        GeneticCode gc(GeneticCode::ELEVEN);
        CharNumConverter cnc(&alphabet);
        
        // create numeric code
        NumGeneticCode gcNum (gc, cnc);
        
        SECTION("Starts") {
            
            // for each start in genetic code, make sure it is recognized by numeric genetic code
            
        }
        
        
        
    }
    
}

