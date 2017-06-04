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
#include "AlphabetDNA.hpp"
#include "NumGeneticCode.hpp"
#include "CharNumConverter.hpp"

#include <set>
#include <string>
#include <iostream>

using namespace std;
using namespace gmsuite;


void testGeneticCode(GeneticCode::gcode_t geneticCode) {
    // create genetic code and numeric converter
    AlphabetDNA alphabet;
    GeneticCode gc(geneticCode);
    CharNumConverter cnc(&alphabet);
    
    // create numeric code
    NumGeneticCode gcNum (gc, cnc);
    
    for (GeneticCode::ttk_const_iterator iter = gc.begin(); iter != gc.end(); iter++) {
        
        // convert codon from string to seq_t
        CharNumConverter::seq_t conv;
        cnc.convert(*iter, conv);
        
        // make sure its status fits with original genetic code
        REQUIRE(gcNum.isStart(conv) == gc.isStart(*iter));
        REQUIRE(gcNum.isStop(conv) == gc.isStop(*iter));
    }

}


TEST_CASE("Testing NumGeneticCode - ") {
    
    GeneticCode::gcode_t geneticCodes [2];
    
    SECTION("Numeric Genetic Code 11") {
        testGeneticCode(GeneticCode::ELEVEN);
    }
    
    SECTION("Numeric Genetic Code 4") {
        testGeneticCode(GeneticCode::FOUR);
    }
    
}

