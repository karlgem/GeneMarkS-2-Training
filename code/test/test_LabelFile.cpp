//
//  test_LabelFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/31/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>

#include "catch.hpp"
#include "LabelFile.hpp"
#include <set>
#include <string>
#include <iostream>

using namespace std;
using namespace gmsuite;

TEST_CASE("Testing LabelFile - ") {
    
    
    
    SECTION("Read LST") {
        LabelFile labfile("/Users/Karl/repos/GeneMarkS-2/code/tmp/ecoli.lst", LabelFile::READ, LabelFile::LST);
        
        vector<Label*> output;
        labfile.read(output);
        
    }
    
}
