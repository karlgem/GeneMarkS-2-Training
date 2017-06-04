//
//  test_ModelFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/9/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include <iostream>

#include "catch.hpp"
#include "ModelFile.hpp"

using namespace std;
using namespace gmsuite;

TEST_CASE("Testing Model File") {
    
    ModelFile mfile("/Users/Karl/repos/GeneMarkS-2/code/tmp/sample.mod", ModelFile::READ);
    
    map<string,string> keyValue;
    mfile.read(keyValue);
    
    for (map<string,string>::const_iterator iter = keyValue.begin(); iter != keyValue.end(); iter++) {
        cout << iter->first << "\n" << iter->second << endl;
    }
    
    ModelFile mfile_write("/Users/Karl/repos/GeneMarkS-2/code/tmp/sample_write.mod", ModelFile::WRITE);
    mfile_write.write(keyValue);

}