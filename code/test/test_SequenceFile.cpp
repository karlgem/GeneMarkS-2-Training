//
//  test_SequenceFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>

#include "catch.hpp"
#include "SequenceFile.hpp"
#include <set>
#include <string>
#include <iostream>

using namespace std;
using namespace gmsuite;

TEST_CASE("Testing SequenceFile - ") {
    
    
    
    SECTION("Read ecoli") {
        SequenceFile seqfile("/Users/Karl/repos/genemark-suite/code/tmp/ecoli.fa", SequenceFile::READ, SequenceFile::FASTA);
        
        vector<Sequence> output;
        seqfile.read(output);
        
    }
    
    SECTION("Read motifs") {
        SequenceFile seqfile("/Users/Karl/repos/genemark-suite/code/tmp/motifs.fa", SequenceFile::READ, SequenceFile::FASTA);
        
        vector<Sequence> output;
        seqfile.read(output);
    }
    
    
}
