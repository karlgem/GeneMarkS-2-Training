//
//  test_PeriodicCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/28/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include "catch.hpp"

#include "NumSequence.hpp"
#include "PeriodicCounts.hpp"

using namespace std;
using namespace gmsuite;


TEST_CASE("Testing PeriodicCounts") {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    PeriodicCounts p (1, 3, numAlph);
    
    vector<Sequence> sequences;
    sequences.push_back(Sequence("ACGT"));
    sequences.push_back(Sequence("AGAG"));
    
    
    
    vector<NumSequence> numSequences;
    for (size_t n = 0; n < sequences.size(); n++)
        numSequences.push_back(NumSequence(sequences[n], cnc));
    
    p.construct(numSequences);
}
