//
//  test_CodingMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include "catch.hpp"

#include "NumSequence.hpp"
#include "PeriodicCounts.hpp"
#include "CodingMarkov.hpp"

using namespace std;
using namespace gmsuite;


TEST_CASE("Testing CodingMarkov") {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    GeneticCode geneticCode(GeneticCode::ELEVEN);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    
    PeriodicCounts p (2, 3, numAlph);
    
    vector<Sequence> sequences;
    sequences.push_back(Sequence("AGG"));
    
    vector<NumSequence> numSequences;
    for (size_t n = 0; n < sequences.size(); n++)
        numSequences.push_back(NumSequence(sequences[n], cnc));
    
    p.construct(numSequences);
    
    CodingMarkov m(2,3,numAlph, numGeneticCode);
    m.construct(&p);
}
