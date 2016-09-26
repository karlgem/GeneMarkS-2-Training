//
//  test_UniformCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include "catch.hpp"
#include "UniformCounts.hpp"
#include "UniformMarkov.hpp"

using namespace std;
using namespace gmsuite;

void testOrder(unsigned order) {
    AlphabetDNA alph;
    
    vector<Sequence> sequences;
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    sequences.push_back(Sequence("AGGAGG"));
    CharNumConverter cnc(&alph);
    
    vector<NumSequence> numSequences;
    for (size_t n = 0; n < sequences.size(); n++)
        numSequences.push_back(NumSequence(sequences[n], cnc));
    
    UniformCounts counts (order, alph, cnc);
    counts.construct(numSequences);
    
    UniformMarkov m(order,alph,cnc);
    m.construct(&counts);
    
    NumSequence s = m.emit(10);
    cout << cnc.convert(s.begin(), s.end()) << endl;
    
}

TEST_CASE("Testing Uniform Markovs") {
    
    SECTION("Order 0") {
        unsigned order = 0;
        testOrder(order);
    }
    
    SECTION("Order 1") {
        unsigned order = 1;
        testOrder(order);
    }
    
    SECTION("Order 2") {
        unsigned order = 2;
        testOrder(order);
    }
}
