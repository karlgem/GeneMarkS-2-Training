//
//  test_MarkovsToString.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/18/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include <iostream>

#include "catch.hpp"
#include "UniformMarkov.hpp"
#include "PeriodicMarkov.hpp"
#include "NonUniformMarkov.hpp"

using namespace std;
using namespace gmsuite;



TEST_CASE("Testing toString of Markovs") {
    
    AlphabetDNA alph;
    
    vector<Sequence> sequences;
    sequences.push_back(Sequence("AAAGGAGGAA"));
    sequences.push_back(Sequence("CCAGGAGGCC"));
    sequences.push_back(Sequence("GGAGGAGGGG"));
    sequences.push_back(Sequence("TTAGGAGGTT"));
    CharNumConverter cnc(&alph);
    
    vector<NumSequence> numSequences;
    for (size_t n = 0; n < sequences.size(); n++)
        numSequences.push_back(NumSequence(sequences[n], cnc));
    
    NumAlphabetDNA numAlph(alph, cnc);
    
    SECTION("Uniform") {
        
        UniformMarkov m0 (0, numAlph);
        m0.construct(numSequences);
        cout << m0.toString() << endl;
        
        UniformMarkov m1 (1, numAlph);
        m1.construct(numSequences);
        cout << m1.toString() << endl;
    }
    
    SECTION("NonUniform") {
        
        size_t length = numSequences[0].size();
        NonUniformMarkov m0 (0, length, numAlph);
        m0.construct(numSequences);
        cout << m0.toString() << endl;
        
        NonUniformMarkov m1 (1, length, numAlph);
        m1.construct(numSequences);
        cout << m1.toString() << endl;
    }
    
    SECTION("Periodic") {
        
        size_t period = 3; //numSequences[0].size();
        PeriodicMarkov m0 (0, period, numAlph);
        m0.construct(numSequences);
        cout << m0.toString() << endl;
        
        PeriodicMarkov m1 (2, period, numAlph);
        m1.construct(numSequences);
        cout << m1.toString() << endl;
    }
    
}
