//
//  test_LongestCommonSubstring.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>
#include <iostream>

#include "catch.hpp"
#include "Sequence.hpp"
#include "NumSequence.hpp"
#include "SequenceAlgorithms.hpp"

using namespace std;
using namespace gmsuite;

TEST_CASE("Testing Longest Commong Substring") {
    
    SECTION("") {
        
        Sequence s1("ATGGAGGAGGCGA");
        Sequence s2("AGGTTCGATAGGAGGAGGTA");
        
        AlphabetDNA alph;
        CharNumConverter cnc(&alph);
        
        NumSequence ns1(s1, cnc);
        NumSequence ns2(s2, cnc);
        
        NumSequence common = SequenceAlgorithms::longestCommonSubstring(ns1, ns2);
        
        cout << "Longest Common Subsequence is: " <<  cnc.convert(common.begin(), common.end()) << endl;   
    }
}
