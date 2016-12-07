//
//  NonCodingCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/7/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NonCodingCounts.hpp"

using namespace gmsuite;

NonCodingCounts::NonCodingCounts(unsigned order, const NumAlphabetDNA &alph) : UniformCounts(order, alph) {
    
}

// count sequence
void NonCodingCounts::count(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement) {
    UniformCounts::count(begin, end, false);
    UniformCounts::count(begin, end, true);
}


// Decount the sequence. This method calls the updateCounts (decrement) method, which
void NonCodingCounts::decount(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement) {
    UniformCounts::decount(begin, end, false);
    UniformCounts::decount(begin, end, true);
}
