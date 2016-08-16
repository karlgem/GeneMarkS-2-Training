//
//  NumSequence.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NumSequence.hpp"

using namespace gmsuite;

// Default constructor: create an empty numeric sequence.
NumSequence::NumSequence() {
    
}

// Constructor: create a numeric sequence from a regular sequence object.
NumSequence::NumSequence(const Sequence &sequence, const CharNumConverter &converter ) {
    
    // TODO: convert letter sequence to numeric
    converter.convert(sequence.begin(), sequence.end(), this->numSeq);
}


// Access an element from a const numeric sequence (i.e. cannot be modified)
const NumSequence::num_t&  NumSequence::operator[](size_type idx) const {
    return numSeq[idx];
}

// access an element
NumSequence::num_t& NumSequence::operator[](size_type idx) {
    return numSeq[idx];
}