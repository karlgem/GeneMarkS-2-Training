//
//  Counts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/26/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Counts.hpp"
#include <stdexcept>

using std::invalid_argument;
using namespace gmsuite;

Counts::Counts(unsigned order, const NumAlphabetDNA &alph) {
    this->order = order;
    this->alphabet = &alph;
}

// Get the model's order
unsigned Counts::getOrder() const {
    return order;
}

// Get the model's Alphabet
const NumAlphabetDNA* Counts::getAlphabet() const {
    return alphabet;
}


// Count the sequence.
void Counts::count(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement) {
    updateCounts(begin, end, "increment", reverseComplement);
}


// Decount the sequence.
void Counts::decount(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement) {
    updateCounts(begin, end, "decrement", reverseComplement);
}
