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

Counts::Counts(unsigned order, const AlphabetDNA* alph) {
    if (alph == NULL)
        throw std::invalid_argument("Alphabet cannot be NULL.");
    
    this->order = order;
    this->alphabet = alph;
}

// Get the model's order
unsigned Counts::getOrder() const {
    return order;
}

// Get the model's Alphabet
const AlphabetDNA* Counts::getAlphabet() const {
    return alphabet;
}


