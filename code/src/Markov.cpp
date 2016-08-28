//
//  Markov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Markov.hpp"
#include <stdexcept>

using std::invalid_argument;
using namespace gmsuite;

// Contructor: Initialize a Markov model with a specific order and alphabet.
Markov::Markov(unsigned order, const AlphabetDNA* alph) {
    
    if (alph == NULL)
        throw invalid_argument("Alphabet cannot be NULL.");
    
    this->order = order;
    this->alphabet = alph;
    
}

// Get the model's order
unsigned Markov::getOrder() const {
    return order;
}

// Get the model's alphabet
const AlphabetDNA* Markov::getAlphabet() const {
    return alphabet;
}



