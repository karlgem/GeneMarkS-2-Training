//
//  NonUniformMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NonUniformMarkov.hpp"
#include "NonUniformCounts.hpp"

#include <stdexcept>

using std::invalid_argument;
using namespace gmsuite;

// Constructor:
NonUniformMarkov::NonUniformMarkov(unsigned order, size_t length, const AlphabetDNA* alph) : Markov(order, alph) {
    this->length = length;
    initialize();
}


// Construct the model probabilities from a list of sequences
void NonUniformMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    
    // get counts
    NonUniformCounts counts (order, length, alphabet);
    counts.construct(sequences);
    
    // construct probabilities from counts
    construct(&counts, pcount);
}

// Construct the model probabilities from existing counts.
void NonUniformMarkov::construct(const Counts* counts, int pcount) {
    
    // counts cannot be NULL
    if (counts == NULL)
        throw invalid_argument("Counts cannot be NULL.");
    
    // counts alphabet must match markov alphabet
    if (counts->getAlphabet() != this->alphabet)
        throw invalid_argument("Counts alphabet must match Markov alphabet.");
    
    // counts order must match Markov order
    if (counts->getOrder() != this->order)
        throw invalid_argument("Counts order must match Markov order.");
    
    // cast counts to PeriodiCounts
    const NonUniformCounts* nonunifCounts = dynamic_cast<const NonUniformCounts*>(counts);
    
    if (nonunifCounts == NULL)
        throw invalid_argument("Counts should have type 'PeriodicCounts'.");
    
    // counts period must match Markov period
    if (nonunifCounts->getLength() != this->length)
        throw invalid_argument("Counts length must match Markov length.");
    
    // start by copying counts
    this->model = nonunifCounts->model;
    
    vector<double> sums (length, 0);            // will contain sum of counts for each period
    
    // add pseudocounts
    for (size_t p = 0; p < length; p++)                         // for each period
        for (size_t n = 0; n < model[p].size(); n++) {          // for each word
            model[p][n] += pcount;                              // add pseudocount
            sums[p] += model[p][n];                             // add to sum
        }
    
    
    // normalize to get joint probabilities (e.g. P(ACG))
    for (size_t p = 0; p < length; p++)                         // for each period
        for (size_t n = 0; n < model[p].size(); n++)            // for each word
            if (sums[p] != 0)                                   // check for division by zero
                model[p][n] /= sums[p];                         // normalize word counts on sum
    
    
    // for each period, convert joint probabilities to Markov (conditional):
    // e.g. P(ACG) -> P(G|AC)
    for (size_t p = 0; p < length; p++)
        jointToMarkov(model[p]);
    
}

// Compute the score of a sequence using the model probabilities
double NonUniformMarkov::evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog) const {
    return 0;
}

// Generate a string representation of the model
string NonUniformMarkov::toString() const {
    return "";
}


// Initialize the model by allocating space, setting the keys, and setting counts to 0
void NonUniformMarkov::initialize() {
    size_t numElements = alphabet->sizeValid();         // the number of elements that can make up valid words (e.g. A,C,G,T)
    
    model.resize(length);
    
    // for each period, allocate space of size 'numWords', and set all elements to zero
    for (size_t p = 0; p < length; p++) {
        
        size_t wordSize = (p+1);                        // the size of a word
        size_t numWords = numElements << wordSize;          // number of possible words of size 'wordSize' with given alphabet
        
        model[p].resize(numWords, 0);
    }
}

// reset counts to zero
void NonUniformMarkov::resetCounts() {
}