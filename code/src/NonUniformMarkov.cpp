//
//  NonUniformMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NonUniformMarkov.hpp"
#include "NonUniformCounts.hpp"

#include <math.h>
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
    
    // if nothing to iterate over, return 0
    if (begin >= end)
        return 0;
    
    double score = 1;       // since we're multiplying (if !useLog)
    if (useLog)
        score = 0;          // since we're summing
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    size_t mask = 0;                        // used to clear wordIndex of junk and capture the relevant bits for the index
    
    size_t position = 0;                       // position of pwm
    
    // To set the mask, we need to set the N lower bits to 1's, where N is the number
    // of bits required to encode a word of 'order+1' elements
    // A word is made up of elementEncondingSize * (order+1) bits. Set them to 1 by simply looping that many times
    
    NumSequence::const_iterator currentElement = begin;
    
    size_t wordIndex = 0;       // contains the index of the current word (made up of order+1 elements)
    
    // loop over first "order" elements to store them as part of the initial word index
    for (size_t i = 0; i < order; i++, currentElement++) {
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        // set the mask to read word of 'i+1' elements
        mask <<= 1;         // shift by one position
        mask |= 1;          // set lowest bit to one
        
        // mask to remove old junk characters (doesn't affect wordIndex here. I do it just for consistency with remaining code)
        wordIndex = wordIndex & mask;
        
        if (useLog)
            score += log2(model[position][wordIndex]);
        else
            score *= model[position][wordIndex];
        
        position++;                                 // increment position
        if (position == length)
            break;                                  // exit loop when reached length of markov model
    }
    
    
    // for every word, add a count in the pwm
    while (currentElement != end) {
        
        // add new element to the wordIndex
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        wordIndex = wordIndex & mask;           // mask to remove old junk characters
        
        if (useLog)
            score += log2(model[position][wordIndex]);
        else
            score *= model[position][wordIndex];
        
        position++;                                // increment frame
        if (position == length)
            break;                              // quit if frame reaches length of model
        
        currentElement++;                       // move to next letter
    }

    return score;
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
