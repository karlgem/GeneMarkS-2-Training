//
//  UniformMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "UniformMarkov.hpp"
#include "UniformCounts.hpp"

#include <math.h>
#include <stdexcept>

using namespace std;
using namespace gmsuite;

// Constructor:
UniformMarkov::UniformMarkov(unsigned order, const AlphabetDNA* alph) : Markov(order, alph) {
    initialize();
}


// Construct the model probabilities from a list of sequences
void UniformMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    // get counts
    UniformCounts counts(order, alphabet);
    counts.construct(sequences);
    
    // construct probabilities from counts
    construct(&counts);
}

// Construct the model probabilities from existing counts.
void UniformMarkov::construct(const Counts* counts, int pcount) {
    
    // counts cannot be NULL
    if (counts == NULL)
        throw invalid_argument("Counts cannot be NULL.");
    
    // counts alphabet must match markov alphabet
    if (counts->getAlphabet() != this->alphabet)
        throw invalid_argument("Counts alphabet must match Markov alphabet.");
    
    // counts order must match Markov order
    if (counts->getOrder() != this->order)
        throw invalid_argument("Counts order must match Markov order.");
    
    // cast counts to UniformCounts
    const UniformCounts* uniformCounts = dynamic_cast<const UniformCounts*>(counts);
    
    if (uniformCounts == NULL)
        throw invalid_argument("Counts should have type 'UniformCounts'.");
    
    // start by copying counts
    this->model = uniformCounts->model;
    
    double sum = 0;         // will contain sum of all elements (for normalization)
    
    // add pseudocounts and gather sum
    for (size_t n = 0; n < model.size(); n++) {                 // for each word
        model[n] += pcount;                                     // add pseudocount
        sum += model[n];                                        // add to sum
    }
    
    // normalize to get joint probabilities (e.g. P(ACG))
    for (size_t n = 0; n < model.size(); n++)                   // for each word
        if (sum != 0)                                           // check for division by zero
            model[n] /= sum;                                    // normalize word counts
    
    // derive joint probabilities for order 'order-1' and lower. This is used to compute words of length shorter than 'order+1'
    
    
    // convert joint probabilities to Markov (conditional)
    // e.g. P(ACG) -> (G|AC)
    jointToMarkov(model);
    

    
}

// Compute the score of a sequence using the model probabilities
double UniformMarkov::evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog) const {
  
    // if nothing to iterate over, return 0
    if (begin >= end) {
        if (useLog)
            return -std::numeric_limits<double>::infinity();
        else
            return 0;
    }
    
    double score = 1;       // since we're multiplying (if !useLog)
    if (useLog)
        score = 0;          // since we're summing
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    size_t mask = 0;                        // used to clear wordIndex of junk and capture the relevant bits for the index
    
    
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

        
        // FIXME: use absolute probabilities of smaller order to account for short words
//        if (useLog)
//            score += log2(model[wordIndex]);
//        else
//            score *= model[wordIndex];
    }
    
    
    // for every word, add a count in the pwm
    while (currentElement != end) {
        
        // add new element to the wordIndex
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        wordIndex = wordIndex & mask;           // mask to remove old junk characters
        
        if (useLog)
            score += log2(model[wordIndex]);
        else
            score *= model[wordIndex];
        
        currentElement++;                       // move to next letter
    }
    
    return score;

}

// Generate a string representation of the model
string UniformMarkov::toString() const {
    return "";
}


// Initialize the model by allocating space, setting the keys, and setting counts to 0
void UniformMarkov::initialize() {
    size_t numElements = alphabet->sizeValid();         // the number of eleents that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                        // size of a word
    size_t numWords = numElements << wordSize;          // number of possible words of size 'wordSize' with the given alphabet
    
    model.resize(numWords, 0);
}








