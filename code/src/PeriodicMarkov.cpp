//
//  PeriodicMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "PeriodicMarkov.hpp"
#include "PeriodicCounts.hpp"

#include <math.h>
#include <stdexcept>

using std::invalid_argument;
using namespace gmsuite;

// Constructor:
PeriodicMarkov::PeriodicMarkov(unsigned order, size_t period, const AlphabetDNA &alph, const CharNumConverter &cn) : Markov(order, alph), cnc(cn) {
    this->period = period;
    initialize();
}


// Construct the model probabilities from a list of sequences
void PeriodicMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    
    // get counts
    PeriodicCounts counts (order, period, *alphabet, cnc);
    counts.construct(sequences);
    
    // construct probabilities from counts
    construct(&counts, pcount);
}

// Construct the model probabilities from existing counts.
void PeriodicMarkov::construct(const Counts* counts, int pcount) {
    
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
    const PeriodicCounts* periodicCounts = dynamic_cast<const PeriodicCounts*>(counts);
    
    if (periodicCounts == NULL)
        throw invalid_argument("Counts should have type 'PeriodicCounts'.");
    
    // counts period must match Markov period
    if (periodicCounts->getPeriod() != this->period)
        throw invalid_argument("Counts period must match Markov period.");
    
    // start by copying counts
    this->model = periodicCounts->model;
    
    vector<double> sums (period, 0);            // will contain sum of counts for each period
    
    // add pseudocounts
    for (size_t p = 0; p < period; p++)                         // for each period
        for (size_t n = 0; n < model[p].size(); n++) {          // for each word
            model[p][n] += pcount;                              // add pseudocount
            sums[p] += model[p][n];                             // add to sum
        }
    
    
    // normalize to get joint probabilities (e.g. P(ACG))
    for (size_t p = 0; p < period; p++)                         // for each period
        for (size_t n = 0; n < model[p].size(); n++)            // for each word
            if (sums[p] != 0)                                   // check for division by zero
                model[p][n] /= sums[p];                         // normalize word counts on sum
    
    
    // create a set of joint distributions for each period
    jointProbs.resize(this->period);
    
    // for each period
    for (size_t p = 0; p < this->period; p++) {
        
        // create a joint distribution for all orders less than or equal to this->order
        jointProbs[p].resize(this->order + 1);      // for order, order-1, order-2, ... 0
        jointProbs[p][this->order] = model[p];         // the joint for 'order' is in 'model'
        
        // derive remaining joint probabilities for 'order-1' and lower. This is used to compute words of length shorter than 'order+1'
        for (unsigned o = this->order; o > 0; o--)
            this->getLowerOrderJoint(o, jointProbs[p][o], jointProbs[p][o-1]);        //  take previous joint probs and reduce (marginalize) it by one
    }
    
    // for each period, convert joint probabilities to Markov (conditional):
    // e.g. P(ACG) -> P(G|AC)
    for (size_t p = 0; p < period; p++)
        jointToMarkov(model[p]);
    
}

// Compute the score of a sequence using the model probabilities
double PeriodicMarkov::evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog) const {
    
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
    
    size_t frame = 0;                       // frame of pwm
    
    // To set the mask, we need to set the N lower bits to 1's, where N is the number
    // of bits required to encode a word of 'order+1' elements
    // A word is made up of elementEncondingSize * (order+1) bits. Set them to 1 by simply looping that many times
    
    NumSequence::const_iterator currentElement = begin;
    
    size_t wordIndex = 0;       // contains the index of the current word (made up of order+1 elements)
    
    // loop over first "order" elements to store them as part of the initial word index
    for (size_t i = 0; i <= order; i++, currentElement++) {
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        // set the mask to read word of 'i+1' elements
        for (size_t n = 0; n < elementEncodingSize; n++) {
            mask <<= 1;         // shift by one position
            mask |= 1;          // set lowest bit to one
        }
        
        // mask to remove old junk characters (doesn't affect wordIndex here. I do it just for consistency with remaining code)
        wordIndex = wordIndex & mask;
        
        // FIXME: use joint probability model (of lower orders, if necessary)
        if (useLog)
            score += log2(model[frame][wordIndex]);
        else
            score *= model[frame][wordIndex];
        
        frame++;                                 // increment frame
        if (frame == period)
            frame = 0;                           // set frame to 0 when period is reached
    }
    
    
    // for every word, add a count in the pwm
    while (currentElement != end) {
        
        // add new element to the wordIndex
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        wordIndex = wordIndex & mask;           // mask to remove old junk characters
        
        if (useLog)
            score += log2(model[frame][wordIndex]);
        else
            score *= model[frame][wordIndex];
        
        frame++;                                 // increment frame
        if (frame == period)
            frame = 0;                           // set frame to 0 when period is reached
        
        currentElement++;                       // move to next letter
    }
    
    return score;
}

// Generate a string representation of the model
string PeriodicMarkov::toString() const {
    return "";
}




// Initialize the model by allocating space, setting the keys, and setting counts to 0
void PeriodicMarkov::initialize() {
    size_t numElements = alphabet->sizeValid();         // the number of elements that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                        // the size of a word
    size_t numWords = pow(numElements, wordSize);       // number of possible words of size 'wordSize' with given alphabet
    
    model.resize(period);
    
    // for each period, allocate space of size 'numWords', and set all elements to zero
    for (size_t p = 0; p < period; p++) {
        model[p].resize(numWords, 0);
    }
}


























