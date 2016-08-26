//
//  GMS2Trainer.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "GMS2Trainer.hpp"

using namespace std;
using namespace gmsuite;

// default constructor
GMS2Trainer::GMS2Trainer() {
    
}


typedef vector<vector<double> > pwm_t;          // positional weight matrix where m[position][element] = probability


void initializePWM(pwm_t &pwm, size_t numElements, size_t order, size_t period) {
    
    if (period == 0)
        throw invalid_argument("Period can't be 0");
    
    size_t wordSize = (order+1);                        // number of elements per word
    size_t numOfWords = numElements << wordSize;        // number of possible words
    
    pwm.resize(period);
    
    // for each period, allocate space of size 'numOfWords', and set all elements to zero
    for (size_t p = 0; p < period; p++) {
        pwm[p].resize(numOfWords, 0);
    }
}

// update counts of coding model from a sequence
void updateCounts(pwm_t &pwm, NumSequence::const_iterator begin, NumSequence::const_iterator end, unsigned order, unsigned period) {
    
    // if sequence does not contain a word of size "order+1", then it doesn't contribute to the counts
    if (distance(begin, end) <= order)
        return;
    
    size_t elementEncodingSize = 2;         // the number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    size_t mask = 0;                        // used to clear wordIndex of junk and capture the relevant bits for the index
    
    size_t frame = 0;                       // frame of pwm (e.g. for 3-periodic models)
    
    // To set the max, we need to set the N lower bits to 1's, where N is the number
    // of bits required to encode a word of 'order+1' elements
    // A word is made up of elementEncondingSize * (order+1) bits. Set them to 1 by simply looping that many times
    for (size_t i = 0; i < elementEncodingSize * (order+1); i++) {
        mask <<= 1;         // shift by one position
        mask |= 1;          // set lowest bit to one
    }
    
    
    NumSequence::const_iterator currentElement = begin;
    
    size_t wordIndex = 0;       // contains the index of the current word (made up of order+1 elements)
    
    // loop over first "order" elements to store them as part of the initial word index
    for (size_t i = 0; i < order; i++, currentElement++) {
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        frame++;                                // increment frame
        if (frame == period)
            frame = 0;                          // reset frame when it reaches "period"
    }
    
    
    // for every word, add a count in the pwm
    while (currentElement != end) {
        
        // add new element to the wordIndex
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        // update count at new word index
        pwm[frame][wordIndex]++;
        
        frame++;                                // increment frame
        if (frame == period)
            frame = 0;                          // reset frame when it reaches "period"
        
        currentElement++;                       // move to next letter
    }
}

void estimateParamtersCoding(const NumSequence &sequence, const vector<gmsuite::Label *> &labels) {
    
    pwm_t pwm;
    initializePWM(pwm, 4, 1, 1);        // TODO
    
    // get counts for 3 period markov model given order
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left;                // get left position of fragment
        size_t right = (*iter)->right;              // get right position of fragment
        size_t length = left + right + 1;           // compute fragment length
        
        updateCounts(pwm, sequence.begin() + left, sequence.begin() + length, 0, 3);
    }
    
}


void GMS2Trainer::estimateParameters(const NumSequence &sequence, const vector<gmsuite::Label *> &labels) {
    
    
    // estimate parameters for coding model
    estimateParamtersCoding(sequence, labels);
    
    // estimate parameters for noncoding model
    
    // estimate parameters for start context
    
    // estimate parameters for motif models
}