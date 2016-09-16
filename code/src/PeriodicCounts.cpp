//
//  PeriodicCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/26/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "PeriodicCounts.hpp"
#include <math.h>
#include <assert.h>
#include <stdexcept>

using namespace std;;
using namespace gmsuite;

// constructor
PeriodicCounts::PeriodicCounts(unsigned order, size_t period, const AlphabetDNA &alph, const CharNumConverter &cnc) : Counts(order, alph, cnc) {
    
    if (period == 0)
        throw std::invalid_argument("Period cannot be 0.");
    
    this->period = period;
    initialize();
}


// construct the model from a list of sequences
void PeriodicCounts::construct(const vector<NumSequence> &sequences) {
    
    // start with zero counts
    resetCounts();
    
    for (vector<NumSequence>::const_iterator iter = sequences.begin(); iter != sequences.end(); iter++)
        count(iter->begin(), iter->end());
}


// generate string representation of the model
string PeriodicCounts::toString() const {
    return "";
}


// Reset counts to zero
void PeriodicCounts::resetCounts() {
    for (size_t p = 0; p < model.size(); p++)
        fill(model[p].begin(), model[p].end(), 0);          // set all values to zero
}

// initialize empty markov model
void PeriodicCounts::initialize() {
    
    size_t numElements = alphabet->sizeValid();         // the number of elements that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                        // the size of a word
    size_t numWords = pow(numElements, wordSize);       // number of possible words of size 'wordSize' with given alphabet
    
    model.resize(period);
    
    // for each period, allocate space of size 'numWords', and set all elements to zero
    for (size_t p = 0; p < period; p++) {
        model[p].resize(numWords, 0);
    }
}


// Update counts for a given sequence, by either incrementing or decrementing them
void PeriodicCounts::updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation, bool reverseComplement) {
    
    if (operation != "increment" && operation != "decrement")
        throw invalid_argument("Operation must either be 'increment' or 'decrement'.");
    
    // if sequence does not contain a word of size "order+1", then it doesn't contribute to the counts
    if (distance(begin, end) <= order)
        return;
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    size_t mask = 0;                        // used to clear wordIndex of junk and capture the relevant bits for the index
    
    size_t frame = 0;                       // frame of pwm (e.g. for 3-periodic models)
    
    // To set the max, we need to set the N lower bits to 1's, where N is the number
    // of bits required to encode a word of 'order+1' elements
    // A word is made up of elementEncondingSize * (order+1) bits. Set them to 1 by simply looping that many times
    for (size_t i = 0; i < elementEncodingSize * (order+1); i++) {
        mask <<= 1;         // shift by one position
        mask |= 1;          // set lowest bit to one
    }
    
    size_t seqSize = distance(begin, end);          // total number of nucleotides in sequence
    size_t remainingSize = seqSize;                 // number of nucleotides still to be counted
    
    NumSequence::const_iterator currentElement = begin;
    if (reverseComplement)
        currentElement = end-1;
    
    size_t wordIndex = 0;       // contains the index of the current word (made up of order+1 elements)
    
    // loop over first "order" elements to store them as part of the initial word index
    for (size_t i = 0; i < order; i++) {
        
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        if (reverseComplement)
            wordIndex += cnc->complement(*currentElement);           // complement element then add to the wordIndex
        else
            wordIndex += *currentElement;                           // add element to the wordIndex
        
        frame++;                                // increment frame
        if (frame == period)
            frame = 0;                          // reset frame when it reaches "period"
        
        // move to next element on negative or positive strand
        (reverseComplement ? currentElement-- : currentElement++);
        
        remainingSize--;
    }
    
    
    // for every word, add a count in the pwm
    while (remainingSize != 0) {
        
        // add new element to the wordIndex
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        
        if (reverseComplement)
            wordIndex += cnc->complement(*currentElement);           // complement element then add to the wordIndex
        else
            wordIndex += *currentElement;                           // add element to the wordIndex
        
        wordIndex = wordIndex & mask;           // mask to remove old junk characters
        
        // increment or decrement
        if (operation == "increment") {
            model[frame][wordIndex]++;              // increment count by 1
        }
        else {
            // update count at new word index
            if (model[frame][wordIndex] == 0)
                throw out_of_range("Cannot decrement sequence counts below 0");
            
            model[frame][wordIndex]--;              // decrement count by 1
        }
        
        frame++;                                // increment frame
        if (frame == period)
            frame = 0;                          // reset frame when it reaches "period"
        
        // move to next element on negative or positive strand
        (reverseComplement ? currentElement-- : currentElement++);
        
        remainingSize--;
    }
}


// Get the model's period
size_t PeriodicCounts::getPeriod() const {
    return period;
}


