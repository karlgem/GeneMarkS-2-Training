//
//  NonUniformCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NonUniformCounts.hpp"
#include <math.h>

using namespace std;;
using namespace gmsuite;

// constructor
NonUniformCounts::NonUniformCounts(unsigned order, size_t length, const AlphabetDNA &alph) : Counts(order, alph) {
    
    if (length == 0)
        throw std::invalid_argument("Length cannot be 0.");
    
    this->length = length;
    initialize();
}


// construct the model from a list of sequences
void NonUniformCounts::construct(const vector<NumSequence> &sequences) {
    
    // start with zero counts
    resetCounts();
    
    for (vector<NumSequence>::const_iterator iter = sequences.begin(); iter != sequences.end(); iter++)
        count(iter->begin(), iter->end());
}


// generate string representation of the model
string NonUniformCounts::toString() const {
    return "";
}


// Reset counts to zero
void NonUniformCounts::resetCounts() {
    for (size_t p = 0; p < model.size(); p++)
        fill(model[p].begin(), model[p].end(), 0);      // set all values to zero
}

// initialize empty markov model
void NonUniformCounts::initialize() {
    
    size_t numElements = alphabet->sizeValid();     // the number of elements that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                    // the size of a word
    size_t numWords = pow(numElements, wordSize);   // number of possible words of size 'wordSize' with given alphabet
    
    model.resize(length);
    
    // for each period, allocate space of size 'numWords', and set all elements to zero
    for (size_t p = 0; p < length; p++) 
        model[p].resize(numWords, 0);
}


// Update counts for a given sequence, by either incrementing or decrementing them
void NonUniformCounts::updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation) {
    
    if (operation != "increment" && operation != "decrement")
        throw invalid_argument("Operation must either be 'increment' or 'decrement'.");
    
    // if nothing to iterate over, return
    if (begin >= end)
        return;
    
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
        
        // increment or decrement
        if (operation == "increment") {
            model[position][wordIndex]++;              // increment count by 1
        }
        else {
            // update count at new word index
            if (model[position][wordIndex] == 0)
                throw out_of_range("Cannot decrement sequence counts below 0");
            
            model[position][wordIndex]--;              // decrement count by 1
        }

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
        
        // increment or decrement
        if (operation == "increment") {
            model[position][wordIndex]++;              // increment count by 1
        }
        else {
            // update count at new word index
            if (model[position][wordIndex] == 0)
                throw out_of_range("Cannot decrement sequence counts below 0");
            
            model[position][wordIndex]--;              // decrement count by 1
        }
        
        position++;                                // increment frame
        if (position == length)
            break;                              // quit if frame reaches length of model
        
        currentElement++;                       // move to next letter
    }
}


// Get the model's period
size_t NonUniformCounts::getLength() const {
    return length;
}


