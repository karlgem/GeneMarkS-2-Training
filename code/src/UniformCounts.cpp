//
//  UniformCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "UniformCounts.hpp"
#include <stdexcept>
#include <math.h>

using namespace std;
using namespace gmsuite;

// Constructor: Create a uniform count model by defining it's order and alphabet
UniformCounts::UniformCounts(unsigned order, const AlphabetDNA* alph) : Counts(order, alph) {
    initialize();
}


// Construct the model counts from a list of sequences
void UniformCounts::construct(const vector<NumSequence> &sequences) {
    // start with zero counts
    resetCounts();
    
    for (vector<NumSequence>::const_iterator iter = sequences.begin(); iter != sequences.end(); iter++)
        count(iter->begin(), iter->end());
}


// Count the sequence.
void UniformCounts::count(NumSequence::const_iterator begin, NumSequence::const_iterator end) {
    updateCounts(begin, end, "increment");
}


// Decount the sequence.
void UniformCounts::decount(NumSequence::const_iterator begin, NumSequence::const_iterator end) {
    updateCounts(begin, end, "decrement");
}


// Generate a string representation of the model
string UniformCounts::toString() const {
    return "";
}


//Reset all counts to zero
void UniformCounts::resetCounts() {
    fill(model.begin(), model.end(), 0);        // set all values to zero
}




// Initialize the model by allocating space and setting counts to 0
void UniformCounts::initialize() {
    size_t numElements = alphabet->sizeValid();         // the number of eleents that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                        // size of a word
    size_t numWords = numElements << wordSize;          // number of possible words of size 'wordSize' with the given alphabet
    
    model.resize(numWords, 0);
}


// Update counts for a given sequence, by either incrementing or decrementing them.
void UniformCounts::updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation) {
    if (operation != "increment" && operation != "decrement")
        throw invalid_argument("Operation must either be 'increment' or 'decrement'.");
    
    // if sequence does not contain a word of size "order+1", then it doesn't contribute to the counts
    if (distance(begin, end) <= order)
        return;
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    size_t mask = 0;                        // used to clear wordIndex of junk and capture the relevant bits for the index
    
    
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
    }
    
    
    // for every word, add a count in the pwm
    while (currentElement != end) {
        
        // add new element to the wordIndex
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        wordIndex = wordIndex & mask;           // mask to remove old junk characters
        
        // increment or decrement
        if (operation == "increment") {
            model[wordIndex]++;              // increment count by 1
        }
        else {
            // update count at new word index
            if (model[wordIndex] == 0)
                throw out_of_range("Cannot decrement sequence counts below 0");
            
            model[wordIndex]--;              // decrement count by 1
        }
        
        currentElement++;                       // move to next letter
    }
}







