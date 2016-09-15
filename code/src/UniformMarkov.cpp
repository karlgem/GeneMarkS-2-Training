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
UniformMarkov::UniformMarkov(unsigned order, const AlphabetDNA &alph, const CharNumConverter &cnc) : Markov(order, alph, cnc) {
    initialize();
}


// Construct the model probabilities from a list of sequences
void UniformMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    // get counts
    UniformCounts counts(order, *alphabet, cnc);
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
    
    // create a joint distribution for all orders less than or equal to this->order
    jointProbs.resize(this->order + 1);      // for order, order-1, order-2, ... 0
    jointProbs[this->order] = model;         // the joint for 'order' is in 'model'
    
    // derive remaining joint probabilities for 'order-1' and lower. This is used to compute words of length shorter than 'order+1'
    for (unsigned o = this->order; o > 0; o--)
        this->getLowerOrderJoint(o, jointProbs[o], jointProbs[o-1]);        //  take previous joint probs and reduce (marginalize) it by one
    
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
    for (size_t i = 0; i <= order; i++) {
        wordIndex <<= elementEncodingSize;      // create space at lower bits for a new element
        wordIndex += *currentElement;           // add new element to the wordIndex
        
        // set the mask to read word of 'i+1' elements
        for (size_t n = 0; n < elementEncodingSize; n++) {
            mask <<= 1;         // shift by one position
            mask |= 1;          // set lowest bit to one
        }
        
        // mask to remove old junk characters (doesn't affect wordIndex here. I do it just for consistency with remaining code)
        wordIndex = wordIndex & mask;

        // FIXME: use absolute probabilities of smaller order to account for short words
//        if (useLog)
//            score += log2(model[wordIndex]);
//        else
//            score *= model[wordIndex];
        
        currentElement++;
        
        // if reached end of sequence, then use joint probabilities to compute sequence prob
        if (currentElement == end) {
            if (useLog)
                return log2(jointProbs[i][wordIndex]);
            else
                return jointProbs[i][wordIndex];
        }
    }
    
    // compute joint probability of first word
    if (useLog)
        score += log2(jointProbs[this->order][wordIndex]);
    else
        score *= jointProbs[this->order][wordIndex];
    
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
    size_t numElements = alphabet->sizeValid();             // the number of elements that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                            // size of a word
    size_t numWords = pow(numElements, wordSize);           // number of possible words of size 'wordSize' with the given alphabet
    
    model.resize(numWords, 0);
}




void UniformMarkov::changeOrder(unsigned newOrder) {
    
    unsigned originalOrder = this->order;
    
    // increase order
    if (newOrder > originalOrder) {
        unsigned orderDifference = newOrder - originalOrder;        // difference in order
        
        // increment order by 1 at a time
        for (unsigned i = 0; i < orderDifference; i++) {
            
            // get joint probs
            jointProbs.push_back(vector<double>());                                                     // add new empty vector to joint probabilities
            incrementOrderByOne(this->order, jointProbs[jointProbs.size()-2], jointProbs[jointProbs.size()-1]);      // fill new vector
            this->order++;
        }
        
        
        // update Markov probabilities from highest joint probability
        this->model = jointProbs[jointProbs.size()-1];
        jointToMarkov(model);
        
    }
    // decrease order
    else if (newOrder < originalOrder) {
        throw logic_error("Decreasing Order functionality not yet implemented.");
    }
    
}



void UniformMarkov::incrementOrderByOne(unsigned currentOrder, const vector<double> &currentProbs, vector<double> &newProbs) const {
    
    size_t numElements = alphabet->sizeValid();                 // number of (valid) elements in alphabet (e.g. A,C,G,T)
    
    // set the size of the new probability space
    size_t newOrder = currentOrder+1;
    size_t wordSize = newOrder+1;                   // number of elements that make up a word
    size_t numWords = pow(numElements, wordSize);   // number of words in the new probability space
    newProbs.resize(numWords, 0);                   // allocate space for set new probability 
    
    // for each key of the current probabilities, the probability value needs to be copied
    // to all keys of "higher" order. Example, since our working representation is in bits,
    // suppose that the element encoding size is 2 (i.e. 2 bits needed to encode a single element).
    // If A = 00, C = 01, G = 10, T = 11, then the key GC=1001
    // In higher order, we have AGC=001001, CGC=011001, GGC=101001, TGC=111001
    // So we bit representations of intergers from 0 up to alphabet size (i.e. 4 in this case),
    // and we add them to GC, but shifted by 2 * element encoding size.
    //
    // Ex: To get TGC:
    // Representation of T:         11
    // Shift T by 2 * 2:        110000
    // Representation of GC:      1001
    // Add shifted T to GC:     111001
    
    size_t elementEncodingSize = ceil(log2(numElements));       // number of bits to encode an element (e.g. 2 bits for A,C,G,T)
    size_t bitsToShift = elementEncodingSize * currentOrder+1;
    
    // for each key in current probability space
    for (size_t key = 0; key < currentProbs.size(); key++) {
        
        // for each element to be added to the current key
        for (size_t e = 0; e < numElements; e++) {
            // shift the element by the number of letters already in the key, times the element encoding size
            size_t newKey = e << bitsToShift;
            
            // add existing key to new partial key; this gives the new full key
            newKey += key;
            
            // set probability value of new key to that of old key
            newProbs[newKey] = currentProbs[key];
        }
    }
}









