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
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace gmsuite;

// Constructor:
NonUniformMarkov::NonUniformMarkov(unsigned order, size_t length, const NumAlphabetDNA &alph) : Markov(order, alph) {
    this->length = length;
    initialize();
}


// Construct the model probabilities from a list of sequences
void NonUniformMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    
    // get counts
    NonUniformCounts counts (order, length, *alphabet);
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
    
    // store joint probabilities
    jointProbs = model;
//    // create a set of joint distributions for each position
//    jointProbs.resize(this->length);
//    
//    // for each position
//    for (size_t p = 0; p < this->length; p++) {
//        
//        // create a joint distribution for all orders less than or equal to this->order
//        jointProbs[p].resize(this->order + 1);          // for order, order-1, order-2, ... 0
//        jointProbs[p][this->order] = model[p];          // the joint for 'order' is in 'model'
//        
//        // derive remaining joint probabilities for 'order-1' and lower. This is used to compute words of length shorter than 'order+1'
//        for (unsigned o = this->order; o > 0; o--)
//            this->getLowerOrderJoint(o, jointProbs[p][o], jointProbs[p][o-1]);        //  take previous joint probs and reduce (marginalize) it by one
//    }
    
    
    
    // for each period, convert joint probabilities to Markov (conditional):
    // e.g. P(ACG) -> P(G|AC)
    for (size_t p = 0; p < length; p++)
        jointToMarkov(model[p]);
    
}

// Compute the score of a sequence using the model probabilities
double NonUniformMarkov::evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog) const {
    
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
    
    size_t position = 0;                       // position of pwm
    
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
    stringstream ssm; ssm << fixed;
    
    // copy all joint probabilities because we may need to modify them by making them all the same order
    nonunif_joint_t joint = this->jointProbs;
    
    // create same set of keys for all positions by increasing the order of the smaller ones
    for (size_t p = 0; p < joint.size(); p++) {
        if (p < order)
            getHigherOrderJoint((unsigned) p, this->jointProbs[p], order, joint[p]);
    }
    
    
    // get length of "word" for current position (hint: they're all the same by now :))
    size_t wordLength = order+1;
    
    // loop over each key first (of any position, since they're all the same)
    for (size_t idx = 0; idx < joint[0].size(); idx++) {
        
        // convert index to numeric sequence
        NumSequence numSeq = this->indexToNumSequence(idx, wordLength);
        
        // convert numeric sequence to string sequence and add to ssm
        ssm << alphabet->getCNC()->convert(numSeq.begin(), numSeq.end());
        
        // loop over all positions, and print probabilities
        for (size_t p = 0; p < joint.size(); p++) {
            ssm << "\t" << joint[p][idx];
        }
        
        ssm << endl;    // goto new line (i.e. new key)
    }

    return ssm.str();
}


// Initialize the model by allocating space, setting the keys, and setting counts to 0
void NonUniformMarkov::initialize() {
    size_t numElements = alphabet->sizeValid();         // the number of elements that can make up valid words (e.g. A,C,G,T)
    
    model.resize(length);
    
    // for each period, allocate space of size 'numWords', and set all elements to zero
    for (size_t p = 0; p < length; p++) {
        
        size_t wordSize = (p<order ? p+1 : order+1);                        // the size of a word
        size_t numWords = pow(numElements, wordSize);   // number of possible words of size 'wordSize' with given alphabet
        
        model[p].resize(numWords, 0);
    }
}


size_t NonUniformMarkov::getLength() const {
    return length;
}




void NonUniformMarkov::changeOrder(unsigned newOrder) {
    
    unsigned originalOrder = this->order;
    
    // increase order
    if (newOrder > originalOrder) {
        unsigned orderDifference = newOrder - originalOrder;        // difference in order
        
        // increment order by 1 at a time
        for (unsigned i = 0; i < orderDifference; i++) {
            
            // increment order at each position
            for (unsigned p = 0; p < jointProbs.size(); p++) {
                if (p <= this->order)
                    continue;
                unsigned orderForPos = (p <= this->order ? p : this->order);
                
                // get joint probs
                vector<double> newProbs;
                incrementOrderByOne(orderForPos, jointProbs[p], newProbs);      // fill new vector
                jointProbs[p] = newProbs;
             
            }
            
           this->order++;           // increase order value
        }
        
        // update Markov probabilities from joint
        this->model = jointProbs;
        
        // for each period, convert joint probabilities to Markov (conditional):
        // e.g. P(ACG) -> P(G|AC)
        for (size_t p = 0; p < length; p++)
            jointToMarkov(model[p]);
        
        
        
        
    }
    // decrease order
    else if (newOrder < originalOrder) {
        throw logic_error("Decreasing Order functionality not yet implemented.");
    }
    
}

















