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
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>


using namespace std;
using namespace gmsuite;

// Constructor:
UniformMarkov::UniformMarkov(unsigned order, const AlphabetDNA &alph, const CharNumConverter &cnc) : Markov(order, alph, cnc) {
    initialize();
}


// Construct the model probabilities from a list of sequences
void UniformMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    // get counts
    UniformCounts counts(order, *alphabet, *cnc);
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
    
    stringstream ssm;
    
    // for each index in join probability distribution
    for (size_t idx = 0; idx < this->model.size(); idx++) {
        // convert index to numeric sequence
        NumSequence numSeq = this->indexToNumSequence(idx, this->order+1);
        
        // convert numeric sequence to string sequence and add to ssm
        ssm << cnc->convert(numSeq.begin(), numSeq.end());
        
        // add probability of current index
        ssm << "\t" << this->jointProbs[order][idx] << endl;
    }
    
    return ssm.str();
    
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



NumSequence UniformMarkov::emit(NumSequence::size_type length) const {
    
    if (length == 0)
        return NumSequence();
    
    // get CDF per conditional
    vector<double> cdf = model;
    Markov::getCDFPerConditional(this->order, cdf);
    
    // to generate first "order+1" elements, or to generate a sequence where length < order+1,
    // get CDF of joint probabilities
    vector<vector<double> > cdfJoint = jointProbs;
    for (unsigned p = 0; p < jointProbs.size(); p++) {
        for (size_t n = 1; n < cdfJoint[p].size(); n++) {
            cdfJoint[p][n] += cdfJoint[p][n-1];
        }
    }
    
    // get random number generator
    typedef boost::mt19937                     ENG;     // Mersenne Twister
    typedef boost::random::uniform_real_distribution<double> DIST;                // Normal Distribution
    typedef boost::variate_generator<ENG,DIST> GEN;     // Variate generator
    
    ENG  eng;
    DIST dist(0,1);
    GEN  gen(eng,dist);
    
    
    // if sequence length <= order+1, simply emit from joint
    if (length <= this->order + 1) {
        double u = gen();
        
        size_t i;       // index of emitted word
        
        // find index
        for (i = 0; i < cdfJoint[length-1].size(); i++) {
            if (u < cdfJoint[length-1][i])
                break;              // sample
        }
        
        if (i == cdfJoint[length-1].size())
            throw logic_error("Could not sample from CDF. Something is wrong.");
        
        // convert word from index
        return Markov::indexToNumSequence(i, length);
    }
    
    
    // otherwise, length > order+1. "order+1" elements from joint distribution,
    // then use conditional to emit one at a time.
    
    // Emit first order+1 elements:
    size_t i;
    double u = gen();
    size_t wordLength = order+1;
    for (i = 0; i < cdfJoint[wordLength-1].size(); i++) {
        if (u < cdfJoint[wordLength-1][i])
            break;              // sample
    }
    
    if (i == cdfJoint[wordLength-1].size())
        throw logic_error("Could not sample from CDF. Something is wrong.");

    // convert word from index
    vector<NumSequence::num_t> numSeq;
    Markov::indexToNumSequence(i, wordLength, numSeq);
    
    // extend numSeq to allocate spacer of size "length" (i.e. for entire sequence)
    numSeq.resize(length);
    
    // Emit remaining sequence. Take last "order" elements from previous word, and sample new element from CDF
    // corresponding to that base
    
    size_t numElements = alphabet->sizeValid();             // number of elements (ACGT)
    size_t elementEncodingSize = ceil(log2(numElements));   // bits required to encode all elements
    size_t maskOne = 0;            // mask single letter
    size_t maskBase = 0;        // masks the last "order" characters of word, to be used as base for new word
    size_t blockSize = alphabet->sizeValid();
    size_t wordIndex = i;           // contains index of the current word (made up of order+1 elements)
    
    
    for (size_t n = 0; n < elementEncodingSize; n++) {
        maskOne <<= 1;          // shift by one
        maskOne |= 1;           // set rightmost bit to 1
    }
    
    // get a mask that captures base
    for (size_t n = 0; n < order; n++) {
        for (size_t n = 0; n < elementEncodingSize; n++) {
            maskBase <<= 1;      // shift by one
            maskBase |= 1;      // set rightmost bit to 1
        }
    }

    // loop over remaining "sequence"
    for (size_t n = wordLength; n < length; n++) {
        
        u = gen();      // generate unif(0,1)
        
        // get base from previous word
        size_t base = wordIndex & maskBase;
        
        // loop over block of that base, and sample from CDF and u
        size_t blockStart = base << elementEncodingSize;
        size_t idx;
        
        for (idx = blockStart; idx < blockStart+blockSize; idx++) {
            if (u < cdf[idx])
                break;
        }
        
        if (idx == blockStart+blockSize)
            throw logic_error("Could not sample from CDF. Something is wrong.");
        
        // add letter to numSeq
        numSeq[n] = (NumSequence::num_t) ( idx & maskOne);
        
        // update wordIndex
        wordIndex = idx;
    }
    
    return NumSequence(numSeq);
}





















