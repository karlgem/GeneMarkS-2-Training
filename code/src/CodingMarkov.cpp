//
//  CodingMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "CodingMarkov.hpp"
#include <math.h>

using namespace std;
using namespace gmsuite;

CodingMarkov::CodingMarkov(unsigned order, size_t period, const NumAlphabetDNA &alph, const NumGeneticCode &geneticCode) : PeriodicMarkov(order, period, alph) {
    this->geneticCode = &geneticCode;
}

CodingMarkov::CodingMarkov(const vector<vector<pair<string, double> > > &keyValue, const NumAlphabetDNA &alph, const CharNumConverter &cnc) : PeriodicMarkov(0, 1, alph) {
    
    unsigned maxOrder = 0;
    
    // loop through keyValue pairs and get the value of the maximum order
    for (size_t pos = 0; pos < keyValue.size(); pos++) {
        for (size_t n = 0; n < keyValue[pos].size(); n++) {
            size_t keySize = keyValue[pos][n].first.size();
            if (keySize > 0)
                if (maxOrder < keySize-1)
                    maxOrder = (unsigned) keySize-1;
        }
    }
    
    this->order = maxOrder;
    this->period = keyValue.size();
    initialize();
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    
    for (size_t frame = 0; frame < keyValue.size(); frame++) {
        for (vector<pair<string, double> >::const_iterator kv = keyValue[frame].begin(); kv != keyValue[frame].end(); kv++) {
            
            vector<NumSequence::num_t> numKey;
            cnc.convert(kv->first, numKey);
            
            // get numeric convert of key
            NumSequence ns (numKey);
            
            size_t wordIndex = 0;       // contains the index of the current key
            
            // get index of key in array
            for (size_t i = 0; i < ns.size(); i++) {
                wordIndex <<= elementEncodingSize;      // create space at lower bits for new element
                wordIndex += ns[i];
            }
            
            this->jointProbs[frame][ns.size()-1][wordIndex] = kv->second;
        }
    }
    
    // conditional probabilities
    for (size_t frame = 0; frame < keyValue.size(); frame++) {
        this->model[frame] = jointProbs[frame][this->order];
        jointToMarkov(this->model[frame]);
    }
    
    // derive remaining join probabilities for 'order-1'
    for (size_t frame = 0; frame < keyValue.size(); frame++) {
        for (unsigned o = this->order; o > 0; o--) {
              this->getLowerOrderJoint(o, jointProbs[frame][o], jointProbs[frame][o-1]);
        }
    }
    
}

void CodingMarkov::addPseudocounts(int pcount) {
    
    size_t wordLength = order+1;
    
    for (size_t p = 0; p < period; p++) {                       // for each period
        for (size_t n = 0; n < model[p].size(); n++) {          // for each word
            
            // get sequence represented by current index
            vector<NumSequence::num_t> numSeq;
            this->Markov::indexToNumSequence(n, wordLength, numSeq);
            
            // if not STOP codon, add pseudocounts
            if (!geneticCode->isStop(numSeq))
                model[p][n] += pcount;                              // add pseudocount
        }
    }
    
}
