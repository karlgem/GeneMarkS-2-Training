//
//  NonCodingMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 11/3/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NonCodingMarkov.hpp"
#include <math.h>

using namespace std;
using namespace gmsuite;


NonCodingMarkov::NonCodingMarkov(const vector<pair<string, double> > &keyValue, const NumAlphabetDNA &alph, const CharNumConverter &cnc) : UniformMarkov(0, alph) {
    
    
    unsigned maxOrder = 0;
    
    // loop through keyValue pairs and get length of maximum order
    for (size_t n = 0; n < keyValue.size(); n++) {
        if (keyValue[n].first.size() > 0)
            if (maxOrder < keyValue[n].first.size() - 1)
                maxOrder = (unsigned) keyValue[n].first.size() - 1;
    }
    
    // reset order and reallocate space
    this->order = maxOrder;
    initialize();
    
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    
    // Joint Probs
    // for each key-value pair, place it in the correct position
    for (vector<pair<string, double> >::const_iterator kv = keyValue.begin(); kv != keyValue.end(); kv++) {
        
        vector<NumSequence::num_t> numKey;
        cnc.convert(kv->first, numKey);
        
        // get numeric convert of key
        NumSequence ns (numKey);
        
        size_t wordIndex = 0;       // contains the index of the current key
        
        // get index of key in array
        for (size_t i = 0; i < ns.size(); i++) {
            wordIndex <<= elementEncodingSize;          // create space at lower bits for new element
            wordIndex += ns[i];
        }
        
        // set the joint probability
        this->jointProbs[ns.size()-1][wordIndex] = kv->second;
    }
    
    // Conditional Probs
    this->model = jointProbs[this->order-1];        // get join of highest order
    jointToMarkov(this->model);                     // convert joint to conditional
}
