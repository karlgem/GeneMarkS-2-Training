//
//  MotifMarkov.cpp
//  Biogem CPP
//
//  Created by Karl Gemayel on 8/20/18.
//  Copyright © 2018 Georgia Institute of Technology. All rights reserved.
//

#include "MotifMarkov.hpp"
#include <math.h>

using namespace std;
using namespace gmsuite;

MotifMarkov::MotifMarkov(const vector<vector<pair<string, double> > > &keyValue, size_t length, const NumAlphabetDNA &alph, const CharNumConverter &cnc) : NonUniformMarkov(0, length, alph) {
    
    unsigned maxOrder = 0;
    
    // loop through key value pairs and get length of maximum order
    for (size_t n = 0; n < keyValue.size(); n++) {
        for (size_t m = 0; m < keyValue[n].size(); m++) {
            if (keyValue[n][m].first.size() > 0) {
                if (maxOrder < keyValue[n][m].first.size() - 1)
                    maxOrder = (unsigned) keyValue[n][m].first.size() - 1;
            }
        }
    }
    
    // reset order and reallocate space
    this->order = maxOrder;
    initialize();
    
    size_t numElements = alphabet->sizeValid();     // number of elements that can make up valid words (e.g. ACGT)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. ACGT requires 2 bits)
    
    
    size_t position = 0;
    // Join Probs
    // for each key-value pair, place it in the correct position
    for (vector<vector<pair<string, double> > >::const_iterator kv = keyValue.begin(); kv != keyValue.end(); kv++) {
        
        // loop over each element in current position
        for (vector<pair<string, double> >::const_iterator kvPos = kv->begin(); kvPos != kv->end(); kvPos++) {
            
            vector<NumSequence::num_t> numKey;
            cnc.convert(kvPos->first, numKey);      // numeric version of key
            
            NumSequence ns (numKey);
            
            size_t wordIndex = 0;       // contains the index of the current key
            
            // get index of key in array
            for (size_t i = 0; i < ns.size(); i++) {
                wordIndex <<= elementEncodingSize;      // create space at lower bits for new element
                wordIndex += ns[i];
            }
            
            // set the joint probability
            this->jointProbs[position][wordIndex] = kvPos->second;
            
        }
        position++;
    }
    
    // conditional probabilities
    this->model = jointProbs;
    
    for (size_t p = 0; p < length; p++)
        jointToMarkov(model[p]);
    
}
