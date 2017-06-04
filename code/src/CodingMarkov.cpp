//
//  CodingMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "CodingMarkov.hpp"

using namespace gmsuite;

CodingMarkov::CodingMarkov(unsigned order, size_t period, const NumAlphabetDNA &alph, const NumGeneticCode &geneticCode) : PeriodicMarkov(order, period, alph) {
    this->geneticCode = &geneticCode;
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
