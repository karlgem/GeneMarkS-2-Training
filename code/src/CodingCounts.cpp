//
//  CodingCounts.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "CodingCounts.hpp"

using namespace gmsuite;

CodingCounts::CodingCounts(unsigned order, size_t period, const NumAlphabetDNA &alph, const NumGeneticCode &geneticCode) : PeriodicCounts(order, period, alph) {
    this->geneticCode = &geneticCode;
}

// update counts
void CodingCounts::updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation, bool reverseComplement) {
    
    size_t seqLen = distance(begin, end);       // sequence length
    
    vector<NumSequence::num_t> codon (3);
    
    // check if contains stop codon
    for (size_t n = 0; n < seqLen-3; n+=3) {
        
        // set new codon
        codon[0] = *(begin+n);
        codon[1] = *(begin+n+1);
        codon[2] = *(begin+n+2);
        
        if (reverseComplement) {
            NumSequence numSeq (codon);
            numSeq.reverseComplement(*this->alphabet->getCNC());
            codon = vector<NumSequence::num_t> (numSeq.begin(), numSeq.end());
        }
        
        if (this->geneticCode->isStop(codon))
            throw std::logic_error("Codons cannot be STOP codons.");
    }
    
    this->PeriodicCounts::updateCounts(begin, end, operation, reverseComplement);
    
}
