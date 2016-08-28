//
//  GMS2Trainer.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "GMS2Trainer.hpp"

#include "PeriodicCounts.hpp"

using namespace std;
using namespace gmsuite;

// default constructor
GMS2Trainer::GMS2Trainer() {
    
}

void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    PeriodicCounts counts (2, 3, &alph);
    
    // get counts for 3 period markov model given order
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left;                // get left position of fragment
        size_t right = (*iter)->right;              // get right position of fragment
        size_t length = left + right + 1;           // compute fragment length
        
        counts.count(sequence.begin()+left, sequence.begin()+length);
    }
    
    // TODO: convert counts to probabilities
    
}

// this function assumes labels are sorted by "left" in increasing order
void GMS2Trainer::estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    PeriodicCounts counts (2, 1, &alph);
    
    size_t leftNoncoding = 0;       // left position of current noncoding region
    
    // get counts for 3 period markov model given order
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left;                // get left position of fragment
        size_t right = (*iter)->right;              // get right position of fragment
        
        if (leftNoncoding < left)
            counts.count(sequence.begin() + leftNoncoding, sequence.begin() + left);
        
        // update left position of (possible) non-coding region after current gene
        leftNoncoding = right+1;
    }
    
    // TODO: convert counts to probabilities
    
}


void GMS2Trainer::estimateParameters(const NumSequence &sequence, const vector<gmsuite::Label *> &labels) {
    
    
    // estimate parameters for coding model
    estimateParamtersCoding(sequence, labels);
    
    // estimate parameters for noncoding model
    estimateParamtersNonCoding(sequence, labels);
    
    // estimate parameters for start context
    
    // estimate parameters for motif models
}