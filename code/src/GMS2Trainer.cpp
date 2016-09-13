//
//  GMS2Trainer.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "GMS2Trainer.hpp"

#include "UniformCounts.hpp"
#include "PeriodicCounts.hpp"
#include "NonUniformCounts.hpp"

using namespace std;
using namespace gmsuite;

// default constructor
GMS2Trainer::GMS2Trainer() {
    
}

void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    PeriodicCounts counts (codingOrder, 3, alph);
    
    // get counts for 3 period markov model given order
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left;                // get left position of fragment
        size_t right = (*iter)->right;              // get right position of fragment
        size_t length = left + right + 1;           // compute fragment length
        
        counts.count(sequence.begin()+left, sequence.begin() + left +length);
    }
    
    // convert counts to probabilities
    coding = new PeriodicMarkov(codingOrder, 3, alph);
    coding->construct(&counts, pcounts);
    
}

// this function assumes labels are sorted by "left" in increasing order
void GMS2Trainer::estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    UniformCounts counts(noncodingOrder, alph);
    
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
    
    // convert counts to probabilities
    noncoding = new UniformMarkov(noncodingOrder, alph);
    noncoding->construct(&counts, pcounts);
}


void GMS2Trainer::estimateParametersStartContext(const NumSequence &sequence, const vector<Label *> &labels) {
    AlphabetDNA alph;
    NonUniformCounts counts(startContextOrder, startContextLength, alph);
    
    // get counts for start context model
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left + 3;                // get left position of start context fragment
        
        // skip sequence if it doesn't have enough nucleotides for start context
        if (left + startContextLength > (*iter)->right)
            continue;
        
        counts.count(sequence.begin() + left, sequence.begin() + left + startContextLength);
    }
    
    // convert counts to probabilities
    startContext = new NonUniformMarkov(startContextOrder, startContextLength, alph);
    startContext->construct(&counts, pcounts);
}


void GMS2Trainer::estimateParameters(const NumSequence &sequence, const vector<gmsuite::Label *> &labels) {
    
    // TODO: sort labels by 'left' in increasing order
    
    // estimate parameters for coding model
    estimateParamtersCoding(sequence, labels);
    
    // estimate parameters for noncoding model
    estimateParamtersNonCoding(sequence, labels);
    
    // estimate parameters for start context
    estimateParametersStartContext(sequence, labels);
    
    // estimate parameters for motif models
}