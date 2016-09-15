//
//  GMS2Trainer.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "GMS2Trainer.hpp"
#include <iostream>
#include "MotifFinder.hpp"
#include "UnivariatePDF.hpp"
#include "UniformCounts.hpp"
#include "PeriodicCounts.hpp"
#include "NonUniformCounts.hpp"
#include "SequenceParser.hpp"

using namespace std;
using namespace gmsuite;

// default constructor
GMS2Trainer::GMS2Trainer() {
    
}

void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    PeriodicCounts counts (codingOrder, 3, alph, *this->cnc);
    
    // get counts for 3 period markov model given order
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left;                // get left position of fragment
        size_t right = (*iter)->right;              // get right position of fragment
        size_t length = right - left + 1;           // compute fragment length
        
        bool reverseComplement = (*iter)->strand == Label::NEG;
        
        counts.count(sequence.begin()+left, sequence.begin() + left +length, reverseComplement);
    }
    
    // convert counts to probabilities
    coding = new PeriodicMarkov(codingOrder, 3, alph, *this->cnc);
    coding->construct(&counts, pcounts);
    
}

// this function assumes labels are sorted by "left" in increasing order
void GMS2Trainer::estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    UniformCounts counts(noncodingOrder, alph, *this->cnc);
    
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
    noncoding = new UniformMarkov(noncodingOrder, alph, *this->cnc);
    noncoding->construct(&counts, pcounts);
}


void GMS2Trainer::estimateParametersStartContext(const NumSequence &sequence, const vector<Label *> &labels) {
    AlphabetDNA alph;
    NonUniformCounts counts(startContextOrder, startContextLength, alph, *this->cnc);
    
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
    startContext = new NonUniformMarkov(startContextOrder, startContextLength, alph, *this->cnc);
    startContext->construct(&counts, pcounts);
}


void GMS2Trainer::estimateParametersMotifModel(const NumSequence &sequence, const vector<Label *> &labels) {
    
    // build motif finder
    MotifFinder::Builder b;
    if (optionsMFinder->align == "left")
        b.setAlign(MFinderModelParams::LEFT);
    else if (optionsMFinder->align == "right")
        b.setAlign(MFinderModelParams::RIGHT);
    else
        b.setAlign(MFinderModelParams::NONE);
    
    b.setWidth(optionsMFinder->width);
    b.setMaxIter(optionsMFinder->maxIter);
    b.setPcounts(optionsMFinder->pcounts);
    b.setNumTries(optionsMFinder->tries);
    b.setBackOrder(optionsMFinder->bkgdOrder);
    b.setMaxEMIter(optionsMFinder->maxEMIter);
    b.setMotifOrder(optionsMFinder->motifOrder);
    b.setShiftEvery(optionsMFinder->shiftEvery);
    
    
    MotifFinder mfinder = b.build();
    
    // if genome is class 1, search for RBS
    if (genomeClass == ProkGeneStartModel::C1) {
        
        // extract upstream of each label
        vector<NumSequence> upstreams;
        SequenceParser::extractUpstreamSequences(sequence, labels, *cnc, upstreamLength, upstreams);
        
        vector<NumSequence::size_type> positions;
        mfinder.findMotifs(upstreams, positions);
        
        // build RBS model
        AlphabetDNA alph;
        NonUniformCounts rbsCounts(optionsMFinder->motifOrder, optionsMFinder->width, alph, *this->cnc);
        for (size_t n = 0; n < upstreams.size(); n++) {
            rbsCounts.count(upstreams[n].begin()+positions[n], upstreams[n].begin()+positions[n]+optionsMFinder->width);
        }
        
        rbs = new NonUniformMarkov(optionsMFinder->motifOrder, optionsMFinder->width, alph, *this->cnc);
        rbs->construct(&rbsCounts, optionsMFinder->pcounts);
        
        // build spacer distribution
        // build histogram from positions
        vector<double> positionCounts (upstreamLength - optionsMFinder->width+1, 0);
        for (size_t n = 0; n < positions.size(); n++) {
            // FIXME account for LEFT alignment
            // below is only for right
            positionCounts[upstreamLength - optionsMFinder->width - positions[n]]++;        // increment position
        }
        
        rbsSpacer = new UnivariatePDF(positionCounts, false, pcounts);
        
    }
    // if genome is class 2, search for weak RBS, and estimate upstream signature pwm
    else if (genomeClass == ProkGeneStartModel::C2) {
        throw logic_error("Code not yet completed.");
    }
    // if genome is class 3, search for RBS and promoter
    else {
        throw logic_error("Code not yet completed.");
    }
    
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
