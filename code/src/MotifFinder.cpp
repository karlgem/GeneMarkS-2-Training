//
//  MotifFinder.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "MotifFinder.hpp"

#include <float.h>              // DBL_MAX
#include <stdlib.h>             // rand
#include <iostream>
#include <algorithm>
#include "CountModels.hpp"
#include "CountModelsV1.hpp"
#include "UnivariatePDF.hpp"
#include "NumAlphabetDNA.hpp"
#include "ProbabilityModels.hpp"
#include "ProbabilityModelsV1.hpp"

using namespace std;
using namespace gmsuite;


// Constructor: create a motif-finder class to find motifs in sequences
MotifFinder::MotifFinder(
                         unsigned width,
                         unsigned motifOrder,
                         unsigned backOrder,
                         double pcounts,
                         MFinderModelParams::align_t align,
                         unsigned tries,
                         unsigned maxIter,
                         unsigned maxEMIter,
                         unsigned shiftEvery) {
    
    this->width = width;
    this->motifOrder = motifOrder;
    this->backOrder = backOrder;
    this->pcounts = pcounts;
    this->align = align;
    this->tries = tries;
    this->maxIter = maxIter;
    this->maxEMIter = maxEMIter;
    this->shiftEvery = shiftEvery;
}


// Find motifs in sequences.
void MotifFinder::findMotifs (const vector<NumSequence> &sequences, vector<NumSequence::size_type> &positions) {
    
    // if no sequences, just return
    if (sequences.size() == 0)
        return;
    
    double maxProbability = -DBL_MAX;                   // maximum probability over all tries
    vector<NumSequence::size_type> maxPositions;        // motif position for max-probability try
    
    // run algorithm several times; choose output with maximum probability
    for (unsigned t = 0; t < tries; t++) {
        
        vector<NumSequence::size_type> tempPosition;
        double tempProb = gibbsFinder(sequences, tempPosition);         // run gibbs finder
        
        // if found better alignment, record it
        if (tempProb > maxProbability) {
            maxProbability = tempProb;
            maxPositions = tempPosition;
        }
    }
    
    // get best positions for output
    positions = maxPositions;
    
}



// Run a single try of gibbs-finder to search for the best motif alignment.
double MotifFinder::gibbsFinder(const vector<NumSequence> &sequences, vector<NumSequence::size_type> &positions) {
    
    vector<NumSequence>::size_type numSeqs = sequences.size();            // number of sequences
    
    /***** Initialize random alignment *****/
    
    // start by randomly selecting motif locations in sequences
    vector<Sequence::size_type> tempPositions (numSeqs);
    
    for (vector<NumSequence>::size_type n = 0; n < numSeqs; n++) {
        // get random position between 0 and number of valid motif positions (i.e. consider motif width)
        tempPositions[n] = rand() % (sequences[n].size() - width + 1);
    }
    
    CharNumConverter cnc(&this->alphabet);
    NumAlphabetDNA numAlphabet(this->alphabet, cnc);
    
    
    /***** Construct initial count models *****/
    CountModels* counts = new CountModelsV1(numAlphabet, width, motifOrder, backOrder, align);
    counts->construct(sequences, tempPositions);
    
    // allocate space for probability models
    ProbabilityModels *probs = new ProbabilityModelsV1(numAlphabet, width, motifOrder, backOrder, pcounts, align);
    
    double maxScore = -DBL_MAX;                 // maximum alignment score
    vector<Sequence::size_type> maxPositions;   // maximum alignment positions
    
    
    vector<size_t> shuffled (numSeqs);
    for (size_t k = 0; k < numSeqs; k++)
        shuffled[k] = k;
    
    // run all iterations until maximum iteration number is reached
    for (size_t iter = 0; iter < maxIter; iter++) {
        
        // shuffle indeces to select sequences in random order
        random_shuffle(shuffled.begin(), shuffled.end());
        
        // 1) select a sequence z
        // 2) remove z from counts
        // 3) build models from remaining sequences
        // 4) find new motif location in z
        // 5) add new z info back to counts
        
        for (vector<NumSequence>::size_type k = 0; k < numSeqs; k++) {
            
            Sequence::size_type zIndex = shuffled[k];                                       // select sequence z
            
            counts->decount(sequences[zIndex], tempPositions[zIndex]);                      // remove z from counts
            
            probs->construct(counts);                                                       // build models from remaining sequences counts
            
            tempPositions[zIndex] = probs->samplePosition(sequences[zIndex]);               // find new motif location in z
            
            counts->count(sequences[zIndex], tempPositions[zIndex]);                        // add new z info back to counts
            
        }
        
        // try shifting motifs left and right to find better locations
        if (iter > 0 && iter % shiftEvery == 0) {
            int amountToShift = attemptShift(sequences, tempPositions);                     // try shifting
            
            // if shift successful
            if (amountToShift != 0) {
                vector<Sequence::size_type> shiftedPositions;
                shiftPositions(tempPositions, amountToShift, sequences, shiftedPositions);  // shift positions by amountToShift
                
                tempPositions = shiftedPositions;                                           // assign new positions
                counts->construct(sequences, tempPositions);                                // construct new counts
            }
        }
        
        // compute probabilities based on latest counts
        probs->construct(counts);
        
        // calculate alignment conditional log-likelihood
        double tempScore = probs->computeCLL();
        
        if (tempScore > maxScore) {
            maxScore = tempScore;
            maxPositions = tempPositions;
        }
    }
    
    
    // get best configuration, and construct new counts
    tempPositions = maxPositions;
    counts->construct(sequences, tempPositions);
    
    // perform EM on best configuration
    for (size_t iter = 0; iter < maxEMIter; iter++) {
        
        // shuffle indeces to select sequences in random order
        random_shuffle(shuffled.begin(), shuffled.end());
        
        for (vector<NumSequence>::size_type k = 0; k < numSeqs; k++) {
            NumSequence::size_type zIndex = shuffled[k];                                    // select sequence z
            counts->decount(sequences[zIndex], tempPositions[zIndex]);                      // remove z from counts
            probs->construct(counts);                                                       // build models from remaining sequences counts
            tempPositions[zIndex] = probs->samplePosition(sequences[zIndex], true);         // find new motif location in z (get max since it's EM)
            counts->count(sequences[zIndex], tempPositions[zIndex]);                        // add new z info back to counts
        }
        
        probs->construct(counts);
        
        // calculate alignment conditional log-likelihood
        double tempScore = probs->computeCLL();
        
        if (tempScore > maxScore) {
            maxScore = tempScore;
            maxPositions = tempPositions;
        }
    }
    
    
    
    // free space
    delete counts;
    delete probs;
    
    positions = maxPositions;       // get best positions for output
    return maxScore;                // return score for best alignment
    
}



void MotifFinder::shiftPositions(const vector<NumSequence::size_type> &original, int shiftAmount, const vector<NumSequence> &sequences, vector<NumSequence::size_type> &result) {
    
    if (shiftAmount == 0) {
        result = original;
        return;
    }
    
    vector<NumSequence>::size_type numSeqs = sequences.size();
    result.resize(numSeqs);
    
    for (size_t i = 0; i < numSeqs; i++) {
        
        Sequence::size_type maxMotifPos = sequences[i].size() - width + 1;     // get maximum valid motif position in sequence
        
        // if shift amount negative
        if (shiftAmount < 0) {
            // check lower boundary
            if (original[i] < abs(shiftAmount))
                result[i] = 0;
            else
                result[i] = original[i] - abs(shiftAmount);
        }
        // if shift amount positive
        else if (shiftAmount > 0) {
            // check upper boundary
            if (maxMotifPos-original[i]-1 < shiftAmount)
                result[i] = maxMotifPos-1;
            else
                result[i] = original[i] + shiftAmount;
        }
        else
            result[i] = original[i];
    }
    
}







int MotifFinder::attemptShift(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions) {
    
    int minShift = -2;      // minimum shift amount
    int maxShift = 2;       // maximum shift amount
    
    int numShifts = abs(minShift) + abs(maxShift) + 1;      // number of shifts to perform
    
    CharNumConverter cnc(&this->alphabet);
    NumAlphabetDNA numAlphabet(this->alphabet, cnc);
    
    // hold shift scores (to sample from)
    vector<double> shiftScores (numShifts, 0);
    
    // for every shift
    for (int shiftIdx = 0; shiftIdx < numShifts; shiftIdx++) {
        
        // get shift amount
        int shift = shiftIdx + minShift;
        
        // shift motif positions
        vector<Sequence::size_type> shiftedPositions;
        shiftPositions(positions, shift, sequences, shiftedPositions);
        
        // construct new probability models from shifts
        ProbabilityModels* probs = new ProbabilityModelsV1(numAlphabet, width, motifOrder, backOrder, pcounts, align);
        probs->construct(sequences, shiftedPositions);
        
        // compute shift score
        shiftScores[shiftIdx] = probs->computeCLL();
        
        // free space for probabilities
        delete probs;
    }
    
    
    // build a distribution from log-space scores
    UnivariatePDF scoreDistribution (shiftScores, true);
    
    int shiftAmount = ((int) scoreDistribution.sample()) + minShift;
    return shiftAmount;
    
}



























