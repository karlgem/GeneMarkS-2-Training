//
//  UniformMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "UniformMarkov.hpp"
#include "UniformCounts.hpp"

#include <stdexcept>

using namespace std;
using namespace gmsuite;

// Constructor:
UniformMarkov::UniformMarkov(unsigned order, const AlphabetDNA* alph) : Markov(order, alph) {
    initialize();
}


// Construct the model probabilities from a list of sequences
void UniformMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    // get counts
    UniformCounts counts(order, alphabet);
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
    
    // convert joint probabilities to Markov (conditional)
    // e.g. P(ACG) -> (G|AC)
    jointToMarkov(model);
    

    
}

// Compute the score of a sequence using the model probabilities
double UniformMarkov::evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog) const {
    return 0;
}

// Generate a string representation of the model
string UniformMarkov::toString() const {
    return "";
}


// Initialize the model by allocating space, setting the keys, and setting counts to 0
void UniformMarkov::initialize() {
    size_t numElements = alphabet->sizeValid();         // the number of eleents that can make up valid words (e.g. A,C,G,T)
    
    size_t wordSize = (order+1);                        // size of a word
    size_t numWords = numElements << wordSize;          // number of possible words of size 'wordSize' with the given alphabet
    
    model.resize(numWords, 0);
}








