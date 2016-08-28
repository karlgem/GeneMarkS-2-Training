//
//  PeriodicMarkov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "PeriodicMarkov.hpp"
#include "PeriodicCounts.hpp"

#include <stdexcept>

using std::invalid_argument;
using namespace gmsuite;

// Constructor:
PeriodicMarkov::PeriodicMarkov(unsigned order, size_t period, const AlphabetDNA* alph) : Markov(order, alph) {
    this->period = period;
    initialize();
}


// Construct the model probabilities from a list of sequences
void PeriodicMarkov::construct(const vector<NumSequence> &sequences, int pcount) {
    
    // get counts
    PeriodicCounts counts (order, period, alphabet);
    counts.construct(sequences);
    
    // construct probabilities from counts
    construct(&counts, pcount);
}

// Construct the model probabilities from existing counts.
void PeriodicMarkov::construct(const Counts* counts, int pcount) {
    
    // counts cannot be NULL
    if (counts == NULL)
        throw invalid_argument("Counts cannot be NULL.");
    
    // counts alphabet must match markov alphabet
    if (counts->getAlphabet() != this->alphabet)
        throw invalid_argument("Counts alphabet must match Markov alphabet.");
    
    // counts order must match Markov order
    if (counts->getOrder() != this->order)
        throw invalid_argument("Counts order must match Markov order.");
    
    // cast counts to PeriodiCounts
    const PeriodicCounts* periodicCounts = dynamic_cast<const PeriodicCounts*>(counts);
    
    if (periodicCounts == NULL)
        throw invalid_argument("Counts should have type 'PeriodicCounts'.");
    
    // counts period must match Markov period
    if (periodicCounts->getPeriod() != this->period)
        throw invalid_argument("Counts period must match Markov period.");
    
    // start by copying counts
    this->model = periodicCounts->model;
    
    vector<double> sums (period, 0);            // will contain sum of counts for each period
    
    // add pseudocounts
    for (size_t p = 0; p < period; p++)                         // for each period
        for (size_t n = 0; n < model[p].size(); n++) {          // for each word
            model[p][n] += pcount;                              // add pseudocount
            sums[p] += model[p][n];                             // add to sum
        }
    
    
    // normalize to get joint probabilities (e.g. P(ACG))
    for (size_t p = 0; p < period; p++)                         // for each period
        for (size_t n = 0; n < model[p].size(); n++)            // for each word
            if (sums[p] != 0)                                   // check for division by zero
                model[p][n] /= sums[p];                         // normalize word counts on sum
    
    
    // for each period, convert joint probabilities to Markov( conditional):
    // e.g. P(ACG) -> P(G|AC)
    for (size_t p = 0; p < period; p++)
        jointToMarkov(model[p]);
    
}

// Compute the score of a sequence using the model probabilities
double PeriodicMarkov::evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog) const {
    return 0;
}

// Generate a string representation of the model
string PeriodicMarkov::toString() const {
    return "";
}


void PeriodicMarkov::convertCountsToProbability() {
    
}

// Initialize the model by allocating space, setting the keys, and setting counts to 0
void PeriodicMarkov::initialize() {
    
}

// reset counts to zero
void PeriodicMarkov::resetCounts() {
    
}


void PeriodicMarkov::jointToMarkov(vector<double> &probs) {
    
    // Problem Context:
    // the probs vector holds element in the following order: e.g.
    //  index   |   element
    //     0    |     A1 A2
    //     1    |     A1 C2
    //     2    |     A1 G2
    //     3    |     A1 T2
    //     4    |     C1 A2
    //     5    |     C1 C2
    //     6    |     C1 G2
    //     7    |     C1 T2
    //     8    |     ...
    //
    // The numbers near the letters indicate the position of the letter in the dinucleotide. E.g.
    // A1C2 is AC, whereas A2C1 == C1A2 == CA
    //
    // To get conditional probabilities, note that
    // P(X2 | X1) = P(X1, X2) / sum_{x} P(X1, x2)           for all values of x (e.g. A,C,G,T)
    //
    // This generalizes to higher orders, where
    // P(Xn | X1,...,Xn-1) = P(X1,...,Xn) / sum_{x} P(X1,...,Xn-1,xn)
    
    // Define a block to be the number of (consecutive) elements contributing
    // to a denominator. E.g. in the above example, AA,AC,AG,AT all contribute
    // to the same denominator sum_{x} P(A1,x2). So the block size is 4
    size_t blockSize = alphabet->sizeValid();
    size_t numOfBlocks = probs.size() / blockSize;
    
    // probs.size() should be divisible by blockSize
    if (probs.size() % blockSize != 0)
        throw std::logic_error("Probs.size() should be divisible by blockSize.");
    
    
    // Start by pre-computing the denominators (i.e. sum_{x} P(X1,...,Xn-1,xn))
    vector<double> denominators (numOfBlocks);
    
    // compute one denominator for every block
    for (size_t b = 0; b < numOfBlocks; b++) {
        
        // get starting index of current block in probs vector
        size_t idx = b * blockSize;
        
        // loop over block and sum probabilities
        for (size_t i = 0; i < blockSize; i++)
            denominators[b] += probs[idx++];
    }
    
    // Now that denominators have been computed, update the probs vector by
    // dividing each element on its corresponding denominator
    for (size_t n = 0; n < probs.size(); n++) {
        
        // get block index that the n'th element belongs to
        size_t blockIdx = (n + n % blockSize) / blockSize;         // derived from: n = blockIdx * blockSize + remainder
        
        if (denominators[blockIdx] != 0)
            probs[n] /= denominators[blockIdx];                     // normalize
    }
    
    

}























