//
//  Markov.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Markov.hpp"
#include <stdexcept>
#include <math.h>

using std::invalid_argument;
using namespace gmsuite;

// Contructor: Initialize a Markov model with a specific order and alphabet.
Markov::Markov(unsigned order, const AlphabetDNA* alph) {
    
    if (alph == NULL)
        throw invalid_argument("Alphabet cannot be NULL.");
    
    this->order = order;
    this->alphabet = alph;
    
}

// Get the model's order
unsigned Markov::getOrder() const {
    return order;
}

// Get the model's alphabet
const AlphabetDNA* Markov::getAlphabet() const {
    return alphabet;
}




void Markov::jointToMarkov(vector<double> &probs) {
    
    // Problem Context:
    // the probs vector holds elements in the following order: e.g.
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
        size_t blockIdx = (n - n % blockSize) / blockSize;         // derived from: n = blockIdx * blockSize + remainder
        
        if (denominators[blockIdx] != 0)
            probs[n] /= denominators[blockIdx];                     // normalize
    }
}



void Markov::getLowerOrderJoint(unsigned currentOrder, const vector<double> &current, vector<double> &result) const {
    
    if (currentOrder == 0)
        throw invalid_argument("Current order cannot be zero when getting lower order probabilities.");

    // Problem Context:
    // the probs vector holds elements in the following order: e.g.
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
    // To get joint probabilities of lower order, note that
    // P(X1) = sum_{x} P(X1,X2)                                 for all values of x (e.g. A,C,G,T)
    //
    // This generalizes to higher orders, where
    // P(X1,...,Xn-1) = sum_{x} P(X1,...,Xn-1,xn)
    
    // Define a block to be the number of (consecutive) elements contributing
    // to the sum. E.g. in the above example, AA,AC,AG,AT all contribute
    // to the same sum_{x} P(A1,x2). So the block size is 4.
    size_t blockSize = alphabet->sizeValid();
    size_t numOfBlocks = current.size() / blockSize;
    
    // current.size() should be divisible by blockSize
    if (current.size() % blockSize != 0)
        throw std::logic_error("Current.size() should be divisible by blockSize.");

    // Compute the marginal probabilities (i.e. P(X1,...,Xn-1) = sum_{x} P(X1,...,Xn-1,xn))
    result.resize(numOfBlocks,0);
    
    // compute one marginal for every block
    for (size_t b = 0; b < numOfBlocks; b++) {
        
        // get starting index of current block in probs vector
        size_t idx = b * blockSize;
        
        // loop over block and sum probabilities
        for (size_t i = 0; i < blockSize; i++)
            result[b] += current[idx++];
    }
}










