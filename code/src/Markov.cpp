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
Markov::Markov(unsigned order, const AlphabetDNA &alph, const CharNumConverter &cnc)  {
    
    this->order = order;
    this->alphabet = &alph;
    this->cnc = &cnc;
    
}

// Get the model's order
unsigned Markov::getOrder() const {
    return order;
}

// Get the model's alphabet
const AlphabetDNA* Markov::getAlphabet() const {
    return alphabet;
}




void Markov::jointToMarkov(vector<double> &probs) const {
    
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


void Markov::getHigherOrderJoint(unsigned currentOrder, const vector<double> &current, unsigned newOrder, vector<double> &result) const {
    if (newOrder < currentOrder)
        throw std::invalid_argument("New order is less than current order.");
    
    if (newOrder == currentOrder) {
        result = current;
        return;
    }
    
    // Problem Context:
    // the probability vector holds elements in the following order:
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
    // To get the joint probabilities of higher order (e.g. from 1 to 3), this means that
    // P(WXYZ) = P(YZ) for all W,X in alphabet
    //
    // Define a block to be the number of  elements that separate the two consecutive P(ABYZ) and P(CDYZ)
    // in the matrix.
    size_t currentWordLength = currentOrder+1;
    size_t blockSize = pow(alphabet->sizeValid(), currentWordLength);
    
    
    
    // assign probability for higher order
    size_t newWordLength = newOrder+1;
    result.resize(pow(alphabet->sizeValid(), newWordLength));
    
    size_t numOfBlocks = result.size() / blockSize;         // number of blocks in result
    // result.size() should be divisible by blockSize
    if (result.size() % blockSize != 0)
        throw std::logic_error("Result.size() should be divisible by blockSize");
    
    // for every word in existing model
    for (size_t currentIdx = 0; currentIdx < current.size(); currentIdx++) {
        
        // assign it to every corresponding higher order word by jumping through blocks
        size_t newIdx = currentIdx;                     // newIdx always starts at currentIdx
        for (size_t b = 0; b < numOfBlocks; b++) {
            
            result[newIdx] = current[currentIdx];
            
            newIdx += blockSize;
        }
    }
}


void Markov::indexToNumSequence(size_t idx, size_t wordLength, vector<NumSequence::num_t> &numSeq) const {
    
    numSeq.resize(wordLength);      // allocate space
    
    size_t numElements = alphabet->sizeValid();             // number of elements that can make up valid words (e.g. A,C,G,T)
    size_t elementEncodingSize = ceil(log2(numElements));   // number of bits required to encode all elements (e.g. A,C,G,T require 2 bits)
    size_t mask = 0;                        // used to clear wordIndex of junk and capture the relevant bits for the index
    
    // get a mask that captures one letter at a time
    for (size_t n = 0; n < elementEncodingSize; n++) {
        mask <<= 1;      // shift by one
        mask |= 1;      // set rightmost bit to 1
    }
    
    // for each character of the word
    for (size_t n = 0; n < wordLength; n++) {
        
        size_t current = idx & mask;            // 'extract' this character from idx
        idx >>= elementEncodingSize;            // remove this character from idx
        
        numSeq[wordLength-n-1] = (NumSequence::num_t) current;       // set current character in numeric sequence
    }
}


NumSequence Markov::indexToNumSequence(size_t idx, size_t wordLength) const {
    vector<NumSequence::num_t>  numSeq;
    
    indexToNumSequence(idx, wordLength, numSeq);
    
    return NumSequence(numSeq);
}




void Markov::incrementOrderByOne(unsigned currentOrder, const vector<double> &currentProbs, vector<double> &newProbs) const {
    
    size_t numElements = alphabet->sizeValid();                 // number of (valid) elements in alphabet (e.g. A,C,G,T)
    
    // set the size of the new probability space
    size_t newOrder = currentOrder+1;
    size_t wordSize = newOrder+1;                   // number of elements that make up a word
    size_t numWords = pow(numElements, wordSize);   // number of words in the new probability space
    newProbs.resize(numWords, 0);                   // allocate space for set new probability
    
    // for each key of the current probabilities, the probability value needs to be copied
    // to all keys of "higher" order. Example, since our working representation is in bits,
    // suppose that the element encoding size is 2 (i.e. 2 bits needed to encode a single element).
    // If A = 00, C = 01, G = 10, T = 11, then the key GC=1001
    // In higher order, we have AGC=001001, CGC=011001, GGC=101001, TGC=111001
    // So we bit representations of intergers from 0 up to alphabet size (i.e. 4 in this case),
    // and we add them to GC, but shifted by 2 * element encoding size.
    //
    // Ex: To get TGC:
    // Representation of T:         11
    // Shift T by 2 * 2:        110000
    // Representation of GC:      1001
    // Add shifted T to GC:     111001
    
    size_t elementEncodingSize = ceil(log2(numElements));       // number of bits to encode an element (e.g. 2 bits for A,C,G,T)
    size_t bitsToShift = elementEncodingSize * currentOrder+1;
    
    // for each key in current probability space
    for (size_t key = 0; key < currentProbs.size(); key++) {
        
        // for each element to be added to the current key
        for (size_t e = 0; e < numElements; e++) {
            // shift the element by the number of letters already in the key, times the element encoding size
            size_t newKey = e << bitsToShift;
            
            // add existing key to new partial key; this gives the new full key
            newKey += key;
            
            // set probability value of new key to that of old key
            newProbs[newKey] = currentProbs[key];
        }
    }
}



void Markov::getCDFPerConditional(unsigned currentOrder, vector<double> &probs) const {
    
    // Problem Context:
    // the probs vector holds elements in the following order: e.g.
    //  index   |   element
    //     0    |     A1 A2 (i.e. A2 | A1) ___
    //     1    |     A1 C2 (i.e. C2 | A1)    |
    //     2    |     A1 G2 (i.e. G2 | A1)    |--> sums to 1
    //     3    |     A1 T2 (i.e. T2 | A1) ___|
    //     4    |     C1 A2 (i.e. A2 | C1) ___
    //     5    |     C1 C2 (i.e. C2 | C1)    |
    //     6    |     C1 G2 (i.e. G2 | C1)    |--> sums to 1
    //     7    |     C1 T2 (i.e. T2 | C1) ___|
    //     8    |     ...
    //
    // The numbers near the letters indicate the position of the letter in the dinucleotide. E.g.
    // A1C2 is AC, whereas A2C1 == C1A2 == CA
    //
    // To get conditional probabilities CDF, note that we just sum consecutive elements
    // that exists with a block (that sums to 1)
    //
    // Define a block to be the number of (consecutive) elements with the same base (that sum to 1)
    // E.g. in the above example, AA,AC,AG,AT all contribute to a single CDF. So the block size is 4.
    size_t blockSize = alphabet->sizeValid();
    size_t numOfBlocks = probs.size() / blockSize;
    
    // probs.size() should be divisible by blockSize
    if (probs.size() % blockSize != 0)
        throw std::logic_error("Probs.size() should be divisible by blockSize.");
    
    
    // sum for every block
    for (size_t b = 0; b < numOfBlocks; b++) {
        
        // get starting index of current block in probs vector
        size_t idx = b * blockSize + 1;         // +1 because first element in a block doesn't change
        
        // loop over block and sum probabilities
        for (size_t i = 1; i < blockSize; i++) {
            probs[idx] += probs[idx-1];         // accumulate sum from previous
            idx++;
        }
    }
}





double KLDivergence(const Markov* P, const Markov* Q) {
    double kl = 0;
    
    
    
    return kl;
}



























