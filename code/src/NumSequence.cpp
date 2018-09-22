//
//  NumSequence.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NumSequence.hpp"
#include <stdexcept>

using std::invalid_argument;
using namespace gmsuite;

const NumSequence::size_type NumSequence::npos;

// Default constructor: create an empty numeric sequence.
NumSequence::NumSequence() {
    
}

// Constructor: create a numeric sequence from a regular sequence object.
NumSequence::NumSequence(const Sequence &sequence, const CharNumConverter &converter ) {
    
    converter.convert(sequence.begin(), sequence.end(), this->numSeq);
}

// Constructor: create a numeric sequence from vector of num_t elements
NumSequence::NumSequence(const vector<num_t> &numSequence) {
    this->numSeq = numSequence;
}

// Access an element from a const numeric sequence (i.e. cannot be modified)
const NumSequence::num_t&  NumSequence::operator[](size_type idx) const {
    return numSeq[idx];
}

// access an element
NumSequence::num_t& NumSequence::operator[](size_type idx) {
    return numSeq[idx];
}

// get subsequence
NumSequence NumSequence::subseq(size_type n, size_type length) const {
    
    if (n >= numSeq.size())
        throw invalid_argument("Input n should be less than sequence length.");
    
    if (n + length > numSeq.size())
        throw invalid_argument("Input n+length should be less than or equal to sequence length.");
    
    NumSequence sub;
    sub.numSeq = vector<num_t>(this->numSeq.begin() + n, this->numSeq.begin() + n + length);
    return sub;
}

// reverse complement in-place
void NumSequence::reverseComplement(const CharNumConverter &cnc) {
    
    // if size zero, nothing to do :)
    if (numSeq.size() == 0)
        return;
    
    size_type left = 0;
    size_type right = numSeq.size()-1;
    
    // for each position < N/2, swap the symmetric 'left' and 'right' elements, and complement each
    // if N is even, the final left will be right before the final right
    // if N is odd, an extra elements remains that needs to be complemented and "swapped" with itself
    for (size_t n = 0; n < numSeq.size()/2; n++) {
        num_t temp = numSeq[left];
        numSeq[left] = cnc.complement(numSeq[right]);   // exchange and complement
        numSeq[right] = cnc.complement(temp);           // exchange and complement
        
        left++;
        right--;
    }
    
    // if N is odd, complement the floor(N/2)'th element
    if (numSeq.size() % 2 != 0)
        numSeq[numSeq.size()/2] = cnc.complement(numSeq[numSeq.size()/2]);
}



// begin iterator
NumSequence::iterator NumSequence::begin() {
    return this->numSeq.begin();
}

// end iterator
NumSequence::iterator NumSequence::end() {
    return this->numSeq.end();
}


// begin const_iterator
NumSequence::const_iterator NumSequence::begin() const {
    return this->numSeq.begin();
}

// end const_iterator
NumSequence::const_iterator NumSequence::end() const {
    return this->numSeq.end();
}

// sequence size
NumSequence::size_type NumSequence::size() const {
    return numSeq.size();
}




bool NumSequence::containsInvalid(const NumAlphabetDNA &alph) const {
    for (size_t n = 0; n < numSeq.size(); n++) {
        if (alph.isAmbiguous(numSeq[n]))
            return true;
    }
    
    return false;
}



NumSequence NumSequence::operator+ (const NumSequence &op2) {
    NumSequence result;
    result.numSeq = this->numSeq;
    result.numSeq.insert(result.numSeq.end(), op2.numSeq.begin(), op2.numSeq.end());     // append op2 vector
    return result;
}
