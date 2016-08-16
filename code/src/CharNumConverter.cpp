//
//  CharNumConverter.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/8/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "CharNumConverter.hpp"

#include <stdexcept>
#include <iterator>     // std::distance


using namespace gmsuite;
using std::invalid_argument;

// Constructor: create a converter using the given alphabet
CharNumConverter::CharNumConverter(const AlphabetDNA* alphabet) {
    
    if (alphabet == NULL) {
        throw invalid_argument("Alphabet cannot be NULL");
    }
    
    this->alphabet = alphabet;
    
    // create converters
    element_t current = 0;
    
    // iterate over each character in the alphabet
    for (AlphabetDNA::const_iterator iter = this->alphabet->begin(); iter != this->alphabet->end(); iter++) {
        this->charToNum[*iter] = current;           // char to number
        this->numToChar[current] = *iter;           // number to char
        
        current++;
    }
    
    // create complement map
    for (AlphabetDNA::const_iterator iter = this->alphabet->beginValid(); iter != this->alphabet->endValid(); iter++) {
        element_t from = charToNum[*iter];
        element_t to = charToNum[this->alphabet->complement(*iter)];
        
        complementDNA[from] = to;
    }
    
}


// Convert a string sequence to its numeric representation.
void CharNumConverter::convert(const string &str, seq_t &result) const {
    
    result.resize(str.length());        // allocate space for result
    
    for (string::size_type i = 0; i < str.size(); i++) {
        result[i] = charToNum.at(str[i]);       // convert each character to number
    }
}


void CharNumConverter::convert(Sequence::const_iterator begin, Sequence::const_iterator end, seq_t &result) const {
    
    Sequence::const_iterator current = begin;
    size_t seqLength = std::distance(begin, end);           // get sequence length
    
    result.resize(seqLength);        // allocate space for result
    
    for (string::size_type i = 0; i < seqLength; i++) {
        result[i] = charToNum.at(*current);          // convert each character to number
        current++;                                  // next character
    }
}


// Convert a single character to its numeric representation.
CharNumConverter::element_t CharNumConverter::convert(char c) const {
    return charToNum.at(c);
}


// Convert a numeric representation sequence back to a string sequence.
string CharNumConverter::convert(seq_t::const_iterator start, seq_t::const_iterator end) const {
    
    string result;
    result.resize(distance(start, end));
    
    seq_t::size_type l = 0;
    
    // iterate from start to end
    for (seq_t::const_iterator curr = start; curr != end; curr++) {
        result[l++] = numToChar.at(*curr);
    }
    
    return result;
    
}


// Convert a numeric element back to a character.
char CharNumConverter::convert(element_t element) const {
    return numToChar.at(element);
}



// Complement DNA sequence
void CharNumConverter::complement(const seq_t &original, seq_t &result) const {
    result.resize(original.size());
    
    for (size_t n = 0; n < original.size(); n++) {
        result[n] = complementDNA.at(original[original.size() - n - 1]);
    }
}


















