//
//  NumAlphabetDNA.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NumAlphabetDNA.hpp"


#include <algorithm>

using namespace gmsuite;


/***** Function Implementations *****/


NumAlphabetDNA::NumAlphabetDNA(const AlphabetDNA &alph, const CharNumConverter &cnc) {
    
    this->cnc = &cnc;

    
    // get alph
    for (AlphabetDNA::const_iterator iter = alph.begin(); iter != alph.end(); iter++)
        this->alph.push_back(cnc.convert(*iter));
    
    // get valid alphabet
    for (AlphabetDNA::const_iterator iter = alph.beginValid(); iter != alph.endValid(); iter++)
        this->alphValid.push_back(cnc.convert(*iter));
    
    // get invalid alphabet
    for (AlphabetDNA::const_iterator iter = alph.beginInvalid(); iter != alph.endInvalid(); iter++)
        this->alphAmbiguous.push_back(cnc.convert(*iter));
    
    // get DNA complement alphabet
    for (AlphabetDNA::const_iterator iter = alph.beginValid(); iter != alph.endValid(); iter++)
        this->complementDNA[cnc.convert(*iter)] = cnc.convert(alph.complement(*iter));
}


// Get the total number of characters in the alphabet.
NumAlphabetDNA::size_type NumAlphabetDNA::size() const {
    return alph.size();
}

// Get the total number of valid characters.
NumAlphabetDNA::size_type NumAlphabetDNA::sizeValid() const {
    return alphValid.size();
}

// Get the total number of ambiguous characters.
NumAlphabetDNA::size_type NumAlphabetDNA::sizeInvalid() const {
    return alphAmbiguous.size();
}

// Check if alphabet contains a character
bool NumAlphabetDNA::contains(element_t c) const {
    return find(alph.begin(), alph.end(), c) != alph.end();
}

// Check if character is valid (i.e. non-ambiguous)
bool NumAlphabetDNA::isValid(element_t c) const {
    return find(alphValid.begin(), alphValid.end(), c) != alphValid.end();
}


// Check if character is ambiguous (i.e. N, R, S, ...)
bool NumAlphabetDNA::isAmbiguous(element_t c) const {
    return find(alphAmbiguous.begin(), alphAmbiguous.end(), c) != alphAmbiguous.end();
}


// complement a character
NumAlphabetDNA::element_t NumAlphabetDNA::complement(element_t c) const {
    map<element_t, element_t>::const_iterator it = complementDNA.find(c);     // find character complement
    
    // if character not found, return it as it is
    if (it == complementDNA.end())
        return c;
    
    // otherwise, return the complement
    else
        return it->second;
}






const CharNumConverter* NumAlphabetDNA::getCNC() const {
    return cnc;
}



















