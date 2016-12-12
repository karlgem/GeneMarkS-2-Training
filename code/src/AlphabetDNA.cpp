//
//  AlphabetDNA.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "AlphabetDNA.hpp"
#include <algorithm>

using namespace gmsuite;

/***** Function Prototypes *****/

vector<char> getValidAlphabet();                // get valid alphabet
map<char, char> getComplementAlphabet();        // get complement of valid alphabet
vector<char> getAmbiguousAlphabet();            // get ambiguous alphabet


/***** Function Implementations *****/


AlphabetDNA::AlphabetDNA() {
    
    this->alphValid = getValidAlphabet();                       // get valid alphabet
    this->alphAmbiguous = getAmbiguousAlphabet();               // get ambiguous alphabet
    
    // concatenate them into full alphabet, starting with valid
    this->alph = this->alphValid;
    this->alph.insert(this->alph.end(), alphAmbiguous.begin(), alphAmbiguous.end());
    
    this->complementDNA = getComplementAlphabet();              // get DNA complement alphabet
}


// Get the total number of characters in the alphabet.
AlphabetDNA::size_type AlphabetDNA::size() const {
    return alph.size();
}

// Get the total number of valid characters.
AlphabetDNA::size_type AlphabetDNA::sizeValid() const {
    return alphValid.size();
}

// Get the total number of ambiguous characters.
AlphabetDNA::size_type AlphabetDNA::sizeInvalid() const {
    return alphAmbiguous.size();
}

// Check if alphabet contains a character
bool AlphabetDNA::contains(char c) const {
    return find(alph.begin(), alph.end(), c) != alph.end();
}

// Check if character is valid (i.e. non-ambiguous)
bool AlphabetDNA::isValid(char c) const {
    return find(alphValid.begin(), alphValid.end(), c) != alphValid.end();
}


// Check if character is ambiguous (i.e. N, R, S, ...)
bool AlphabetDNA::isAmbiguous(char c) const {
    return find(alphAmbiguous.begin(), alphAmbiguous.end(), c) != alphAmbiguous.end();
}


// complement a character
char AlphabetDNA::complement(char c) const {
    map<char, char>::const_iterator it = complementDNA.find(c);     // find character complement
    
    // if character not found, return it as it is
    if (it == complementDNA.end())
        return c;
    
    // otherwise, return the complement
    else
        return it->second;
}

// complement a character
string AlphabetDNA::reverseComplement(const string &original) const {
    string result;
    result.resize(original.size());
    
    for (size_t n = 0; n < original.size(); n++) {
        
        size_t posInResult = result.size() - n - 1;
        
        char c = original[n];
        
        map<char, char>::const_iterator it = complementDNA.find(c);     // find character complement
        
        // if character not found, return it as it is
        if (it == complementDNA.end())
            result[posInResult] = c;
        
        // otherwise, return the complement
        else
            result[posInResult] = it->second;
    }
    
    return result;
}





// get valid alphabet
vector<char> getValidAlphabet() {
    
    string str = "ACGT";
    
    vector<char> alph (str.size());
    for (size_t i = 0; i < str.size(); i++)
        alph[i] = str[i];
    
    sort(alph.begin(), alph.end());         // make sure it's sorted
    
    return alph;
}

// get complement of valid alphabet
map<char, char> getComplementAlphabet() {
    map<char, char> comp;
    comp['A'] = 'T';
    comp['T'] = 'A';
    comp['C'] = 'G';
    comp['G'] = 'C';
    return comp;
}

// get ambiguous alphabet
vector<char> getAmbiguousAlphabet() {
    string str = "WSMKRYBDHVNZ";
    
    vector<char> alph (str.size());
    for (size_t i = 0; i < str.size(); i++)
        alph[i] = str[i];
    
    sort(alph.begin(), alph.end());         // make sure it's sorted
    
    return alph;
}


























