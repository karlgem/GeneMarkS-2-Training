//
//  AlphabetDNA.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef AlphabetDNA_hpp
#define AlphabetDNA_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#include <set>

using std::vector;
using std::string;
using std::set;
using std::map;

namespace gmsuite {
    
    
    /**
     * @class AlphabetDNA
     * @brief AlphabetDNA is the DNA alphabet used in sequences
     *
     * This class holds the alphabet used by DNA sequences, as well as the ambiguous
     * characters that might be found, according to IUPAC notation.
     *
     * The alphabet is composed of two classes of characters:
     *  - Valid: standard characters making up a sequence (A,C,G,T)
     *  - Ambiguous: characters making up ambiguous elements (W,S,M,K,R,Y,B,D,H,V,N,Z)
     */
    class AlphabetDNA {
        
    public:
        
        typedef vector<char>::size_type size_type;              /**< Type for number of characters */
        typedef vector<char>::const_iterator const_iterator;    /**< Const iterator over characters */
        
        
        /**
         * Constructor: initialize DNA alphabet
         *
         * The DNA alphabet is made up of two classes of characters
         *   - Valid: standard characters making up a sequence (A,C,G,T)
         *   - Ambiguous: characters making up ambiguous elements (W,S,M,K,R,Y,B,D,H,V,N,Z)
         */
        AlphabetDNA();
        
        /**
         * Get the total number of characters in the alphabet.
         *
         * @return the number of characters
         */
        size_type size() const;
        
        /**
         * Get the total number of valid characters.
         *
         * @return the number of valid characters.
         */
        size_type sizeValid() const;
        
        /**
         * Get the total number of ambiguous characters.
         *
         * @return the number of ambiguous characters.
         */
        size_type sizeInvalid() const;
        
        /**
         * Check if alphabet contains a character
         *
         * @param c the character to find
         * @return true if c is found; false otherwise.
         */
        bool contains(char c) const;
        
        /**
         * Check if character is valid (i.e. non-ambiguous)
         *
         * @param c the character
         * @return true if c is valid; false otherwise.
         */
        bool isValid(char c) const;
        
        
        /**
         * Check if character is ambiguous (i.e. N, R, S, ...)
         *
         * @param c the character
         * @return true if c is valid; false otherwise.
         */
        bool isAmbiguous(char c) const;
        
        /**
         * Complement a character
         */
        char complement(char c) const;
        
        /**
         * Complement a string
         */
        string reverseComplement(const string &s) const;
        
        
        // Iterators:
        
        const_iterator begin() const { return alph.begin(); }               /**< Begin iterator over all characters */
        const_iterator end() const { return alph.end(); }                   /**< Begin const iterator over all characters */
        
        const_iterator beginValid() const { return alphValid.begin(); }          /**< Begin iterator over valid characters */
        const_iterator endValid() const { return alphValid.end(); }              /**< Begin const iterator over valid characters */
        
        const_iterator beginInvalid() const { return alphAmbiguous.begin(); }        /**< Begin iterator over ambiguous characters */
        const_iterator endInvalid() const { return alphAmbiguous.end(); }            /**< Begin const iterator over ambiguous characters */
        
        
    private:
        
        vector<char> alph;                  /**< All characters that make up the alphabet (including ambiguities) */
        vector<char> alphValid;             /**< Characters that make up the valid alphabet (e.g. A, C, G, T) */
        vector<char> alphAmbiguous;         /**< Characters that make up the invalid alphabet (e.g. N, R, S, ...) */
        
        map<char, char> complementDNA;      /**< complement nucleotides for DNA alphabet */
        
        
    };
}


#endif /* AlphabetDNA_hpp */
