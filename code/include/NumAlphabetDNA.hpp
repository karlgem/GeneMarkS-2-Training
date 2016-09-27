//
//  NumAlphabetDNA.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NumAlphabetDNA_hpp
#define NumAlphabetDNA_hpp

#include <stdio.h>

#include "AlphabetDNA.hpp"
#include "CharNumConverter.hpp"

namespace gmsuite {
    
    /**
     * @class NumAlphabetDNA
     * @brief A class that represents alphabet DNA in numeric form
     */
    class NumAlphabetDNA {
        
    public:
        
        typedef CharNumConverter::element_t element_t;
        typedef vector<element_t>::size_type size_type;              /**< Type for number of characters */
        typedef vector<element_t>::const_iterator const_iterator;    /**< Const iterator over characters */
        
        
        /**
         * Constructor: initialize DNA alphabet
         *
         * The DNA alphabet is made up of two classes of characters
         *   - Valid: standard characters making up a sequence (A,C,G,T)
         *   - Ambiguous: characters making up ambiguous elements (W,S,M,K,R,Y,B,D,H,V,N,Z)
         */
        NumAlphabetDNA(const AlphabetDNA &alph, const CharNumConverter &cnc);
        
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
        bool contains(element_t c) const;
        
        /**
         * Check if character is valid (i.e. non-ambiguous)
         *
         * @param c the character
         * @return true if c is valid; false otherwise.
         */
        bool isValid(element_t c) const;
        
        
        /**
         * Check if character is ambiguous (i.e. N, R, S, ...)
         *
         * @param c the character
         * @return true if c is valid; false otherwise.
         */
        bool isAmbiguous(element_t c) const;
        
        /**
         * Complement a character
         */
        element_t complement(element_t c) const;
        
        
        // Iterators:
        
        const_iterator begin() const { return alph.begin(); }               /**< Begin iterator over all characters */
        const_iterator end() const { return alph.end(); }                   /**< Begin const iterator over all characters */
        
        const_iterator beginValid() const { return alphValid.begin(); }          /**< Begin iterator over valid characters */
        const_iterator endValid() const { return alphValid.end(); }              /**< Begin const iterator over valid characters */
        
        const_iterator beginInvalid() const { return alphAmbiguous.begin(); }        /**< Begin iterator over ambiguous characters */
        const_iterator endInvalid() const { return alphAmbiguous.end(); }            /**< Begin const iterator over ambiguous characters */
        
        const CharNumConverter* getCNC() const;
        
    private:
        const CharNumConverter *cnc;
        vector<element_t> alph;                  /**< All characters that make up the alphabet (including ambiguities) */
        vector<element_t> alphValid;             /**< Characters that make up the valid alphabet (e.g. A, C, G, T) */
        vector<element_t> alphAmbiguous;         /**< Characters that make up the invalid alphabet (e.g. N, R, S, ...) */
        
        map<element_t, element_t> complementDNA;      /**< complement nucleotides for DNA alphabet */
        
    };
    
    
}

#endif /* NumAlphabetDNA_hpp */
