//
//  CharNumConverter.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/8/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef CharNumConverter_hpp
#define CharNumConverter_hpp

#include <stdio.h>
#include <vector>

#include "Sequence.hpp"
#include "AlphabetDNA.hpp"

using std::vector;

namespace gmsuite {
    
    /**
     * @class CharNumConverter
     * @brief Manage conversion between character and numeric sequences
     *
     * For efficiency reasons, the underlying representation of sequences in many cases 
     * is not a string of characters, but rather a sequence of numbers, each number representing
     * a character in the alphabet.
     *
     * Accordingly, this class handles the conversion of characters to and from numbers
     * based on the supplied alphabet.
     */
    class CharNumConverter {
        
    public:
        
        typedef int element_t;                  /**< type of the element */
        typedef vector<element_t> seq_t;        /**< type of a sequence of elements */
        
        /**
         * Constructor: create a converter using the given alphabet
         *
         * @param alphabet the alphabet used for conversion
         *
         * @exception invalid_argument if alphabet is NULL
         */
        CharNumConverter(const AlphabetDNA* alphabet);
        
        
        /**
         * Convert a string sequence to its numeric representation.
         *
         * @param str the string sequence
         * @param result the output numeric representation of the string
         *
         * @exception out_of_range if character is not in the converter's alphabet
         */
        void convert(const string &str, seq_t &result) const;
        
        
        /**
         * Convert a Sequence to its numeric representation.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         * @param result the output numeric representation of the string
         *
         * @exception out_of_range if character is not in the converter's alphabet
         */
        void convert(Sequence::const_iterator begin, Sequence::const_iterator end, seq_t &result) const;
        
        
        /**
         * Convert a single character to its numeric representation.
         *
         * @param c the character
         * @return the numeric representation
         *
         * @exception out_of_range if character is not in the converter's alphabet
         */
        element_t convert(char c) const;
        
        
        /**
         * Convert a numeric representation sequence back to a string sequence.
         *
         * @param start an interator pointing to the first element for conversion
         * @param end an interator pointing to the end of the conversion (i.e. not included in conversion)
         * @return a string representation of the numeric sequence
         *
         * @exception out_of_range if numeric element is not in the valid range
         */
        string convert(seq_t::const_iterator start, seq_t::const_iterator end) const;
        
        
        /**
         * Convert a numeric element back to a character.
         *
         * @param element the numeric element
         * @return a character representation of the numeric sequence
         *
         * @exception out_of_range if numeric element is not in the valid range
         */
        char convert(element_t element) const;
        
        
        /**
         * Complement DNA sequence
         *
         * @param original the original sequence
         * @param result the output
         * @exception logic_error DNA complements can only be done for DNA alphabets
         * @exception out_of_range if numeric element is not in the valid conversion range
         */
        void complement(const seq_t &original, seq_t &result) const;
        
        /**
         * Complement DNA element
         *
         * @param original the original element
         * @exception logic_error DNA complements can only be done for DNA alphabets
         * @exception out_of_range if numeric element is not in the valid conversion range
         */
        element_t complement(element_t original) const;
        
        
    private:
        
        const AlphabetDNA* alphabet;                    /**< the alphabet used by the converter */
        map<char, element_t> charToNum;                 /**< convert from character to element */
        map<element_t, char> numToChar;                 /**< convert from element to character */
        map<element_t, element_t> complementDNA;        /**< complement DNA "elements" */
        
    };

    
}

#endif /* CharNumConverter_hpp */
