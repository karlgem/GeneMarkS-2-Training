//
//  NumSequence.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NumSequence_hpp
#define NumSequence_hpp

#include <stdio.h>
#include <vector>

#include "Sequence.hpp"
#include "NumAlphabetDNA.hpp"
#include "CharNumConverter.hpp"

using std::vector;
using std::ostream;

namespace gmsuite {
    
    /**
     * @class NumSequence
     * @brief A numeric representation of a DNA sequence
     *
     * Numeric representation of DNA sequences allow for quick computations
     * e.g. in HMM
     *
     * Note to Alex: For now, I'm keeping the variable numSeq as public. This variable contains
     * the actual numeric sequence. That said, you should initialize the variable using the constructor
     * that takes a Sequence and CharNumConverter (which is used for converting the characters to numbers).
     *
     * @see Sequence
     * @see CharNumConverter
     */
    class NumSequence {
        
        
    public:
        
        typedef vector<int>::size_type size_type;                       // type for numeric sequence size
        static const size_type npos = Sequence::npos;                   /**< Returned to indicate no matches */
        typedef CharNumConverter::element_t num_t;                      /**< define generic type for number @see CharNumConverter */
        
        /**
         * Default constructor: create an empty numeric sequence.
         */
        NumSequence();
        
        /**
         * Constructor: create a numeric sequence from a regular sequence object.
         *
         * @param sequence the original letter-based sequence
         * @param converter the object used to convert character sequences to numeric sequences
         */
        NumSequence(const Sequence &sequence, const CharNumConverter &converter);
        
        
        /**
         * Constructor: create a numeric sequence a vector of num_t elements. This simply
         * copies the vector into the NumSequence class.
         *
         * @param numSequence the vector of num_t elements constituting the sequence
         */
        NumSequence(const vector<num_t> &numSequence);
        
        
        /**
         * Access an element from a const numeric sequence (i.e. cannot be modified)
         *
         * @param idx the index of the element
         */
        const int&  operator[](size_type idx) const;
        
        /**
         * Access an element from a numeric sequence
         *
         * @param idx the index of the element
         */
        int& operator[](size_type idx);
        
        /**
         * Get the size of the sequence (equivalent to sequence length()).
         *
         * @return the size of the sequence.
         */
        virtual size_type size() const;
        
        
        /**
         * Get a subsequence
         *
         * @param n the start index of the subsequence (inclusive)
         * @param length the length of the subsequence
         *
         * @exception std::invalid_argument thrown if n is larger than the sequence's
         * length, or if n + length is larger than the sequence's length.
         */
        NumSequence subseq(size_type n, size_type length) const;
        
        
        /**
         * Reverse complement the numeric sequence
         *
         * @param cnc the char-num converter that holds the DNA-complement information for the 
         * alphabet used by this sequence.
         */
        void reverseComplement(const CharNumConverter &cnc);
        
        
        bool containsInvalid(const NumAlphabetDNA &alph) const;
        
        
        NumSequence operator+ (const NumSequence &op2);
        
        
        /*************** Sequence Iterators *******************/
        
        // Iterators
        typedef vector<num_t>::iterator iterator;                   /**< Sequence iterator */
        
        virtual iterator begin();                                   /**< Start of iterator */
        virtual iterator end();                                     /**< End of iterator   */
        
        
        // Const Iterators
        typedef vector<num_t>::const_iterator const_iterator;       /**< Const sequence iterator */
        
        virtual const_iterator begin() const;                       /**< Start of const iterator */
        virtual const_iterator end() const;                         /**< End of const iterator   */
    
    private:
        
        vector<num_t> numSeq;                                           /**< Numeric sequence */
    };
    
}

#endif /* NumSequence_hpp */
