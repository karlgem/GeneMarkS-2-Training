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
        typedef CharNumConverter::element_t num_t;                      /**< define generic type for number @see CharNumConverter */
        
        vector<num_t> numSeq;                                           /**< Numeric sequence */
        
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
    
    private:

    };
    
}

#endif /* NumSequence_hpp */
