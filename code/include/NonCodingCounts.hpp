//
//  NonCodingCounts.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/7/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NonCodingCounts_hpp
#define NonCodingCounts_hpp

#include <stdio.h>
#include "UniformCounts.hpp"

namespace gmsuite {
    
    
    class NonCodingCounts : public UniformCounts {
        
    public:
        
        /**
         * Constructor: Create a uniform count model by defining it's order and alphabet
         *
         * @param order the model's order
         * @param alph the alphabet used by the model
         */
        NonCodingCounts(unsigned order, const NumAlphabetDNA &alph);
        
        /**
         * Count the sequence. This method calls the updateCounts (increment) method, which
         * should be implemented by each derived class.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         */
        void count(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement=false);
        
        
        /**
         * Decount the sequence. This method calls the updateCounts (decrement) method, which
         * should be implemented by each derived class.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         */
        void decount(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement=false);
        
    };
    
}


#endif /* NonCodingCounts_hpp */
