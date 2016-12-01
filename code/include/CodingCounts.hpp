//
//  CodingCounts.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef CodingCounts_hpp
#define CodingCounts_hpp

#include <stdio.h>
#include "PeriodicCounts.hpp"
#include "NumGeneticCode.hpp"

namespace gmsuite {

    /**
     * @class CodingCounts
     * @brief Get frequency counts for coding region
     */
    class CodingCounts : public PeriodicCounts {
        
    public:
        
        /**
         * Constructor: Create a periodic count model by defining it's order,
         * period, and alphabet.
         *
         * @param order the model's order (at every period)
         * @param period the model's period; i.e. the number of different count models
         * @param alph the alphabet used by the model
         */
        CodingCounts(unsigned order, size_t period, const NumAlphabetDNA &alph, const NumGeneticCode &geneticCode);

        
        
        /**
         * Update counts for a given sequence, by either incrementing or decrementing them. This provides
         * a common implementation for count/decount methods. Each derived class of Counts should implement
         * this method, as it is called by the count/decount methods.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         * @param operation what to do: either "increment" or "decrement"
         *
         * @throw invalid_argument if operation is neither "increment" or "decrement"
         */
        void updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation, bool reverseComplement = false);
        
    protected:
        const NumGeneticCode* geneticCode;      /**< Numeric version of the genetic code */
    };
}

#endif /* CodingCounts_hpp */
