//
//  CountModels.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef CountModels_hpp
#define CountModels_hpp

#include <stdio.h>

#include "NumSequence.hpp"

namespace gmsuite {
    
    /**
     * @class CountModels
     * @brief A class to hold counts (used in MotifFinder)
     */
    class CountModels {
        
    public:
        
        /**
         * Destructor
         */
        virtual ~CountModels();
        
        /**
         * Construct counts from a set of motifs.
         *
         * @param sequences the sequences
         * @param positions the positions of motifs in these sequences
         */
        virtual void construct(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions) = 0;
        
        
        /**
         * Decount a motif
         *
         * @param sequence the sequence to decount
         * @param pos the position of the motif in the sequence
         *
         * @exception invalid_argument if pos is not a valid motif location in the sequence
         */
        virtual void decount(const NumSequence &sequence, NumSequence::size_type pos) = 0;
        
        /**
         * Count a motif
         *
         * @param sequence the sequence to count
         * @param pos the position of the motif in the sequence
         *
         * @exception invalid_argument if pos is not a valid motif location in the sequence
         */
        virtual void count(const NumSequence &sequence, NumSequence::size_type pos) = 0;
        
        /**
         * Get string representation of counting models
         */
        virtual string toString() const = 0;
        
    };
}


#endif /* CountModels_hpp */
