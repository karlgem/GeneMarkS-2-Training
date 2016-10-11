//
//  CountModelsV1.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef CountModelsV1_hpp
#define CountModelsV1_hpp

#include <stdio.h>

#include "NumSequence.hpp"
#include "CountModels.hpp"
#include "UniformCounts.hpp"
#include "NonUniformCounts.hpp"
#include "MFinderModelParams.hpp"
#include "ProbabilityModelsV1.hpp"


namespace gmsuite {
    
    class CountModelsV1 : public CountModels {
        
        friend class ProbabilityModelsV1;
        
    public:
        
        /**
         * Constructor
         */
        CountModelsV1(const NumAlphabetDNA &alphabet, Sequence::size_type width, unsigned motifOrder, unsigned backOrder, MFinderModelParams::align_t align);
        
        /**
         * Copy constructor
         */
        CountModelsV1(const CountModelsV1 &obj);
        
        /**
         * Copy assignment operator
         */
        CountModelsV1& operator=(const CountModelsV1& other);
        
        /**
         * Destructor
         */
        ~CountModelsV1();
        
        
        
        /**
         * Construct counts from a set of motifs.
         *
         * @param sequences the sequences
         * @param positions the positions of motifs in these sequences
         */
        void construct(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions);
        
        
        /**
         * Decount a motif
         *
         * @param sequence the sequence to decount
         * @param pos the position of the motif in the sequence
         *
         * @exception invalid_argument if pos is not a valid motif location in the sequence
         */
        void decount(const NumSequence &sequence, NumSequence::size_type pos);
        
        /**
         * Count a motif
         *
         * @param sequence the sequence to count
         * @param pos the position of the motif in the sequence
         *
         * @exception invalid_argument if pos is not a valid motif location in the sequence
         */
        void count(const NumSequence &sequence, NumSequence::size_type pos);
        
        
        /**
         * Get string representation of counting models
         */
        string toString() const;
        
        
    private:
        
        const NumAlphabetDNA* alphabet;                /**< alphabet */
        NumSequence::size_type width;               /**< motif width */
        unsigned motifOrder;                        /**< order for motif model */
        unsigned backOrder;                         /**< order for background model */
        MFinderModelParams::align_t align;          /**< align */
        
        NonUniformCounts* mMotif;                   /**< motif counts */
        UniformCounts* mBack;                       /**< background counts */
        vector<double> positionCounts;              /**< counts of motif positions */
    };
}


#endif /* CountModelsV1_hpp */
