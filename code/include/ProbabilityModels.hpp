//
//  ProbabilityModels.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ProbabilityModels_hpp
#define ProbabilityModels_hpp

#include <stdio.h>

#include "Sequence.hpp"
#include "CountModels.hpp"

namespace gmsuite {
    
    /**
     * @class ProbabilityModels
     * @brief A class to hold probabilities (used in MotifFinder)
     */
    class ProbabilityModels {
        
    public:
        
        /**
         * Destructor
         */
        virtual ~ProbabilityModels();
        
        /**
         * Construct counts from a set of motifs.
         *
         * @param sequences the sequences
         * @param positions the positions of motifs in these sequences
         */
        virtual void construct(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions) = 0;
        
        
        /**
         * Construct counts from a set of motifs.
         *
         * @param counts sequences
         */
        virtual void construct(const CountModels* counts) = 0;
        
        
        /**
         * Compute conditional log-likelihood
         */
        virtual double computeCLL() = 0;
        
        
        /**
         * Compute the score for a given motif position
         *
         * @param sequence the sequence
         * @param pos the position of the motif in the sequence
         */
        virtual double computePositionScore(const NumSequence &sequence, NumSequence::size_type pos) = 0;
        
        
        /**
         * Compute the scores for each valid motif position in the sequence
         *
         * @param sequence the sequence
         * @param scores the output scores of all valid positions in the sequence
         */
        virtual void computePositionScores(const NumSequence &sequence, vector<double> &scores) = 0;
        
        /**
         * Sample the position for motif in sequences
         *
         * @param sequence the sequence
         * @param getMax if set, the position with the highest probability is returned; otherwise it is sampled.
         *
         * @return the position of a motif
         */
        virtual NumSequence::size_type samplePosition(const NumSequence &sequence, bool getMax = false) = 0;
        
        /**
         * Get string representation of counting models
         */
        virtual string toString() const = 0;
        
    };
}



#endif /* ProbabilityModels_hpp */
