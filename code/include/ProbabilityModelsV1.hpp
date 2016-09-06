//
//  ProbabilityModelsV1.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ProbabilityModelsV1_hpp
#define ProbabilityModelsV1_hpp

#include <stdio.h>

#include "UniformMarkov.hpp"
#include "NonUniformMarkov.hpp"
#include "ProbabilityModels.hpp"
#include "MFinderModelParams.hpp"
#include "CountModels.hpp"
#include "NonUniformCounts.hpp"
#include "UnivariatePDF.hpp"
//#include "Distribution.hpp"

namespace gmsuite {
    
    /**
     * @class ProbabilityModels
     * @brief A class to hold probabilities (used in MotifFinder)
     */
    class ProbabilityModelsV1 : public ProbabilityModels {
        
    public:
        
        /**
         * Constructor
         */
        ProbabilityModelsV1(const AlphabetDNA &alphabet, NumSequence::size_type width, unsigned motifOrder, unsigned backOrder, double pcounts, MFinderModelParams::align_t align);
        
        /**
         * Copy constructor
         */
        ProbabilityModelsV1(const ProbabilityModelsV1 &obj);
        
        /**
         * Copy assignment operator
         */
        ProbabilityModelsV1& operator=(const ProbabilityModelsV1& other);
        
        /**
         * Destructor
         */
        ~ProbabilityModelsV1();
        
        
        /**
         * Construct probabilities from a set of motifs.
         *
         * @param sequences the sequences
         * @param positions the positions of motifs in these sequences
         */
        void construct(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions);
        
        
        /**
         * Construct probabilities from counts
         *
         * @param sequences the sequences
         * @param positions the positions of motifs in these sequences
         */
        void construct(const CountModels* counts);
        
        
        /**
         * Compute conditional log-likelihood
         */
        double computeCLL();
        
        
        /**
         * Compute the score for a given motif position
         *
         * @param sequence the sequence
         * @param pos the position of the motif in the sequence
         */
        double computePositionScore(const NumSequence &sequence, NumSequence::size_type pos);
        
        
        /**
         * Compute the scores for each valid motif position in the sequence
         *
         * @param sequence the sequence
         * @param scores the output scores of all valid positions in the sequence
         */
        void computePositionScores(const NumSequence &sequence, vector<double> &scores);
        
        
        /**
         * Sample the position for motif in sequences
         *
         * @param sequence the sequence
         * @param getMax if set, the position with the highest probability is returned; otherwise it is sampled.
         *
         * @return the position of a motif
         */
        NumSequence::size_type samplePosition(const NumSequence &sequence, bool getMax = false);
        
        
        /**
         * Get string representation of counting models
         */
        string toString() const;
        
        
    private:
        
        const AlphabetDNA *alphabet;               /**< alphabet */
        NumSequence::size_type width;              /**< motif width */
        unsigned motifOrder;                    /**< order for motif model */
        unsigned backOrder;                     /**< order for background model */
        double pcounts;                         /**< pseudocounts */
        MFinderModelParams::align_t align;         /**< align */
        
        NonUniformMarkov* mMotif;                  /**< motif probabilities */
        NonUniformCounts* mMotifCounts;            /**< motif counts */
        UniformMarkov* mBack;                      /**< background probabilities */
        UnivariatePDF *positionDistribution;     /**< position distribution; defined by alignment */
        vector<double> positionCounts;          /**< position counts; defined by alignment */
        
        
    };
}



#endif /* ProbabilityModelsV1_hpp */
