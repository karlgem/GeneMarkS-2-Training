//
//  NonUniformMarkov.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NonUniformMarkov_hpp
#define NonUniformMarkov_hpp

#include <stdio.h>


#include "Markov.hpp"

namespace gmsuite {
    
    
    /**
     * @class NonUniformMarkov
     * @brief A Non-uniform Markov model
     *
     * This class represents a non-uniform Markov model of any order.
     */
    class NonUniformMarkov : public Markov {
        
    public:
        
        /**
         * Constructor:
         *
         * @param order the model's order
         * @param length the model's length
         * @param alphabet the alphabet used by the model
         */
        NonUniformMarkov(unsigned order, size_t length, const AlphabetDNA &alph, const CharNumConverter &cnc);
        
        
        /**
         * Construct the model probabilities from a list of sequences
         *
         * @param sequences the list of sequences
         * @param pcount the pseudocounts
         */
        void construct(const vector<NumSequence> &sequences, int pcount = 0);
        
        /**
         * Construct the model probabilities from existing counts.
         *
         * Note, the derived Counts subclass should be the complementary version
         * of the derived Markov subclass (e.g. UnifCounts and UnifMarkov).
         * Otherwise, an exception is thrown.
         *
         * @param counts the counts model
         * @param pcount the pseudocounts
         */
        void construct(const Counts* counts, int pcount = 0);
        
        /**
         * Compute the score of a sequence using the model probabilities
         *
         * @param begin where to begin evaluation in the sequence
         * @param end where to end evaluation in the sequences (exclusive)
         * @param useLog whether log-form should be used.
         */
        double evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog = false) const;
        
        /**
         * Generate a string representation of the model
         *
         * @return a string representation of the model
         */
        string toString() const;
        
        
    protected:
        
        /**
         * Initialize the model by allocating space, setting the keys, and setting counts to 0
         */
        void initialize();
        
        
        // The structure of the model 'm' is a vector of vectors, where m[p] holds
        // the probabilities for position 'p' of the model. If the model order is 2, the length
        // is 4, and the alphabet is made up of 2 letters A,B, then the model structure will
        // look like:
        //    m[0]  m[1]  m[2]  m[3]
        //      A   AA    AAA   AAA
        //      B   AB    AAB   AAB
        //          BA    ABA   ABA
        //          BB    ABB   ABB
        //                BAA   BAA
        //                BAB   BAB
        //                BBA   BBA
        //                BBB   BBB
        typedef vector<vector<double> > nonunif_markov_t;              // for probabilities
        
        
        // For order 2, to compute P(ACGT), we can break it down into
        // P(ACGT) = P(AC) * P(G|AC) * P(T|CG)
        // These joint probabilities P(AC) are stored in the below model,
        // since they don't fit easily in the above structure. Since this is a nonuniform model,
        // each position has its own set of joint probability distributions
        typedef vector<vector<double> > nonunif_joint_t;   // for joint probabilities of each position
        
        nonunif_markov_t model;             // to store probabilities
        nonunif_joint_t jointProbs;         // to store joint probabilities
        size_t length;                      // model length
        
    };
    
}



#endif /* NonUniformMarkov_hpp */
