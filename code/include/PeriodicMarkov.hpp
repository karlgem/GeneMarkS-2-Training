//
//  PeriodicMarkov.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef PeriodicMarkov_hpp
#define PeriodicMarkov_hpp

#include <stdio.h>

#include "Markov.hpp"

namespace gmsuite {
    
    
    /**
     * @class PeriodicMarkov
     * @brief A Periodic Markov model
     *
     * This class represents a periodic Markov model of any order.
     */
    class PeriodicMarkov : public Markov {
        
    public:
        
        /**
         * Constructor:
         *
         * @param order the model's order
         * @param period the model's period
         * @param alphabet the alphabet used by the model
         */
        PeriodicMarkov(unsigned order, size_t period, const AlphabetDNA* alph);
        
        
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
        
        
    private:
        
        /**
         * Initialize the model by allocating space, setting the keys, and setting counts to 0
         */
        void initialize();
        
        
        // The structure of the model 'm' is a vector of vectors, where m[p] holds
        // the probabilities for frame 'p' of the model. If the model order is 2, the period
        // is 3, and the alphabet is made up of 2 letters A,B, then the model structure will
        // look like:
        //    m[0]  m[1]  m[2]
        //    AAA   AAA   AAA
        //    AAB   AAB   AAB
        //    ABA   ABA   ABA
        //    ABB   ABB   ABB
        //    BAA   BAA   BAA
        //    BAB   BAB   BAB
        //    BBA   BBA   BBA
        //    BBB   BBB   BBB
        typedef vector<vector<double> > period_markov_t;              // for markov probabilities
        
        
        // For order 2, to compute P(ACGT), we can break it down into
        // P(ACGT) = P(AC) * P(G|AC) * P(T|CG)
        // These joint probabilities P(AC) are stored in the below model,
        // since they don't fit easily in the above structure. Since this is a periodic model,
        // each frame has its own set of joint probability distributions
        typedef vector<vector<vector<double > > > period_joint_t;   // for joint probabilities of lower orders
        
        period_markov_t model;          // to store probabilities
        period_joint_t jointProbs;      // to store joint probabilities
        size_t period;                  // model length
    };
    
}

#endif /* PeriodicMarkov_hpp */
