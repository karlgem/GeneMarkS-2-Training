//
//  UniformMarkov.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef UniformMarkov_hpp
#define UniformMarkov_hpp

#include <stdio.h>
#include "Markov.hpp"

namespace gmsuite {
    
    /**
     * @class UniformMarkov
     * @brief A uniform Markov model
     *
     * This class represents a uniform Markov model of any order.
     */
    class UniformMarkov : public Markov {
        
    public:
        
        /**
         * Constructor:
         *
         * @param order the model's order
         * @param alph the alphabet used by the model
         */
        UniformMarkov(unsigned order, const NumAlphabetDNA &alph);
        
        
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
        
        
        /**
         * Change the model's order. Note: changing to a lower order model will
         * incur and irreversable loss of information; i.e. lowering an order then 
         * raising it again will not return the original model, as the distribution
         * has already been marginalized.
         *
         * @param newOrder the new order
         */
        void changeOrder(unsigned newOrder);
        
        
        /**
         * Emit a sequence from the markov model
         *
         * @param length the length of the sequence
         */
        NumSequence emit(NumSequence::size_type length) const;
        
        
    protected:

        /**
         * Initialize the model by allocating space, setting the keys, and setting counts to 0
         */
        void initialize();
        
        
        
        // The structure of the model 'm' is a vector of doubles, where m holds
        // the counts of the model. If the model order is 2, and the alphabet is
        // made up of 2 letters A,B, then the model structure will look like:
        //     m
        //    AAA
        //    AAB
        //    ABA
        //    ABB
        //    BAA
        //    BAB
        //    BBA
        //    BBB
        typedef vector<double> uniform_markov_t;            // for markov probabilities
        
        // For order 2, to compute P(ACGT), we can break it down into
        // P(ACGT) = P(AC) * P(G|AC) * P(T|CG)
        // These joint probabilities P(AC) are stored in the below model,
        // since they don't fit easily in the above structure. 
        typedef vector<vector<double > > uniform_joint_t;   // for joint probabilities of lower orders
        
        uniform_markov_t model;                         // to store probabilities
        uniform_joint_t jointProbs;                     // to store joint probabilities of lower orders
    };
    
}

#endif /* UniformMarkov_hpp */
