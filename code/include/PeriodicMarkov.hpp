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
     * @class UnifMarkov
     * @brief A Uniform Markov model
     *
     * This class represents a uniform Markov model of any order.
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
        
        // Define: type to store probabilities of non-uniform Markov model.
        typedef vector<vector<double> > period_markov_t;              // for probabilities
        typedef vector<vector<string> > keys_t;                 // for keys (i.e. terms)
        
        
        period_markov_t model;                // to store probabilities
        size_t period;                  // model length
        
        /**
         * Convert counts in model to probability.
         */
        void convertCountsToProbability();
        
        /**
         * Initialize the model by allocating space, setting the keys, and setting counts to 0
         */
        void initialize();
        
        /**
         * Reset all counts to zero
         */
        void resetCounts();
        
        /**
         * Convert joing probabilities to Markov (conditional) probabilities.
         * E.g. P(ACG) -> P(G|CA)
         */
        void jointToMarkov(vector<double> &probs);
        
        
    };
    
}

#endif /* PeriodicMarkov_hpp */
