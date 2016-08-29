//
//  Markov.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/27/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Markov_hpp
#define Markov_hpp

#include <stdio.h>

#include "Counts.hpp"
#include "NumSequence.hpp"
#include "AlphabetDNA.hpp"

namespace gmsuite {
    
    /**
     * @class Markov
     * @brief An abstract class for Markov models
     *
     * A Markov class representing generic Markov model probabilities.
     */
    class Markov {
        
    public:
        
        /*
         * Contructor: Initialize a Markov model with a specific order and alphabet.
         *
         * @param order the order of the model
         * @param alph the alphabet used by the sequences handled by this model
         *
         * @throw invalid_argument if alph is NULL
         */
        Markov(unsigned order, const AlphabetDNA* alph);
        
        
        /**
         * Construct the model probabilities from a list of sequences
         *
         * @param sequences the list of sequences
         * @param pcount the pseudocounts
         */
        virtual void construct(const vector<NumSequence> &sequences, int pcount = 0) = 0;
        
        
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
        virtual void construct(const Counts* counts, int pcount = 0) = 0;
        
        
        /**
         * Compute the score of a sequence using the model probabilities
         *
         * @param begin where to begin evaluation in the sequence
         * @param end where to end evaluation in the sequences (exclusive)
         * @param useLog whether log-form should be used.
         */
        virtual double evaluate(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool useLog = false) const = 0;
        

        /**
         * Get the model's order
         *
         * @return the model's order
         */
        unsigned getOrder() const;
        
        
        /**
         * Get the model's alphabet
         *
         * @return the model's alphabet
         */
        const AlphabetDNA* getAlphabet() const;
        
        
        /**
         * Generate a string representation of the model
         *
         * @return a string representation of the model
         */
        virtual string toString() const = 0;
        
    protected:
        
        unsigned order;                             /**< the model's order */
        const AlphabetDNA* alphabet;                /**< the alphabet */
        
        
        /**
         * Convert joing probabilities to Markov (conditional) probabilities.
         * E.g. P(ACG) -> P(G|CA)
         *
         * @param probs vector of joint probabilities
         */
        void jointToMarkov(vector<double> &probs);
        
    };
    
}

#endif /* Markov_hpp */
