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
#include "NumAlphabetDNA.hpp"

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
        Markov(unsigned order, const NumAlphabetDNA &alph);
        
        
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
        const NumAlphabetDNA* getAlphabet() const;
        
        
        /**
         * Generate a string representation of the model
         *
         * @return a string representation of the model
         */
        virtual string toString() const = 0;
        
        
        /**
         * Compute the KL-divergence of two Markov models
         *
         * @param P first markov model
         * @param Q second markov model
         * @return KL(P || Q)
         */
        friend double KLDivergence(const Markov* P, const Markov* Q);
        
        
    protected:
        
        unsigned order;                             /**< the model's order */
        const NumAlphabetDNA* alphabet;                /**< the alphabet */
        const CharNumConverter *cnc;                /**< the char-number converter; for reverse complementation */
        
        /**
         * Convert joing probabilities to Markov (conditional) probabilities.
         * E.g. P(ACG) -> P(G|CA)
         *
         * @param probs vector of joint probabilities
         */
        void jointToMarkov(vector<double> &probs) const;
        
        
        /**
         * Compute the 'direct' (i.e. currentOrder-1) lower order distribution of joint probabilities.
         *
         * @param currentOrder the current order of the joint probabilities
         * @param current a vector holding the joint probabilities for order <currentOrder>
         * @param result the vector where the result will be stored. That is, the joint probabilities
         * for order = currentOrder-1
         *
         * @throw invalid_argument if currentOrder is equal to zero
         */
        void getLowerOrderJoint(unsigned currentOrder, const vector<double> &current, vector<double> &result) const;
        
        
        /**
         * Compute higher order joint probabilities.
         *
         * @param currentOrder the current order of the joint probabilities
         * @param current a vector holding the joint probabilities for order <currentOrder>
         * @param newOrder the new order of joint probabilities
         * @param result the vector where the result will be stored. That is, the joint probabilities
         * for order=newOrder
         *
         * @exception invalid_argument if newOrder < currentOrder
         */
        void getHigherOrderJoint(unsigned currentOrder, const vector<double> &current, unsigned newOrder, vector<double> &result) const;
        
        /**
         * Increment order of joint probabilities by one.
         *
         * @param currentOrder the order of the "original" probabilities
         * @param currentProbs the "original" probabilities
         * @param newProbs a vector where the new probabilities (of order currentOrder+1) will be stored
         */
        void incrementOrderByOne(unsigned currentOrder, const vector<double> &currentProbs, vector<double> &newProbs) const;
        
        /**
         * Get cumulative distribution frequency for each conditional probability.
         *
         * @param currentOrder the order of the "original" probabilities
         * @param probs the "original" probabilities
         */
        void getCDFPerConditional(unsigned currentOrder, vector<double> &probs) const;
        
        /**
         * Convert the array index to the numeric sequence 'located' at that index.
         * For example, for index 7 and wordLength 3 (i..e 000111 -> ACT), the numeric
         * sequence will be 013
         */
        NumSequence indexToNumSequence(size_t idx, size_t wordLength) const;
        
        
        /**
         * Convert the array index to the numeric sequence 'located' at that index.
         * For example, for index 7 and wordLength 3 (i..e 000111 -> ACT), the numeric
         * sequence will be 013
         */
        void indexToNumSequence(size_t idx, size_t wordLength, vector<NumSequence::num_t> &numSeq) const;
        
    };
    
}

#endif /* Markov_hpp */
