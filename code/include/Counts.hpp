//
//  Counts.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/26/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Counts_hpp
#define Counts_hpp

#include <stdio.h>

#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"

namespace gmsuite {
    
    /**
     * @class Counts
     * @brief An abstract class that represents Markov counts
     *
     * A counts class representing generic Markov model counts. Given a set of sequences,
     * this class can extract the element counts of that sequence based on the underlying
     * Markov properties.
     */
    class Counts {
        
    public:
        
        /**
         * Constructor: Initialize a Counts model with a specific order and alphabet.
         *
         * @param order the order of the model
         * @param alph the alphabet used by the sequences handled by this model
         *
         * @throw invalid_argument if alph is NULL
         */
        Counts(unsigned order, const AlphabetDNA &alph, const CharNumConverter &cnc);
        
        
        /**
         * Count the sequence. This method calls the updateCounts (increment) method, which
         * should be implemented by each derived class.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         */
        void count(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement=false);
        
        
        /**
         * Decount the sequence. This method calls the updateCounts (decrement) method, which
         * should be implemented by each derived class.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         */
        void decount(NumSequence::const_iterator begin, NumSequence::const_iterator end, bool reverseComplement=false);
        
        
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
         * Construct the model counts from a list of sequences
         *
         * @param sequences the list of sequences
         */
        virtual void construct(const vector<NumSequence> &sequences) = 0;
        
        
        /**
         * Generate a string representation of the model
         *
         * @return a string representation of the model
         */
        virtual string toString() const = 0;
        
        
        /**
         * Reset all counts to zero
         */
        virtual void resetCounts() = 0;
        
    protected:
        
        unsigned order;                     /**< The model's order */
        const AlphabetDNA *alphabet;        /**< The alphabet */
        const CharNumConverter *cnc;        /**< The character-number converter. Contains info on complementation */
        
        
        /**
         * Update counts for a given sequence, by either incrementing or decrementing them. This provides
         * a common implementation for count/decount methods. Each derived class of Counts should implement
         * this method, as it is called by the count/decount methods.
         *
         * @param begin the start of the sequence
         * @param end the end of the sequence
         * @param operation what to do: either "increment" or "decrement"
         *
         * @throw invalid_argument if operation is neither "increment" or "decrement"
         */
        virtual void updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation, bool reverseComplement = false) = 0;
        
    };
    
}

#endif /* Counts_hpp */
