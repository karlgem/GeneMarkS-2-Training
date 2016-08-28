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
     * A counts class representing generic Markov model counts
     */
    class Counts {
        
    public:
        
        /**
         * Constructor: Initialize a Counts model witha specific order and alphabet
         *
         * @param order the order of the model
         * @param alph the alphabet used by the sequences handled by this model
         *
         * @throw invalid_argument if alph is NULL
         */
        Counts(unsigned order, const AlphabetDNA* alph);
        
        /**
         * Copy constructor: handle data management during copies
         *
         * @param other the model to be copied
         */
        Counts(const Counts &other);
        
        
        /**
         * Copy assignment operator: handle data management during copies
         *
         * @param other the model to be copied
         */
        Counts& operator=(const Counts &other);
        
        
        /**
         * Destructor
         */
        virtual ~Counts();
        
        
        /**
         * Construct the model counts from a list of sequences
         *
         * @param sequences the list of sequences
         */
        virtual void construct(const vector<NumSequence> &sequences) = 0;
        
        
        /**
         * Count the sequence.
         *
         * @param sequence the sequence
         */
        virtual void count(NumSequence::const_iterator begin, NumSequence::const_iterator end) = 0;
        
        
        /**
         * Decount the sequence.
         *
         * @param sequence the sequence
         */
        virtual void decount(NumSequence::const_iterator begin, NumSequence::const_iterator end) = 0;
        
        
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
        
        
        /**
         * Reset all counts to zero
         */
        virtual void resetCounts() = 0;
        
    protected:
        
        unsigned order;                     /**< The model's order */
        const AlphabetDNA *alphabet;        /**< The alphabet */
        
    };
    
}

#endif /* Counts_hpp */
