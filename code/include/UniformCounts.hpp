//
//  UniformCounts.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef UniformCounts_hpp
#define UniformCounts_hpp

#include <stdio.h>
#include "Counts.hpp"

namespace gmsuite {
    
    /**
     * @class UniformCounts
     * @brief A class to handle counts through the uniform (i.e. standard) Markov setting
     */
    class UniformCounts : public Counts {
        
        friend class UniformMarkov;
       
    public:
        
        /**
         * Constructor: Create a uniform count model by defining it's order and alphabet
         *
         * @param order the model's order
         * @param alph the alphabet used by the model
         */
        UniformCounts(unsigned order, const AlphabetDNA &alph, const CharNumConverter &cnc);
        
        
        /**
         * Construct the model counts from a list of sequences
         *
         * @param sequences the list of sequences
         */
        void construct(const vector<NumSequence> &sequences);
        
        
        /**
         * Generate a string representation of the model
         *
         * @return a string representation of the model
         */
        string toString() const;
        
        
        /**
         * Reset all counts to zero
         */
        virtual void resetCounts();

    
    protected:
        
        /**
         * Initialize the model by allocating space and setting counts to 0
         */
        void initialize();
        
        
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
        void updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation, bool reverseComplement = false);
        
        
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
        typedef vector<double> uniform_counts_t;        // to store counts
        
        uniform_counts_t model;                         // to store counts
        
    };
}

#endif /* UniformCounts_hpp */
