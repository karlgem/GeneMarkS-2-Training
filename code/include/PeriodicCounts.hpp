//
//  PeriodicCounts.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/26/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef PeriodicCounts_hpp
#define PeriodicCounts_hpp

#include <stdio.h>

#include "Counts.hpp"

namespace gmsuite {
    
    /**
     * @class PeriodicCounts
     * @brief A class to handle counts through a periodic setting
     *
     * Periodicity here implies that, given a sequence, every i < P element
     * of that sequence contributes to the same i'th count model, where
     * P is called the period of the model.
     */
    class PeriodicCounts : public Counts {
        
        friend class PeriodicMarkov;
        
    public:
        
        /**
         * Constructor: Create a periodic count model by defining it's order,
         * period, and alphabet.
         *
         * @param order the model's order (at every period)
         * @param period the model's period; i.e. the number of different count models
         * @param alph the alphabet used by the model
         */
        PeriodicCounts(unsigned order, size_t period, const AlphabetDNA *alph);
        
        
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
         * Get the model's period/
         * 
         * @return the model's period.
         */
        size_t getPeriod() const;
        
        
        /**
         * Reset all counts to zero
         */
        void resetCounts();

        
    private:
        
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
        void updateCounts(NumSequence::const_iterator begin, NumSequence::const_iterator end, string operation);
        
        
        // The structure of the model 'm' is a vector of vectors, where m[p] holds
        // the counts for frame 'p' of the model. If the model order is 2, the period
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
        typedef vector<vector<double> > period_counts_t;          // to store counts
        
        period_counts_t model;          // to store counts
        size_t period;                  // model's period
        
    };
    
}

#endif /* PeriodicCounts_hpp */
