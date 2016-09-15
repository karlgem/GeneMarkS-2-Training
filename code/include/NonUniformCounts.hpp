//
//  NonUniformCounts.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/29/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NonUniformCounts_hpp
#define NonUniformCounts_hpp

#include <stdio.h>
#include "Counts.hpp"

namespace gmsuite {
    
    /**
     * @class NonUniformCounts
     * @brief A class to handle counts through a nonuniform setting
     *
     * Nonuniformity here implies that, given a set of sequences, the i'th
     * element of each sequence contributes to the frequencies of the i'th
     * count model. Parts of sequences beyond the model's length are
     * ignored.
     */
    class NonUniformCounts : public Counts {
        
        friend class NonUniformMarkov;
        
    public:
        
        /**
         * Constructor: Create a nonuniform count model by defining it's order,
         * length, and alphabet.
         *
         * @param order the model's order (at every 'length')
         * @param length the model's length; i.e. the number of different count models
         * @param alph the alphabet used by the model
         * @param cnc the character-number converter; for reverse complementation
         */
        NonUniformCounts(unsigned order, size_t length, const AlphabetDNA &alph, const CharNumConverter &cnc);
        
        
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
         * Get the model's length
         *
         * @return the model's length.
         */
        size_t getLength() const;
        
        
        /**
         * Reset all counts to zero
         */
        void resetCounts();
        
        
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
        
        
        // The structure of the model 'm' is a vector of vectors, where m[p] holds
        // the counts for position 'p' of the model. If the model order is 2, the length
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
        typedef vector<vector<double> > nonunif_counts_t;          // to store counts
        
        nonunif_counts_t model;             // to store counts
        size_t length;                      // model's length
        
    };

}

#endif /* NonUniformCounts_hpp */
