//
//  MotifFinder.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef MotifFinder_hpp
#define MotifFinder_hpp

#include <stdio.h>
#include <limits.h>
#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"
#include "MFinderModelParams.hpp"
#include "OptionsMFinder.hpp"
#include "ProbabilityModels.hpp"

namespace gmsuite {
    
    /**
     * @class MotifFinder
     * @brief This class finds motifs in sequences
     *
     * Given a set of sequences, find a motif of length k that is
     * common amongst all sequences, under a certain probabilistic model.
     */
    class MotifFinder {
        
    public:
        
        class Builder;          // used for easy building of the many MotifFinder parameters
        
        /**
         * Constructor: create a motif-finder class to find motifs in sequences
         *
         * @param width the motif's length
         * @param motifOrder the order for the motif Markov model
         * @param backOrder the order for the background model
         * @param pcounts the pseudocounts
         * @param align whether to align sequences and use a length distribution (see class def)
         * @param tries the number of times with different starting points the algorithm is executed
         * @param maxIter the maximum number of gibbs iterations per try
         * @param maxEMIter the maximum number of iterations for the EM phase
         * @param shiftEvery the number of iterations between every shift-operation
         */
        MotifFinder(unsigned width,
                    unsigned motifOrder     =   0,
                    unsigned backOrder      =   0,
                    double pcounts          =   1,
                    MFinderModelParams::align_t align           =   MFinderModelParams::NONE,
                    unsigned tries          =   10,
                    unsigned maxIter        =   100,
                    unsigned maxEMIter      =   10,
                    unsigned shiftEvery     =   20,
                    double filterThresh     =   -std::numeric_limits<double>::infinity());
        
        
        /**
         * Find motifs in sequences.
         *
         * @param sequences the sequences to be searched
         * @param positions the motif positions in each sequence
         */
        void findMotifs (const vector<NumSequence> &sequences, vector<NumSequence::size_type> &positions);
        void findMotifs (const vector<NumSequence> &sequences, vector<NumSequence::size_type> &positions, ProbabilityModels *&probs);
        
        
    private:
        
        
        /**
         * Run a single try of gibbs-finder to search for the best motif alignment.
         *
         * @param sequences the sequences to be searched
         * @param positions the motif positions in each sequence
         * @return the alignment's probability
         */
        double gibbsFinder(const vector<NumSequence> &sequences, vector<NumSequence::size_type> &positions, ProbabilityModels *&probs);
        
        
        /**
         * Attempt to shift positions, and return the sampled amount to shift
         */
        int attemptShift(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions);
        
        /**
         * Shift positions by a certain amount
         */
        void shiftPositions(const vector<NumSequence::size_type> &original, int shiftAmount, const vector<NumSequence> &sequences, vector<NumSequence::size_type> &result);
        
        
        
        // member variables
        unsigned width                              ;       /**< motif length */
        unsigned motifOrder                         ;       /**< motif model order */
        unsigned backOrder                          ;       /**< background model order */
        double   pcounts                            ;       /**< pseudocounts */
        unsigned tries                              ;       /**< number of tries (starting points) */
        unsigned maxIter                            ;       /**< maximum number of iterations per try */
        unsigned maxEMIter                          ;       /**< maximum number of EM iterations per try */
        unsigned shiftEvery                         ;       /**< number of iterations between every shift-operation */
        AlphabetDNA alphabet                        ;       /**< the alphabet used by the sampler; at the moment, only DNA is available */
        MFinderModelParams::align_t align           ;       /**< whether to align sequences first and use length distribution */
        double filterThresh                         ;       /**< allows filtering of sequences with low scores */
        
    };
    
    
    
    // Class for building a MotifFinder (for advanced users)
    class MotifFinder::Builder {
        
    private:
        // variables needed for MotifFinder class
        
        unsigned    tries;
        unsigned    MAX_ITER;
        unsigned    MAX_EM_ITER;
        unsigned    shiftEvery;
        unsigned    motifOrder;
        unsigned    backOrder;
        unsigned    width;
        double      pcounts;
        bool        allSeqsPerIter;
        MFinderModelParams::align_t align;
        double        filterThresh;
        
        
    public:
        
        Builder() {
            tries           =   10;
            MAX_ITER        =   60;
            shiftEvery      =   10;
            MAX_EM_ITER     =   10;
            motifOrder      =   0;
            backOrder       =   0;
            width           =   6;
            pcounts         =   1;
            align           =   MFinderModelParams::NONE;
            allSeqsPerIter  =   true;
            filterThresh    =   -std::numeric_limits<double>::infinity();
        }
        
        
        // set custom values for MotifFinder creation
        // returns Builder for shorthand inline usage
        Builder& setNumTries    (const unsigned v)                          { this->tries = v;      return *this;   }
        Builder& setMaxIter     (const unsigned v)                          { this->MAX_ITER = v;   return *this;   }
        Builder& setMaxEMIter   (const unsigned v)                          { this->MAX_EM_ITER = v;return *this;   }
        Builder& setShiftEvery  (const unsigned v)                          { this->shiftEvery = v; return *this;   }
        Builder& setMotifOrder  (const unsigned v)                          { this->motifOrder = v; return *this;   }
        Builder& setBackOrder   (const unsigned v)                          { this->backOrder = v;  return *this;   }
        Builder& setWidth       (const unsigned v)                          { this->width = v;      return *this;   }
        Builder& setPcounts     (const unsigned v)                          { this->pcounts = v;    return *this;   }
        Builder& setAlign       (const MFinderModelParams::align_t v)       { this->align = v;      return *this;   }
        Builder& fullLoopPerIter(const bool v)                              { this->allSeqsPerIter = v; return *this; }
        Builder& setFilterThresh(const double v)                            { this->filterThresh = v; return *this; }
        
        // build motif finder with set parameters
        MotifFinder build() {
            return MotifFinder(width, motifOrder, backOrder, pcounts, align, tries, MAX_ITER, MAX_EM_ITER, shiftEvery, filterThresh);
        }
        
        MotifFinder build(const OptionsMFinder &options) {
            this->setNumTries(options.tries).setMaxIter(options.maxIter).setMaxEMIter(options.maxEMIter).setShiftEvery(options.shiftEvery);
            this->setMotifOrder(options.motifOrder).setBackOrder(options.bkgdOrder).setWidth(options.width).setPcounts(options.pcounts);
            this->setAlign(options.align);
            this->setFilterThresh(options.filterThresh);
            return build();
        }

    };


}

#endif /* MotifFinder_hpp */
