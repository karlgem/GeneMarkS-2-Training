//
//  GMS2Trainer.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef GMS2Trainer_hpp
#define GMS2Trainer_hpp

#include <stdio.h>
#include "Label.hpp"
#include "NumSequence.hpp"
#include "UniformMarkov.hpp"
#include "PeriodicMarkov.hpp"
#include "NonUniformMarkov.hpp"


namespace gmsuite {
    
    
    /**
     * @class GMS2Trainer
     * @brief Train model parameters for GMS2
     *
     * This class trains model parameters used by GMS2, such as coding model, non-coding model, etc...
     */
    class GMS2Trainer {
        
    public:
        
        /**
         * Default constructor:
         */
        GMS2Trainer();
        
        /**
         * Train
         */
        void estimateParameters(const NumSequence &sequence, const vector<Label*> &labels);
        
        void estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersStartContext(const NumSequence &sequence, const vector<Label *> &labels);
        
        unsigned pcounts;
        unsigned codingOrder;
        unsigned noncodingOrder;
        unsigned startContextOrder;
        NumSequence::size_type startContextLength;
        
        // public variables for models
        NonUniformMarkov *motif;
        UniformMarkov *noncoding;
        PeriodicMarkov *coding;
        NonUniformMarkov *startContext;
        
        
    };
}

#endif /* GMS2Trainer_hpp */
