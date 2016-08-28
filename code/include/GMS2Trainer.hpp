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
        
        
        
    };
}

#endif /* GMS2Trainer_hpp */
