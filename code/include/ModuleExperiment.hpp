//
//  ModuleExperiment.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/14/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ModuleExperiment_hpp
#define ModuleExperiment_hpp

#include <stdio.h>

#include "Module.hpp"
#include "OptionsExperiment.hpp"

namespace gmsuite {
    
    /**
     * @class ModuleExperiment
     * @brief A class that runs some experiments
     */
    class ModuleExperiment : public Module {
        
    public:
        
        /**
         * Constructor: initialize the ModuleExperiment module with optino parameters
         *
         * @param options the parameters defining how experiments are to be run
         */
        ModuleExperiment(const OptionsExperiment& options);
        
        /**
         * Run module according to the provided options
         */
        void run();
        
    private:
        
        const OptionsExperiment& options;       /**< Module options */
        
        // Each experiment requires a separate run command
        void runMatchSeqToUpstream();
        void runMatchSeqToNoncoding();
        void runBuildStartModels();
        void runBuildStartModels2();
        void runBuildStartModels3();
        void runMatchRBSTo16S();
        void runPromoterIsValidForAchaea();
        void runPromoterIsValidForBacteria();
        
        void runStartModelStrategy2();
        
        void runScoreStarts();
    };
}

#endif /* ModuleExperiment_hpp */
