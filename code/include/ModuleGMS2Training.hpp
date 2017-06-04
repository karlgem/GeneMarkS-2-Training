//
//  ModuleGMS2Training.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ModuleGMS2Training_hpp
#define ModuleGMS2Training_hpp

#include <stdio.h>
#include "Module.hpp"
#include "Sequence.hpp"
#include "OptionsGMS2Training.hpp"

namespace gmsuite {
    
    /**
     * @class ModuleGMS2Training
     * @brief A class than runs the GMS2Training module for estimating parameters for GMS2
     */
    class ModuleGMS2Training : public Module {
        
    public:
        
        /**
         * Constructor: initialize the MFinder module with option parameters
         *
         * @param options the parameters defining how the module is to be run
         */
        ModuleGMS2Training(const OptionsGMS2Training& options);
        
        /**
         * Run the GMS2 module according to the provided options.
         */
        void run();
        
    private:
        
        const OptionsGMS2Training& options;         /**< Module option */
        
        
        
        
        
    };
    
}


#endif /* ModuleGMS2Training_hpp */
