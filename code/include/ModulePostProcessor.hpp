//
//  ModulePostProcessor.hpp
//  Biogem CPP
//
//  Created by Karl Gemayel on 9/18/18.
//  Copyright © 2018 Georgia Institute of Technology. All rights reserved.
//

#ifndef ModulePostProcessor_hpp
#define ModulePostProcessor_hpp

#include <stdio.h>

#include "Module.hpp"
#include "OptionsPostProcessor.hpp"

namespace gmsuite {
    
    /**
     * @class ModulePostProcessor
     * @brief A class than runs the Post-Processor module for updating predictions
     */
    class ModulePostProcessor : public Module {
        
    public:
        
        /**
         * Constructor: initialize the Post-processor module with option parameters
         *
         * @param options the parameters defining how the module is to be run
         */
        ModulePostProcessor(const OptionsPostProcessor& options);
        
        /**
         * Run the Post-processor module according to the provided options.
         */
        void run();
        
    private:
        
        const OptionsPostProcessor& options;         /**< Module option */
        
        
        
        
        
    };
    
}


#endif /* ModulePostProcessor_hpp */
