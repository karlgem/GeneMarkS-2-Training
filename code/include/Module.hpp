//
//  Module.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/17/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Module_hpp
#define Module_hpp

#include <stdio.h>


namespace gmsuite {
    
    /**
     * @class Module
     * @brief An (abstract) class than defines a generic module
     */
    class Module {
        
    public:
        
        /**
         * Run the module.
         */
        virtual void run() = 0;
        
    };
    
}


#endif /* Module_hpp */
