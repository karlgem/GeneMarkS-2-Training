//
//  ModuleUtilities.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ModuleUtilities_hpp
#define ModuleUtilities_hpp

#include <stdio.h>

#include "Module.hpp"
#include "Sequence.hpp"
#include "OptionsUtilities.hpp"

namespace gmsuite {
    
    /**
     * @class ModuleUtilities
     * @brief A class than runs some utility modules
     */
    class ModuleUtilities : public Module {
        
    public:
        
        /**
         * Constructor: initialize the ModuleUtilities module with option parameters
         *
         * @param options the parameters defining how the module is to be run
         */
        ModuleUtilities(const OptionsUtilities& options);
        
        /**
         * Run module according to the provided options.
         */
        void run();
        
    private:
        
        const OptionsUtilities& options;         /**< Module option */
        
        // Each submodule (see OptionsUtilities) requires a separate run command
        void runExtractUpstream();
        void runStartModelInfo();
        void runEmitNonCoding();
        void runCountNumORF();
        void runMatchSeqToUpstream();
        void runMatchSeqToNoncoding();
        void runLabelsSimilarityCheck();
        void runExtractStartContextPerOperonStatus();
        void runExtractStartContextPerMotifStatus();
        void runComputeGC();
        void runSeparateFGIOAndIG();
        void runExtractStartContext();
        void runDNAToAA();
        void runChangeOrderNonCoding();
    };
    
}


#endif /* ModuleUtilities_hpp */
