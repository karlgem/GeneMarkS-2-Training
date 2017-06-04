//
//  ModuleMFinder.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ModuleMFinder_hpp
#define ModuleMFinder_hpp

#include <stdio.h>

#include "Module.hpp"
#include "Sequence.hpp"
#include "OptionsMFinder.hpp"

namespace gmsuite {
    
    /**
     * @class ModuleMFinder
     * @brief A class than runs the MFinder module for motif search
     */
    class ModuleMFinder : public Module {
        
    public:
        
        /**
         * Constructor: initialize the MFinder module with option parameters
         *
         * @param options the parameters defining hwo the module is to be run
         */
        ModuleMFinder(const OptionsMFinder& options);
        
        /**
         * Run the GMS2 module according to the provided options.
         */
        void run();
        
    private:
        
        const OptionsMFinder& options;         /**< Module option */
        
        
        /**
         * Read input sequence from file.
         * TODO: For now, this assumes a single FASTA sequence in the file. This should
         * be updated to handle genomes (i.e. multiple chromosomes, plasmids, etc...)
         *
         * @param filename the name of the file containing the sequence
         */
        Sequence readInputSequence(string filename) const;
        
        
        
    };
    
}

#endif /* ModuleMFinder_hpp */
