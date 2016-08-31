//
//  ModuleGMS2.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/16/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ModuleGMS2_hpp
#define ModuleGMS2_hpp

#include <stdio.h>

#include "Module.hpp"
#include "Sequence.hpp"
#include "OptionsGMS2.hpp"

namespace gmsuite {
    
    /**
     * @class ModuleGMS2
     * @brief A class than runs the GeneMarkS2 module
     */
    class ModuleGMS2 : public Module {
        
    public:
        
        /**
         * Constructor: initialize the GMS2 module with option parameters
         *
         * @param options the parameters defining hwo the module is to be run
         */
        ModuleGMS2(const OptionsGMS2& options);
        
        /**
         * Run the GMS2 module according to the provided options.
         */
        void run();
        
    private:
        
        const OptionsGMS2& options;         /**< Module option */
        
        
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

#endif /* ModuleGMS2_hpp */
