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

#include "Label.hpp"
#include "Module.hpp"
#include "Sequence.hpp"
#include "GMS2Trainer.hpp"
#include "NumSequence.hpp"
#include "OptionsGMS2.hpp"
#include "CharNumConverter.hpp"
#include "ProkGeneStartModel.hpp"

namespace gmsuite {
    
    /**
     * @class ModuleGMS2
     * @brief A class than runs the GeneMarkS2 module
     */
    class ModuleGMS2 : public Module {
        
    public:
        
        typedef ProkGeneStartModel::genome_group_t genome_group_t;
        
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
        
        
        /**
         * Estimate parameters for GMS2 where, given parsing (labels) of the sequence and the genome class,
         * parameters are estimated for the corresponding models
         */
        void estimateModelParameters(const NumSequence &sequence, const vector<Label*> &labels, genome_group_t genomeClass, GMS2Trainer &trainer);
        
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
        
        
        
        genome_group_t classifyGenome(const NumSequence &numSeq, const CharNumConverter &cnc, const vector<Label*> labels, NumSequence::size_type upstrLength) const;
        
    };
    
}

#endif /* ModuleGMS2_hpp */
