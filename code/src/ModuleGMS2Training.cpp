//
//  ModuleGMS2Training.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleGMS2Training.hpp"
#include "GMS2Trainer.hpp"

#include "ModelFile.hpp"
#include "LabelFile.hpp"
#include "SequenceFile.hpp"

using namespace gmsuite;

// constructor initializing the module's options
ModuleGMS2Training::ModuleGMS2Training(const OptionsGMS2Training &opt) : options(opt) {
    
}


void ModuleGMS2Training::run() {
    
    AlphabetDNA alph;
    GeneticCode geneticCode (GeneticCode::ELEVEN);
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    // read sequence from file
    SequenceFile seqFile(options.fn_sequence, SequenceFile::READ);
    Sequence strSeq = seqFile.read();
    
    // get numeric sequence
    NumSequence sequence(strSeq, cnc);
    
    // read labels from file
    LabelFile labFile(options.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labFile.read(labels);
    
    
    // set up trainer
    GMS2Trainer trainer (options.pcounts, options.codingOrder, options.noncodingOrder, options.startContextOrder, options.upstreamLength, options.startContextLength, options.genomeClass, options.optionsMFinder, numAlph, options.MIN_GENE_LEN, numGeneticCode, options.startContextMargin);
    
    trainer.estimateParameters(sequence, labels);
    
    // write parameters to file
    map<string, string> toMod;
    trainer.toModFile(toMod);
    
    ModelFile modFile(options.fn_outmod, ModelFile::WRITE);
    modFile.write(toMod, "NATIVE");
    
    
    // free up memory for labels
    for (size_t n = 0; n < labels.size(); n++)
        if (labels[n] != NULL)
            delete labels[n];
}

