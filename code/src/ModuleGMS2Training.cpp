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
#include <iostream>

using namespace std;
using namespace gmsuite;

// constructor initializing the module's options
ModuleGMS2Training::ModuleGMS2Training(const OptionsGMS2Training &opt) : options(opt) {
    
}


void ModuleGMS2Training::run() {
    
    AlphabetDNA alph;
    GeneticCode geneticCode (options.geneticCode);
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
    GMS2Trainer trainer (options.pcounts, options.codingOrder, options.noncodingOrder, options.startContextOrder, options.upstreamLength, options.startContextLength, options.genomeClass, options.optionsMFinder, numAlph, options.MIN_GENE_LEN, numGeneticCode, options.startContextMargin, options.trainOnNativeOnly, options.runMotifSearch, options.upstrLenFGIO, options.widthArchaeaPromoter);
    
    trainer.estimateParameters(sequence, labels);
    
    // get parameters from training
    vector<pair<string, string> > toMod;
    trainer.toModFile(toMod, options);
    
    // write settings that were given in from the settings file
    if (!options.fn_settings.empty()) {
        ModelFile settingsMFile(options.fn_settings, ModelFile::READ);
        
        map<string, string> settings;
        settingsMFile.read(settings);
        
        for (map<string, string>::const_iterator iter = settings.begin(); iter != settings.end(); iter++) {
            // search for key
            vector<pair<string, string> >::iterator key;
            for (key = toMod.begin(); key != toMod.end(); key++) {
                if (key->first == iter->first)      // if key found
                    break;
            }
            
            // if key not found, add it
            if (key == toMod.end()) {
                toMod.insert(toMod.begin()+1, pair<string, string>(iter->first, iter->second));       // add element to top of list
            }
            // if key found, replace value
            else {
                key->second = iter->second;
            }
        }
    }
    
    // write parameters to file
    ModelFile modFile(options.fn_outmod, ModelFile::WRITE);
    modFile.write(toMod, "NATIVE");
    
    
    // free up memory for labels
    for (size_t n = 0; n < labels.size(); n++)
        if (labels[n] != NULL)
            delete labels[n];
}

