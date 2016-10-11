//
//  test_GMS2Trainer.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/14/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include <stdio.h>

#include "catch.hpp"
#include "Label.hpp"
#include "LabelFile.hpp"
#include "GMS2Trainer.hpp"
#include "SequenceFile.hpp"

using namespace gmsuite;

TEST_CASE("Test GMS2 Trainer") {
    
    // setup mfinder options
    OptionsMFinder optionsMFinder("mfinder");
    optionsMFinder.align = "right";
    optionsMFinder.width = 6;
    optionsMFinder.motifOrder = 0;
    optionsMFinder.bkgdOrder = 0;
    optionsMFinder.pcounts = 1;
    optionsMFinder.tries = 10;
    optionsMFinder.maxIter = 60;
    optionsMFinder.shiftEvery = 20;
    optionsMFinder.maxEMIter = 10;
    
    
    // setup trainer options
    GMS2Trainer trainer;
    trainer.pcounts = 1;;
    trainer.codingOrder = 4;
    trainer.noncodingOrder = 2;
    trainer.startContextOrder = 0;
    trainer.upstreamLength = 40;
    trainer.startContextLength = 12;
    trainer.optionsMFinder = &optionsMFinder;
    
    trainer.genomeClass = ProkGeneStartModel::C1;       // genome class 1 for ecoli
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    trainer.alphabet = &numAlph;
    
    // read sequence from file
    SequenceFile seqFile ("/Users/Karl/repos/GeneMarkS-2/code/tmp/ecoli.fa", SequenceFile::READ);        // open sequence file for reading
    Sequence strSequence = seqFile.read();                                                          // read single sequence from file.
    
    NumSequence sequence(strSequence, cnc);
    
    // read labels from file
    vector<Label*> labels;
    LabelFile labelFile("/Users/Karl/repos/GeneMarkS-2/code/tmp/ecoli.lst", LabelFile::READ);
    labelFile.read(labels);
    
    trainer.estimateParameters(sequence, labels);
    
    
}



























