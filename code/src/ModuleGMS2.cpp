//
//  ModuleGMS2.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/16/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleGMS2.hpp"
#include <iostream>

#include "Label.hpp"
#include "Sequence.hpp"
#include "SequenceFile.hpp"

#include "GMS2Trainer.hpp"
#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"
#include "GeneticCode.hpp"
#include "NumGeneticCode.hpp"
#include "CharNumConverter.hpp"


using namespace std;
using namespace gmsuite;

// constructor initializing the module's options
ModuleGMS2::ModuleGMS2(const OptionsGMS2& opt) : options(opt) {
}


// run GMS2 module
void ModuleGMS2::run() {
    
    const unsigned MAX_ITER = 10;
    
    // read single sequence from file
    Sequence seq;
    try {
        seq = readInputSequence(options.fname_in);
    }
    // if file could not be opened
    catch (ios_base::failure &ex) {
        cerr << "Error: Could not open file: " << options.fname_in << endl;
        return;
    }
    
    // read labels from file
    vector<Label*> labels;
    
    
    AlphabetDNA alphabet;                       // sequence alphabet
    GeneticCode gc (GeneticCode::ELEVEN);       // Create genetic code 11
    CharNumConverter cnc (&alphabet);           // converts characters to numbers
    
    // numeric versions of sequence and genetic code
    NumSequence numSeq (seq, cnc);              // convert character sequence to numeric form
    NumGeneticCode gcNum (gc, cnc);             // create numeric genetic code 11
    
    
    /*******************************\
     *      Step 1: Initiation     *
    \*******************************/
    
    // (1) Read heuristic parameters
    // (2) Use heuristic parameters to get initial parse (prediction) of the sequence
    // (3) Classific genome into one of 3 types
    
    
    
    
    /*******************************\
     *      Step 2: Main Cycle     *
    \*******************************/
    
    // Perform MAX_ITER iterations, where
    //  (a) Perform motif search based on latest predicted starts
    //  (b) Estimate model parameters for coding, non-coding, motif models, ...
    //  (c) Pass parameters to HMM predictor, and get back new genome predictions
    //  (d) If convergence reached, end loop
    for (size_t iter = 0; iter < MAX_ITER; iter++) {
        GMS2Trainer trainer;
        trainer.estimateParameters(numSeq, labels);
    }
    
    
    
    
    /*******************************\
     *  Step 3: Adaptive Training  *
    \*******************************/
    
    // Repeat:
    //  (a) Generate non-coding sequences from non-coding model
    //  (b) Perform gene prediction and compute number of false positives
    //  (c) Update parameters based on threshold
    
    
    
    
    /*******************************\
     *      Step 4: Output         *
    \*******************************/
    
    // print output
    
    
    
    
}




Sequence ModuleGMS2::readInputSequence(string filename) const {
    Sequence seq;
    
    // read single sequence from file
    SequenceFile seqFile (options.fname_in, SequenceFile::READ);        // open sequence file for reading
    return seqFile.read();                                              // read single sequence from file.
}





