//
//  ModuleGMS2.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/16/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleGMS2.hpp"
#include <iostream>

#include "Sequence.hpp"
#include "SequenceFile.hpp"

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
    
    Sequence seq;
    
    // read single sequence from file
    try {
        SequenceFile seqFile (options.fname_in, SequenceFile::READ);        // open sequence file for reading
    
        seq = seqFile.read();          // read single sequence from file.
    
        cout << seq.toString() << endl;
    }
    // if file could not be opened
    catch (ios_base::failure &ex) {
        cerr << "Error: Could not open file: " << options.fname_in << endl;
        return;
    }
    
    AlphabetDNA alphabet;
    
    CharNumConverter cnc (&alphabet);
    
    // convert sequence to numeric form
    NumSequence numSeq (seq, cnc);
    
    for (NumSequence::size_type i = 0; i < numSeq.numSeq.size(); i++) {
        cout << numSeq[i];
    }
    
    cout << endl;
    
    // Create genetic code 11
    GeneticCode gc (GeneticCode::ELEVEN);
    
    cout << "ATG: " << gc.isStart("ATG") << endl;
    cout << "GTG: " << gc.isStart("GTG") << endl;
    cout << "TTG: " << gc.isStart("TTG") << endl;
    
    cout << "TAA: " << gc.isStop("TAA") << endl;
    cout << "TGA: " << gc.isStop("TGA") << endl;
    cout << "TAG: " << gc.isStop("TAG") << endl;
    
    // Create numeric genetic code 11
    NumGeneticCode gcNum (gc, cnc);
    
    
}










