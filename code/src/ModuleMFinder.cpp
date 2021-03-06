//
//  ModuleMFinder.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleMFinder.hpp"

#include "MotifFinder.hpp"
#include "SequenceFile.hpp"

#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"
#include "GeneticCode.hpp"
#include "NumGeneticCode.hpp"
#include "CharNumConverter.hpp"

#include <iostream>

using namespace std;
using namespace gmsuite;

ModuleMFinder::ModuleMFinder(const OptionsMFinder& opt) : options(opt) {
    
}

void ModuleMFinder::run() {
    
    
    // set motif finder options
    MotifFinder::Builder b;
    b.setAlign(options.align).setWidth(options.width).setMaxIter(options.maxIter).setMaxEMIter(options.maxEMIter).setNumTries(options.tries);
    b.setPcounts(options.pcounts).setMotifOrder(options.motifOrder).setBackOrder(options.bkgdOrder).setShiftEvery(options.shiftEvery);
    
    // build motif finder from above options
    MotifFinder mfinder = b.build();
    
    
    // read sequences from file
    vector<Sequence> sequences;
    try {
        // open sequence file for reading
        SequenceFile seqFile (options.fname_in, SequenceFile::READ);
        seqFile.read(sequences);                  // read single sequence from file.
        
    }
    // if file could not be opened
    catch (ios_base::failure &ex) {
        cerr << "Error: Could not open file: " << options.fname_in << endl;
        return;
    }
    
    
    AlphabetDNA alphabet;                       // sequence alphabet
    GeneticCode gc (GeneticCode::ELEVEN);       // Create genetic code 11
    CharNumConverter cnc (&alphabet);           // converts characters to numbers
    
    // numeric versions of sequences and genetic code
    vector<NumSequence> numSequences (sequences.size());
    for (size_t n = 0; n < sequences.size(); n++)
        numSequences[n] = NumSequence(sequences[n], cnc);
    
    NumGeneticCode gcNum (gc, cnc);             // create numeric genetic code 11
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(numSequences, positions);
    
    printf("NAME MotifFinder\n");
    printf("VERSION 1.0\n");
    
    // print positions
    for (size_t n = 0; n < numSequences.size(); n++) {
        cout << cnc.convert(numSequences[n].begin() + positions[n], numSequences[n].begin() + positions[n] + options.width);
        cout << "\t" << positions[n] + 1 << "\t" << numSequences[n].size() << endl;
    }
}
