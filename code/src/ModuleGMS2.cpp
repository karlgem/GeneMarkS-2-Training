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

#include "UnivariatePDF.hpp"
#include "MotifFinder.hpp"
#include "GMS2Trainer.hpp"
#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"
#include "GeneticCode.hpp"
#include "NumGeneticCode.hpp"
#include "CharNumConverter.hpp"


using namespace std;
using namespace gmsuite;


// Extract upstream sequence: throws exception when 'not enough' sequence to extract
NumSequence extractUpstreamSequence(const NumSequence& sequence, const Label &label, const CharNumConverter &cnc, NumSequence::size_type upstrLength) {
    
    if (upstrLength == 0)
        throw out_of_range("Cannot extract upstream sequence of length 0");
    
    if (label.strand == Label::POS) {           // positive strand
        
        if (label.left < upstrLength)
            //throw out_of_range("Not enough sequence for position " + to_string(label.left) + " to extract upstream of length " + to_string(upstrLength));
            throw out_of_range("Not enough upstream sequence.");

        size_t left = label.left - upstrLength;                     // left idx of upstream sequence
        size_t right = label.left-1;                                // right idx of upstream sequence
        size_t length = right - left + 1;                           // length of upstream sequence
        
        NumSequence upstream = sequence.subseq(left, length);       // get upstream subsequence
        return upstream;
    }
    else {                                      // negative strand; reverse complement
        
        if (label.right + upstrLength >= sequence.size())
            //throw out_of_range("Not enought sequence for position " + to_string(label.right) + " on negative strand to extract upstream of length " + to_string(upstrLength));
            throw out_of_range("Not enough upstream sequence.");
        
        size_t left = label.right + 1;                              // left idx of upstream sequence
        size_t right = label.right + upstrLength;                   // right idx of upstream sequence
        size_t length = right - left + 1;                           // length of upstream sequence
        
        NumSequence upstream = sequence.subseq(left, length);       // get upstream subsequence
        upstream.reverseComplement(cnc);                            // reverse complement
        return upstream;
    }
}

// Classify genome
ModuleGMS2::genome_class_t ModuleGMS2::classifyGenome(const NumSequence &numSeq, const CharNumConverter &cnc, const vector<Label*> labels, NumSequence::size_type upstrLength) const {
    
    // extract upstream region, and reverse complement when on negative strand
    vector<NumSequence> upstreamRegions (labels.size());
    
    size_t numSkipped = 0;          // number of labels skipped: i.e. won't be used in genome classification (because of no upstream sequence)
    
    for (size_t n = 0; n < labels.size(); n++) {
        try {
            upstreamRegions[n-numSkipped] = extractUpstreamSequence(numSeq, *labels[n], cnc, upstrLength);
        }
        catch (out_of_range) {
            numSkipped++;
        }
    }
    
    // resize vector to remove skipped (empty) numsequences
    upstreamRegions.resize(upstreamRegions.size()-numSkipped);
    
    
    // Now that we have the upstream regions, perform motif search
    
    // start by building motif finder
    // convert align option from string format to align_t format
    MFinderModelParams::align_t align = MFinderModelParams::NONE;
    if (options.optionsMFinder.align == "left")
        align = MFinderModelParams::LEFT;
    else if (options.optionsMFinder.align == "right")
        align = MFinderModelParams::RIGHT;
    else if (options.optionsMFinder.align == "none")
        align = MFinderModelParams::NONE;
    else
        throw invalid_argument("Align option must be one of: left, right, none");
    
    MotifFinder::Builder b;
    b.setAlign(align).setWidth(options.optionsMFinder.width).setMaxIter(options.optionsMFinder.maxIter).setMaxEMIter(options.optionsMFinder.maxEMIter).setNumTries(options.optionsMFinder.tries);
    b.setPcounts(options.optionsMFinder.pcounts).setMotifOrder(options.optionsMFinder.motifOrder).setBackOrder(options.optionsMFinder.bkgdOrder).setShiftEvery(options.optionsMFinder.shiftEvery);

    MotifFinder mfinder = b.build();        // build mfinder
    
    // perform motif search
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(upstreamRegions, positions);
    
    // build histogram from positions
    vector<double> positionCounts (upstrLength - options.optionsMFinder.width+1, 0);
    for (size_t n = 0; n < positions.size(); n++) {
        // FIXME account for LEFT alignment
        // below is only for right
        positionCounts[upstrLength - options.optionsMFinder.width - positions[n]]++;        // increment position
    }
    
    UnivariatePDF spacerDistribution(positionCounts);
    
    genome_class_t genomeClass;
    
    // TODO: analyze spacer and get class
    double maxProb = 0;
    size_t posOfMax = 0;
    bool maxFound = false;
    
    // find max prob
    for (size_t n = 0; n < spacerDistribution.size(); n++) {
        if (spacerDistribution[n] > maxProb) {
            maxProb = spacerDistribution[n];
            posOfMax = n;
            maxFound = true;
        }
    }
    
    // default is class 1
    if (!maxFound)
        return ProkGeneStartModel::C1;
    
    // decide start class
    if (maxProb > options.CLASS_PROB_THRESHOLD) {
        
        // Class 1: rbs
        if (posOfMax < options.CLASS_DIST_THRESHOLD)
            genomeClass = ProkGeneStartModel::C1;
        
        // Class 2: promoter
        else
            genomeClass = ProkGeneStartModel::C3;
        
    }
    // Class 3
    else
        genomeClass = ProkGeneStartModel::C2;
        
    
    
    
    return genomeClass;
}


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
    // (3) Classify genome into one of 3 types
    
    // Step 2: Get initial parse from heuristic parameters          // TODO: set initial
    vector<Label*> initialParse;
    
    // Step 3: Classify genome into type
    genome_class_t genomeClass = classifyGenome(numSeq, cnc, initialParse, options.upstrLength);
    
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





