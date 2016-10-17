//
//  ModuleExperiment.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/14/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleExperiment.hpp"

#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"
#include "CharNumConverter.hpp"
#include "SequenceFile.hpp"
#include "LabelFile.hpp"
#include "SequenceParser.hpp"
#include "SequenceAlgorithms.hpp"
#include "MotifFinder.hpp"
#include "GMS2Trainer.hpp"

#include <iostream>

using namespace std;
using namespace gmsuite;


ModuleExperiment::ModuleExperiment(const OptionsExperiment& opt) : options(opt) {
    
}



void ModuleExperiment::run() {
    
    if (options.experiment == OptionsExperiment::MATCH_SEQ_TO_UPSTREAM)
        runMatchSeqToUpstream();
    else if (options.experiment == OptionsExperiment::MATCH_SEQ_TO_NONCODING)
        runMatchSeqToNoncoding();
    else if (options.experiment == OptionsExperiment::BuildStartModelsOptions)
        runBuildStartModels();
    
}

// match query to upstream regions
void ModuleExperiment::runMatchSeqToUpstream() {
    
    OptionsExperiment::MatchSeqToUpstreamOptions expOptions = options.matchSeqToUpstream;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extractlk upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreams, expOptions.allowOverlaps, expOptions.minGeneLength);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    unsigned matchThresh = 4;            // threshold for nonmatches
    vector<NumSequence> nonMatch;           // keep track of 'non matching'
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        NumSequence sub = upstreams[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, sub);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(upstreams[i]);
    }
    
    // run motif search on non-match
    MotifFinder::Builder b;
    b.setAlign(MFinderModelParams::RIGHT);
    
    MotifFinder mfinder = b.build();
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonMatch, positions);
    
    cout << "The number of remaining sequences is: " << nonMatch.size() << endl;
    for (size_t n = 0; n < nonMatch.size(); n++) {
        cout << cnc.convert(nonMatch[n].begin() + positions[n], nonMatch[n].begin() + positions[n] + 6);
        cout << "\t" << positions[n] + 1 << "\t" << nonMatch[n].size() << endl;
    }
}

// match query to simulated noncoding sequences
void ModuleExperiment::runMatchSeqToNoncoding() {
    
    OptionsExperiment::MatchSeqToNoncodingOptions expOptions = options.matchSeqToNoncoding;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    NumSequence numSequence (strSequence, cnc);
    
    OptionsMFinder optionsMFinder;
    
    // run training step
    GMS2Trainer trainer((unsigned) expOptions.pcounts, 0, expOptions.order, 0, 40, 0, ProkGeneStartModel::C1, optionsMFinder, numAlph, 300);
    
    trainer.estimateParameters(numSequence, labels);
    
    // Generate non-coding sequences
    vector<NumSequence> simNonCoding (expOptions.numNoncoding);
    
    for (size_t n = 0; n < simNonCoding.size(); n++)
        simNonCoding[n] = trainer.noncoding->emit(expOptions.length);
    
    // get sequence to match with
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    // for each noncoding sequence, match it
    for (size_t i = 0; i < simNonCoding.size(); i++) {
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, simNonCoding[i]);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
    }
        
}





// match query to upstream regions
void ModuleExperiment::runBuildStartModels() {
    
    OptionsExperiment::BuildStartModelsOptions expOptions = options.buildStartModels;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extractlk upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreams, expOptions.allowOverlaps, expOptions.minGeneLength);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    unsigned matchThresh = expOptions.min16SMatch;          // threshold for nonmatches
    vector<NumSequence> nonMatch;                           // keep track of 'non matching'
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        NumSequence sub = upstreams[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, sub);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(upstreams[i]);
    }
    
    // run motif search on non-match
    MotifFinder::Builder b;
    b.build(expOptions.mfinderOptions);
    
    MotifFinder mfinder = b.build();
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonMatch, positions);
    
    cout << "The number of remaining sequences is: " << nonMatch.size() << endl;
    for (size_t n = 0; n < nonMatch.size(); n++) {
        cout << cnc.convert(nonMatch[n].begin() + positions[n], nonMatch[n].begin() + positions[n] + 6);
        cout << "\t" << positions[n] + 1 << "\t" << nonMatch[n].size() << endl;
    }
}




































