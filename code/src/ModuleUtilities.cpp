//
//  ModuleUtilities.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleUtilities.hpp"

#include <iostream>

#include "SequenceAlgorithms.hpp"
#include "UnivariatePDF.hpp"
#include "NonUniformMarkov.hpp"
#include "NonUniformCounts.hpp"
#include "UniformCounts.hpp"
#include "UniformMarkov.hpp"
#include "MotifFinder.hpp"
#include "LabelFile.hpp"
#include "SequenceFile.hpp"
#include "SequenceParser.hpp"
#include "GMS2Trainer.hpp"

using namespace std;
using namespace gmsuite;


// constructor
ModuleUtilities::ModuleUtilities(const OptionsUtilities& opt) : options(opt) {
    
}

// Run module according to the provided options.
void ModuleUtilities::run() {
    
    if (options.utility == EXTRACT_UPSTR) {
        runExtractUpstream();
    }
    else if (options.utility == START_MODEL_INFO)
        runStartModelInfo();
    else if (options.utility == MATCH_SEQ_TO_UPSTREAM)
        runMatchSeqToUpstream();
    else if (options.utility == MATCH_SEQ_TO_NONCODING)
        runMatchSeqToNoncoding();
    else            // unrecognized utility to run
        throw invalid_argument("Unknown utility function " + options.utility);
    
}



void ModuleUtilities::runExtractUpstream() {
    
    // read sequence file
    SequenceFile sequenceFile (options.extractUpstreamUtility.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (options.extractUpstreamUtility.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extract upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, options.extractUpstreamUtility.length, upstreams, options.extractUpstreamUtility.allowOverlaps, options.extractUpstreamUtility.minimumGeneLength);
    
    // convert upstream numeric sequences back to string
    vector<Sequence> strUpstreams (upstreams.size());
    
    for (size_t n = 0; n < upstreams.size(); n++) {
        strUpstreams[n] = Sequence(cnc.convert(upstreams[n].begin(), upstreams[n].end()));
    }
    
    
    // write sequences to file
    SequenceFile outputFile(options.extractUpstreamUtility.fn_output, SequenceFile::WRITE);
    outputFile.write(strUpstreams);
    
    
}











// a class to compute KL divergence (until I can figure out a better method of doing it)
class KLDivergence : public UniformMarkov, public NonUniformMarkov {
    
public:
    KLDivergence(const NonUniformMarkov *P, const UniformMarkov *Q) : NonUniformMarkov(*P), UniformMarkov(*Q)  {
        
        unsigned maxOrder = max(P->getOrder(), Q->getOrder());
        
        // raise order of background model to that of motif model
        UniformMarkov::changeOrder(maxOrder);
        NonUniformMarkov::changeOrder(maxOrder);
        
        if (this->UniformMarkov::order != this->NonUniformMarkov::order)
            throw logic_error("ScoreComputer: Orders of models should be the same.");
    }
    
    double computeKL() const {
        double score = 0;
        unsigned order = this->NonUniformMarkov::order;
        
        // for uniform model, get conditional pdf for every order <= "motif model order"
        vector<vector<double> > mBackConditionals (order+1);
        for (size_t n = 0; n < mBackConditionals.size()-1; n++) {
            mBackConditionals[n] = this->UniformMarkov::jointProbs[n];      // copy joint
            UniformMarkov::jointToMarkov(mBackConditionals[n]);             // convert joint to markov
        }
        
        // for last one, just copy it
        mBackConditionals[mBackConditionals.size()-1] = this->UniformMarkov::model;
        
        
        
        // for each position in motif model
        for (NonUniformMarkov::nonunif_markov_t::size_type pos = 0; pos < this->NonUniformMarkov::model.size(); pos++) {
            
            // get size of word (i.e. depending on order and position)
            size_t wordSize = order+1;
            if (pos < order)
                wordSize = pos+1;
            
            
            
            // for every word in that position
            for (size_t word = 0; word < this->NonUniformMarkov::model[pos].size(); word++) {
                double ratio = 0;
                if (mBackConditionals[wordSize-1][word] != 0) {
                    ratio = this->NonUniformMarkov::model[pos][word] / mBackConditionals[wordSize-1][word];
                }
                
                if (ratio != 0)
                    score += this->NonUniformMarkov::model[pos][word] * log(ratio);
            }
        }
        
        return score;
    }
    
};





void ModuleUtilities::runStartModelInfo() {
    
    // read sequence file
    SequenceFile sequenceFile (options.startModelInfoUtility.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (options.startModelInfoUtility.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    NumSequence numSequence (strSequence, cnc);
    
    // run training step
    const OptionsGMS2Training* optTrain = &options.startModelInfoUtility.optionsGMS2Training;
    GMS2Trainer trainer (optTrain->pcounts, optTrain->codingOrder, optTrain->noncodingOrder, optTrain->startContextOrder, optTrain->upstreamLength, optTrain->startContextLength, optTrain->genomeClass, optTrain->optionsMFinder, numAlph, optTrain->MIN_GENE_LEN);
    
    trainer.estimateParameters(numSequence, labels);
    
    
    
    
    // Generate non-coding sequences
    vector<NumSequence> simNonCoding (options.startModelInfoUtility.numOfSimNonCoding);
    
    for (size_t n = 0; n < simNonCoding.size(); n++) {
        simNonCoding[n] = trainer.noncoding->emit(optTrain->upstreamLength);
//        cout << cnc.convert(simNonCoding[n].begin(), simNonCoding[n].end()) << endl;
    }
    
    // convert align option from string format to align_t format
    MFinderModelParams::align_t align = MFinderModelParams::NONE;
    const OptionsMFinder *optionsMFinder = &options.startModelInfoUtility.optionsGMS2Training.optionsMFinder;
    if (optionsMFinder->align == "left")
        align = MFinderModelParams::LEFT;
    else if (optionsMFinder->align == "right")
        align = MFinderModelParams::RIGHT;
    else if (optionsMFinder->align == "none")
        align = MFinderModelParams::NONE;
    else
        throw invalid_argument("Align option must be one of: left, right, none");
    
    // set motif finder options
    MotifFinder::Builder b;
    b.setAlign(align).setWidth(optionsMFinder->width).setMaxIter(optionsMFinder->maxIter).setMaxEMIter(optionsMFinder->maxEMIter).setNumTries(optionsMFinder->tries);
    b.setPcounts(optionsMFinder->pcounts).setMotifOrder(optionsMFinder->motifOrder).setBackOrder(optionsMFinder->bkgdOrder).setShiftEvery(optionsMFinder->shiftEvery);
    
    // build motif finder from above options
    MotifFinder mfinder = b.build();
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(simNonCoding, positions);
    
    // build RBS model
    NonUniformCounts rbsCounts(optionsMFinder->motifOrder, optionsMFinder->width, numAlph);
    for (size_t n = 0; n < simNonCoding.size(); n++) {
        rbsCounts.count(simNonCoding[n].begin()+positions[n], simNonCoding[n].begin()+positions[n]+optionsMFinder->width);
    }
    
    NonUniformMarkov rbsSim(optionsMFinder->motifOrder, optionsMFinder->width, numAlph);
    rbsSim.construct(&rbsCounts, optionsMFinder->pcounts);
    
    // build spacer distribution
    // build histogram from positions
    vector<double> positionCounts (optTrain->upstreamLength - optionsMFinder->width+1, 0);
    for (size_t n = 0; n < positions.size(); n++) {
        // FIXME account for LEFT alignment
        // below is only for right
        positionCounts[optTrain->upstreamLength - optionsMFinder->width - positions[n]]++;        // increment position
    }
    
    UnivariatePDF rbsSpacerSim(positionCounts, false, optionsMFinder->pcounts);

    
    
    // compute KL of motif versus noncoding, and spacer versus uniform
    KLDivergence klDivergence(trainer.rbs, trainer.noncoding);
    double klMotif = klDivergence.computeKL();
    
    // compute kl of spacer vs uniform
    double klSpacer = 0;
    for (size_t n = 0; n < trainer.rbsSpacer->size(); n++) {
        double ratio = (*trainer.rbsSpacer)[n] / (1.0/optTrain->upstreamLength);
        
        if (ratio != 0)
            klSpacer += (*trainer.rbsSpacer)[n] * log2(ratio);
    }
    
    // compute KL of motif versus noncoding, and spacer versus uniform
    KLDivergence klDivergenceSim(&rbsSim, trainer.noncoding);
    double klMotifSim = klDivergenceSim.computeKL();
    
    // compute kl of spacer vs uniform
    double klSpacerSim = 0;
    for (size_t n = 0; n < rbsSpacerSim.size(); n++) {
        double ratio = rbsSpacerSim[n] / (1.0/optTrain->upstreamLength);
        
        if (ratio != 0)
            klSpacerSim += rbsSpacerSim[n] * log2(ratio);
    }
    
    
    
    
    cout << klMotif << "\t" << klSpacer << "\t" << klMotifSim << "\t" << klSpacerSim << endl;
    
//    cout << trainer.rbs->toString() << endl;
    

    
}




void ModuleUtilities::runMatchSeqToUpstream() {
    // read sequence file
    SequenceFile sequenceFile (options.matchSeqWithUpstream.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (options.matchSeqWithUpstream.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extract upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, options.matchSeqWithUpstream.length, upstreams, options.matchSeqWithUpstream.allowOverlaps, options.matchSeqWithUpstream.minimumGeneLength);
    
    // get sequence to match with
    Sequence strMatchSeq (options.matchSeqWithUpstream.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, upstreams[i]);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
    }
    

}



void ModuleUtilities::runMatchSeqToNoncoding() {
    
    // read sequence file
    SequenceFile sequenceFile (options.matchSeqWithNoncoding.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (options.matchSeqWithNoncoding.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    NumSequence numSequence (strSequence, cnc);
    
    // run training step
    const OptionsGMS2Training* optTrain = &options.matchSeqWithNoncoding.optionsGMS2Training;
    GMS2Trainer trainer (optTrain->pcounts, optTrain->codingOrder, optTrain->noncodingOrder, optTrain->startContextOrder, optTrain->upstreamLength, optTrain->startContextLength, optTrain->genomeClass, optTrain->optionsMFinder, numAlph, optTrain->MIN_GENE_LEN);
    
    trainer.estimateParameters(numSequence, labels);
    
    // Generate non-coding sequences
    vector<NumSequence> simNonCoding (options.matchSeqWithNoncoding.numOfSimNonCoding);
    
    for (size_t n = 0; n < simNonCoding.size(); n++) {
        simNonCoding[n] = trainer.noncoding->emit(optTrain->upstreamLength);
    }
    
    
    // get sequence to match with
    Sequence strMatchSeq (options.matchSeqWithNoncoding.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < simNonCoding.size(); i++) {
        
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, simNonCoding[i]);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
    }
    
}









