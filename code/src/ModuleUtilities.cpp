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
#include "ModelFile.hpp"
#include "NonCodingMarkov.hpp"
#include "MotifMarkov.hpp"
#include "LabelsParser.hpp"
#include <time.h>

using namespace std;
using namespace gmsuite;


// constructor
ModuleUtilities::ModuleUtilities(const OptionsUtilities& opt) : options(opt) {
    
}

// Run module according to the provided options.
void ModuleUtilities::run() {
    
    if (options.utility == OptionsUtilities::EXTRACT_UPSTR) {
        runExtractUpstream();
    }
    else if (options.utility == OptionsUtilities::START_MODEL_INFO)
        runStartModelInfo();
    else if (options.utility == OptionsUtilities::MATCH_SEQ_TO_UPSTREAM)
        runMatchSeqToUpstream();
    else if (options.utility == OptionsUtilities::MATCH_SEQ_TO_NONCODING)
        runMatchSeqToNoncoding();
    else if (options.utility == OptionsUtilities::LABELS_SIMILARITY_CHECK)
        runLabelsSimilarityCheck();
    else if (options.utility == OptionsUtilities::EMIT_NON_CODING)
        runEmitNonCoding();
    else if (options.utility == OptionsUtilities::COUNT_NUM_ORF)
        runCountNumORF();
    else if (options.utility == OptionsUtilities::EXTRACT_SC_PER_OPERON_STATUS)
        runExtractStartContextPerOperonStatus();
    else if (options.utility == OptionsUtilities::EXTRACT_SC_PER_MOTIF_STATUS)
        runExtractStartContextPerMotifStatus();
    else if (options.utility == OptionsUtilities::COMPUTE_GC)
        runComputeGC();
    else if (options.utility == OptionsUtilities::SEPARATE_FGIO_AND_IG)
        runSeparateFGIOAndIG();
    else if (options.utility == OptionsUtilities::EXTRACT_START_CONTEXT)
        runExtractStartContext();
    else if (options.utility == OptionsUtilities::DNA_TO_AA)
        runDNAToAA();
    else if (options.utility == OptionsUtilities::CHANGE_ORDER_NONCODING)
        runChangeOrderNonCoding();
    else if (options.utility == OptionsUtilities::COMPUTE_KL)
        runComputeKL();
    else if (options.utility == OptionsUtilities::AB_FILTER)
        runABFilter();
    else if (options.utility == OptionsUtilities::EXTRACT_SPACER_NT_MODEL)
        runExtractSpacerNTModel();
    else if (options.utility == OptionsUtilities::EXTRACT_LORF)
        runExtractLORF();
    
//    else            // unrecognized utility to run
//        throw invalid_argument("Unknown utility function " + options.utility);
    
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




void ModuleUtilities::runLabelsSimilarityCheck() {
    
    // open label files
    LabelFile labelFileA (options.labelsSimilarityCheck.fn_labelsA, LabelFile::READ);
    LabelFile labelFileB (options.labelsSimilarityCheck.fn_labelsB, LabelFile::READ);
    
    // read labels from files
    vector<Label*> labelsA, labelsB;
    labelFileA.read(labelsA);
    labelFileB.read(labelsB);
    
    // count number of similar stops
    size_t identicalStops = 0;
    size_t identicalStarts = 0;
    
    // loop over first set: A
    for (vector<Label*>::const_iterator a = labelsA.begin(); a != labelsA.end(); a++) {
        
        Label::strand_t strandA = (*a)->strand;
        size_t startA = (strandA == Label::POS ? (*a)->left  : (*a)->right);
        size_t stopA  = (strandA == Label::POS ? (*a)->right : (*a)->left);
        
        // loop over second set: B
        for (vector<Label*>::const_iterator b = labelsB.begin(); b != labelsB.end(); b++) {
            
            Label::strand_t strandB = (*b)->strand;
            size_t startB = (strandB == Label::POS ? (*b)->left  : (*b)->right);
            size_t stopB  = (strandB == Label::POS ? (*b)->right : (*b)->left);
            
            if (strandA == strandB) {
                if (startA == startB)
                    identicalStarts++;
                if (stopA == stopB)
                    identicalStops++;
            }
        }
    }
    
    double similarity = 0;
    if (labelsA.size() + labelsB.size() != 0)
        similarity = (2.0*identicalStarts) / (labelsA.size() + labelsB.size());
    
    cout << similarity << endl;
    
    
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
                    score += this->NonUniformMarkov::model[pos][word] * log2(ratio);            
            }
//            cout << endl;
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
    GeneticCode geneticCode(GeneticCode::ELEVEN);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    NumAlphabetDNA numAlph(alph, cnc);
    
    NumSequence numSequence (strSequence, cnc);
    
    // run training step
    const OptionsGMS2Training* optTrain = &options.startModelInfoUtility.optionsGMS2Training;
    GMS2Trainer trainer (optTrain->pcounts, optTrain->codingOrder, optTrain->noncodingOrder, optTrain->startContextOrder, optTrain->upstreamLength, optTrain->startContextLength, optTrain->genomeGroup, optTrain->optionsMFinder, numAlph, optTrain->MIN_GENE_LEN, numGeneticCode, optTrain->startContextMargin, optTrain->trainOnNativeOnly);
    
    trainer.estimateParameters(numSequence, labels);
    
    
    
    
    // Generate non-coding sequences
    vector<NumSequence> simNonCoding (options.startModelInfoUtility.numOfSimNonCoding);
    
    for (size_t n = 0; n < simNonCoding.size(); n++) {
        simNonCoding[n] = trainer.noncoding->emit(optTrain->upstreamLength);
//        cout << cnc.convert(simNonCoding[n].begin(), simNonCoding[n].end()) << endl;
    }
    
    // convert align option from string format to align_t format
    
    const OptionsMFinder *optionsMFinder = &options.startModelInfoUtility.optionsGMS2Training.optionsMFinder;
    
    
    // set motif finder options
    MotifFinder::Builder b;
    b.setAlign(optionsMFinder->align).setWidth(optionsMFinder->width).setMaxIter(optionsMFinder->maxIter).setMaxEMIter(optionsMFinder->maxEMIter).setNumTries(optionsMFinder->tries);
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
    
    vector<NumSequence> nonRBS;
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, upstreams[i]);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
        
        if (match.size() <= 2)
            nonRBS.push_back(upstreams[i]);
    }
    
    
    MotifFinder::Builder b;
    b.setAlign(MFinderModelParams::RIGHT);
    
    MotifFinder mfinder = b.build();
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonRBS, positions);
    
    cout << "The number of remaining sequences is: " << nonRBS.size() << endl;
    // print positions
    for (size_t n = 0; n < nonRBS.size(); n++) {
        cout << cnc.convert(nonRBS[n].begin() + positions[n], nonRBS[n].begin() + positions[n] + 6);
        cout << "\t" << positions[n] + 1 << "\t" << nonRBS[n].size() << endl;
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
    GeneticCode geneticCode(GeneticCode::ELEVEN);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    
    NumSequence numSequence (strSequence, cnc);
    
    // run training step
    const OptionsGMS2Training* optTrain = &options.matchSeqWithNoncoding.optionsGMS2Training;
    GMS2Trainer trainer (optTrain->pcounts, optTrain->codingOrder, optTrain->noncodingOrder, optTrain->startContextOrder, optTrain->upstreamLength, optTrain->startContextLength, optTrain->genomeGroup, optTrain->optionsMFinder, numAlph, optTrain->MIN_GENE_LEN, numGeneticCode, optTrain->startContextMargin, optTrain->trainOnNativeOnly);
    
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


void ModuleUtilities::runEmitNonCoding() {
    
    srand(time(NULL));
    
    // read model file
    ModelFile mfile (options.emitNonCoding.fn_mod, ModelFile::READ);
    
    // set keys for non-coding model
    vector<string> keys;
    keys.push_back("NON_MAT");
    keys.push_back("NON_ORDER");
    
    // get key-value pair for non-coding model
    map<string, string> mKeyValuePair;
    mfile.read(mKeyValuePair, keys);
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // break non-coding string into vector of codon-probability pairs
    vector<pair<string, double> > nonCodingProbs;
    istringstream ssm(mKeyValuePair["NON_MAT"]);
    string line;
    while (std::getline(ssm, line)) {
        
        stringstream lineStream (line);
        string codon;
        double prob;
        
//        std::getline(lineStream, line, '\t');
        
        lineStream >> codon >> prob;
        
        nonCodingProbs.push_back(pair<string, double> (codon, prob));
    }
    
    // build non-coding model
    NonCodingMarkov nonCodingMarkov(nonCodingProbs, numAlph, cnc);
    
    // if new-order option set, change model order to given value
    if (options.emitNonCoding.order != numeric_limits<int>::infinity()) {
        nonCodingMarkov.changeOrder(options.emitNonCoding.order);
    }
    
    // emit noncoding sequence
    NumSequence emittedNumSeq = nonCodingMarkov.emit(options.emitNonCoding.length);
    
    // get string form of sequence
    Sequence emittedSeq = Sequence(cnc.convert(emittedNumSeq.begin(), emittedNumSeq.end()));
    
    // print emitted sequence to file
    SequenceFile sfile (options.emitNonCoding.fn_out, SequenceFile::WRITE);
    vector<Sequence> v;
    v.push_back(emittedSeq);
    sfile.write(v);
}


void ModuleUtilities::runCountNumORF() {
    
    OptionsUtilities::CountNumORF utilOpt = options.countNumORF;
    
    // read model file
    ModelFile mfile (utilOpt.fn_mod, ModelFile::READ);
    
    // set keys for useful parameters
    vector<string> keys;
    keys.push_back("GCODE");
    keys.push_back("GENE_MIN_LENGTH");
    
    // get key-value pair for parameters
    map<string, string> mKeyValuePair;
    mfile.read(mKeyValuePair, keys);
    
    
    // read sequence file
    SequenceFile sfile(utilOpt.fn_sequence, SequenceFile::READ);
    Sequence seq = sfile.read();
    
    GeneticCode::gcode_t genCodeValue;
    if (mKeyValuePair["GCODE"] == "11")
        genCodeValue = GeneticCode::ELEVEN;
    else if (mKeyValuePair["GCODE"] == "4")
        genCodeValue = GeneticCode::FOUR;
    else
        throw logic_error("Genetic code invalid: " + mKeyValuePair["GCODE"]);
    
    GeneticCode geneticCode(genCodeValue);
    AlphabetDNA alph;
    
    // count number of ORFs
    size_t numORF = 0;
    size_t minORFLength = boost::lexical_cast<size_t> (mKeyValuePair["GENE_MIN_LENGTH"]);
    
    // min ORF length can't be zero or larger than sequence size
    if (minORFLength == 0 || minORFLength > seq.size()) {
        cout << 0 << endl;
        return;
    }
    
    // if minORFLength is not divisible by 3, "round" it up
    if (minORFLength % 3 != 0)
        minORFLength = (((size_t) (minORFLength/3)) + 1) * 3;
    
    size_t codonLen = 3;
    
    for (size_t n = 0; n < seq.size(); n++) {
        
        // INV: n in [0, N-4]  =>   n+3 in [3, N-1]
        
        // check for ORF on positive strand
        if (n > codonLen) {
            
            string candStop = seq.toString(n-2,codonLen);       // get stop codon
            
            if (geneticCode.isStop(candStop)) {                 // is stop on positive strand?
                
                // check if there's enough nt for ORF of min length
                if (n >= minORFLength-1) {
                    
                    // search for start longer than minORFLength
                    size_t m = n - 5;
                    
                    while (true) {
                        
                        string cand = seq.toString(m,codonLen);         // get candidate start
                        size_t currLen = n-m+1;
                        // found start?
                        if (currLen >= minORFLength  && geneticCode.isStart(cand)) {
                            numORF++;
                            if (utilOpt.printSeq)
                                cout << m << "\t" << n << "\t" << "+" << "\t" << seq.toString(m, n-m+1) << endl;
                            break;
                        }
                        // found in-frame stop on same strand?
                        else if (geneticCode.isStop(cand))
                            break;
                        
                        // check if no more codons exist
                        if (m < codonLen)
                            break;
                        
                        // otherwise move to previous codon
                        m -= codonLen;
                    }
                }
            }
            
            
        }
        
        
        // check for ORF on negative strand
        if (n <= seq.size() - minORFLength) {
            
            string candStop = alph.reverseComplement(seq.toString(n,codonLen));       // get stop codon on negative strand
            
            // is stop on negative strand
            if (geneticCode.isStop(candStop)) {
                
                // check if there's enought nt for ORF of min length
                if (n < seq.size() - (minORFLength-1)) {
                    
                    // search for start longer than minORFLength
                    size_t m = n + 5;
                    
                    while (true) {
                        string cand = alph.reverseComplement(seq.toString(m-2,codonLen));
                        
                        size_t currLength = m - n + 1;
                        // found start?
                        if (currLength >= minORFLength && geneticCode.isStart(cand)) {
                            numORF++;
                            if (utilOpt.printSeq)
                                cout << n << "\t" << m << "\t" << currLength << "\t" << "-" << "\t" <<  alph.reverseComplement(seq.toString(n, m - n + 1)) << endl;
                            break;
                        }
                        else if (geneticCode.isStop(cand))
                            break;
                        
                        // if no more codons exist
                        if (m >= seq.size() - 3)
                            break;
                        
                        // otherwise move to previous codon
                        m += 3;
                    }
                }
            }
        }
    }
    
    if (!utilOpt.printSeq)
        cout << numORF << endl;
    
}


void ModuleUtilities::runExtractStartContext() {
    
    OptionsUtilities::ExtractStartContext utilOpt = options.extractStartContext;
    
    // read sequence file
    SequenceFile sequenceFile(utilOpt.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile(utilOpt.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    GeneticCode geneticCode(GeneticCode::ELEVEN);           // FIXME: allow genetic code 4
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    
    NumSequence sequence(strSequence, cnc);
    vector<NumSequence> startContexts;
    size_t scLength = utilOpt.rightRelativeToStart - utilOpt.leftRelativeToStart + 1;
    SequenceParser::extractStartContextSequences(sequence, labels, cnc, utilOpt.leftRelativeToStart, scLength, startContexts);
    
    if (startContexts.size() != labels.size())
        throw logic_error("Number of start contexts differs from number of labels - will cause misalignment in fasta defs");
    
    // FIXME: Allow printing to file
    for (size_t n = 0; n < startContexts.size(); n++) {
        if (startContexts[n].size() == 0)
            continue;
        
        if (utilOpt.outputFastaDefs)
            cout << ">" << labels[n]->toString(true) << endl;
        
        string strSC = cnc.convert(startContexts[n].begin(), startContexts[n].end());
        cout << strSC << endl;
    }
}


void ModuleUtilities::runExtractStartContextPerOperonStatus() {
    OptionsUtilities::ExtractStartContextPerOperonStatus utilOpt = options.extractStartContextPerOperonStatus;
    
    // read sequence file
    SequenceFile sequenceFile (utilOpt.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (utilOpt.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    GeneticCode geneticCode(GeneticCode::ELEVEN);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    
    NumSequence sequence (strSequence, cnc);
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, utilOpt.distThreshFGIO, utilOpt.distThreshIG, operonStatuses);
    
    vector<NumSequence> startContexts;
    SequenceParser::extractStartContextSequences(sequence, labels, cnc, -3, 18, startContexts);
    
    for (size_t n = 0; n < startContexts.size(); n++) {
        if (startContexts[n].size() == 0)
            continue;
        
        string label = "AMBIG";
        if (operonStatuses[n] == LabelsParser::FGIO) label = "FGIO";
        if (operonStatuses[n] == LabelsParser::NFGIO) label = "IG";
        
        string strSC = cnc.convert(startContexts[n].begin(), startContexts[n].end());
        cout << label << "\t" << strSC << endl;
    }
    
}

void ModuleUtilities::runExtractStartContextPerMotifStatus() {
    
    OptionsUtilities::ExtractStartContextPerMotifStatus utilOpt = options.extractStartContextPerMotifStatus;
    
    cout << utilOpt.fn_sequence << endl;
    cout << utilOpt.fn_label << endl;
    
    // read sequence file
    SequenceFile sequenceFile (utilOpt.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (utilOpt.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    GeneticCode geneticCode(GeneticCode::ELEVEN);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    
    NumSequence sequence (strSequence, cnc);
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, utilOpt.distThreshFGIO, utilOpt.distThreshIG, operonStatuses);
    
    vector<Label*> labelsFGIO;
    vector<Label*> labelsIG;
    for (size_t n = 0; n < labels.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)    labelsFGIO.push_back(labels[n]);        // fgio
        if (operonStatuses[n] == LabelsParser::NFGIO)   labelsIG.push_back(labels[n]);          // ig
    }
    
    vector<NumSequence> upstreamsFGIO, upstreamsIG;
    SequenceParser::extractUpstreamSequences(sequence, labelsFGIO, cnc, 20, upstreamsFGIO, true, 0);
    SequenceParser::extractUpstreamSequences(sequence, labelsIG, cnc, 20, upstreamsIG, true, 0);
    
    vector<NumSequence> startContextsFGIO, startContextsIG;
    SequenceParser::extractStartContextSequences(sequence, labelsFGIO, cnc, -3, 18, startContextsFGIO);
    SequenceParser::extractStartContextSequences(sequence, labelsIG, cnc, -3, 18, startContextsIG);
    
    Sequence matchToStr (utilOpt.matchTo);
    NumSequence matchSeq (matchToStr, cnc);
    
    vector<NumSequence> startContextsRBS;
    vector<NumSequence> startContextsPromoter;
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (utilOpt.allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    for (size_t n = 0; n < startContextsFGIO.size(); n++) {
        if (startContextsFGIO[n].size() == 0)
            continue;
        
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreamsFGIO[n], positionsOfMatches, substitutions);
        
        // keep track of nonmatches
        if (match.size() < utilOpt.matchThresh)
            startContextsPromoter.push_back(startContextsFGIO[n]);
//        else
//            startContextsRBS.push_back(startContextsFGIO[n]);
    }
    
    startContextsRBS.insert(startContextsRBS.end(), startContextsIG.begin(), startContextsIG.end());
    
    // print promoter-based start-contexts
    for (size_t n = 0; n < startContextsPromoter.size(); n++) {
        cout << "PROM" << "\t" << cnc.convert(startContextsPromoter[n].begin(), startContextsPromoter[n].end()) << endl;
    }
    
    // print RBS-based
    for (size_t n = 0; n < startContextsRBS.size(); n++) {
        cout << "RBS" << "\t" << cnc.convert(startContextsRBS[n].begin(), startContextsRBS[n].end()) << endl;
    }
    
}







void ModuleUtilities::runComputeGC() {
    OptionsUtilities::ComputeGC utilOpt = options.computeGC;
    
    // read sequence file
    SequenceFile sequenceFile (utilOpt.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file (if given)
    vector<Label*> labels;
    
    if (!utilOpt.fn_label.empty()) {
        LabelFile labelFile (utilOpt.fn_label, LabelFile::READ);
        labelFile.read(labels);
        
        vector<double> gcs;
        SequenceAlgorithms::computeGC(strSequence, labels, gcs);
        for (size_t n = 0; n < gcs.size(); n++)
            cout << gcs[n] << endl;
        
    }
    // for entire sequence
    else {
        double gc = SequenceAlgorithms::computeGC(strSequence);
        cout << gc << endl;
    }
}


void ModuleUtilities::runComputeKL() {
    OptionsUtilities::ComputeKL utilOpt = options.computeKL;
    
    // read model file
    ModelFile modFile (utilOpt.fn_mod, ModelFile::READ);
    
    map<string, string> mKeyValuePair;
    modFile.read(mKeyValuePair);
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // break non-coding string into vector of codon-probability pairs
    vector<pair<string, double> > nonCodingProbs;
    istringstream ssm(mKeyValuePair["NON_MAT"]);
    string line;
    while (std::getline(ssm, line)) {
        
        stringstream lineStream (line);
        string codon;
        double prob;
        
        lineStream >> codon >> prob;
        
        nonCodingProbs.push_back(pair<string, double> (codon, prob));
    }
    
    // build non-coding model
    NonCodingMarkov nonCodingMarkov(nonCodingProbs, numAlph, cnc);
    
    
    // Read RBS model
    
    
    istringstream ssmMotifWidth(mKeyValuePair["RBS_WIDTH"]);    // FIXME: any motif type
    std::getline(ssmMotifWidth, line);
    int motifWidth = (int) strtol(line.c_str(), NULL, 10);
    
    vector<vector<pair<string, double> > > motifProbs (motifWidth);
    istringstream ssmMotif(mKeyValuePair[utilOpt.motifLabel]);
    
    line = "";
    while (std::getline(ssmMotif, line)) {
        stringstream lineStream (line);
        string nucleotides;
        vector<double> prob (motifWidth);
        
        lineStream >> nucleotides;
        
        for (int i = 0; i < motifWidth; i++) {
            lineStream >> prob[i];
            motifProbs[i].push_back(pair<string, double> (nucleotides, prob[i]));
        }
        
        
    }
    
    MotifMarkov motifMarkov (motifProbs, (size_t) motifWidth, numAlph, cnc);
    
    
    KLDivergence div (&motifMarkov, &nonCodingMarkov);
    
    
    
    
    
    
    // spacer KL
    istringstream motifDurMax(mKeyValuePair["RBS_MAX_DUR"]);        // FIXME: any motif spacer
    std::getline(motifDurMax, line);
    int maxDur = (int) strtol(line.c_str(), NULL, 10);
    int spacerLen = maxDur+1;
    
    istringstream motifDur(mKeyValuePair["RBS_POS_DISTR"]);
    vector<double> positionProbs (spacerLen, 0);
    line = "";
    while (std::getline(motifDur, line)) {
        stringstream lineStream (line);
        int pos;
        double prob;
        
        lineStream >> pos >> prob;
        
        positionProbs[pos] = prob;
    }
    
    
    UnivariatePDF motifSpacerPDF(positionProbs, false);
    
    
    // compute kl of spacer vs uniform
    double klSpacerUnif = 0;
    for (size_t n = 0; n < spacerLen; n++) {
        double ratio = motifSpacerPDF[n] / (1.0/spacerLen);
        
        if (ratio != 0)
            klSpacerUnif += motifSpacerPDF[n] * log2(ratio);
    }
    
    double motifKL = div.computeKL();
    
    cout << motifKL << "\t" << klSpacerUnif << "\t" << motifKL + klSpacerUnif <<  endl;;
}


void ModuleUtilities::runSeparateFGIOAndIG() {
    OptionsUtilities::SeparateFGIOAndIG utilOpt = options.separateFGIOAndIG;
    
    // read labels file
    LabelFile labelFile (utilOpt.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, utilOpt.distThreshFGIO, utilOpt.distThreshIG, operonStatuses);
    
    // get stats of operons
    size_t numFGIO = 0, numIG = 0, numUNK = 0;
    for (size_t n = 0; n < operonStatuses.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)        numFGIO++;
        else if (operonStatuses[n] == LabelsParser::NFGIO)  numIG++;
        else
            numUNK++;
    }
    
    // get FGIO and IG labels
    vector<Label*> labelsFGIO (numFGIO);
    vector<Label*> labelsIG (numIG);
    size_t currFGIO = 0, currIG = 0;
    
    for (size_t n = 0; n < operonStatuses.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)        labelsFGIO[currFGIO++] = labels[n];
        else if (operonStatuses[n] == LabelsParser::NFGIO)  labelsIG[currIG++] = labels[n];
    }
    
    // write labels to file
    if (!utilOpt.fnout_fgio.empty()) {
        LabelFile fgioFile (utilOpt.fnout_fgio, LabelFile::WRITE);
        fgioFile.write(labelsFGIO);
    }
    
    if (!utilOpt.fnout_ig.empty()) {
        LabelFile igFile (utilOpt.fnout_ig, LabelFile::WRITE);
        igFile.write(labelsIG);
    }
}







void ModuleUtilities::runDNAToAA() {
    OptionsUtilities::DNAToAA utilOpt = options.dnaToAA;
    
    // read sequences
    SequenceFile sequenceFile (utilOpt.fnseqs, SequenceFile::READ);
    
    vector<gmsuite::Sequence> sequences;
    sequenceFile.read(sequences);
    
    AlphabetDNA alph;
    GeneticCode gcode(GeneticCode::ELEVEN);           // FIXME: allow genetic code 4
    
    for (vector<Sequence>::iterator iter = sequences.begin(); iter != sequences.end(); iter++) {
        string seqStr = iter->toString();
        
        if (seqStr.length() == 0)
            continue;
        
        // assert length is multiple of 3
        if (seqStr.length() % 3 != 0) {
            throw invalid_argument("Sequence length isn't multiple of 3");
        }
        
        string AA = "";
        
        // loop over codons
        for (int i = 0; i < seqStr.length()-3; i += 3) {
            AA += gcode.translateCodon(seqStr.substr(i, 3));
        }
        
        if (utilOpt.outputFastaDefs) {
            cout << ">" << iter->getMetaData() << endl;
        }
        
        cout << AA << endl;
    }
    
    
}


void ModuleUtilities::runABFilter() {
    OptionsUtilities::ABFilter utilOpt = options.abFilter;
    
    
    // Read two label files
    vector<Label*> labelsA, labelsB, labelsOut;        // where labels will be stored and placed at output
    
    LabelFile flabA (utilOpt.fnA, LabelFile::READ);
    LabelFile flabB (utilOpt.fnB, LabelFile::READ);
    
    flabA.read(labelsA);
    flabB.read(labelsB);
    
    // make sure they're sorted
    sort(labelsA.begin(), labelsA.end(), Label::compareByLeftAndStrand);
    sort(labelsB.begin(), labelsB.end(), Label::compareByLeftAndStrand);
    
    // create a hash map for labelsB, where keys are of the form "3prime;strand", and the value is the
    // index of that gene in labelsB
    stringstream ssm;
    
    map<string, size_t> idToIdxB;
    for (size_t idx = 0; idx < labelsB.size(); idx++) {
        
        Label* label = labelsB[idx];
        
        ssm.str(std::string());        // clear the stream
        
        // determine 3-prime end
        size_t loc3Prime = label->get3Prime();
        ssm << loc3Prime << ";" << label->strandToStr();
        
        // add key
        idToIdxB[ssm.str()] = idx;
    }
    
    
    // For each label lA in A, find it's 3-prime end counterpart in array B (using hash)
    // Then, find if there is a label in B close by to the start of lA (i.e. close to
    // 5prime of lA).
    for (size_t idx = 0; idx < labelsA.size(); idx++) {
        Label* labA = labelsA[idx];
        
        // construct id from A
        ssm.str(std::string());
        
        size_t loc3Prime = labA->get3Prime();
        ssm << loc3Prime << ";" << labA->strandToStr();
    
        size_t idxInB = idToIdxB[ssm.str()];
        
        // get distance to upstream gene from 5prime of labA
        int distanceToUpstreamLabel = numeric_limits<int>::infinity();          // negative values indicate overlap
        
        bool keepLabel = false;
        
        if (labA->strand == Label::POS) {
            // if no label upstream of current, add to output list
            if (idxInB == 0)
                keepLabel = true;
            // otherwise, check distance to upstream
            else {
                size_t prevIdxInB = idxInB-1;
                Label* prevLabB = labelsB[prevIdxInB];
                
                distanceToUpstreamLabel = (int)labA->left - (int)prevLabB->right + 1;
            }
        }
        // negative strand
        else {
            // if no label upstream of current, add to output list
            if (idxInB == labelsB.size())
                keepLabel = true;
            // otherwise, check distance to upstream
            else {
                size_t prevIdxInB = idxInB+1;
                Label* prevLabB = labelsB[prevIdxInB];
                
                distanceToUpstreamLabel = (int) prevLabB->left - (int)labA->right + 1;
            }
        }
        
        // if distance less than threshold, place in output vector. Otherwise, discard
        if (distanceToUpstreamLabel > utilOpt.threshDistToUpstream)
            keepLabel = true;
        
        // reverse label if reverse option is set
        if (utilOpt.selectNoSatisfy)
            keepLabel = !keepLabel;
        
        if (keepLabel)
            labelsOut.push_back(labA);
    }
    
    
    // print output labels
    LabelFile fout (utilOpt.fnout, LabelFile::WRITE);
    fout.write(labelsOut);
}

void ModuleUtilities::runChangeOrderNonCoding() {
    OptionsUtilities::ChangeOrderNonCoding utilOpt = options.changeOrderNonCoding;
    
    // read input model file
    ModelFile fin(utilOpt.fnin, ModelFile::READ);
    
    map<string,string> modelIn;
    fin.read(modelIn);
    
    
    // set keys for non-coding model
    vector<string> keys;
    keys.push_back("NON_MAT");
    keys.push_back("NON_ORDER");
    
    // get key-value pair for non-coding model
    map<string, string> mKeyValuePair;
    fin.read(mKeyValuePair, keys);
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // break non-coding string into vector of codon-probability pairs
    vector<pair<string, double> > nonCodingProbs;
    istringstream ssm(mKeyValuePair["NON_MAT"]);
    string line;
    while (std::getline(ssm, line)) {
        
        stringstream lineStream (line);
        string codon;
        double prob;
        
        lineStream >> codon >> prob;
        nonCodingProbs.push_back(pair<string, double> (codon, prob));
    }
    
    // build non-coding model
    NonCodingMarkov nonCodingMarkov(nonCodingProbs, numAlph, cnc);
    
    nonCodingMarkov.changeOrder(utilOpt.newOrder);
    
    modelIn["NON_MAT"] = nonCodingMarkov.toString();
    
    stringstream ssm_order;
    ssm_order << utilOpt.newOrder;
    modelIn["NON_ORDER"] = (ssm_order.str());
    
    // write to output
    ModelFile fout (utilOpt.fnout, ModelFile::WRITE);
    
    vector<pair<string, string> > modelOut;
    
    for (map<string, string>::iterator iter = modelIn.begin(); iter != modelIn.end(); iter++) {
        pair<string, string> p;
        p.first = iter->first;
        p.second = iter->second;
        modelOut.push_back(p);
    }
    
    fout.write(modelOut);
}





pair<NumSequence::size_type, double> findBestMotifLocation2(NumSequence::const_iterator begin, NumSequence::const_iterator end, boost::shared_ptr<const MotifMarkov> motifMarkov, boost::shared_ptr<const UnivariatePDF> motifSpacer) {
    
    NumSequence::size_type bestLoc = 0;
    double bestScore = -numeric_limits<double>::infinity();
    NumSequence::size_type motifWidth = motifMarkov->getLength();
    NumSequence::size_type fragLen = std::distance(begin, end);
    
    // if motif can't fit in fragment, return -infinity score
    if (fragLen < motifWidth)
        return pair<NumSequence::size_type, double> (NumSequence::npos, bestScore);
    
    // loop over all valid positions from motif
    NumSequence::size_type currLoc = 0;
    for (NumSequence::const_iterator currElement = begin; currElement < end-motifWidth; currElement++, currLoc++) {
        
        double score = motifMarkov->evaluate(currElement, currElement+motifWidth, true);
        
        NumSequence::size_type posFromRight = fragLen - motifWidth - currLoc;
        score += log2(motifSpacer->operator[](posFromRight));
        
        if (score > bestScore) {
            bestScore = score;
            bestLoc = currLoc;
        }
        
//        AlphabetDNA alph;
//        CharNumConverter cnc(&alph);
//        NumAlphabetDNA numAlph(alph, cnc);
//
//        cout << cnc.convert(currElement, currElement+motifWidth) << "\t" << score << endl;
    }
    
    return pair<NumSequence::size_type, double> (bestLoc, bestScore);
}



void ModuleUtilities::runExtractSpacerNTModel() {
    OptionsUtilities::ExtractSpacerNTModel utilOpt = options.extractSpacerNTModel;
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);

    
    // read sequence from file
    SequenceFile seqFile(utilOpt.fnsequences, SequenceFile::READ);
    Sequence strSeq = seqFile.read();
    NumSequence numSeq (strSeq, cnc);
    
    // read labels from file
    LabelFile labelFile (utilOpt.fnlabels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // read model file
    ModelFile modFile(utilOpt.fnmod, ModelFile::READ);
    map<string, string> keyValPair;
    modFile.read(keyValPair);
    
    // check if models enabled
    bool rbsEnabled = false, promoterEnabled = false;
    
    boost::shared_ptr<MotifMarkov> rbsMarkov, promoterMarkov;
    boost::shared_ptr<UnivariatePDF> rbsSpacerPDF, promoterSpacerPDF;
    
    // RBS enabled?
    if (keyValPair.count("RBS") > 0 && keyValPair["RBS"] == "1")
        rbsEnabled = true;
    
    if (keyValPair.count("PROMOTER") > 0 && keyValPair["PROMOTER"] == "1")
        promoterEnabled = true;
    
    // if neither promoter nor RBS enabled, then nothing to extract
    if (!rbsEnabled && !promoterEnabled)
        return;
    
    istringstream motifEnabled(keyValPair["RBS"]);
    
    string line;
    
    // RBS model
    if (rbsEnabled) {
        
        istringstream ssmMotifWidth(keyValPair["RBS_WIDTH"]);    // FIXME: any motif type
        std::getline(ssmMotifWidth, line);
        int motifWidth = (int) strtol(line.c_str(), NULL, 10);
        
        
        vector<vector<pair<string, double> > > rbsProbs (motifWidth);
        istringstream ssmRBS(keyValPair["RBS_MAT"]);
        line = "";
        while (std::getline(ssmRBS, line)) {
            stringstream lineStream (line);
            string nt;
            vector<double> prob (motifWidth);
            
            lineStream >> nt;
            
            for (int i = 0; i < motifWidth; i++) {
                lineStream >> prob[i];
                rbsProbs[i].push_back(pair<string, double> (nt, prob[i]));
            }
        }
        
        
        rbsMarkov = boost::shared_ptr<MotifMarkov> (new MotifMarkov(rbsProbs, (size_t) motifWidth, numAlph, cnc));
        
        // rbs spacer
        istringstream rbsDurMax(keyValPair["RBS_MAX_DUR"]);        // FIXME: any motif spacer
        std::getline(rbsDurMax, line);
        int maxRBSDur = (int) strtol(line.c_str(), NULL, 10);
        int rbsSpacerLen = maxRBSDur+1;
        
        istringstream rbsDur(keyValPair["RBS_POS_DISTR"]);
        vector<double> rbsPositionProbs (rbsSpacerLen, 0);
        line = "";
        while (std::getline(rbsDur, line)) {
            stringstream lineStream (line);
            int pos;
            double prob;
            
            lineStream >> pos >> prob;
            
            rbsPositionProbs[pos] = prob;
        }
        rbsSpacerPDF = boost::shared_ptr<UnivariatePDF> (new UnivariatePDF(rbsPositionProbs, false));
    }
    
    if (promoterEnabled) {
        
        if (keyValPair.count("PROMOTER_MAX_DUR") > 0) {
            istringstream promDurMax(keyValPair["PROMOTER_MAX_DUR"]);
            std::getline(promDurMax, line);
            int maxPromDur = (int) strtol(line.c_str(), NULL, 10);
            int promSpacerLen = maxPromDur + 1;
            
            istringstream promDur (keyValPair["PROMOTER_POS_DISTR"]);
            vector<double> promPositionsProbs (promSpacerLen, 0);
            line = "";
            while(std::getline(promDur, line)) {
                stringstream lineStream (line);
                int pos;
                double prob;
                
                lineStream >> pos >> prob;
                promPositionsProbs[pos] = prob;
            }
            promoterSpacerPDF = boost::shared_ptr<UnivariatePDF>  (new UnivariatePDF(promPositionsProbs, false));
        }
        
        // Promoter Model
        if (keyValPair.count("PROMOTER_WIDTH") < 0) {
            istringstream ssmPromWidth(keyValPair["PROMOTER_WIDTH"]);
            std::getline(ssmPromWidth, line);
            int promWidth = (int) strtol(line.c_str(), NULL, 10);
            
            vector<vector<pair<string, double> > > promProbs (promWidth);
            istringstream ssmProm(keyValPair["PROMOTER_MAT"]);
            line = "";
            while (std::getline(ssmProm, line)) {
                stringstream lineStream (line);
                string nt;
                vector<double> prob (promWidth);
                
                lineStream >> nt;
                
                for (int i = 0; i < promWidth; i++) {
                    lineStream >> prob[i];
                    promProbs[i].push_back(pair<string, double> (nt, prob[i]));
                }
            }
            
            promoterMarkov = boost::shared_ptr<MotifMarkov>  (new MotifMarkov(promProbs, (size_t) promWidth, numAlph, cnc));
        }
    }
    
    
    UniformCounts spacerCounts(2, numAlph);
    // for each label
    for (vector<Label*>::const_iterator iter = labels.begin(); iter!= labels.end(); iter++) {
        // get search window around labelled start
        size_t labelLeft = (*iter)->left;
        size_t labelRight = (*iter)->right;
        Label::strand_t strand = (*iter)->strand;
            
        if (numSeq.size() > 0) {
            
            // extract upstream region
            NumSequence upstrSeq;
            
            if (strand == Label::POS) {
                if (labelLeft < utilOpt.upstreamLength)
                    continue;
                
                size_t upstrLeft = labelLeft - utilOpt.upstreamLength;
                upstrSeq = numSeq.subseq(upstrLeft, utilOpt.upstreamLength);
            }
            else {
                if (labelRight + utilOpt.upstreamLength >= numSeq.size())
                    continue;
                
                upstrSeq = numSeq.subseq(labelRight+1, utilOpt.upstreamLength);
                upstrSeq.reverseComplement(cnc);
            }
            
            pair<NumSequence::size_type, double> bestRBS (NumSequence::npos, -numeric_limits<double>::infinity());
            pair<NumSequence::size_type, double> bestProm (NumSequence::npos, -numeric_limits<double>::infinity());

//            if (utilOpt.debug)
//                cout << "#" << cnc.convert(upstrSeq.begin(), upstrSeq.end()) << endl;
            
            // compute best locations of motifs
            if (rbsMarkov)
                bestRBS = findBestMotifLocation2(upstrSeq.begin(), upstrSeq.end(), rbsMarkov, rbsSpacerPDF);
            if (promoterMarkov)
                bestProm = findBestMotifLocation2(upstrSeq.begin(), upstrSeq.end(), promoterMarkov, promoterSpacerPDF);
            
            // if RBS only needed and promoter has larger value, skip
            if (utilOpt.RBSOnly && bestProm.second > bestRBS.second) {
                continue;
            }
            
            pair<NumSequence::size_type, double> bestMotif = bestRBS;
            size_t bestMotifWidth = 0;
            
            if (rbsMarkov)
                bestMotifWidth = rbsMarkov->getLength();
            
            if (bestProm.second > bestRBS.second) {
                bestMotif = bestProm;
                bestMotifWidth = promoterMarkov->getLength();
            }
            
            if (utilOpt.debug) {
                cout  << "#" << cnc.convert(upstrSeq.begin() + bestMotif.first, upstrSeq.begin()+bestMotif.first + bestMotifWidth) << endl;
            }
            
            NumSequence::const_iterator spacerBegin = upstrSeq.begin() + bestMotif.first + bestMotifWidth;
            NumSequence::const_iterator spacerEnd = upstrSeq.end();
            
            if (utilOpt.debug) {
                cout << "#SP:\t" << cnc.convert(spacerBegin, spacerEnd) << endl;
            }

            // extract region between RBS and end of upstream region (near start)
            spacerCounts.count(spacerBegin, spacerEnd);
        }
    }
    
    UniformMarkov spacerMarkov(2, numAlph);
    spacerMarkov.construct(&spacerCounts);
    
    cout << spacerMarkov.toString() << endl;
    
    
}

void reverseComplementInPlace(string &nt) {
    map<char, char> cmap;
    cmap.insert(std::pair<char,char>('A', 'T'));
    cmap.insert(std::pair<char,char>('T', 'A'));
    cmap.insert(std::pair<char,char>('G', 'C'));
    cmap.insert(std::pair<char,char>('C', 'G'));
    
    size_t fromLeft = 0;
    size_t fromRight = nt.size()-1;
    
    while (fromLeft <= fromRight) {
        
        char cFromLeft =  nt[fromLeft];
        char cFromRight = nt[fromRight];
        
        // complement the characters
        cFromLeft  =    cmap.count(cFromLeft)  > 0 ? cmap.at(cFromLeft)  : cFromLeft;
        cFromRight =    cmap.count(cFromRight) > 0 ? cmap.at(cFromRight) : cFromRight;
        
        // swap them
        nt[fromLeft]  = cFromRight;
        nt[fromRight] = cFromLeft;
        
        fromLeft++;
        if (fromRight == 0)
            break;
        fromRight--;
    }
}

void ModuleUtilities::runExtractLORF() {
    OptionsUtilities::ExtractLORF utilOpt = options.extractLORF;

    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    GeneticCode gcode(utilOpt.gcode);
    NumGeneticCode numGcode(gcode, cnc);

    // get sequences
    SequenceFile fseq (utilOpt.fnsequences, SequenceFile::READ);
    vector<Sequence> sequences;     fseq.read(sequences);               // read sequences

    // get labels
    LabelFile flabel (utilOpt.fnlabels, LabelFile::READ);
    vector<Label*> labels;  flabel.read(labels);

    // for each label, find longest orf
    for (size_t n = 0; n < labels.size(); n++) {
        Label* currLabel = labels[n];
        
        
        size_t currSeqIdx = 0;
        
        // find sequence with seqname same as current label
        while (currSeqIdx < sequences.size()) {
            if (currLabel->meta == sequences[currSeqIdx].getMetaData()) {
                break;
            }
            currSeqIdx++;
        }
        
        if (currSeqIdx == sequences.size()) {
            cout << "Warning: could not find sequence name " << currLabel->meta << endl;
            continue;
        }
        
        string frag = "";
        string fastaHeader = "";
        
        // find longest ORF
        if (currLabel->strand == Label::POS) {
            
            size_t currPos = currLabel->left;
            size_t lorfLoc = currPos;
            while (currPos >= 3) {
                currPos -= 3;        // go back one codon
                
                string codon = sequences[currSeqIdx].toString(currPos, 3);
                
                if (gcode.isStart(codon))
                    lorfLoc = currPos;
                if (gcode.isStop(codon))
                    break;
            }
            
            size_t newLength = currLabel->right - lorfLoc + 1;
            
            // extract lorf
            frag = sequences[currSeqIdx].toString(lorfLoc, newLength);
            
            // create new fasta def
            stringstream ssm;
            ssm << currLabel->meta << ";";
            ssm << lorfLoc+1            << ";" << currLabel->right+1 << ";" << currLabel->strandToStr() << ";";
            ssm << currLabel->left+1    << ";" << currLabel->right+1 << ";" << currLabel->strandToStr();
            fastaHeader = ssm.str();
            
            
        }
        // on negative strand
        else {
            // extract
            
            size_t currPos = currLabel->right;
            size_t lorfLoc = currPos;
            
            while (currPos < sequences[currSeqIdx].size()-3) {
                currPos += 3;            // go "back" one codon
                
                string codon = sequences[currSeqIdx].toString(currPos-2, 3);
                reverseComplementInPlace(codon);
                
                if (gcode.isStart(codon))
                    lorfLoc = currPos;
                if (gcode.isStop(codon))
                    break;
            }
            
            size_t newLength = lorfLoc - currLabel->left + 1;
            
            frag = sequences[currSeqIdx].toString(currLabel->left, newLength);
            reverseComplementInPlace(frag);
            
            // create new fasta def
            stringstream ssm;
            ssm << currLabel->meta << ";";
            ssm << currLabel->left+1    << ";" << lorfLoc+1             << ";" << currLabel->strandToStr() << ";";
            ssm << currLabel->left+1    << ";" << currLabel->right+1    << ";" << currLabel->strandToStr();
            fastaHeader = ssm.str();
        }
        
        cout << fastaHeader << endl << frag << endl;
    }
}

