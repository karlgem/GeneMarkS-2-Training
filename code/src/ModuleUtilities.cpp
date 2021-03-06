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
#include "LabelsParser.hpp"

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












