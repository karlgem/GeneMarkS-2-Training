//
//  GMS2Trainer.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "GMS2Trainer.hpp"
#include <iostream>
#include "MotifFinder.hpp"
#include "UnivariatePDF.hpp"
#include "UniformCounts.hpp"
#include "PeriodicCounts.hpp"
#include "CodingMarkov.hpp"
#include "NonUniformCounts.hpp"
#include "LabelsParser.hpp"
#include "NonCodingCounts.hpp"
#include "SequenceParser.hpp"
#include "CodingCounts.hpp"
#include <boost/lexical_cast.hpp>
#include "OptionsGMS2Training.hpp"
#include "SequenceAlgorithms.hpp"

using namespace std;
using namespace gmsuite;


typedef enum {FIRST_OP, NOT_FIRST_OP, IGNORE} GeneStat;         // gene status

void filterNs (const vector<NumSequence> &original, vector<NumSequence> &filtered) {
    
    
    if (original.size() == 0)
        return;
    
    AlphabetDNA alphabet;
    CharNumConverter cnc(&alphabet);
    NumAlphabetDNA numAlph(alphabet, cnc);
    
    if (alphabet.sizeInvalid() > 0) {
        
        // remove all sequences containing N's
        for (size_t n = 0; n < original.size(); n++) {
            
            // if sequence does not contain N
            if (!original[n].containsInvalid(numAlph)) {
                filtered.push_back(original[n]);
            }
        }
    }
    else {
        filtered = original;
    }
}



// default constructor
GMS2Trainer::GMS2Trainer() {
    // DEPRECATED
    // public variables for models
    noncoding = NULL;
    coding = NULL;
    startContext = NULL;
    
    // start models
    rbs = NULL;
    promoter = NULL;
    upstreamSignature = NULL;
    rbsSpacer = NULL;
    promoterSpacer = NULL;
    
}

GMS2Trainer::GMS2Trainer(unsigned pcounts,
                         unsigned codingOrder,
                         unsigned noncodingOrder,
                         unsigned startContextOrder,
                         NumSequence::size_type upstreamLength,
                         NumSequence::size_type startContextLength,
                         genome_class_t genomeClass,
                         const OptionsMFinder &optionsMFinder,
                         const NumAlphabetDNA &alph,
                         const NumSequence::size_type MIN_GENE_LEN,
                         const NumGeneticCode &numGenCode,
                         int scMargin,
                         bool trainOnNative,
                         bool runMotifSearch,
                         NumSequence::size_type upstrFGIO,
                         unsigned widthArchaeaPromoter,
                         string matchTo,
                         bool allowAGSubstitution,
                         unsigned matchThresh,
                         NumSequence::size_type upstreamSignatureLength,
                         unsigned upstreamSignatureOrder,
                         bool trainNonCodingOnFullGenome,
                         unsigned FGIO_DIST_THRESH) {
    
    this->pcounts = pcounts;
    this->codingOrder = codingOrder;
    this->noncodingOrder = noncodingOrder;
    this->startContextOrder = startContextOrder;
    this->upstreamLength = upstreamLength;
    this->startContextLength = startContextLength;
    this->genomeClass = genomeClass;
    this->optionsMFinder = &optionsMFinder;
    this->alphabet = &alph;
    this->MIN_GENE_LEN = MIN_GENE_LEN;
    this->numGeneticCode = &numGenCode;
    this->scMargin = scMargin;
    this->trainOnNative = trainOnNative;
    this->runMotifSearch = runMotifSearch;
    this->UPSTR_LEN_FGIO = upstrFGIO;
    this->widthArchaeaPromoter = widthArchaeaPromoter;
    this->upstreamSignatureLength = upstreamSignatureLength;
    this->upstreamSignatureOrder = upstreamSignatureOrder;
    this->trainNonCodingOnFullGenome = trainNonCodingOnFullGenome;
    
    this->matchTo = matchTo;
    this->allowAGSubstitution = allowAGSubstitution;
    this->matchThresh = matchThresh;
    
    this->FGIO_DIST_THRESH = FGIO_DIST_THRESH;
    this->NFGIO_DIST_THRES = 22;
    
    if (genomeClass == ProkGeneStartModel::C2 || genomeClass == ProkGeneStartModel::C3 || genomeClass == ProkGeneStartModel::C5) {
        this->UPSTR_LEN_IG = upstreamLength;
    }
    
    // public variables for models
    noncoding = NULL;
    coding = NULL;
    startContext = NULL;
    startContextRBS = NULL;
    startContextPromoter = NULL;
    
    // start models
    rbs = NULL;
    promoter = NULL;
    upstreamSignature = NULL;
    rbsSpacer = NULL;
    promoterSpacer = NULL;

    this->numLeaderless = 0;
    this->numFGIO = 0;
    this->genomeType = "no-motif";
}


void GMS2Trainer::estimateParametersStartStopCodons(const NumSequence &sequence, const vector<Label*> &labels, const vector<bool> &use) {
    
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
    // get starts and stops from genetic code definition
    vector<CharNumConverter::seq_t> starts = this->numGeneticCode->getStarts();
    vector<CharNumConverter::seq_t> stops = this->numGeneticCode->getStops();
    
    // fill in maps
    for (size_t n = 0; n < starts.size(); n++)
        this->startProbs[starts[n]] = 0;
    for (size_t n = 0; n < stops.size(); n++)
        this->stopProbs[stops[n]] = 0;              // stop probabilities to 1
    
    size_t totalStarts = 0;
    size_t totalStops = 0;
    
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        CharNumConverter::seq_t start;
        CharNumConverter::seq_t stop;
        
        // get start from positive strand
        if ((*iter)->strand == Label::POS) {
            start = CharNumConverter::seq_t(sequence.begin() + (*iter)->left, sequence.begin() + (*iter)->left + 3);
            stop = CharNumConverter::seq_t(sequence.begin() + (*iter)->right-2, sequence.begin() + (*iter)->right + 1);
        }
        // get start from negative strand and reverse complement
        else {
            start = CharNumConverter::seq_t (sequence.begin() + (*iter)->right-2, sequence.begin()+(*iter)->right+1);
            NumSequence s(start); s.reverseComplement(*this->alphabet->getCNC());   // reverse complement
            start = CharNumConverter::seq_t (s.begin(), s.end());
            
            stop = CharNumConverter::seq_t (sequence.begin() + (*iter)->left, sequence.begin() + (*iter)->left + 3);
            NumSequence ss(stop); ss.reverseComplement(*this->alphabet->getCNC());  // reverse complement
            stop = CharNumConverter::seq_t (ss.begin(), ss.end());
        }
        
        // if it is a start, add it to counts
        if (this->numGeneticCode->isStart(start)) {
            this->startProbs[start] += 1;
            totalStarts += 1;
        }
        
        if (this->numGeneticCode->isStop(stop)) {
            this->stopProbs[stop] += 1;
            totalStops += 1;
        }
    }
    
    // add pseudocounts
    for (map<CharNumConverter::seq_t, double>::iterator  iter = this->startProbs.begin(); iter != this->startProbs.end(); iter++) {
        iter->second += pcounts;
        totalStarts += pcounts;
    }
    
    for (map<CharNumConverter::seq_t, double>::iterator iter = this->stopProbs.begin(); iter != this->stopProbs.end(); iter++) {
        iter->second += pcounts;
        totalStops += pcounts;
    }
    
    // normalize start counts
    if (totalStarts > 0)
        for (map<CharNumConverter::seq_t, double>::iterator  iter = this->startProbs.begin(); iter != this->startProbs.end(); iter++)
            iter->second /= totalStarts;
    
    if (totalStops > 0)
        for (map<CharNumConverter::seq_t, double>::iterator iter = this->stopProbs.begin(); iter != this->stopProbs.end(); iter++)
            iter->second /= totalStops;
    
}

void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels, NumSequence::size_type scSize, const vector<bool> &use) {
    
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
//    PeriodicCounts counts (codingOrder, 3, *this->alphabet);
    CodingCounts counts (codingOrder, 3, *this->alphabet, *this->numGeneticCode);
    
    // get counts for 3 period markov model given order
    size_t n = 0;
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        if (!useAll && !use[n++])
            continue;       // skip unwanted genes
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        size_t left = (*iter)->left+3;                // get left position of fragment
        size_t right = (*iter)->right-3;              // get right position of fragment
        
        if (left > right)
            continue;
        
        size_t length = right - left + 1;           // compute fragment length
        
        bool reverseComplement = (*iter)->strand == Label::NEG;
        
        counts.count(sequence.begin()+left, sequence.begin() + left +length, reverseComplement);
        
        // FIXME:
        // if start context > 0, remove segement near start
        if (scSize > 0 && scMargin < -3) {
            size_t scSizeInCoding = abs(scMargin) - 3;          // length of start context that overlaps with CDS
            if (reverseComplement)
                counts.decount(sequence.begin()+right-scSizeInCoding+1, sequence.begin()+right+1, reverseComplement);
            else
                counts.decount(sequence.begin()+left, sequence.begin()+left+scSizeInCoding);
        }
        
    }
    
    // convert counts to probabilities
//    coding = new PeriodicMarkov(codingOrder, 3, *this->alphabet);
    coding = new CodingMarkov(codingOrder, 3, *this->alphabet, *this->numGeneticCode);
    coding->construct(&counts, pcounts);
    
}

// this function assumes labels are sorted by "left" in increasing order
void GMS2Trainer::estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
//    UniformCounts counts(noncodingOrder, *this->alphabet);
    NonCodingCounts counts(noncodingOrder, *this->alphabet);
    
    // train on full genome (i.e. consider it as a single non-coding sequence)
    if (this->trainNonCodingOnFullGenome) {
        counts.count(sequence.begin(), sequence.end());
    }
    // train on labels
    else {
        size_t leftNoncoding = 0;       // left position of current noncoding region
        
        // get counts for 3 period markov model given order
        size_t n = 0;
        for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
            if (!useAll && !use[n++])
                continue;       // skip unwanted genes
            
            if (*iter == NULL)
                throw invalid_argument("Label can't be NULL");
            
            size_t left = (*iter)->left;                // get left position of fragment
            size_t right = (*iter)->right;              // get right position of fragment
            
            if (leftNoncoding < left) {
                counts.count(sequence.begin() + leftNoncoding, sequence.begin() + left);
            }
            
            // update left position of (possible) non-coding region after current gene
            leftNoncoding = right+1;
            
        }
        
        // add last non-coding sequence
        counts.count(sequence.begin() + leftNoncoding, sequence.begin() + sequence.size());
    }
    
    // convert counts to probabilities
    noncoding = new UniformMarkov(noncodingOrder, *this->alphabet);
    noncoding->construct(&counts, pcounts);
}


void GMS2Trainer::estimateParametersStartContext(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    
    NonUniformCounts counts(startContextOrder, startContextLength, *this->alphabet);
    
    // get counts for start context model
    size_t n = 0;
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        if (!useAll && !use[n++])
            continue;       // skip unwanted genes
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        // skip sequence if it doesn't have enough nucleotides for start context
        size_t ntBeforeStart = (startContextLength + scMargin);
        if ((*iter)->strand == Label::POS && (*iter)->left < ntBeforeStart) {
            continue;
        }
        else if ((*iter)->strand == Label::NEG && (*iter)->right > sequence.size() - ntBeforeStart) {
            continue;
        }
        if ((*iter)->right - (*iter)->left + 1 - 6 < startContextLength)        // get length of sequence, and -6 to account for start/stop codons
            continue;
        
        size_t left;
        
        if ((*iter)->strand == Label::POS)
//            left = (*iter)->left + 3;
//            left = (*iter)->left-3;    // left = 3;   (3 + -(18-15)) = 0
            left = (*iter)->left - ((int)startContextLength + scMargin);
        else
//            left = (*iter)->right - 2 - startContextLength;
//            left = (*iter)->right+3 - startContextLength + 1;      // right = 20:    20 - 15
            left = (*iter)->right + scMargin + 1;
        
        bool reverseComplement = (*iter)->strand == Label::NEG;
        counts.count(sequence.begin() + left, sequence.begin() + left + startContextLength, reverseComplement);
        
//        NumSequence s = sequence.subseq(left, startContextLength);
//        if (reverseComplement)
//            s.reverseComplement(cnc);
//        
//        cout << cnc.convert(s.begin(), s.end()) << endl;
        
    }
    
    // convert counts to probabilities
    if (!runMotifSearch) {
        startContext = new NonUniformMarkov(startContextOrder, startContextLength, *this->alphabet);
        startContext->construct(&counts, pcounts);
    }
    else {
        startContextRBS = new NonUniformMarkov(startContextOrder, startContextLength, *this->alphabet);
        startContextRBS->construct(&counts, pcounts);
    }
}


// Output: should set "rbsSpacer" and "rbs" variables
void GMS2Trainer::estimateParametersMotifModel(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    // check if all labels should be used
    if (use.size() > 0) {
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    MotifFinder::Builder b;
    b.setAlign(optionsMFinder->align);
    b.setWidth(optionsMFinder->width);
    b.setMaxIter(optionsMFinder->maxIter);
    b.setPcounts(optionsMFinder->pcounts);
    b.setNumTries(optionsMFinder->tries);
    b.setBackOrder(optionsMFinder->bkgdOrder);
    b.setMaxEMIter(optionsMFinder->maxEMIter);
    b.setMotifOrder(optionsMFinder->motifOrder);
    b.setShiftEvery(optionsMFinder->shiftEvery);
    
    
    MotifFinder mfinder = b.build();
    
    // if genome is class 1, search for RBS
    if (genomeClass == ProkGeneStartModel::C1 || genomeClass == ProkGeneStartModel::C5) {
        this->genomeType = "pure-rbs";
        if (genomeClass == ProkGeneStartModel::C5)
            this->genomeType = "class-c";
        
        // START: REMOVE THIS
        // copy only usable labels
        size_t numUse = labels.size();
        if (use.size() > 0) {
            numUse = 0;
            for (size_t n = 0; n < labels.size(); n++)
                if (use[n])
                    numUse++;
        }
        
        vector<Label*> useLabels (numUse);
        size_t pos = 0;
        for (size_t n = 0; n < labels.size(); n++) {
            if (use[n])
                useLabels[pos++] = labels[n];
        }
        
        // split labels into sets based on operon status
        vector<LabelsParser::operon_status_t> operonStatuses;
        LabelsParser::partitionBasedOnOperonStatus(useLabels, FGIO_DIST_THRESH, NFGIO_DIST_THRES, operonStatuses);
        
        size_t numFGIO = 0, numIG = 0, numUNK = 0;
        for (size_t n = 0; n < operonStatuses.size(); n++) {
            if (operonStatuses[n] == LabelsParser::FGIO)        numFGIO++;
            else if (operonStatuses[n] == LabelsParser::NFGIO)  numIG++;
            else
                numUNK++;
        }
        this->numFGIO = numFGIO;
        // END: REMOVE THIS
        
        // extract upstream of each label
        vector<NumSequence> upstreamsRaw;
        SequenceParser::extractUpstreamSequences(sequence, labels, *alphabet->getCNC(), upstreamLength, upstreamsRaw, false, MIN_GENE_LEN, use);
        
        vector<NumSequence> upstreams;
        for (size_t n = 0; n < upstreamsRaw.size(); n++) {
            if (!upstreamsRaw[n].containsInvalid(numAlph))
                upstreams.push_back(upstreamsRaw[n]);
        }
        
        vector<NumSequence::size_type> positions;
        mfinder.findMotifs(upstreams, positions);
        
        // build RBS model
        NonUniformCounts rbsCounts(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet);
        for (size_t n = 0; n < upstreams.size(); n++) {
            rbsCounts.count(upstreams[n].begin()+positions[n], upstreams[n].begin()+positions[n]+optionsMFinder->width);
        }
        
        rbs = new NonUniformMarkov(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet);
        rbs->construct(&rbsCounts, optionsMFinder->pcounts);
        
        // build spacer distribution
        // build histogram from positions
        vector<double> positionCounts (upstreamLength - optionsMFinder->width+1, 0);
        for (size_t n = 0; n < positions.size(); n++) {
            // FIXME account for LEFT alignment
            // below is only for right
            positionCounts[upstreamLength - optionsMFinder->width - positions[n]]++;        // increment position
        }
        
        rbsSpacer = new UnivariatePDF(positionCounts, false, pcounts);
        
    }
    // if genome is class 2: promoter and RBS in archaea
    else if (genomeClass == ProkGeneStartModel::C2) {
        this->genomeType = "archaea-promoter";
        estimateParametersMotifModel_Promoter(sequence, labels, use);
    }
    // if genome is class 3: promoter and RBS in bacteria
    else if (genomeClass == ProkGeneStartModel::C3) {
        this->genomeType = "bacteria-promoter";
        estimateParametersMotifModel_Tuberculosis(sequence, labels, use);
    }
    // if genome is class 4: RBS and upstream signature
    else {
        this->genomeType = "upstream-signature";
        estimateParametersMotifModel_Synechocystis(sequence, labels, use);
    }
    
}


void runMotifFinder(const vector<NumSequence> &sequencesRaw, const OptionsMFinder &optionsMFinder, const NumAlphabetDNA  &numAlph, size_t upstreamLength, NonUniformMarkov* &motifMarkov, UnivariatePDF* &motifSpacer) {
    
//    AlphabetDNA alph;
//    CharNumConverter cnc(&alph);
    
    vector<NumSequence> upstreams;
    for (size_t n = 0; n < sequencesRaw.size(); n++) {
        if (!sequencesRaw[n].containsInvalid(numAlph))
            upstreams.push_back(sequencesRaw[n]);
    }
    
    MotifFinder::Builder b;
    MotifFinder mfinder = b.build(optionsMFinder);
    
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(upstreams, positions);
    
    // build RBS model
    NonUniformCounts motifCounts(optionsMFinder.motifOrder, optionsMFinder.width, numAlph);
    for (size_t n = 0; n < upstreams.size(); n++) {
        motifCounts.count(upstreams[n].begin()+positions[n], upstreams[n].begin()+positions[n]+optionsMFinder.width);
    }
    
    motifMarkov = new NonUniformMarkov(optionsMFinder.motifOrder, optionsMFinder.width, numAlph);
    motifMarkov->construct(&motifCounts, optionsMFinder.pcounts);
    
    // build spacer distribution
    // build histogram from positions
    vector<double> positionCounts (upstreamLength - optionsMFinder.width+1, 0);
    for (size_t n = 0; n < positions.size(); n++) {
        // FIXME account for LEFT alignment
        // below is only for right
        positionCounts[upstreamLength - optionsMFinder.width - positions[n]]++;        // increment position
    }
    
    motifSpacer = new UnivariatePDF(positionCounts, false, optionsMFinder.pcounts);
    
    
}

vector<GeneStat> separateLabelsViaOperonStatus(const vector<Label*> &labels, unsigned thresh, unsigned threshNFIO) {
    
    
    vector<GeneStat> result (labels.size());
    
    // for every label
    for (size_t n = 0; n < labels.size(); n++) {
        
        // get label
        const Label* lab = labels[n];
        
        // if on positive strand
        if (lab->strand == Label::POS) {
            
            size_t start = lab->left;          // start location
            
            // if no gene before it, set as first in operon
            if (n == 0) {
                result[n] = FIRST_OP;
            }
            // otherwise, check to see if stop of previous gene is nearby
            else {
                
                const Label* prevLab = labels[n-1];         // previous label
                
                if (prevLab->strand == Label::POS) {   // should be on positive strand
                    
                    size_t prevStop = prevLab->right;
                    
                    // if too far away, then current gene is first in op
                    if (prevStop < start - thresh)
                        result[n] = FIRST_OP;
                    // if threshNFIO is set and gene is too close
                    else if (threshNFIO != 0) {
                        if (prevStop > start - threshNFIO)
                            result[n] = NOT_FIRST_OP;
                        else
                            result[n] = IGNORE;
                    }
                    // otherwise, current gene belongs to same operon as previous
                    else
                        result[n] = NOT_FIRST_OP;
                    
                }
                // otherwise it's on negative strand, so first in operon
                else
                    result[n] = FIRST_OP;
            }
        }
        
        // if on negative strand
        else if (lab->strand == Label::NEG) {
            
            size_t start = lab->right;        // start location
            
            // if no gene after it, then set as first in operon
            if (n == labels.size()-1) {
                result[n] = FIRST_OP;
            }
            // otherwise, check to see if stop of previous gene (i.e. to the right) is nearby
            else {
                const Label* prevLab = labels[n+1];     // "next" label
                
                if (prevLab->strand == Label::NEG) {      // should be on negative strand
                    
                    size_t prevStop = prevLab->left;
                    
                    // if too far away, then current gene is first in op
                    if (prevStop > start + thresh)
                        result[n] = FIRST_OP;
                    // if threshNFIO is set and gene is too close
                    else if (threshNFIO != 0) {
                        if (prevStop < start + threshNFIO)
                            result[n] = NOT_FIRST_OP;
                        else
                            result[n] = IGNORE;
                    }
                    // otherwise, current gene belongs to same operon as previous
                    else
                        result[n] = NOT_FIRST_OP;
                }
                // oterhwise, previous gene is on different strand, so current is first in operon
                else
                    result[n] = FIRST_OP;
            }
        }
        
    }
    
    
    return result;
    
}


void GMS2Trainer::estimateParametersMotifModel_Promoter(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    
    // copy only usable labels
    size_t numUse = labels.size();
    if (use.size() > 0) {
        numUse = 0;
        for (size_t n = 0; n < labels.size(); n++)
            if (use[n])
                numUse++;
    }
    
    vector<Label*> useLabels (numUse);
    size_t pos = 0;
    for (size_t n = 0; n < labels.size(); n++) {
        if (use[n])
            useLabels[pos++] = labels[n];
    }
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(useLabels, FGIO_DIST_THRESH, NFGIO_DIST_THRES, operonStatuses);
    
    size_t numFGIO = 0, numIG = 0, numUNK = 0;
    for (size_t n = 0; n < operonStatuses.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)        numFGIO++;
        else if (operonStatuses[n] == LabelsParser::NFGIO)  numIG++;
        else
            numUNK++;
    }
    
    // get FGIO and IG upstreams and run motif search
    vector<Label*> labelsFGIO (numFGIO);
    vector<Label*> labelsIG (numIG);
    size_t currFGIO = 0, currIG = 0;        // indices
    
    for (size_t n = 0; n < operonStatuses.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)        labelsFGIO[currFGIO++] = useLabels[n];
        else if (operonStatuses[n] == LabelsParser::NFGIO)  labelsIG[currIG++] = useLabels[n];
    }
    
    vector<NumSequence> upstreamsFGIO;
    SequenceParser::extractUpstreamSequences(sequence, labelsFGIO, cnc, this->UPSTR_LEN_FGIO, upstreamsFGIO);
    
    vector<NumSequence> upstreamsIG;
    SequenceParser::extractUpstreamSequences(sequence, labelsIG, cnc, this->UPSTR_LEN_IG, upstreamsIG);
    
//    void runMotifFinder(const vector<NumSequence> &sequencesRaw, OptionsMFinder &optionsMFinder, size_t upstreamLength, NonUniformMarkov* motifMarkov, UnivariatePDF* motifSpacer) {
//
    this->numLeaderless = upstreamsFGIO.size();
    this->numFGIO = upstreamsFGIO.size();

    OptionsMFinder optionMFinderFGIO (*this->optionsMFinder);
    optionMFinderFGIO.width = widthArchaeaPromoter;
    
    runMotifFinder(upstreamsFGIO, optionMFinderFGIO, *this->alphabet, this->UPSTR_LEN_FGIO, this->promoter, this->promoterSpacer);
    runMotifFinder(upstreamsIG, *this->optionsMFinder, *this->alphabet, this->UPSTR_LEN_IG, this->rbs, this->rbsSpacer);
    
}



void GMS2Trainer::estimateParametersMotifModel_Tuberculosis(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    
    // copy only usable labels
    size_t numUse = labels.size();
    if (use.size() > 0) {
        numUse = 0;
        for (size_t n = 0; n < labels.size(); n++)
            if (use[n])
                numUse++;
    }
    
    vector<Label*> useLabels (numUse);
    size_t pos = 0;
    for (size_t n = 0; n < labels.size(); n++) {
        if (use[n])
            useLabels[pos++] = labels[n];
    }
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(useLabels, FGIO_DIST_THRESH, NFGIO_DIST_THRES, operonStatuses);
    
    size_t numFGIO = 0, numIG = 0, numUNK = 0;
    for (size_t n = 0; n < operonStatuses.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)        numFGIO++;
        else if (operonStatuses[n] == LabelsParser::NFGIO)  numIG++;
        else
            numUNK++;
    }
    
    // get FGIO and IG upstreams and run motif search
    vector<Label*> labelsFGIO (numFGIO);
    vector<Label*> labelsIG (numIG);
    size_t currFGIO = 0, currIG = 0;        // indices
    
    for (size_t n = 0; n < operonStatuses.size(); n++) {
        if (operonStatuses[n] == LabelsParser::FGIO)        labelsFGIO[currFGIO++] = useLabels[n];
        else if (operonStatuses[n] == LabelsParser::NFGIO)  labelsIG[currIG++] = useLabels[n];
    }
    
    vector<NumSequence> upstreamsRBS;
    vector<NumSequence> upstreamsPromoter;
    
    size_t upstrLen = 20;
    
    // match FGIO to 16S tail
    vector<NumSequence> upstreamsFGIO;
    SequenceParser::extractUpstreamSequences(sequence, labelsFGIO, cnc, upstrLen, upstreamsFGIO);
    
    Sequence strMatchSeq (matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));

    size_t skipFromStart = 3;
    for (size_t n = 0; n < upstreamsFGIO.size(); n++) {
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreamsFGIO[n], positionsOfMatches, substitutions);

        // keep track of nonmatches
        if (match.size() < matchThresh)
            upstreamsPromoter.push_back(upstreamsFGIO[n].subseq(0, upstreamsFGIO[n].size() - skipFromStart));
        else
            upstreamsRBS.push_back(upstreamsFGIO[n]);
    }
    
    
    vector<NumSequence> upstreamsIG;
    SequenceParser::extractUpstreamSequences(sequence, labelsIG, cnc, upstrLen, upstreamsIG);
    for (size_t n = 0; n < upstreamsIG.size(); n++) {
        upstreamsRBS.push_back(upstreamsIG[n]);
    }
    
    this->numLeaderless = upstreamsPromoter.size();
    this->numFGIO = upstreamsFGIO.size();
    
    
    runMotifFinder(upstreamsPromoter, *this->optionsMFinder, *this->alphabet, upstrLen-skipFromStart, this->promoter, this->promoterSpacer);
    runMotifFinder(upstreamsRBS, *this->optionsMFinder, *this->alphabet, upstrLen, this->rbs, this->rbsSpacer);
    
    // shift probabilities
    vector<double> extendedProbs (promoterSpacer->size()+skipFromStart, 0);
    for (size_t n = 0; n < promoterSpacer->size(); n++) {
        extendedProbs[n+skipFromStart] = (*promoterSpacer)[n];
    }
    
    delete promoterSpacer;
    promoterSpacer = new UnivariatePDF(extendedProbs);
    
    
}


void GMS2Trainer::estimateParametersMotifModel_Synechocystis(const NumSequence &sequence, const vector<Label*> &labels, const vector<bool> &use) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // copy only usable labels
    size_t numUse = labels.size();
    if (use.size() > 0) {
        numUse = 0;
        for (size_t n = 0; n < labels.size(); n++)
            if (use[n])
                numUse++;
    }
    
    vector<Label*> useLabels (numUse);
    size_t pos = 0;
    for (size_t n = 0; n < labels.size(); n++) {
        if (use[n])
            useLabels[pos++] = labels[n];
    }
    
    size_t upstrLen = 20;
    Sequence strMatchSeq (matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    
    // extract upstream for every sequence and match it to 16S tail
    vector<NumSequence> upstreams (useLabels.size());
    SequenceParser::extractUpstreamSequences(sequence, useLabels, cnc, upstrLen, upstreams);
    size_t skipFromStart = 0;
    
    vector<Label*> labelsSig;
    vector<Label*> labelsRBS;
    
    for (size_t n = 0; n < upstreams.size(); n++) {
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreams[n], positionsOfMatches, substitutions);
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            labelsSig.push_back(useLabels[n]);
        else
            labelsRBS.push_back(useLabels[n]);
    }
    
    // for all non-Sig sequences, append "N" to mask Sig sequences
    string Ns = "";
    size_t numNs = 0;
    if ( this->upstreamSignatureLength >= this->startContextLength &&
         this->upstreamSignatureLength - this->startContextLength >= this->upstreamSignatureOrder)
        numNs = (this->upstreamSignatureLength - this->startContextLength) - this->upstreamSignatureOrder;
    
    for (size_t n = 0; n < numNs; n++)
        Ns += "N";
    
    NumSequence numSeqNs (Sequence(Ns), cnc);       // numeric sequence of N's
    
    NonUniformCounts counts(upstreamSignatureOrder, upstreamSignatureLength, *this->alphabet);
    
    vector<NumSequence> contextsRBS;
    long long posRelToStart = - (startContextLength + scMargin + upstreamSignatureOrder);
    SequenceParser::extractStartContextSequences(sequence, labelsRBS, cnc, posRelToStart, startContextLength + this->upstreamSignatureOrder, contextsRBS);

    for (size_t n = 0; n < contextsRBS.size(); n++) {
        NumSequence withNs = numSeqNs + contextsRBS[n];     // append N's
        counts.count(withNs.begin(), withNs.end());
//        cout << cnc.convert(withNs.begin(), withNs.end()) << endl;
    }
    
    // add Sig sequences
    vector<NumSequence> contextsSig;
    SequenceParser::extractStartContextSequences(sequence, labelsSig, cnc, -( (int)upstreamSignatureLength + scMargin), upstreamSignatureLength, contextsSig);
    
    for (size_t n = 0; n < contextsSig.size(); n++) {
        counts.count(contextsSig[n].begin(), contextsSig[n].end());
//        cout << cnc.convert(contextsSig[n].begin(), contextsSig[n].end()) << endl;
    }
    
    // for sequences without motifs
    startContext = new NonUniformMarkov(upstreamSignatureOrder, upstreamSignatureLength, *this->alphabet);
    startContext->construct(&counts, pcounts);
    
    // run motif search for RBS
    vector<NumSequence> upstreamsRBS;
    SequenceParser::extractUpstreamSequences(sequence, labelsRBS, cnc, this->upstreamLength, upstreamsRBS);
    runMotifFinder(upstreamsRBS, *this->optionsMFinder, *this->alphabet, this->upstreamLength, this->rbs, this->rbsSpacer);
    
    
    
}


// DEPRECATED
void GMS2Trainer::estimateParametersMotifModel_Promoter_DEPRECATED(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    // split sequences into first in operon and the rest
    
    // FIXME: make it more efficient
    
    // copy only usable labels
    size_t numUse = labels.size();
    if (use.size() > 0) {
        numUse = 0;
        for (size_t n = 0; n < labels.size(); n++)
            if (use[n])
                numUse++;
    }
    
    vector<Label*> useLabels (numUse);
    size_t pos = 0;
    for (size_t n = 0; n < labels.size(); n++) {
        if (use[n])
            useLabels[pos++] = labels[n];
    }
    
    
    // get upstream regions
    vector<NumSequence> upstreamsFGIO_style;
    vector<NumSequence> upstreamsNFGIO_style;
    
    long long upstr_min_value = MIN_UPSTR_LEN_FGIO;
    NumSequence::size_type upstrLength = this->upstreamLength - upstr_min_value;
    
    SequenceParser::extractStartContextSequences(sequence, useLabels, *alphabet->getCNC(), -this->upstreamLength, upstrLength, upstreamsFGIO_style);
    
    if (UPSTR_LEN_NFGIO == 0)
        upstreamsNFGIO_style = upstreamsFGIO_style;
    else
        SequenceParser::extractStartContextSequences(sequence, useLabels, *alphabet->getCNC(), -this->UPSTR_LEN_NFGIO, this->UPSTR_LEN_NFGIO, upstreamsNFGIO_style);
    
    
    vector<GeneStat> status = separateLabelsViaOperonStatus(useLabels, FGIO_DIST_THRESH, NFGIO_DIST_THRES);
    
    vector<NumSequence> upstreamsFGIO;
    vector<NumSequence> upstreamsNFGIO;
    
    // split sequences
    for (size_t n = 0; n < status.size(); n++) {
        if (status[n] == IGNORE)
            continue;
        
        if (status[n] == FIRST_OP)
            upstreamsFGIO.push_back(upstreamsFGIO_style[n]);
        else
            upstreamsNFGIO.push_back(upstreamsNFGIO_style[n]);
    }
    
    
    // build motif finder
    MotifFinder::Builder b;
    
    b.setAlign(optionsMFinder->align);
    b.setWidth(optionsMFinder->width);
    b.setMaxIter(optionsMFinder->maxIter);
    b.setPcounts(optionsMFinder->pcounts);
    b.setNumTries(optionsMFinder->tries);
    b.setBackOrder(optionsMFinder->bkgdOrder);
    b.setMaxEMIter(optionsMFinder->maxEMIter);
    b.setMotifOrder(optionsMFinder->motifOrder);
    b.setShiftEvery(optionsMFinder->shiftEvery);
    
    MotifFinder mfinder = b.build();
    
    // remove N's
    vector<NumSequence> upstreamsFGIO_noN;          // upstreams of first genes in operon that don't contain N
    vector<NumSequence> upstreamsNFGIO_noN;         // upstreams of non-first genes in operon that don't contain N
    
    filterNs(upstreamsFGIO, upstreamsFGIO_noN);
    filterNs(upstreamsNFGIO, upstreamsNFGIO_noN);
    
    vector<NumSequence::size_type> positionsFGIO, positionsNFGIO;
    
    // run MFinder in operons
    mfinder.findMotifs(upstreamsFGIO_noN, positionsFGIO);
    
    // build Promoter model
    NonUniformCounts promoterCounts(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet);
    for (size_t n = 0; n < upstreamsFGIO_noN.size(); n++) {
        promoterCounts.count(upstreamsFGIO_noN[n].begin()+positionsFGIO[n], upstreamsFGIO_noN[n].begin()+positionsFGIO[n]+optionsMFinder->width);
    }
    
    promoter = new NonUniformMarkov(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet);
    promoter->construct(&promoterCounts, optionsMFinder->pcounts);
    
    // build spacer distribution
    // build histogram from positions
    vector<double> positionCounts_promoter (this->upstreamLength - optionsMFinder->width+1, 0);
    for (size_t n = 0; n < positionsFGIO.size(); n++) {
        // FIXME account for LEFT alignment
        // below is only for right
        positionCounts_promoter[upstreamLength - optionsMFinder->width - positionsFGIO[n]]++;        // increment position
    }
    
    promoterSpacer = new UnivariatePDF(positionCounts_promoter, false, pcounts);
    
    
    // run MFinder on RBS
    mfinder.findMotifs(upstreamsNFGIO_noN, positionsNFGIO);
    
    
    // build RBS model
    NonUniformCounts rbsCounts(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet);
    for (size_t n = 0; n < upstreamsNFGIO_noN.size(); n++) {
        rbsCounts.count(upstreamsNFGIO_noN[n].begin()+positionsNFGIO[n], upstreamsNFGIO_noN[n].begin()+positionsNFGIO[n]+optionsMFinder->width);
    }
    
    rbs = new NonUniformMarkov(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet);
    rbs->construct(&rbsCounts, optionsMFinder->pcounts);
    
    // build spacer distribution
    // build histogram from positions
    size_t upstrLenNFGIO = (UPSTR_LEN_NFGIO == 0 ? upstreamLength : UPSTR_LEN_NFGIO);
    vector<double> positionCounts_rbs (upstrLenNFGIO - optionsMFinder->width+1, 0);
    for (size_t n = 0; n < positionsNFGIO.size(); n++) {
        // FIXME account for LEFT alignment
        // below is only for right
        positionCounts_rbs[upstrLenNFGIO - optionsMFinder->width - positionsFGIO[n]]++;        // increment position
    }
    
    rbsSpacer = new UnivariatePDF(positionCounts_rbs, false, pcounts);
    
    
    
}

void GMS2Trainer::estimateParameters(const NumSequence &sequence, const vector<gmsuite::Label *> &labels) {
    
    // reset all models
    deallocAllModels();
    
    // TODO: sort labels by 'left' in increasing order
    vector<bool> useCoding (labels.size(), true);               // assume all labels should be used for coding model
    vector<bool> useNonCoding (labels.size(), true);            // assume all labels should be used for noncoding model
    vector<bool> useStartContext (labels.size(), true);         // assume all labels should be used for startcontext model
    vector<bool> useMotif (labels.size(), true);                // assume all labels should be used for motif search model
    
    
    // estimate parameters for coding model
    selectLabelsForCodingParameters(labels, useCoding);
    estimateParamtersCoding(sequence, labels, startContextLength, useCoding);
    
    estimateParametersStartStopCodons(sequence, labels);
    
    // estimate parameters for noncoding model
    useNonCoding = useCoding;                           // for now, assume
    estimateParamtersNonCoding(sequence, labels);
    
    // estimate parameters for start context
    useStartContext = useCoding;
    estimateParametersStartContext(sequence, labels);
    
    // estimate parameters for motif models
    if (runMotifSearch) {
        useMotif = useCoding;
        estimateParametersMotifModel(sequence, labels, useMotif);
    }
}




void GMS2Trainer::deallocAllModels() {
    
    // dealloc public variables for models
    if (noncoding != NULL) delete noncoding;
    if (coding != NULL) delete coding;
    if (startContext != NULL) delete startContext;
    if (rbs != NULL) delete rbs;
    if (promoter != NULL) delete promoter;
    if (upstreamSignature != NULL) delete upstreamSignature;
    if (rbsSpacer != NULL) delete rbsSpacer;
    if (promoterSpacer != NULL) delete promoterSpacer;
}

// destructor
GMS2Trainer::~GMS2Trainer() {
    deallocAllModels();
}


void codonFrequencyToMod(const map<CharNumConverter::seq_t, double> &codons, const CharNumConverter &cnc, vector<pair<string, string> > &toMod) {
    
    for (map<CharNumConverter::seq_t, double>::const_iterator iter = codons.begin(); iter != codons.end(); iter++) {
        stringstream ssm; ssm << fixed;
        
        string cod = cnc.convert(iter->first.begin(), iter->first.end());
        ssm << iter->second;
        
        toMod.push_back(pair<string,string>(cod, ssm.str()));
    }
    
}


void GMS2Trainer::toModFile(vector<pair<string, string> > &toMod, const OptionsGMS2Training &options) const {
    
    typedef pair<string, string> mpair;
    
    // name and genetic code
    toMod.push_back(mpair("NAME", "gms2-training"));
    toMod.push_back(mpair("GCODE", this->numGeneticCode->getName()));
    toMod.push_back(mpair("NON_DURATION_DECAY", boost::lexical_cast<string>(options.nonDurationDecay)));
    toMod.push_back(mpair("COD_DURATION_DECAY", boost::lexical_cast<string>(options.codDurationDecay)));
    toMod.push_back(mpair("COD_P_N", boost::lexical_cast<string>(options.codProbN)));
    toMod.push_back(mpair("NON_P_N", boost::lexical_cast<string>(options.nonProbN)));
    toMod.push_back(mpair("GENE_MIN_LENGTH", boost::lexical_cast<string>(options.geneMinLengthPrediction)));
    
    
    // add start/stop codon probabilities
    codonFrequencyToMod(startProbs, *this->alphabet->getCNC(), toMod);
    codonFrequencyToMod(stopProbs,  *this->alphabet->getCNC(), toMod);
    
    // add description to mod file
    if (coding != NULL) {
        toMod.push_back(mpair("COD_ORDER", boost::lexical_cast<string>(coding->getOrder())));
        toMod.push_back(mpair("COD_MAT", coding->toString()));
    }
    
    if (noncoding != NULL) {
        toMod.push_back(mpair("NON_ORDER", boost::lexical_cast<string>(noncoding->getOrder())));
        toMod.push_back(mpair("NON_MAT",  noncoding->toString()));
    }
    
    if (startContext != NULL) {
        toMod.push_back(mpair("SC", "1"));
        toMod.push_back(mpair("SC_ORDER", boost::lexical_cast<string>(startContext->getOrder())));
        toMod.push_back(mpair("SC_WIDTH", boost::lexical_cast<string>(startContext->getLength())));
        toMod.push_back(mpair("SC_MARGIN", boost::lexical_cast<string>(scMargin)));
        toMod.push_back(mpair("SC_MAT", startContext->toString()));
    }
    
    if (rbs != NULL) {
        toMod.push_back(mpair("RBS", "1"));
        toMod.push_back(mpair("RBS_ORDER", boost::lexical_cast<string>(rbs->getOrder())));
        toMod.push_back(mpair("RBS_WIDTH", boost::lexical_cast<string>(rbs->getLength())));
        toMod.push_back(mpair("RBS_MARGIN", "0"));
        toMod.push_back(mpair("RBS_MAT", rbs->toString()));
        
        if (startContextRBS != NULL) {
            toMod.push_back(mpair("SC_RBS", "1"));
            toMod.push_back(mpair("SC_RBS_ORDER", boost::lexical_cast<string>(startContextRBS->getOrder())));
            toMod.push_back(mpair("SC_RBS_WIDTH", boost::lexical_cast<string>(startContextRBS->getLength())));
            toMod.push_back(mpair("SC_RBS_MARGIN", boost::lexical_cast<string>(scMargin)));
            toMod.push_back(mpair("SC_RBS_MAT", startContextRBS->toString()));
        }
    }
    
    if (rbsSpacer != NULL && rbsSpacer->size() > 0) {
        toMod.push_back(mpair("RBS_MAX_DUR", boost::lexical_cast<string>(rbsSpacer->size() - 1)));
        toMod.push_back(mpair("RBS_POS_DISTR", rbsSpacer->toString()));
    }
    
    toMod.push_back(mpair("PROMOTER_NUM_FGIO", boost::lexical_cast<string>(this->numFGIO)));
    
    if (promoter != NULL) {
        toMod.push_back(mpair("PROMOTER", "1"));
        toMod.push_back(mpair("PROMOTER_ORDER", boost::lexical_cast<string>(promoter->getOrder())));
        toMod.push_back(mpair("PROMOTER_WIDTH", boost::lexical_cast<string>(promoter->getLength())));
        toMod.push_back(mpair("PROMOTER_MARGIN", "0"));
        toMod.push_back(mpair("PROMOTER_MAT", promoter->toString()));
        toMod.push_back(mpair("PROMOTER_NUM_LEADERLESS", boost::lexical_cast<string>(this->numLeaderless)));
        
        
        if (startContextPromoter != NULL) {
            toMod.push_back(mpair("SC_PROMOTER", "1"));
            toMod.push_back(mpair("SC_PROMOTER_ORDER", boost::lexical_cast<string>(startContextPromoter->getOrder())));
            toMod.push_back(mpair("SC_PROMOTER_WIDTH", boost::lexical_cast<string>(startContextPromoter->getLength())));
            toMod.push_back(mpair("SC_PROMOTER_MARGIN", boost::lexical_cast<string>(scMargin)));
            toMod.push_back(mpair("SC_PROMOTER_MAT", startContextPromoter->toString()));
        }
        else if (startContextRBS != NULL) {                 // copy promoter start context from RBS start context
            toMod.push_back(mpair("SC_PROMOTER", "1"));
            toMod.push_back(mpair("SC_PROMOTER_ORDER", boost::lexical_cast<string>(startContextRBS->getOrder())));
            toMod.push_back(mpair("SC_PROMOTER_WIDTH", boost::lexical_cast<string>(startContextRBS->getLength())));
            toMod.push_back(mpair("SC_PROMOTER_MARGIN", boost::lexical_cast<string>(scMargin)));
            toMod.push_back(mpair("SC_PROMOTER_MAT", startContextRBS->toString()));
        }
    }
    
    if (promoterSpacer != NULL && promoterSpacer->size() > 0) {
        toMod.push_back(mpair("PROMOTER_MAX_DUR", boost::lexical_cast<string>(promoterSpacer->size() - 1)));
        toMod.push_back(mpair("PROMOTER_POS_DISTR", promoterSpacer->toString()));
    }

    toMod.push_back(mpair("GENOME_TYPE", this->genomeType));
    
//    if (upstreamSignature != NULL) {
//        int upstrSigMargin = 0;
//        if (startContext != NULL)
//            upstrSigMargin= max(0, (int)startContext->getLength() - scMargin);
//            
//        toMod.push_back(mpair("UPSTR_SIG_ORDER", "1"));
//        toMod.push_back(mpair("UPSTR_SIG_ORDER", boost::lexical_cast<string>(upstreamSignature->getOrder())));
//        toMod.push_back(mpair("UPSTR_SIG_WIDTH", boost::lexical_cast<string>(upstreamSignature->getLength())));
//        toMod.push_back(mpair("UPSTR_SIGN_MARGIN", boost::lexical_cast<string>(scMargin)));
//        toMod.push_back(mpair("UPSTR_SIG_MAT", upstreamSignature->toString()));
//    }
}




// assume useCoding is same size as labels vector
void GMS2Trainer::selectLabelsForCodingParameters(const vector<Label*> &labels, vector<bool> &useCoding) const {
    
    if (useCoding.size() != labels.size())
        throw invalid_argument("Labels and useCoding vectors should have the same size");
    
    
    for (size_t n = 0; n < labels.size(); n++) {
        if (labels[n] == NULL)
            throw invalid_argument("Label cannot be null");
        
        // "remove" short genes
        if (labels[n]->right - labels[n]->left + 1 < MIN_GENE_LEN)
            useCoding[n] = false;
        
        // train on native
        if (trainOnNative) {
            useCoding[n] = (labels[n]->geneClass.find("native") != string::npos);
        }
    }
}













