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
#include "NonUniformCounts.hpp"
#include "SequenceParser.hpp"
#include <boost/lexical_cast.hpp>

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
                         const NumGeneticCode &numGenCode) {
    
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
        this->stopProbs[stops[n]] = 1;              // stop probabilities to 1
    
    size_t totalStarts = 0;
    
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        
        CharNumConverter::seq_t start;
        
        // get start from positive strand
        if ((*iter)->strand == Label::POS) {
            start = CharNumConverter::seq_t(sequence.begin() + (*iter)->left, sequence.begin() + (*iter)->left + 3);
        }
        // get start from negative strand and reverse complement
        else {
            start = CharNumConverter::seq_t (sequence.begin() + (*iter)->right-2, sequence.begin()+(*iter)->right+1);
            NumSequence s(start); s.reverseComplement(*this->alphabet->getCNC());   // reverse complement
            start = CharNumConverter::seq_t (s.begin(), s.end());
        }
        
        // if it is a start, add it to counts
        if (this->numGeneticCode->isStart(start)) {
            this->startProbs[start] += 1;
            totalStarts += 1;
        }
    }
    
    // normalize start counts
    if (totalStarts > 0)
        for (map<CharNumConverter::seq_t, double>::iterator  iter = this->startProbs.begin(); iter != this->startProbs.end(); iter++)
            iter->second /= totalStarts;
    
    
}

void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels, NumSequence::size_type scSize, const vector<bool> &use) {
    
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
    PeriodicCounts counts (codingOrder, 3, *this->alphabet);
    
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
        
        // if start context > 0, remove segement near start
        if (scSize > 0) {
            if (reverseComplement)
                counts.decount(sequence.begin()+right-scSize+1, sequence.begin()+right+1, reverseComplement);
            else
                counts.decount(sequence.begin()+left, sequence.begin()+left+scSize);
        }
        
    }
    
    // convert counts to probabilities
    coding = new PeriodicMarkov(codingOrder, 3, *this->alphabet);
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
    
    UniformCounts counts(noncodingOrder, *this->alphabet);
    
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
    
    NonUniformCounts counts(startContextOrder, startContextLength, *this->alphabet);
    
    // get counts for start context model
    size_t n = 0;
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        if (!useAll && !use[n++])
            continue;       // skip unwanted genes
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        // skip sequence if it doesn't have enough nucleotides for start context
        if ((*iter)->right - (*iter)->left + 1 - 6 < startContextLength)        // get lenght of sequence, and -6 to account for start/stop codons
            continue;
        
        size_t left;
        
        if ((*iter)->strand == Label::POS)
            left = (*iter)->left + 3;
        else
            left = (*iter)->right - 2 - startContextLength;
        
        bool reverseComplement = (*iter)->strand == Label::NEG;
        counts.count(sequence.begin() + left, sequence.begin() + left + startContextLength, reverseComplement);
    }
    
    // convert counts to probabilities
    startContext = new NonUniformMarkov(startContextOrder, startContextLength, *this->alphabet);
    startContext->construct(&counts, pcounts);
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
    if (genomeClass == ProkGeneStartModel::C1) {
        
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
    // if genome is class 2, search for weak RBS, and estimate upstream signature pwm
    else if (genomeClass == ProkGeneStartModel::C2) {
        throw logic_error("Code not yet completed.");
    }
    // if genome is class 3, search for RBS and promoter
    else {
        throw logic_error("Code not yet completed.");
        estimateParametersMotifModel_Promoter(sequence, labels, use);
    }
    
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
    useMotif = useCoding;
    estimateParametersMotifModel(sequence, labels);
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


void codonFrequencyToMod(const map<CharNumConverter::seq_t, double> &codons, const CharNumConverter &cnc, map<string, string> &toMod) {
    
    for (map<CharNumConverter::seq_t, double>::const_iterator iter = codons.begin(); iter != codons.end(); iter++) {
        stringstream ssm;
        ssm << cnc.convert(iter->first.begin(), iter->first.end()); // << "\t" << iter->second << endl;
        
        toMod[ssm.str()] = boost::lexical_cast<string>(iter->second);
    }
    
}


void GMS2Trainer::toModFile(map<string, string> &toMod) const {
    
    
    
    // name and genetic code
    toMod["NAME"] = "gms2-training";
    toMod["GCODE"] = "11";
    toMod["NON_DURATION_DECAY"] = "150";
    toMod["COD_DURATION_DECAY"] = "300";
    toMod["COD_P_N"] = "0.4";
    toMod["NON_P_N"] = "0.6";
    toMod["GENE_MIN_LENGTH"] = "89";
    
    
    // add start/stop codon probabilities
    codonFrequencyToMod(startProbs, *this->alphabet->getCNC(), toMod);
    codonFrequencyToMod(stopProbs,  *this->alphabet->getCNC(), toMod);
    
    // add description to mod file
    if (noncoding != NULL) {
        toMod["NON_ORDER"] = boost::lexical_cast<string>(noncoding->getOrder());
        toMod["NON_MAT"] =  noncoding->toString();
    }
    
    if (coding != NULL) {
        toMod["COD_ORDER"] = boost::lexical_cast<string>(coding->getOrder());
        toMod["COD_MAT"] = coding->toString();
    }
    
    if (startContext != NULL) {
        toMod["SC"] = "1";
        toMod["SC_ORDER"] = boost::lexical_cast<string>(startContext->getOrder());
        toMod["SC_WIDTH"] = boost::lexical_cast<string>(startContext->getLength());
        toMod["SC_MAT"] = startContext->toString();
        toMod["SC_POS_DISTR"] = "\n" + boost::lexical_cast<string>(-3 - (int)startContext->getLength()) + "\t1.0\n";
    }
    
    if (rbs != NULL) {
        toMod["RBS"] = "1";
        toMod["RBS_ORDER"] = boost::lexical_cast<string>(rbs->getOrder());
        toMod["RBS_WIDTH"] = boost::lexical_cast<string>(rbs->getLength());
        toMod["RBS_MAT"] = rbs->toString();
    }
    
    if (promoter != NULL) {
        toMod["PROMOTER"] = "1";
        toMod["PROMOTER_ORDER"] = boost::lexical_cast<string>(promoter->getOrder());
        toMod["PROMOTER_WIDTH"] = boost::lexical_cast<string>(promoter->getLength());
        toMod["PROMOTER_MAT"] = promoter->toString();
    }
    
    if (upstreamSignature != NULL) {
        toMod["UPSTR_SIG_ORDER"] = boost::lexical_cast<string>(upstreamSignature->getOrder());
        toMod["UPSTR_SIG_WIDTH"] = boost::lexical_cast<string>(upstreamSignature->getLength());
        toMod["UPSTR_SIG_MAT"] = upstreamSignature->toString();
    }
    
    if (rbsSpacer != NULL) {
        toMod["RBS_POS_DISTR"] = rbsSpacer->toString();
    }
    
    if (promoterSpacer != NULL) {
        toMod["PROMOTER_POS_DISTR"] = promoterSpacer->toString();
    }
    
}




// assume useCoding is same size as labels vector
void GMS2Trainer::selectLabelsForCodingParameters(const vector<Label*> &labels, vector<bool> &useCoding) const {
    
    if (useCoding.size() != labels.size())
        throw invalid_argument("Labels and useCoding vectors should have the same size");
    
    // "remove" short genes
    for (size_t n = 0; n < labels.size(); n++) {
        if (labels[n] == NULL)
            throw invalid_argument("Label cannot be null");
        
        if (labels[n]->right - labels[n]->left + 1 < MIN_GENE_LEN)
            useCoding[n] = false;
    }
    
    
    
}













