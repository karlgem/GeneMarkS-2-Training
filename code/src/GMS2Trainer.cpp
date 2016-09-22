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
                         const CharNumConverter &cnc,
                         const AlphabetDNA &alph,
                         const NumSequence::size_type MIN_GENE_LEN) {
    
    this->pcounts = pcounts;
    this->codingOrder = codingOrder;
    this->noncodingOrder = noncodingOrder;
    this->startContextOrder = startContextOrder;
    this->upstreamLength = upstreamLength;
    this->startContextLength = startContextLength;
    this->genomeClass = genomeClass;
    this->optionsMFinder = &optionsMFinder;
    this->cnc = &cnc;
    this->alphabet = &alph;
    this->MIN_GENE_LEN = MIN_GENE_LEN;
    
    
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

void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels, NumSequence::size_type scSize, const vector<bool> &use) {
    
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
    PeriodicCounts counts (codingOrder, 3, *this->alphabet, *this->cnc);
    
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
    coding = new PeriodicMarkov(codingOrder, 3, *this->alphabet, *this->cnc);
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
    
    UniformCounts counts(noncodingOrder, *this->alphabet, *this->cnc);
    
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
    noncoding = new UniformMarkov(noncodingOrder, *this->alphabet, *this->cnc);
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
    
    NonUniformCounts counts(startContextOrder, startContextLength, *this->alphabet, *this->cnc);
    
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
    startContext = new NonUniformMarkov(startContextOrder, startContextLength, *this->alphabet, *this->cnc);
    startContext->construct(&counts, pcounts);
}


void GMS2Trainer::estimateParametersMotifModel(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    // check if all labels should be used
    if (use.size() > 0) {
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
    // build motif finder
    MotifFinder::Builder b;
    if (optionsMFinder->align == "left")
        b.setAlign(MFinderModelParams::LEFT);
    else if (optionsMFinder->align == "right")
        b.setAlign(MFinderModelParams::RIGHT);
    else
        b.setAlign(MFinderModelParams::NONE);
    
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
        vector<NumSequence> upstreams; 
        SequenceParser::extractUpstreamSequences(sequence, labels, *cnc, upstreamLength, upstreams, false, MIN_GENE_LEN, use);
        
        vector<NumSequence::size_type> positions;
        mfinder.findMotifs(upstreams, positions);
        
        // build RBS model
        NonUniformCounts rbsCounts(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet, *this->cnc);
        for (size_t n = 0; n < upstreams.size(); n++) {
            rbsCounts.count(upstreams[n].begin()+positions[n], upstreams[n].begin()+positions[n]+optionsMFinder->width);
        }
        
        rbs = new NonUniformMarkov(optionsMFinder->motifOrder, optionsMFinder->width, *this->alphabet, *this->cnc);
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
    }
    
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



void GMS2Trainer::toModFile(map<string, string> &toMod) const {
    
    // dealloc public variables for models
    if (noncoding != NULL) {
        toMod["NONC_ORDER"] = boost::lexical_cast<string>(noncoding->getOrder());
        toMod["NONC_MAT"] =  noncoding->toString();
    }
    
    if (coding != NULL) {
        toMod["COD_ORDER"] = boost::lexical_cast<string>(coding->getOrder());
        toMod["COD_MAT"] = coding->toString();
    }
    
    if (startContext != NULL) {
        toMod["SC_ORDER"] = boost::lexical_cast<string>(startContext->getOrder());
        toMod["SC_MAT"] = startContext->toString();
    }
    
    if (rbs != NULL) {
        toMod["RBS_ORDER"] = boost::lexical_cast<string>(rbs->getOrder());
        toMod["RBS_WIDTH"] = boost::lexical_cast<string>(rbs->getLength());
        toMod["RBS_MAT"] = rbs->toString();
    }
    
    if (promoter != NULL) {
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













