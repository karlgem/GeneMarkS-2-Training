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

// Function Prototypes
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
                         genome_group_t genomeClass,
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
                         unsigned FGIO_DIST_THRESH,
                         bool cutPromTrainSeqs) {
    
//    this->params.pcounts = pcounts;
//    this->params.orderCoding = codingOrder;
//    this->params.orderNonCoding = noncodingOrder;
//    this->params.orderStartContext = startContextOrder;
//    this->upstreamLength = upstreamLength;
//    this->startContextLength = startContextLength;
//    this->genomeClass = genomeClass;
//    this->optionsMFinder = &optionsMFinder;
//    this->alphabet = &alph;
//    this->MIN_GENE_LEN = MIN_GENE_LEN;
//    this->numGeneticCode = &numGenCode;
//    this->scMargin = scMargin;
//    this->trainOnNative = trainOnNative;
//    this->runMotifSearch = runMotifSearch;
//    this->UPSTR_LEN_FGIO = upstrFGIO;
//    this->widthArchaeaPromoter = widthArchaeaPromoter;
//    this->upstreamSignatureLength = upstreamSignatureLength;
//    this->upstreamSignatureOrder = upstreamSignatureOrder;
//    this->trainNonCodingOnFullGenome = trainNonCodingOnFullGenome;
//
//    this->matchTo = matchTo;
//    this->allowAGSubstitution = allowAGSubstitution;
//    this->matchThresh = matchThresh;
//
//    this->FGIO_DIST_THRESH = FGIO_DIST_THRESH;
//    this->NFGIO_DIST_THRES = 22;
//    this->cutPromTrainSeqs = cutPromTrainSeqs;
//
//    if (genomeClass == ProkGeneStartModel::A || genomeClass == ProkGeneStartModel::B || genomeClass == ProkGeneStartModel::C || genomeClass == ProkGeneStartModel::A2) {
//        this->UPSTR_LEN_IG = upstreamLength;
//    }
//
//    // public variables for models
//    noncoding = NULL;
//    coding = NULL;
//    startContext = NULL;
//    startContextRBS = NULL;
//    startContextPromoter = NULL;
//
//    // start models
//    rbs = NULL;
//    promoter = NULL;
//    upstreamSignature = NULL;
//    rbsSpacer = NULL;
//    promoterSpacer = NULL;
//
//    this->numLeaderless = 0;
//    this->numFGIO = 0;
//    this->genomeType = "no-motif";
}


GMS2Trainer::GMS2Trainer(
            // Coding and Noncoding Models
            unsigned                orderCoding                         ,
            unsigned                orderNonCoding                      ,
            unsigned                orderStartContext                   ,
            NumSequence::size_type  lengthStartContext                  ,
            int                     marginStartContext                  ,
            // Misc Variables
            NumSequence::size_type  fgioDistanceThresh                  ,
            NumSequence::size_type  igioDistanceThresh                  ,
            unsigned                pcounts                             ,
            genome_group_t          genomeGroup                         ,
            GeneticCode::gcode_t    gcode                               ,
            NumSequence::size_type  minimumGeneLengthTraining           ,
            bool                    onlyTrainOnNativeGenes              ,
            bool                    runMotifSearch                      ,
            const OptionsMFinder&   optionsMFinder                      ,
            // Group-A
            unsigned                groupA_widthPromoter                ,
            unsigned                groupA_widthRBS                     ,
            NumSequence::size_type  groupA_upstreamLengthPromoter       ,
            NumSequence::size_type  groupA_upstreamLengthRBS            ,
            double                  groupA_spacerScoreThresh            ,
            NumSequence::size_type  groupA_spacerDistThresh             ,
            NumSequence::size_type  groupA_spacerWindowSize             ,
             string                 groupA_extendedSD                   ,
             unsigned               groupA_minMatchToExtendedSD         ,
             bool                   groupA_allowAGSubstitution          ,
            // Group B
            unsigned                groupB_widthPromoter                ,
            unsigned                groupB_widthRBS                     ,
            NumSequence::size_type  groupB_upstreamLengthPromoter       ,
            NumSequence::size_type  groupB_upstreamLengthRBS            ,
            double                  groupB_spacerScoreThresh            ,
            NumSequence::size_type  groupB_spacerDistThresh             ,
            NumSequence::size_type  groupB_spacerWindowSize             ,
            string                  groupB_extendedSD                   ,
            unsigned                groupB_minMatchToExtendedSD         ,
            bool                    groupB_allowAGSubstitution          ,
            // Group-C
            unsigned                groupC_widthRBS                     ,
            NumSequence::size_type  groupC_upstreamLengthRBS            ,
            NumSequence::size_type  groupC_upstreamRegion3Prime         ,
            unsigned                groupC_minMatchRBSPromoter          ,
            unsigned                groupC_minMatchToExtendedSD         ,
            string                  groupC_extendedSD                   ,
            // Group-C2
            unsigned                groupC2_widthSDRBS                  ,
            unsigned                groupC2_widthNonSDRBS               ,
            NumSequence::size_type  groupC2_upstreamLengthSDRBS         ,
            NumSequence::size_type  groupC2_upstreamLengthNonSDRBS      ,
            NumSequence::size_type  groupC2_upstreamRegion3Prime        ,
            unsigned                groupC2_minMatchToExtendedSD        ,
            string                  groupC2_extendedSD                  ,
            // Group-D
            unsigned                groupD_widthRBS                     ,
            NumSequence::size_type  groupD_upstreamLengthRBS            ,
            double                  groupD_percentMatchRBS              ,
            string                  groupD_extendedSD                   ,
            unsigned                groupD_minMatchToExtendedSD         ,
            bool                    groupD_allowAGSubstitution          ,
            // Group-E
            unsigned                groupE_widthRBS                     ,
            NumSequence::size_type  groupE_upstreamLengthRBS            ,
            NumSequence::size_type  groupE_lengthUpstreamSignature      ,
            unsigned                groupE_orderUpstreamSignature       ,
            string                  groupE_extendedSD                   ,
            unsigned                groupE_minMatchToExtendedSD         ,
            bool                    groupE_allowAGSubstitution
                              ) {
    
    
    this->params.orderCoding                     =  orderCoding                           ;
    this->params.orderNonCoding                  =  orderNonCoding                        ;
    this->params.orderStartContext               =  orderStartContext                     ;
    this->params.lengthStartContext              =  lengthStartContext                    ;
    this->params.marginStartContext              =  marginStartContext                    ;
    this->params.fgioDistanceThresh              =  fgioDistanceThresh                    ;
    this->params.igioDistanceThresh              =  igioDistanceThresh                    ;
    this->params.pcounts                         =  pcounts                               ;
    this->params.genomeGroup                     =  genomeGroup                           ;
    this->params.gcode                           =  gcode                                 ;
    this->params.minimumGeneLengthTraining       =  minimumGeneLengthTraining             ;
    this->params.onlyTrainOnNativeGenes          =  onlyTrainOnNativeGenes                ;
    this->params.runMotifSearch                  =  runMotifSearch                        ;
    this->params.optionsMFinder                  =  &optionsMFinder                       ;
    this->params.groupA_widthPromoter            =  groupA_widthPromoter                  ;
    this->params.groupA_widthRBS                 =  groupA_widthRBS                       ;
    this->params.groupA_upstreamLengthPromoter   =  groupA_upstreamLengthPromoter         ;
    this->params.groupA_upstreamLengthRBS        =  groupA_upstreamLengthRBS              ;
    this->params.groupA_spacerScoreThresh        =  groupA_spacerScoreThresh              ;
    this->params.groupA_spacerDistThresh         =  groupA_spacerDistThresh               ;
    this->params.groupA_spacerWindowSize         =  groupA_spacerWindowSize               ;
    this->params.groupA_extendedSD               =  groupA_extendedSD                     ;
    this->params.groupA_minMatchToExtendedSD     =  groupA_minMatchToExtendedSD           ;
    this->params.groupA_allowAGSubstitution      =  groupA_allowAGSubstitution            ;
    this->params.groupB_widthPromoter            =  groupB_widthPromoter                  ;
    this->params.groupB_widthRBS                 =  groupB_widthRBS                       ;
    this->params.groupB_upstreamLengthPromoter   =  groupB_upstreamLengthPromoter         ;
    this->params.groupB_upstreamLengthRBS        =  groupB_upstreamLengthRBS              ;
    this->params.groupB_spacerScoreThresh        =  groupB_spacerScoreThresh              ;
    this->params.groupB_spacerDistThresh         =  groupB_spacerDistThresh               ;
    this->params.groupB_spacerWindowSize         =  groupB_spacerWindowSize               ;
    this->params.groupB_extendedSD               =  groupB_extendedSD                     ;
    this->params.groupB_minMatchToExtendedSD     =  groupB_minMatchToExtendedSD           ;
    this->params.groupB_allowAGSubstitution      =  groupB_allowAGSubstitution            ;
    this->params.groupC_widthRBS                 =  groupC_widthRBS                       ;
    this->params.groupC_upstreamLengthRBS        =  groupC_upstreamLengthRBS              ;
    this->params.groupC_upstreamRegion3Prime     =  groupC_upstreamRegion3Prime           ;
    this->params.groupC_minMatchRBSPromoter      =  groupC_minMatchRBSPromoter            ;
    this->params.groupC_minMatchToExtendedSD     =  groupC_minMatchToExtendedSD           ;
    this->params.groupC_extendedSD               =  groupC_extendedSD                     ;
    this->params.groupC2_widthSDRBS                  = groupC2_widthSDRBS                    ;
    this->params.groupC2_widthNonSDRBS               = groupC2_widthNonSDRBS                 ;
    this->params.groupC2_upstreamLengthSDRBS         = groupC2_upstreamLengthSDRBS           ;
    this->params.groupC2_upstreamLengthNonSDRBS      = groupC2_upstreamLengthNonSDRBS        ;
    this->params.groupC2_upstreamRegion3Prime        = groupC2_upstreamRegion3Prime          ;
    this->params.groupC2_minMatchToExtendedSD        = groupC2_minMatchToExtendedSD          ;
    this->params.groupC2_extendedSD                  = groupC2_extendedSD                    ;
    this->params.groupD_widthRBS                 =  groupD_widthRBS                       ;
    this->params.groupD_upstreamLengthRBS        =  groupD_upstreamLengthRBS              ;
    this->params.groupD_percentMatchRBS          =  groupD_percentMatchRBS                ;
    this->params.groupD_extendedSD               =  groupD_extendedSD                     ;
    this->params.groupD_minMatchToExtendedSD     =  groupD_minMatchToExtendedSD           ;
    this->params.groupD_allowAGSubstitution      =  groupD_allowAGSubstitution            ;
    this->params.groupE_widthRBS                 =  groupE_widthRBS                       ;
    this->params.groupE_upstreamLengthRBS        =  groupE_upstreamLengthRBS              ;
    this->params.groupE_lengthUpstreamSignature  =  groupE_lengthUpstreamSignature        ;
    this->params.groupE_orderUpstreamSignature   =  groupE_orderUpstreamSignature         ;
    this->params.groupE_extendedSD               =  groupE_extendedSD                     ;
    this->params.groupE_minMatchToExtendedSD     =  groupE_minMatchToExtendedSD           ;
    this->params.groupE_allowAGSubstitution      =  groupE_allowAGSubstitution            ;
    
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
    
    charAlph = new AlphabetDNA();
    cnc = new CharNumConverter(charAlph);
    alphabet = new NumAlphabetDNA(*charAlph, *cnc);
    geneticCode = new GeneticCode(params.gcode);
    numGeneticCode = new NumGeneticCode(*geneticCode, *cnc);
    
    cutPromTrainSeqs = false;
    
}

// Esimate parameters for start/stop codons
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
        iter->second += params.pcounts;
        totalStarts += params.pcounts;
    }
    
    for (map<CharNumConverter::seq_t, double>::iterator iter = this->stopProbs.begin(); iter != this->stopProbs.end(); iter++) {
        iter->second += params.pcounts;
        totalStops += params.pcounts;
    }
    
    // normalize start counts
    if (totalStarts > 0)
        for (map<CharNumConverter::seq_t, double>::iterator  iter = this->startProbs.begin(); iter != this->startProbs.end(); iter++)
            iter->second /= totalStarts;
    
    if (totalStops > 0)
        for (map<CharNumConverter::seq_t, double>::iterator iter = this->stopProbs.begin(); iter != this->stopProbs.end(); iter++)
            iter->second /= totalStops;
    
}

// Estimate parameters for gene coding model
void GMS2Trainer::estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels, NumSequence::size_type scSize, const vector<bool> &use) {
    
    // check if all labels should be used
    bool useAll = true;
    if (use.size() > 0) {
        useAll = false;
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
//    PeriodicCounts counts (codingOrder, 3, *this->alphabet);
    CodingCounts counts (params.orderCoding, 3, *this->alphabet, *this->numGeneticCode);
    
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
        if (scSize > 0 && params.marginStartContext < -3) {
            size_t scSizeInCoding = abs(params.marginStartContext) - 3;          // length of start context that overlaps with CDS
            if (reverseComplement)
                counts.decount(sequence.begin()+right-scSizeInCoding+1, sequence.begin()+right+1, reverseComplement);
            else
                counts.decount(sequence.begin()+left, sequence.begin()+left+scSizeInCoding);
        }
        
    }
    
    // convert counts to probabilities
//    coding = new PeriodicMarkov(codingOrder, 3, *this->alphabet);
    coding = new CodingMarkov(params.orderCoding, 3, *this->alphabet, *this->numGeneticCode);
    coding->construct(&counts, params.pcounts);
    
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
    NonCodingCounts counts(params.orderNonCoding, *this->alphabet);
    
    // train non-coding on labels
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
    noncoding = new UniformMarkov(params.orderNonCoding, *this->alphabet);
    noncoding->construct(&counts, params.pcounts);
}

// Estimate parameters for start-context model
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
    
    NonUniformCounts counts(params.orderStartContext, params.lengthStartContext, *this->alphabet);
    
    // get counts for start context model
    size_t n = 0;
    for (vector<Label*>::const_iterator iter = labels.begin(); iter != labels.end(); iter++) {
        if (!useAll && !use[n++])
            continue;       // skip unwanted genes
        
        if (*iter == NULL)
            throw invalid_argument("Label can't be NULL");
        
        // skip sequence if it doesn't have enough nucleotides for start context
        size_t ntBeforeStart = (params.lengthStartContext + params.marginStartContext);
        if ((*iter)->strand == Label::POS && (*iter)->left < ntBeforeStart) {
            continue;
        }
        else if ((*iter)->strand == Label::NEG && (*iter)->right > sequence.size() - ntBeforeStart) {
            continue;
        }
        if ((*iter)->right - (*iter)->left + 1 - 6 < params.lengthStartContext)        // get length of sequence, and -6 to account for start/stop codons
            continue;
        
        size_t left;
        
        if ((*iter)->strand == Label::POS)
//            left = (*iter)->left + 3;
//            left = (*iter)->left-3;    // left = 3;   (3 + -(18-15)) = 0
            left = (*iter)->left - ((int)params.lengthStartContext + params.marginStartContext);
        else
//            left = (*iter)->right - 2 - startContextLength;
//            left = (*iter)->right+3 - startContextLength + 1;      // right = 20:    20 - 15
            left = (*iter)->right + params.marginStartContext + 1;
        
        bool reverseComplement = (*iter)->strand == Label::NEG;
        counts.count(sequence.begin() + left, sequence.begin() + left + params.lengthStartContext, reverseComplement);
        
//        NumSequence s = sequence.subseq(left, startContextLength);
//        if (reverseComplement)
//            s.reverseComplement(cnc);
//        
//        cout << cnc.convert(s.begin(), s.end()) << endl;
        
    }
    
    // convert counts to probabilities
    if (!params.runMotifSearch) {
        startContext = new NonUniformMarkov(params.orderStartContext, params.lengthStartContext, *this->alphabet);
        startContext->construct(&counts, params.pcounts);
    }
    else {
        startContextRBS = new NonUniformMarkov(params.orderStartContext, params.lengthStartContext, *this->alphabet);
        startContextRBS->construct(&counts, params.pcounts);
    }
}


// Estimate parameters for motif models (based on genome group)
void GMS2Trainer::estimateParametersMotifModel(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use) {
    
    // check if all labels should be used
    if (use.size() > 0) {
        if (use.size() != labels.size())
            throw invalid_argument("Labels and Use vector should have the same length");
    }
    
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

    
    if (params.genomeGroup == ProkGeneStartModel::A) {
        this->genomeType = "group-a";
        estimateParametersMotifModel_GroupA(sequence, useLabels);
    }
    else if (params.genomeGroup == ProkGeneStartModel::A2) {
        this->genomeType = "group-a2";
        estimateParametersMotifModel_groupA2(sequence, useLabels);
    }
    else if (params.genomeGroup == ProkGeneStartModel::B) {
        this->genomeType = "group-b";
        estimateParametersMotifModel_GroupB(sequence, useLabels);
    }
    else if (params.genomeGroup == ProkGeneStartModel::C) {
        this->genomeType = "group-c";
        estimateParametersMotifModel_GroupC(sequence, useLabels);
    }
    else if (params.genomeGroup == ProkGeneStartModel::C2) {
        this->genomeType = "group-c2";
        estimateParametersMotifModel_GroupC2(sequence,useLabels);
    }
    else if (params.genomeGroup == ProkGeneStartModel::D) {
        this->genomeType = "group-d";
        estimateParametersMotifModel_GroupD(sequence, useLabels);
    }
    else {
        this->genomeType = "group-e";
        estimateParametersMotifModel_GroupE(sequence, useLabels);
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
    estimateParamtersCoding(sequence, labels, params.lengthStartContext, useCoding);
    
    estimateParametersStartStopCodons(sequence, labels);
    
    // estimate parameters for noncoding model
    useNonCoding = useCoding;                           // for now, assume
    estimateParamtersNonCoding(sequence, labels);
    
    // estimate parameters for start context
    useStartContext = useCoding;
    estimateParametersStartContext(sequence, labels);
    
    // estimate parameters for motif models
    if (params.runMotifSearch) {
        useMotif = useCoding;
        estimateParametersMotifModel(sequence, labels, useMotif);
    }
}



// Deallocate memory for all models
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

// Destructor
GMS2Trainer::~GMS2Trainer() {
    deallocAllModels();
}

// convert codon frequency models to model file output
void codonFrequencyToMod(const map<CharNumConverter::seq_t, double> &codons, const CharNumConverter &cnc, vector<pair<string, string> > &toMod) {
    
    for (map<CharNumConverter::seq_t, double>::const_iterator iter = codons.begin(); iter != codons.end(); iter++) {
        stringstream ssm; ssm << fixed;
        
        string cod = cnc.convert(iter->first.begin(), iter->first.end());
        ssm << iter->second;
        
        toMod.push_back(pair<string,string>(cod, ssm.str()));
    }
    
}

// Convert models and parameters to string in model file format
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
        toMod.push_back(mpair("SC_MARGIN", boost::lexical_cast<string>(params.marginStartContext)));
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
            toMod.push_back(mpair("SC_RBS_MARGIN", boost::lexical_cast<string>(params.marginStartContext)));
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
            toMod.push_back(mpair("SC_PROMOTER_MARGIN", boost::lexical_cast<string>(params.marginStartContext)));
            toMod.push_back(mpair("SC_PROMOTER_MAT", startContextPromoter->toString()));
        }
        else if (startContextRBS != NULL) {                 // copy promoter start context from RBS start context
            toMod.push_back(mpair("SC_PROMOTER", "1"));
            toMod.push_back(mpair("SC_PROMOTER_ORDER", boost::lexical_cast<string>(startContextRBS->getOrder())));
            toMod.push_back(mpair("SC_PROMOTER_WIDTH", boost::lexical_cast<string>(startContextRBS->getLength())));
            toMod.push_back(mpair("SC_PROMOTER_MARGIN", boost::lexical_cast<string>(params.marginStartContext)));
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




// Select labels to be used for estimation of coding model parameters
void GMS2Trainer::selectLabelsForCodingParameters(const vector<Label*> &labels, vector<bool> &useCoding) const {
    
    if (useCoding.size() != labels.size())
        throw invalid_argument("Labels and useCoding vectors should have the same size");
    
    
    for (size_t n = 0; n < labels.size(); n++) {
        if (labels[n] == NULL)
            throw invalid_argument("Label cannot be null");
        
        // "remove" short genes
        if (labels[n]->right - labels[n]->left + 1 <= params.minimumGeneLengthTraining)
            useCoding[n] = false;
        
        // remove all atypical genes, and keep all native (including short)
        if (params.onlyTrainOnNativeGenes) {
            useCoding[n] = (labels[n]->geneClass.find("native") != string::npos);
        }
    }
}






/*************************\
 *      Start Models     *
\*************************/




void GMS2Trainer::estimateParametersMotifModel_groupA2(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, params.fgioDistanceThresh, params.igioDistanceThresh, operonStatuses);
    
    vector<Label*> labelsFGIO, labelsIG, labelsUNK;
    LabelsParser::splitBasedOnPartition(labels, operonStatuses, labelsFGIO, labelsIG, labelsUNK);
    
    vector<NumSequence> upstreamsRBS;
    vector<NumSequence> upstreamsPromoter;
    
    
    
    // match FGIO to 16S tail
    vector<NumSequence> upstreamsFGIOForMatching, upstreamsFGIOForPromoter;
    SequenceParser::extractUpstreamSequences(sequence, labelsFGIO, cnc, params.groupA_upstreamLengthPromoter, upstreamsFGIOForPromoter, true);
    
    upstreamsFGIOForMatching.resize(upstreamsFGIOForPromoter.size());
    for (size_t n = 0; n < upstreamsFGIOForPromoter.size(); n++) {
        upstreamsFGIOForMatching[n] = upstreamsFGIOForPromoter[n].subseq(params.groupA_upstreamLengthPromoter - params.groupA_upstreamLengthRBS, params.groupA_upstreamLengthRBS);
        assert(upstreamsFGIOForMatching[n].size() == params.groupA_upstreamLengthRBS);
    }
    
    assert(upstreamsFGIOForPromoter.size() == upstreamsFGIOForMatching.size());
    
    Sequence strMatchSeq (params.groupA_extendedSD);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (params.groupA_allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    size_t skipFromStart = 3;
    for (size_t n = 0; n < upstreamsFGIOForMatching.size(); n++) {
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreamsFGIOForMatching[n], positionsOfMatches, substitutions);
        
        // keep track of nonmatches
        if (match.size() < params.groupA_minMatchToExtendedSD)
            upstreamsPromoter.push_back(upstreamsFGIOForPromoter[n]);
        else
            upstreamsRBS.push_back(upstreamsFGIOForMatching[n]);
    }
    
    
    vector<NumSequence> upstreamsIG;
    SequenceParser::extractUpstreamSequences(sequence, labelsIG, cnc, params.groupA_upstreamLengthRBS, upstreamsIG);
    for (size_t n = 0; n < upstreamsIG.size(); n++) {
        upstreamsRBS.push_back(upstreamsIG[n]);
    }
    
    this->numLeaderless = upstreamsPromoter.size();
    this->numFGIO = upstreamsFGIOForMatching.size();
    
    OptionsMFinder optionMFinderPromoter (*this->params.optionsMFinder);
    optionMFinderPromoter.width = params.groupA_widthPromoter;
    
    // take first
    if (cutPromTrainSeqs) {
        for (size_t n = 0; n < upstreamsPromoter.size(); n++) {
            upstreamsPromoter[n] = upstreamsPromoter[n].subseq(0, params.groupA_upstreamLengthPromoter - 15);
        }
    }
    
    runMotifFinder(upstreamsPromoter, optionMFinderPromoter, *this->alphabet, params.groupA_upstreamLengthPromoter, this->promoter, this->promoterSpacer);
    runMotifFinder(upstreamsRBS, *this->params.optionsMFinder, *this->alphabet, params.groupA_upstreamLengthRBS, this->rbs, this->rbsSpacer);
    
    //    // shift probabilities
    //    vector<double> extendedProbs (promoterSpacer->size()+skipFromStart, 0);
    //    for (size_t n = 0; n < promoterSpacer->size(); n++) {
    //        extendedProbs[n+skipFromStart] = (*promoterSpacer)[n];
    //    }
    //
    //    delete promoterSpacer;
    //    promoterSpacer = new UnivariatePDF(extendedProbs);
    
    // shift probabilities
    if (cutPromTrainSeqs) {
        vector<double> extendedProbs (promoterSpacer->size()+15, 0);
        for (size_t n = 0; n < promoterSpacer->size(); n++) {
            extendedProbs[n+15] = (*promoterSpacer)[n];
        }
        
        delete promoterSpacer;
        promoterSpacer = new UnivariatePDF(extendedProbs);
    }
    
}


void GMS2Trainer::estimateParametersMotifModel_GroupA(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, params.fgioDistanceThresh, params.igioDistanceThresh, operonStatuses);
    
    vector<Label*> labelsFGIO, labelsIG, labelsUNK;
    LabelsParser::splitBasedOnPartition(labels, operonStatuses, labelsFGIO, labelsIG, labelsUNK);
    
    vector<NumSequence> upstreamsFGIO;
    SequenceParser::extractUpstreamSequences(sequence, labelsFGIO, cnc, this->params.groupA_upstreamLengthPromoter, upstreamsFGIO, true);
    
    // take first
    if (cutPromTrainSeqs) {
        for (size_t n = 0; n < upstreamsFGIO.size(); n++) {
            upstreamsFGIO[n] = upstreamsFGIO[n].subseq(0, this->params.groupA_upstreamLengthPromoter - 15);
        }
    }
    
    vector<NumSequence> upstreamsIG;
    SequenceParser::extractUpstreamSequences(sequence, labelsIG, cnc, this->params.groupA_upstreamLengthRBS, upstreamsIG);
    
    //    void runMotifFinder(const vector<NumSequence> &sequencesRaw, OptionsMFinder &optionsMFinder, size_t upstreamLength, NonUniformMarkov* motifMarkov, UnivariatePDF* motifSpacer) {
    //
    this->numLeaderless = upstreamsFGIO.size();
    this->numFGIO = upstreamsFGIO.size();
    
    OptionsMFinder optionMFinderFGIO (*this->params.optionsMFinder);
    optionMFinderFGIO.width = params.groupA_widthPromoter;
    
    runMotifFinder(upstreamsFGIO, optionMFinderFGIO, *this->alphabet, this->params.groupA_upstreamLengthPromoter, this->promoter, this->promoterSpacer);
    runMotifFinder(upstreamsIG, *this->params.optionsMFinder, *this->alphabet, this->params.groupA_upstreamLengthRBS, this->rbs, this->rbsSpacer);
    
    
    // shift probabilities
    if (cutPromTrainSeqs) {
        vector<double> extendedProbs (promoterSpacer->size()+15, 0);
        for (size_t n = 0; n < promoterSpacer->size(); n++) {
            extendedProbs[n+15] = (*promoterSpacer)[n];
        }
        
        delete promoterSpacer;
        promoterSpacer = new UnivariatePDF(extendedProbs);
    }
    
}


void GMS2Trainer::estimateParametersMotifModel_GroupB(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, params.fgioDistanceThresh, params.igioDistanceThresh, operonStatuses);
    
    vector<Label*> labelsFGIO, labelsIG, labelsUNK;
    LabelsParser::splitBasedOnPartition(labels, operonStatuses, labelsFGIO, labelsIG, labelsUNK);
    
    vector<NumSequence> upstreamsRBS;
    vector<NumSequence> upstreamsPromoter;
    
    // match FGIO to 16S tail
    vector<NumSequence> upstreamsFGIO;
    SequenceParser::extractUpstreamSequences(sequence, labelsFGIO, cnc, params.groupB_upstreamLengthPromoter, upstreamsFGIO);
    
    Sequence strMatchSeq (params.groupB_extendedSD);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (params.groupB_allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    size_t skipFromStart = 3;
    for (size_t n = 0; n < upstreamsFGIO.size(); n++) {
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreamsFGIO[n], positionsOfMatches, substitutions);
        
        // keep track of nonmatches
        if (match.size() < params.groupB_minMatchToExtendedSD)
            upstreamsPromoter.push_back(upstreamsFGIO[n].subseq(0, upstreamsFGIO[n].size() - skipFromStart));
        else
            upstreamsRBS.push_back(upstreamsFGIO[n]);
    }
    
    
    vector<NumSequence> upstreamsIG;
    SequenceParser::extractUpstreamSequences(sequence, labelsIG, cnc, params.groupB_upstreamLengthRBS, upstreamsIG);
    for (size_t n = 0; n < upstreamsIG.size(); n++) {
        upstreamsRBS.push_back(upstreamsIG[n]);
    }
    
    this->numLeaderless = upstreamsPromoter.size();
    this->numFGIO = upstreamsFGIO.size();
    
    
    runMotifFinder(upstreamsPromoter, *this->params.optionsMFinder, *this->alphabet, params.groupB_upstreamLengthPromoter-skipFromStart, this->promoter, this->promoterSpacer);
    runMotifFinder(upstreamsRBS, *this->params.optionsMFinder, *this->alphabet, params.groupB_upstreamLengthRBS, this->rbs, this->rbsSpacer);
    
    // shift probabilities
    vector<double> extendedProbs (promoterSpacer->size()+skipFromStart, 0);
    for (size_t n = 0; n < promoterSpacer->size(); n++) {
        extendedProbs[n+skipFromStart] = (*promoterSpacer)[n];
    }
    
    delete promoterSpacer;
    promoterSpacer = new UnivariatePDF(extendedProbs);
    
    
}


void GMS2Trainer::estimateParametersMotifModel_GroupC(const NumSequence &sequence, const vector<Label *> &labels) {
    this->estimateParametersMotifModel_GroupD(sequence, labels);
}


void GMS2Trainer::estimateParametersMotifModel_GroupC2(const NumSequence &sequence, const vector<Label *> &labels) {
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, params.fgioDistanceThresh, params.igioDistanceThresh, operonStatuses);
    
    // scope since variables won't (shouldn't) be used outside
    {
    vector<Label*> labelsFGIO, labelsIGIO, labelsAMBIG;
    LabelsParser::splitBasedOnPartition(labels, operonStatuses, labelsFGIO, labelsIGIO, labelsAMBIG);
    
    this->numFGIO = labelsFGIO.size();
    }
    
    // extract upstream of each label
    vector<NumSequence> upstreamsRaw;
    SequenceParser::extractUpstreamSequences(sequence, labels, *alphabet->getCNC(), params.groupC2_upstreamLengthSDRBS, upstreamsRaw, false, params.minimumGeneLengthTraining);
    
    vector<NumSequence> upstreams;
    for (size_t n = 0; n < upstreamsRaw.size(); n++) {
        if (!upstreamsRaw[n].containsInvalid(*this->alphabet))
            upstreams.push_back(upstreamsRaw[n].subseq(0, upstreamsRaw[n].size() - params.groupC2_upstreamRegion3Prime));
    }
    
    
    vector<NumSequence> upstreamsSD, upstreamsNonSD;
    
    // match against SD
    Sequence strMatchSeq (params.groupC2_extendedSD);
    NumSequence matchSeq (strMatchSeq, *this->cnc);
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (params.groupB_allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (this->cnc->convert('A'), this->cnc->convert('G')));
    
    size_t skipFromStart = 3;
    for (size_t n = 0; n < upstreams.size(); n++) {
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreams[n], positionsOfMatches, substitutions);
        
        // keep track of nonmatches
        if (match.size() < params.groupC2_minMatchToExtendedSD)
            upstreamsNonSD.push_back(upstreams[n].subseq(0, upstreams[n].size() - skipFromStart));
        else
            upstreamsSD.push_back(upstreams[n]);
    }
    
    
    runMotifFinder(upstreamsSD, *this->params.optionsMFinder, *this->alphabet, params.groupC2_upstreamLengthSDRBS-skipFromStart, this->promoter, this->promoterSpacer);
    runMotifFinder(upstreamsNonSD, *this->params.optionsMFinder, *this->alphabet, params.groupC2_upstreamLengthNonSDRBS, this->rbs, this->rbsSpacer);
    
    // shift probabilities
    vector<double> extendedProbs (promoterSpacer->size()+skipFromStart, 0);
    for (size_t n = 0; n < promoterSpacer->size(); n++) {
        extendedProbs[n+skipFromStart] = (*promoterSpacer)[n];
    }
    
    delete promoterSpacer;
    promoterSpacer = new UnivariatePDF(extendedProbs);
    
    this->numLeaderless = upstreamsSD.size();
    
}


void GMS2Trainer::estimateParametersMotifModel_GroupD(const NumSequence &sequence, const vector<Label *> &labels) {
    
    MotifFinder::Builder b;
    MotifFinder mfinder = b.build(*this->params.optionsMFinder);
    
    
    // split labels into sets based on operon status
    vector<LabelsParser::operon_status_t> operonStatuses;
    LabelsParser::partitionBasedOnOperonStatus(labels, params.fgioDistanceThresh, params.igioDistanceThresh, operonStatuses);
    
    vector<Label*> labelsFGIO, labelsIGIO, labelsAMBIG;
    LabelsParser::splitBasedOnPartition(labels, operonStatuses, labelsFGIO, labelsIGIO, labelsAMBIG);
    
    this->numFGIO = labelsFGIO.size();
    
    // extract upstream of each label
    vector<NumSequence> upstreamsRaw;
    SequenceParser::extractUpstreamSequences(sequence, labels, *alphabet->getCNC(), params.groupD_upstreamLengthRBS, upstreamsRaw, false, params.minimumGeneLengthTraining);
    
    vector<NumSequence> upstreams;
    for (size_t n = 0; n < upstreamsRaw.size(); n++) {
        if (!upstreamsRaw[n].containsInvalid(*this->alphabet))
            upstreams.push_back(upstreamsRaw[n].subseq(0, upstreamsRaw[n].size() - params.groupC_upstreamRegion3Prime));
    }
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(upstreams, positions);
    
    // build RBS model
    NonUniformCounts rbsCounts(params.optionsMFinder->motifOrder, params.optionsMFinder->width, *this->alphabet);
    for (size_t n = 0; n < upstreams.size(); n++) {
        rbsCounts.count(upstreams[n].begin()+positions[n], upstreams[n].begin()+positions[n]+params.optionsMFinder->width);
    }
    
    rbs = new NonUniformMarkov(params.optionsMFinder->motifOrder, params.optionsMFinder->width, *this->alphabet);
    rbs->construct(&rbsCounts, params.optionsMFinder->pcounts);
    
    // build spacer distribution
    // build histogram from positions
    vector<double> positionCounts (params.groupD_upstreamLengthRBS - params.optionsMFinder->width+1, 0);
    for (size_t n = 0; n < positions.size(); n++) {
        // FIXME account for LEFT alignment
        // below is only for right
        positionCounts[params.groupD_upstreamLengthRBS - params.optionsMFinder->width - positions[n]]++;        // increment position
    }
    
    rbsSpacer = new UnivariatePDF(positionCounts, false, params.pcounts);
    
}


void GMS2Trainer::estimateParametersMotifModel_GroupE(const NumSequence &sequence, const vector<Label *> &labels) {
    
    AlphabetDNA alph;
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    
    Sequence strMatchSeq (params.groupE_extendedSD);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    pair<NumSequence::size_type, NumSequence::size_type> positionsOfMatches;
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (params.groupE_allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    
    // extract upstream for every sequence and match it to 16S tail
    vector<NumSequence> upstreams (labels.size());
    SequenceParser::extractUpstreamSequences(sequence, labels, cnc, params.groupE_upstreamLengthRBS, upstreams);
    size_t skipFromStart = 0;
    
    vector<Label*> labelsSig;
    vector<Label*> labelsRBS;
    
    for (size_t n = 0; n < upstreams.size(); n++) {
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, upstreams[n], positionsOfMatches, substitutions);
        
        // keep track of nonmatches
        if (match.size() < params.groupE_minMatchToExtendedSD)
            labelsSig.push_back(labels[n]);
        else
            labelsRBS.push_back(labels[n]);
    }
    
    // for all non-Sig sequences, append "N" to mask Sig sequences
    string Ns = "";
    size_t numNs = 0;
    if ( this->params.groupE_lengthUpstreamSignature >= this->params.lengthStartContext &&
        this->params.groupE_lengthUpstreamSignature - this->params.lengthStartContext >= this->params.groupE_orderUpstreamSignature)
        numNs = (this->params.groupE_lengthUpstreamSignature - this->params.lengthStartContext) - this->params.groupE_orderUpstreamSignature;
    
    for (size_t n = 0; n < numNs; n++)
        Ns += "N";
    
    NumSequence numSeqNs (Sequence(Ns), cnc);       // numeric sequence of N's
    
    NonUniformCounts counts(params.groupE_orderUpstreamSignature, params.groupE_lengthUpstreamSignature, *this->alphabet);
    
    vector<NumSequence> contextsRBS;
    long long posRelToStart = - (params.lengthStartContext + params.marginStartContext + params.groupE_orderUpstreamSignature);
    SequenceParser::extractStartContextSequences(sequence, labelsRBS, cnc, posRelToStart, params.lengthStartContext + this->params.groupE_orderUpstreamSignature, contextsRBS);
    
    for (size_t n = 0; n < contextsRBS.size(); n++) {
        NumSequence withNs = numSeqNs + contextsRBS[n];     // append N's
        counts.count(withNs.begin(), withNs.end());
        //        cout << cnc.convert(withNs.begin(), withNs.end()) << endl;
    }
    
    // add Sig sequences
    vector<NumSequence> contextsSig;
    SequenceParser::extractStartContextSequences(sequence, labelsSig, cnc, -( (int)params.groupE_lengthUpstreamSignature + params.marginStartContext), params.groupE_lengthUpstreamSignature, contextsSig);
    
    for (size_t n = 0; n < contextsSig.size(); n++) {
        counts.count(contextsSig[n].begin(), contextsSig[n].end());
        //        cout << cnc.convert(contextsSig[n].begin(), contextsSig[n].end()) << endl;
    }
    
    // for sequences without motifs
    startContext = new NonUniformMarkov(params.groupE_orderUpstreamSignature, params.groupE_lengthUpstreamSignature, *this->alphabet);
    startContext->construct(&counts, params.pcounts);
    
    // run motif search for RBS
    vector<NumSequence> upstreamsRBS;
    SequenceParser::extractUpstreamSequences(sequence, labelsRBS, cnc, this->params.groupE_upstreamLengthRBS, upstreamsRBS);
    runMotifFinder(upstreamsRBS, *this->params.optionsMFinder, *this->alphabet, this->params.groupE_upstreamLengthRBS, this->rbs, this->rbsSpacer);
    
    
    
}




















// REMOVE
//vector<GeneStat> separateLabelsViaOperonStatus(const vector<Label*> &labels, unsigned thresh, unsigned threshNFIO) {
//
//
//    vector<GeneStat> result (labels.size());
//
//    // for every label
//    for (size_t n = 0; n < labels.size(); n++) {
//
//        // get label
//        const Label* lab = labels[n];
//
//        // if on positive strand
//        if (lab->strand == Label::POS) {
//
//            size_t start = lab->left;          // start location
//
//            // if no gene before it, set as first in operon
//            if (n == 0) {
//                result[n] = FIRST_OP;
//            }
//            // otherwise, check to see if stop of previous gene is nearby
//            else {
//
//                const Label* prevLab = labels[n-1];         // previous label
//
//                if (prevLab->strand == Label::POS) {   // should be on positive strand
//
//                    size_t prevStop = prevLab->right;
//
//                    // if too far away, then current gene is first in op
//                    if (prevStop < start - thresh)
//                        result[n] = FIRST_OP;
//                    // if threshNFIO is set and gene is too close
//                    else if (threshNFIO != 0) {
//                        if (prevStop > start - threshNFIO)
//                            result[n] = NOT_FIRST_OP;
//                        else
//                            result[n] = IGNORE;
//                    }
//                    // otherwise, current gene belongs to same operon as previous
//                    else
//                        result[n] = NOT_FIRST_OP;
//
//                }
//                // otherwise it's on negative strand, so first in operon
//                else
//                    result[n] = FIRST_OP;
//            }
//        }
//
//        // if on negative strand
//        else if (lab->strand == Label::NEG) {
//
//            size_t start = lab->right;        // start location
//
//            // if no gene after it, then set as first in operon
//            if (n == labels.size()-1) {
//                result[n] = FIRST_OP;
//            }
//            // otherwise, check to see if stop of previous gene (i.e. to the right) is nearby
//            else {
//                const Label* prevLab = labels[n+1];     // "next" label
//
//                if (prevLab->strand == Label::NEG) {      // should be on negative strand
//
//                    size_t prevStop = prevLab->left;
//
//                    // if too far away, then current gene is first in op
//                    if (prevStop > start + thresh)
//                        result[n] = FIRST_OP;
//                    // if threshNFIO is set and gene is too close
//                    else if (threshNFIO != 0) {
//                        if (prevStop < start + threshNFIO)
//                            result[n] = NOT_FIRST_OP;
//                        else
//                            result[n] = IGNORE;
//                    }
//                    // otherwise, current gene belongs to same operon as previous
//                    else
//                        result[n] = NOT_FIRST_OP;
//                }
//                // oterhwise, previous gene is on different strand, so current is first in operon
//                else
//                    result[n] = FIRST_OP;
//            }
//        }
//
//    }
//
//
//    return result;
//
//}

