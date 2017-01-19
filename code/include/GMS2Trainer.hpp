//
//  GMS2Trainer.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef GMS2Trainer_hpp
#define GMS2Trainer_hpp

#include <stdio.h>
#include <string>
#include <map>

#include "Label.hpp"

#include "UnivariatePDF.hpp"
#include "NumSequence.hpp"
#include "NumGeneticCode.hpp"
#include "UniformMarkov.hpp"
#include "OptionsMFinder.hpp"
#include "PeriodicMarkov.hpp"
#include "NonUniformMarkov.hpp"
#include "ProkGeneStartModel.hpp"
#include "OptionsGMS2Training.hpp"

using std::map;
using std::string;
using std::pair;

namespace gmsuite {
    
    
    /**
     * @class GMS2Trainer
     * @brief Train model parameters for GMS2
     *
     * This class trains model parameters used by GMS2, such as coding model, non-coding model, etc...
     */
    class GMS2Trainer {
        
    public:
        
        typedef ProkGeneStartModel::genome_class_t genome_class_t;
        
        /**
         * Default constructor:
         */
        GMS2Trainer();
        
        GMS2Trainer(unsigned pcounts,
                    unsigned codingOrder,
                    unsigned noncodingOrder,
                    unsigned startContextOrder,
                    NumSequence::size_type upstreamLength,
                    NumSequence::size_type startContextLength,
                    genome_class_t genomeClass,
                    const OptionsMFinder &optionsMFinder,
                    const NumAlphabetDNA &alph,
                    const NumSequence::size_type MIN_GENE_LEN,
                    const NumGeneticCode &numGeneticCode,
                    int scMargin,
                    bool trainOnNative,
                    bool runMotifSearch=true,
                    NumSequence::size_type upstrFGIO = 40,
                    unsigned widthArchaeaPromoter = 12,
                    string matchTo = "TAAGGAGGTGA",
                    bool allowAGSubstitution = true,
                    unsigned matchThresh = 4);
        
        ~GMS2Trainer();
        
        /**
         * Train
         */
        void estimateParameters(const NumSequence &sequence, const vector<Label*> &labels);
        
        void estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels, NumSequence::size_type scSize = 0, const vector<bool> &use = vector<bool>());
        void estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        void estimateParametersStartContext(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        void estimateParametersMotifModel(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        void estimateParametersStartStopCodons(const NumSequence &sequence, const vector<Label*> &labels, const vector<bool> &use = vector<bool>());
        
        void estimateParametersMotifModel_Promoter(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        
        void estimateParametersMotifModel_Tuberculosis(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        void estimateParametersMotifModel_Synechocystis(const NumSequence &sequence, const vector<Label*> &labels, const vector<bool> &use = vector<bool>());
        
        void estimateParametersMotifModel_Promoter_DEPRECATED(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        
        // parameters
        unsigned pcounts;
        unsigned codingOrder;
        unsigned noncodingOrder;
        unsigned startContextOrder;
        NumSequence::size_type upstreamLength;          // for genome
        NumSequence::size_type startContextLength;
        genome_class_t genomeClass;
        const OptionsMFinder* optionsMFinder;
        const NumAlphabetDNA *alphabet;
        const NumGeneticCode *numGeneticCode;
        NumSequence::size_type MIN_GENE_LEN;            // minimum gene length
        NumSequence::size_type MIN_UPSTR_LEN_FGIO;           // minimum upstream length for first-genes-in-operon
        NumSequence::size_type UPSTR_LEN_NFGIO;         // deprecated
        NumSequence::size_type UPSTR_LEN_IG;    
        NumSequence::size_type UPSTR_LEN_FGIO;
        unsigned FGIO_DIST_THRESH;
        unsigned NFGIO_DIST_THRES;
        int scMargin;
        bool trainOnNative;
        bool runMotifSearch;
        unsigned widthArchaeaPromoter;
        
        string matchTo;
        bool allowAGSubstitution;
        unsigned matchThresh;
        
        NumSequence::size_type upstreamSignatureLength;
        unsigned upstreamSignatureOrder;
        
        // public variables for models
//        NonUniformMarkov *motif;
        UniformMarkov *noncoding;
        PeriodicMarkov *coding;
        NonUniformMarkov *startContext;
        
        // start models
        NonUniformMarkov *rbs;
        NonUniformMarkov *promoter;
        NonUniformMarkov *upstreamSignature;
        UnivariatePDF *rbsSpacer;
        UnivariatePDF *promoterSpacer;
        
        map<CharNumConverter::seq_t, double> startProbs;
        map<CharNumConverter::seq_t, double> stopProbs;
        
        
        
        void toModFile(vector<pair<string, string> > &toMod, const OptionsGMS2Training &options) const;
        
    private:
        void deallocAllModels();
        
        
        void selectLabelsForCodingParameters(const vector<Label*> &labels, vector<bool> &useCoding) const;
    };
}

#endif /* GMS2Trainer_hpp */
