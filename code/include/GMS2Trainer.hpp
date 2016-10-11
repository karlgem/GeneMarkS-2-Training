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
#include "UniformMarkov.hpp"
#include "OptionsMFinder.hpp"
#include "PeriodicMarkov.hpp"
#include "NonUniformMarkov.hpp"
#include "ProkGeneStartModel.hpp"

using std::map;
using std::string;

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
                    const NumSequence::size_type MIN_GENE_LEN);
        
        ~GMS2Trainer();
        
        /**
         * Train
         */
        void estimateParameters(const NumSequence &sequence, const vector<Label*> &labels);
        
        void estimateParamtersCoding(const NumSequence &sequence, const vector<Label *> &labels, NumSequence::size_type scSize = 0, const vector<bool> &use = vector<bool>());
        void estimateParamtersNonCoding(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        void estimateParametersStartContext(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        void estimateParametersMotifModel(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        
        
        void estimateParametersMotifModel_Promoter(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use = vector<bool>());
        
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
        NumSequence::size_type MIN_GENE_LEN;            // minimum gene length
        NumSequence::size_type MIN_UPSTR_LEN_FGIO;           // minimum upstream length for first-genes-in-operon
        NumSequence::size_type UPSTR_LEN_NFGIO;
        unsigned FGIO_DIST_THRESH;
        unsigned NFGIO_DIST_THRES;
        
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
        
        
        
        void toModFile(map<string, string> &toMod) const;
        
    private:
        void deallocAllModels();
        
        
        void selectLabelsForCodingParameters(const vector<Label*> &labels, vector<bool> &useCoding) const;
    };
}

#endif /* GMS2Trainer_hpp */
