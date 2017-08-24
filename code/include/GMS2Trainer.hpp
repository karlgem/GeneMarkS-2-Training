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
        
        typedef ProkGeneStartModel::genome_group_t genome_group_t;
        
        class Builder;                  // used for easy building of the many GMS2Trainer parameters
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
                    genome_group_t genomeClass,
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
                    unsigned matchThresh = 4,
                    NumSequence::size_type upstreamSignatureLength = 35,
                    unsigned upstreamSignatureOrder = 2,
                    bool trainNonCodingOnFullGenome=false,
                    unsigned FGIO_DIST_THRESH = 25,
                    bool cutPromTrainSeqs = false);
        
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
        void estimateParametersMotifModel_groupA2(const NumSequence &sequence, const vector<Label *> &labels, const vector<bool> &use);
            
        // parameters
        unsigned pcounts;
        unsigned codingOrder;
        unsigned noncodingOrder;
        unsigned startContextOrder;
        NumSequence::size_type upstreamLength;          // for genome
        NumSequence::size_type startContextLength;
        genome_group_t genomeClass;
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
        bool trainNonCodingOnFullGenome;
        size_t numLeaderless;
        size_t numFGIO;
        string genomeType;
        
        string matchTo;
        bool allowAGSubstitution;
        unsigned matchThresh;
        
        NumSequence::size_type upstreamSignatureLength;
        unsigned upstreamSignatureOrder;
        
        // public variables for models
//        NonUniformMarkov *motif;
        UniformMarkov *noncoding;
        PeriodicMarkov *coding;
        NonUniformMarkov *startContextRBS;
        NonUniformMarkov *startContextPromoter;
        NonUniformMarkov *startContext;             // used for synechocystis-type or when no motif is allowed
        
        // start models
        NonUniformMarkov *rbs;
        NonUniformMarkov *promoter;
        NonUniformMarkov *upstreamSignature;
        UnivariatePDF *rbsSpacer;
        UnivariatePDF *promoterSpacer;
        
        map<CharNumConverter::seq_t, double> startProbs;
        map<CharNumConverter::seq_t, double> stopProbs;
        
        bool cutPromTrainSeqs;      // when set, promoters are trained on fragment of total sequence
        
        
        
        void toModFile(vector<pair<string, string> > &toMod, const OptionsGMS2Training &options) const;
        
    private:
        void deallocAllModels();
        
        
        void selectLabelsForCodingParameters(const vector<Label*> &labels, vector<bool> &useCoding) const;
    };
    
    
    // Class for building a GMS2Trainer (for advanced users)
    class GMS2Trainer::Builder {
        
        typedef ProkGeneStartModel::genome_group_t genome_group_t;
        
    private:
        // variables needed for the GMS2Trainer class
        
        // Group-A
        unsigned                groupA_widthPromoter                ;
        unsigned                groupA_widthRBS                     ;
        NumSequence::size_type  groupA_upstreamLengthPromoter       ;
        NumSequence::size_type  groupA_upstreamLengthRBS            ;
        double                  groupA_spacerScoreThresh            ;
        NumSequence::size_type  groupA_spacerDistThresh             ;
        NumSequence::size_type  groupA_spacerWindowSize             ;
        
        // Group B
        unsigned                groupB_widthPromoter                ;
        unsigned                groupB_widthRBS                     ;
        NumSequence::size_type  groupB_upstreamLengthPromoter       ;
        NumSequence::size_type  groupB_upstreamLengthRBS            ;
        double                  groupB_spacerScoreThresh            ;
        NumSequence::size_type  groupB_spacerDistThresh             ;
        NumSequence::size_type  groupB_spacerWindowSize             ;
        string                  groupB_extendedSD                   ;
        unsigned                groupB_minMatchToExtendedSD         ;
        bool                    groupB_allowAGSubstitution          ;
        
        // Group-C
        unsigned                groupC_widthRBS                     ;
        NumSequence::size_type  groupC_upstreamLengthRBS            ;
        unsigned                groupC_minMatchRBSPromoter          ;
        unsigned                groupC_minMatchToExtendedSD         ;
        string                  groupC_extendedSD                   ;
        
        // Group-D
        unsigned                groupD_widthRBS                     ;
        NumSequence::size_type  groupD_upstreamLengthRBS            ;
        double                  groupD_percentMatchRBS              ;
        string                  groupD_extendedSD                   ;
        unsigned                groupD_minMatchToExtendedSD         ;
        bool                    groupD_allowAGSubstitution          ;
        
        // Group-E
        unsigned                groupE_widthRBS                     ;
        NumSequence::size_type  groupE_upstreamLengthRBS            ;
        NumSequence::size_type  groupE_lengthUpstreamSignature      ;
        unsigned                groupE_orderUpstreamSignature       ;
        string                  groupE_extendedSD                   ;
        bool                    groupE_allowAGSubstitution          ;
        
        // Coding and Noncoding Models
        unsigned                orderCoding                         ;
        unsigned                orderNonCoding                      ;
        unsigned                orderStartContext                   ;
        NumSequence::size_type  lengthStartContext                  ;
        int                     marginStartContext                  ;
        
        // Misc Variables
        NumSequence::size_type  fgioDistanceThresh                  ;
        NumSequence::size_type  igioDistanceThresh                  ;
        unsigned                pcounts                             ;
        genome_group_t          genomeGroup                         ;
        GeneticCode::gcode_t    gcode                               ;
        NumSequence::size_type  minimumGeneLengthTraining           ;
        bool                    onlyTrainOnNativeGenes              ;
        bool                    runMotifSearch                      ;
        
        
    public:
        
        Builder() {
            // Group-A
            groupA_widthPromoter                = 12                          ;
            groupA_widthRBS                     = 6                           ;
            groupA_upstreamLengthPromoter       = 20                          ;
            groupA_upstreamLengthRBS            = 20                          ;
            groupA_spacerScoreThresh            = 0.1                         ;
            groupA_spacerDistThresh             = 14                          ;
            groupA_spacerWindowSize             = 1                           ;
            
            // Group-B
            groupB_widthPromoter                = 6                           ;
            groupB_widthRBS                     = 6                           ;
            groupB_upstreamLengthPromoter       = 20                          ;
            groupB_upstreamLengthRBS            = 20                          ;
            groupB_spacerScoreThresh            = 0.25                        ;
            groupB_spacerDistThresh             = 14                          ;
            groupB_spacerWindowSize             = 1                           ;
            groupB_extendedSD                   = "TAAGGAGGTGA"               ;
            groupB_minMatchToExtendedSD         = 4                           ;
            groupB_allowAGSubstitution          = true                        ;
            
            // Group-C
            groupC_widthRBS                     = 6                           ;
            groupC_upstreamLengthRBS            = 20                          ;
            groupC_minMatchRBSPromoter          = 3                           ;
            groupC_minMatchToExtendedSD         = 4                           ;
            groupC_extendedSD                   = "TAAGGAGGTGA"               ;
            
            // Group-D
            groupD_widthRBS                     = 6                           ;
            groupD_upstreamLengthRBS            = 20                          ;
            groupD_percentMatchRBS              = 0.5                         ;
            groupD_extendedSD                   = "TAAGGAGGTGA"               ;
            groupD_minMatchToExtendedSD         = 4                           ;
            groupD_allowAGSubstitution          = true                        ;
            
            // Group-E
            groupE_widthRBS                     = 6                           ;
            groupE_upstreamLengthRBS            = 20                          ;
            groupE_lengthUpstreamSignature      = 35                          ;
            groupE_orderUpstreamSignature       = 2                           ;
            groupE_extendedSD                   = "TAAGGAGGTGA"               ;
            groupE_allowAGSubstitution          = true                        ;
            
            // Models
            orderCoding                         = 5                           ;
            orderNonCoding                      = 2                           ;
            orderStartContext                   = 2                           ;
            lengthStartContext                  = 18                          ;
            marginStartContext                  = 15                          ;
            
            // Misc
            fgioDistanceThresh                  = 25                          ;
            igioDistanceThresh                  = 22                          ;
            pcounts                             = 1                           ;
            genomeGroup                         = ProkGeneStartModel::D       ;
            gcode                               = GeneticCode::ELEVEN         ;
            minimumGeneLengthTraining           = 300                         ;
            onlyTrainOnNativeGenes              = false                       ;
            runMotifSearch                      = true                        ;
        }
        
        GMS2Trainer build() {
            return GMS2Trainer ();
        }
        
        GMS2Trainer build(const OptionsGMS2Training &options) {
            // set all parameters
            return build();
        }
        
        // set custom values for MotifFinder creation
        // returns Builder for shorthand inline usage
        Builder& setGroupA_widthPromoter            (const unsigned v)                  { this->groupA_widthPromoter           = v;    return *this; }
        Builder& setGroupA_widthRBS                 (const unsigned v)                  { this->groupA_widthRBS                = v;    return *this; }
        Builder& setGroupA_upstreamLengthPromoter   (const NumSequence::size_type v)    { this->groupA_upstreamLengthPromoter  = v;    return *this; }
        Builder& setGroupA_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupA_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupA_spacerScoreThresh        (const double v)                    { this->groupA_spacerScoreThresh       = v;    return *this; }
        Builder& setGroupA_spacerDistThresh         (const NumSequence::size_type v)    { this->groupA_spacerDistThresh        = v;    return *this; }
        Builder& setGroupA_spacerWindowSize         (const NumSequence::size_type v)    { this->groupA_spacerWindowSize        = v;    return *this; }
        Builder& setGroupB_widthPromoter            (const unsigned v)                  { this->groupB_widthPromoter           = v;    return *this; }
        Builder& setGroupB_widthRBS                 (const unsigned v)                  { this->groupB_widthRBS                = v;    return *this; }
        Builder& setGroupB_upstreamLengthPromoter   (const NumSequence::size_type v)    { this->groupB_upstreamLengthPromoter  = v;    return *this; }
        Builder& setGroupB_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupB_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupB_spacerScoreThresh        (const double v)                    { this->groupB_spacerScoreThresh       = v;    return *this; }
        Builder& setGroupB_spacerDistThresh         (const NumSequence::size_type v)    { this->groupB_spacerDistThresh        = v;    return *this; }
        Builder& setGroupB_spacerWindowSize         (const NumSequence::size_type v)    { this->groupB_spacerWindowSize        = v;    return *this; }
        Builder& setGroupB_extendedSD               (const string v)                    { this->groupB_extendedSD              = v;    return *this; }
        Builder& setGroupB_minMatchToExtendedSD     (const unsigned v)                  { this->groupB_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupB_allowAGSubstitution      (const bool v)                      { this->groupB_allowAGSubstitution     = v;    return *this; }
        Builder& setGroupC_widthRBS                 (const unsigned v)                  { this->groupC_widthRBS                = v;    return *this; }
        Builder& setGroupC_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupC_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupC_minMatchRBSPromoter      (const unsigned v)                  { this->groupC_minMatchRBSPromoter     = v;    return *this; }
        Builder& setGroupC_minMatchToExtendedSD     (const unsigned v)                  { this->groupC_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupC_extendedSD               (const string v)                    { this->groupC_extendedSD              = v;    return *this; }
        Builder& setGroupD_widthRBS                 (const unsigned v)                  { this->groupD_widthRBS                = v;    return *this; }
        Builder& setGroupD_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupD_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupD_percentMatchRBS          (const double v)                    { this->groupD_percentMatchRBS         = v;    return *this; }
        Builder& setGroupD_extendedSD               (const string v)                    { this->groupD_extendedSD              = v;    return *this; }
        Builder& setGroupD_minMatchToExtendedSD     (const unsigned v)                  { this->groupD_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupD_allowAGSubstitution      (const bool v)                      { this->groupD_allowAGSubstitution     = v;    return *this; }
        Builder& setGroupE_widthRBS                 (const unsigned v)                  { this->groupE_widthRBS                = v;    return *this; }
        Builder& setGroupE_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupE_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupE_lengthUpstreamSignature  (const NumSequence::size_type v)    { this->groupE_lengthUpstreamSignature = v;    return *this; }
        Builder& setGroupE_orderUpstreamSignature   (const unsigned v)                  { this->groupE_orderUpstreamSignature  = v;    return *this; }
        Builder& setGroupE_extendedSD               (const string v)                    { this->groupE_extendedSD              = v;    return *this; }
        Builder& setGroupE_allowAGSubstitution      (const bool v)                      { this->groupE_allowAGSubstitution     = v;    return *this; }
        
        
        
        Builder& setOrderCoding                     (const unsigned v)                  {  orderCoding               = v;     return *this; }
        Builder& setOrderNonCoding                  (const unsigned v)                  {  orderNonCoding            = v;     return *this; }
        Builder& setOrderStartContext               (const unsigned v)                  {  orderStartContext         = v;     return *this; }
        Builder& setLengthStartContext              (const NumSequence::size_type v)    {  lengthStartContext        = v;     return *this; }
        Builder& setMarginStartContext              (const int v)                       {  marginStartContext        = v;     return *this; }
        Builder& setFgioDistanceThresh              (const NumSequence::size_type v)    {  fgioDistanceThresh        = v;     return *this; }
        Builder& setIgioDistanceThresh              (const NumSequence::size_type v)    {  igioDistanceThresh        = v;     return *this; }
        Builder& setPcounts                         (const unsigned v)                  {  pcounts                   = v;     return *this; }
        Builder& setGenomeGroup                     (const genome_group_t v)            {  genomeGroup               = v;     return *this; }
        Builder& setGcode                           (const GeneticCode::gcode_t v)      {  gcode                     = v;     return *this; }
        Builder& setMinimumGeneLengthTraining       (const NumSequence::size_type v)    {  minimumGeneLengthTraining = v;     return *this; }
        Builder& setOnlyTrainOnNativeGenes          (const bool v)                      {  onlyTrainOnNativeGenes    = v;     return *this; }
        Builder& setRunMotifSearch                  (const bool v)                      {  runMotifSearch            = v;     return *this; }
    };
    
    
    
    
}


#endif /* GMS2Trainer_hpp */











































