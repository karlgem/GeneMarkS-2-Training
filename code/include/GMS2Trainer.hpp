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
    
    
    class GMS2TrainerParameters {
        
        typedef ProkGeneStartModel::genome_group_t genome_group_t;

    public:
        
        // Group-A
        unsigned                groupD_widthPromoter                ;
        unsigned                groupD_widthRBS                     ;
        NumSequence::size_type  groupD_upstreamLengthPromoter       ;
        NumSequence::size_type  groupD_upstreamLengthRBS            ;
        double                  groupD_spacerScoreThresh            ;
        NumSequence::size_type  groupD_spacerDistThresh             ;
        NumSequence::size_type  groupD_spacerWindowSize             ;
        string                  groupD_extendedSD                   ;
        unsigned                groupD_minMatchToExtendedSD         ;
        bool                    groupD_allowAGSubstitution          ;
        
        // Group B
        unsigned                groupC_widthPromoter                ;
        unsigned                groupC_widthRBS                     ;
        NumSequence::size_type  groupC_upstreamLengthPromoter       ;
        NumSequence::size_type  groupC_upstreamLengthRBS            ;
        double                  groupC_spacerScoreThresh            ;
        NumSequence::size_type  groupC_spacerDistThresh             ;
        NumSequence::size_type  groupC_spacerWindowSize             ;
        string                  groupC_extendedSD                   ;
        unsigned                groupC_minMatchToExtendedSD         ;
        bool                    groupC_allowAGSubstitution          ;
        
        // Group-C
        unsigned                groupB_widthRBS                     ;
        NumSequence::size_type  groupB_upstreamLengthRBS            ;
        NumSequence::size_type  groupB_upstreamRegion3Prime         ;
        unsigned                groupB_minMatchRBSPromoter          ;
        unsigned                groupB_minMatchToExtendedSD         ;
        string                  groupB_extendedSD                   ;
        
        // Group-C-2
        unsigned                groupB2_widthSDRBS                  ;
        unsigned                groupB2_widthNonSDRBS               ;
        NumSequence::size_type  groupB2_upstreamLengthSDRBS         ;
        NumSequence::size_type  groupB2_upstreamLengthNonSDRBS      ;
        NumSequence::size_type  groupB2_upstreamRegion3Prime        ;
        unsigned                groupB2_minMatchToExtendedSD        ;
        string                  groupB2_extendedSD                  ;
        
        // Group-D
        unsigned                groupA_widthRBS                     ;
        NumSequence::size_type  groupA_upstreamLengthRBS            ;
        double                  groupA_percentMatchRBS              ;
        string                  groupA_extendedSD                   ;
        unsigned                groupA_minMatchToExtendedSD         ;
        bool                    groupA_allowAGSubstitution          ;
        
        // Group-E
        unsigned                groupX_widthRBS                     ;
        NumSequence::size_type  groupX_upstreamLengthRBS            ;
        NumSequence::size_type  groupX_lengthUpstreamSignature      ;
        unsigned                groupX_orderUpstreamSignature       ;
        string                  groupX_extendedSD                   ;
        unsigned                groupX_minMatchToExtendedSD         ;
        bool                    groupX_allowAGSubstitution          ;
        
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
        const OptionsMFinder*   optionsMFinder                      ;
    };
    
    
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
                    bool cutPromTrainSeqs = false);             // Deprecated
        
        GMS2Trainer(// Coding and Noncoding Models
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
                    unsigned                groupD_widthPromoter                ,
                    unsigned                groupD_widthRBS                     ,
                    NumSequence::size_type  groupD_upstreamLengthPromoter       ,
                    NumSequence::size_type  groupD_upstreamLengthRBS            ,
                    double                  groupD_spacerScoreThresh            ,
                    NumSequence::size_type  groupD_spacerDistThresh             ,
                    NumSequence::size_type  groupD_spacerWindowSize             ,
                    string                  groupD_extendedSD                   ,
                    unsigned                groupD_minMatchToExtendedSD         ,
                    bool                    groupD_allowAGSubstitution          ,
                    // Group B
                    unsigned                groupC_widthPromoter                ,
                    unsigned                groupC_widthRBS                     ,
                    NumSequence::size_type  groupC_upstreamLengthPromoter       ,
                    NumSequence::size_type  groupC_upstreamLengthRBS            ,
                    double                  groupC_spacerScoreThresh            ,
                    NumSequence::size_type  groupC_spacerDistThresh             ,
                    NumSequence::size_type  groupC_spacerWindowSize             ,
                    string                  groupC_extendedSD                   ,
                    unsigned                groupC_minMatchToExtendedSD         ,
                    bool                    groupC_allowAGSubstitution          ,
                    // Group-C
                    unsigned                groupB_widthRBS                     ,
                    NumSequence::size_type  groupB_upstreamLengthRBS            ,
                    NumSequence::size_type  groupB_upstreamRegion3Prime         ,
                    unsigned                groupB_minMatchRBSPromoter          ,
                    unsigned                groupB_minMatchToExtendedSD         ,
                    string                  groupB_extendedSD                   ,
                    // Group-C-2
                    unsigned                groupB2_widthSDRBS                  ,
                    unsigned                groupB2_widthNonSDRBS               ,
                    NumSequence::size_type  groupB2_upstreamLengthSDRBS         ,
                    NumSequence::size_type  groupB2_upstreamLengthNonSDRBS      ,
                    NumSequence::size_type  groupB2_upstreamRegion3Prime        ,
                    unsigned                groupB2_minMatchToExtendedSD        ,
                    string                  groupB2_extendedSD                  ,
                    // Group-D
                    unsigned                groupA_widthRBS                     ,
                    NumSequence::size_type  groupA_upstreamLengthRBS            ,
                    double                  groupA_percentMatchRBS              ,
                    string                  groupA_extendedSD                   ,
                    unsigned                groupA_minMatchToExtendedSD         ,
                    bool                    groupA_allowAGSubstitution          ,
                    // Group-E
                    unsigned                groupX_widthRBS                     ,
                    NumSequence::size_type  groupX_upstreamLengthRBS            ,
                    NumSequence::size_type  groupX_lengthUpstreamSignature      ,
                    unsigned                groupX_orderUpstreamSignature       ,
                    string                  groupX_extendedSD                   ,
                    unsigned                groupX_minMatchToExtendedSD         ,
                    bool                    groupX_allowAGSubstitution
        );
        
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
        
        void estimateParametersMotifModel_GroupD(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersMotifModel_groupD2(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersMotifModel_GroupC(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersMotifModel_GroupB(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersMotifModel_GroupB2(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersMotifModel_GroupA(const NumSequence &sequence, const vector<Label *> &labels);
        void estimateParametersMotifModel_GroupX(const NumSequence &sequence, const vector<Label*> &labels);
        
        
        // public variables for models
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
        
        
    public:                 // parameters
        
        GMS2TrainerParameters params;
        size_t numLeaderless;
        size_t numFGIO;
        string genomeType;

        const AlphabetDNA *charAlph;
        const NumAlphabetDNA *alphabet;
        const GeneticCode *geneticCode;
        const NumGeneticCode *numGeneticCode;
        const CharNumConverter *cnc;
        
        
        
    };
    
    
    // Class for building a GMS2Trainer (for advanced users)
    class GMS2Trainer::Builder {
        
        typedef ProkGeneStartModel::genome_group_t genome_group_t;
        
    private:
        // variables needed for the GMS2Trainer class
        
        // Group-A
        unsigned                groupD_widthPromoter                ;
        unsigned                groupD_widthRBS                     ;
        NumSequence::size_type  groupD_upstreamLengthPromoter       ;
        NumSequence::size_type  groupD_upstreamLengthRBS            ;
        double                  groupD_spacerScoreThresh            ;
        NumSequence::size_type  groupD_spacerDistThresh             ;
        NumSequence::size_type  groupD_spacerWindowSize             ;
        string                  groupD_extendedSD                   ;
        unsigned                groupD_minMatchToExtendedSD         ;
        bool                    groupD_allowAGSubstitution          ;
        
        // Group B
        unsigned                groupC_widthPromoter                ;
        unsigned                groupC_widthRBS                     ;
        NumSequence::size_type  groupC_upstreamLengthPromoter       ;
        NumSequence::size_type  groupC_upstreamLengthRBS            ;
        double                  groupC_spacerScoreThresh            ;
        NumSequence::size_type  groupC_spacerDistThresh             ;
        NumSequence::size_type  groupC_spacerWindowSize             ;
        string                  groupC_extendedSD                   ;
        unsigned                groupC_minMatchToExtendedSD         ;
        bool                    groupC_allowAGSubstitution          ;
        
        // Group-C
        unsigned                groupB_widthRBS                     ;
        NumSequence::size_type  groupB_upstreamLengthRBS            ;
        NumSequence::size_type  groupB_upstreamRegion3Prime         ;
        unsigned                groupB_minMatchRBSPromoter          ;
        unsigned                groupB_minMatchToExtendedSD         ;
        string                  groupB_extendedSD                   ;
        
        // Group-C-2
        unsigned                groupB2_widthSDRBS                  ;
        unsigned                groupB2_widthNonSDRBS               ;
        NumSequence::size_type  groupB2_upstreamLengthSDRBS         ;
        NumSequence::size_type  groupB2_upstreamLengthNonSDRBS      ;
        NumSequence::size_type  groupB2_upstreamRegion3Prime        ;
        unsigned                groupB2_minMatchToExtendedSD        ;
        string                  groupB2_extendedSD                  ;
        
        // Group-D
        unsigned                groupA_widthRBS                     ;
        NumSequence::size_type  groupA_upstreamLengthRBS            ;
        double                  groupA_percentMatchRBS              ;
        string                  groupA_extendedSD                   ;
        unsigned                groupA_minMatchToExtendedSD         ;
        bool                    groupA_allowAGSubstitution          ;
        
        // Group-E
        unsigned                groupX_widthRBS                     ;
        NumSequence::size_type  groupX_upstreamLengthRBS            ;
        NumSequence::size_type  groupX_lengthUpstreamSignature      ;
        unsigned                groupX_orderUpstreamSignature       ;
        string                  groupX_extendedSD                   ;
        unsigned                groupX_minMatchToExtendedSD         ;
        bool                    groupX_allowAGSubstitution          ;
        
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
        const OptionsMFinder*   optionsMFinder                      ;
        
        
    public:
        
        Builder() {
            // Group-A
            groupD_widthPromoter                = 12                          ;
            groupD_widthRBS                     = 6                           ;
            groupD_upstreamLengthPromoter       = 20                          ;
            groupD_upstreamLengthRBS            = 20                          ;
            groupD_spacerScoreThresh            = 0.1                         ;
            groupD_spacerDistThresh             = 14                          ;
            groupD_spacerWindowSize             = 1                           ;
            groupD_extendedSD                   = "TAAGGAGGTGA"               ;
            groupD_minMatchToExtendedSD         = 4                           ;
            groupD_allowAGSubstitution          = true                        ;
            
            // Group-B
            groupC_widthPromoter                = 6                           ;
            groupC_widthRBS                     = 6                           ;
            groupC_upstreamLengthPromoter       = 20                          ;
            groupC_upstreamLengthRBS            = 20                          ;
            groupC_spacerScoreThresh            = 0.25                        ;
            groupC_spacerDistThresh             = 14                          ;
            groupC_spacerWindowSize             = 1                           ;
            groupC_extendedSD                   = "TAAGGAGGTGA"               ;
            groupC_minMatchToExtendedSD         = 4                           ;
            groupC_allowAGSubstitution          = true                        ;
            
            // Group-C
            groupB_widthRBS                     = 6                           ;
            groupB_upstreamLengthRBS            = 20                          ;
            groupB_upstreamRegion3Prime         = 0                           ;
            groupB_minMatchRBSPromoter          = 3                           ;
            groupB_minMatchToExtendedSD         = 4                           ;
            groupB_extendedSD                   = "TAAGGAGGTGA"               ;
            
            // Group-C-2
            groupB2_widthSDRBS                  = 6                           ;
            groupB2_upstreamLengthSDRBS         = 20                          ;
            groupB2_widthNonSDRBS               = 6                           ;
            groupB2_upstreamLengthNonSDRBS      = 20                          ;
            groupB2_upstreamRegion3Prime        = 0                           ;
            groupB2_minMatchToExtendedSD        = 4                           ;
            groupB2_extendedSD                  = "TAAGGAGGTGA"               ;
            
            // Group-D
            groupA_widthRBS                     = 6                           ;
            groupA_upstreamLengthRBS            = 20                          ;
            groupA_percentMatchRBS              = 0.5                         ;
            groupA_extendedSD                   = "TAAGGAGGTGA"               ;
            groupA_minMatchToExtendedSD         = 4                           ;
            groupA_allowAGSubstitution          = true                        ;
            
            // Group-E
            groupX_widthRBS                     = 6                           ;
            groupX_upstreamLengthRBS            = 20                          ;
            groupX_lengthUpstreamSignature      = 35                          ;
            groupX_orderUpstreamSignature       = 2                           ;
            groupX_extendedSD                   = "TAAGGAGGTGA"               ;
            groupX_minMatchToExtendedSD         = 4                           ;
            groupX_allowAGSubstitution          = true                        ;
            
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
            genomeGroup                         = ProkGeneStartModel::A       ;
            gcode                               = GeneticCode::ELEVEN         ;
            minimumGeneLengthTraining           = 300                         ;
            onlyTrainOnNativeGenes              = false                       ;
            runMotifSearch                      = true                        ;
            optionsMFinder                      = NULL                        ;     // FIXME: figure out how to set default MFinder options
        }
        
        GMS2Trainer build() {
            return GMS2Trainer (orderCoding, orderNonCoding, orderStartContext, lengthStartContext, marginStartContext, fgioDistanceThresh, igioDistanceThresh, pcounts, genomeGroup, gcode, minimumGeneLengthTraining, onlyTrainOnNativeGenes, runMotifSearch, *optionsMFinder, groupD_widthPromoter, groupD_widthRBS, groupD_upstreamLengthPromoter, groupD_upstreamLengthRBS, groupD_spacerScoreThresh, groupD_spacerDistThresh, groupD_spacerWindowSize, groupD_extendedSD, groupD_minMatchToExtendedSD, groupD_allowAGSubstitution, groupC_widthPromoter, groupC_widthRBS, groupC_upstreamLengthPromoter, groupC_upstreamLengthRBS, groupC_spacerScoreThresh, groupC_spacerDistThresh, groupC_spacerWindowSize, groupC_extendedSD, groupC_minMatchToExtendedSD, groupC_allowAGSubstitution, groupB_widthRBS, groupB_upstreamLengthRBS, groupB_upstreamRegion3Prime, groupB_minMatchRBSPromoter, groupB_minMatchToExtendedSD, groupB_extendedSD, groupB2_widthSDRBS, groupB2_widthNonSDRBS, groupB2_upstreamLengthSDRBS, groupB2_upstreamLengthNonSDRBS, groupB2_upstreamRegion3Prime, groupB2_minMatchToExtendedSD, groupB2_extendedSD, groupA_widthRBS, groupA_upstreamLengthRBS, groupA_percentMatchRBS, groupA_extendedSD, groupA_minMatchToExtendedSD, groupA_allowAGSubstitution, groupX_widthRBS, groupX_upstreamLengthRBS, groupX_lengthUpstreamSignature, groupX_orderUpstreamSignature, groupX_extendedSD, groupX_minMatchToExtendedSD, groupX_allowAGSubstitution);
        }
        
        GMS2Trainer build(const OptionsGMS2Training &options) {
            
            setGroupD_widthPromoter          (options.groupD_widthPromoter          );
            setGroupD_widthRBS               (options.groupD_widthRBS               );
            setGroupD_upstreamLengthPromoter (options.groupD_upstreamLengthPromoter );
            setGroupD_upstreamLengthRBS      (options.groupD_upstreamLengthRBS      );
            setGroupD_spacerScoreThresh      (options.groupD_spacerScoreThresh      );
            setGroupD_spacerDistThresh       (options.groupD_spacerDistThresh       );
            setGroupD_spacerWindowSize       (options.groupD_spacerWindowSize       );
            setGroupD_extendedSD             (options.groupD_extendedSD             );
            setGroupD_minMatchToExtendedSD   (options.groupD_minMatchToExtendedSD   );
            setGroupD_allowAGSubstitution    (options.groupD_allowAGSubstitution    );
            setGroupC_widthPromoter          (options.groupC_widthPromoter          );
            setGroupC_widthRBS               (options.groupC_widthRBS               );
            setGroupC_upstreamLengthPromoter (options.groupC_upstreamLengthPromoter );
            setGroupC_upstreamLengthRBS      (options.groupC_upstreamLengthRBS      );
            setGroupC_spacerScoreThresh      (options.groupC_spacerScoreThresh      );
            setGroupC_spacerDistThresh       (options.groupC_spacerDistThresh       );
            setGroupC_spacerWindowSize       (options.groupC_spacerWindowSize       );
            setGroupC_extendedSD             (options.groupC_extendedSD             );
            setGroupC_minMatchToExtendedSD   (options.groupC_minMatchToExtendedSD   );
            setGroupC_allowAGSubstitution    (options.groupC_allowAGSubstitution    );
            setGroupB_widthRBS               (options.groupB_widthRBS               );
            setGroupB_upstreamLengthRBS      (options.groupB_upstreamLengthRBS      );
            setGroupB_upstreamRegion3Prime   (options.groupB_upstreamRegion3Prime   );
            setGroupB_minMatchRBSPromoter    (options.groupB_minMatchRBSPromoter    );
            setGroupB_minMatchToExtendedSD   (options.groupB_minMatchToExtendedSD   );
            setGroupB_extendedSD             (options.groupB_extendedSD             );
            
            setGroupB2_widthSDRBS              (options.groupB2_widthSDRBS             );
            setGroupB2_widthNonSDRBS           (options.groupB2_widthNonSDRBS          );
            setGroupB2_upstreamLengthSDRBS     (options.groupB2_upstreamLengthSDRBS    );
            setGroupB2_upstreamLengthNonSDRBS  (options.groupB2_upstreamLengthNonSDRBS );
            setGroupB2_upstreamRegion3Prime    (options.groupB2_upstreamRegion3Prime   );
            setGroupB2_minMatchToExtendedSD    (options.groupB2_minMatchToExtendedSD   );
            setGroupB2_extendedSD              (options.groupB2_extendedSD             );
            
            
            setGroupA_widthRBS               (options.groupA_widthRBS               );
            setGroupA_upstreamLengthRBS      (options.groupA_upstreamLengthRBS      );
            setGroupA_percentMatchRBS        (options.groupA_percentMatchRBS        );
            setGroupA_extendedSD             (options.groupA_extendedSD             );
            setGroupA_minMatchToExtendedSD   (options.groupA_minMatchToExtendedSD   );
            setGroupA_allowAGSubstitution    (options.groupA_allowAGSubstitution    );
            setGroupX_widthRBS               (options.groupX_widthRBS               );
            setGroupX_upstreamLengthRBS      (options.groupX_upstreamLengthRBS      );
            setGroupX_lengthUpstreamSignature(options.groupX_lengthUpstreamSignature);
            setGroupX_orderUpstreamSignature (options.groupX_orderUpstreamSignature );
            setGroupX_extendedSD             (options.groupX_extendedSD             );
            setGroupX_minMatchToExtendedSD   (options.groupX_minMatchToExtendedSD   );
            setGroupX_allowAGSubstitution    (options.groupX_allowAGSubstitution    );
            setOrderCoding                   (options.orderCoding                   );
            setOrderNonCoding                (options.orderNonCoding                );
            setOrderStartContext             (options.orderStartContext             );
            setLengthStartContext            (options.lengthStartContext            );
            setMarginStartContext            (options.marginStartContext            );
            setFgioDistanceThresh            (options.fgioDistanceThresh            );
            setIgioDistanceThresh            (options.igioDistanceThresh            );
            setPcounts                       (options.pcounts                       );
            setGenomeGroup                   (options.genomeGroup                   );
            setGcode                         (options.gcode                         );
            setMinimumGeneLengthTraining     (options.minimumGeneLengthTraining     );
            setOnlyTrainOnNativeGenes        (options.onlyTrainOnNativeGenes        );
            setRunMotifSearch                (options.runMotifSearch                );
            setOptionsMFinder                (options.optionsMFinder                );
            
            return build();
        }
        
        // set custom values for MotifFinder creation
        // returns Builder for shorthand inline usage
        Builder& setGroupD_widthPromoter            (const unsigned v)                  { this->groupD_widthPromoter           = v;    return *this; }
        Builder& setGroupD_widthRBS                 (const unsigned v)                  { this->groupD_widthRBS                = v;    return *this; }
        Builder& setGroupD_upstreamLengthPromoter   (const NumSequence::size_type v)    { this->groupD_upstreamLengthPromoter  = v;    return *this; }
        Builder& setGroupD_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupD_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupD_spacerScoreThresh        (const double v)                    { this->groupD_spacerScoreThresh       = v;    return *this; }
        Builder& setGroupD_spacerDistThresh         (const NumSequence::size_type v)    { this->groupD_spacerDistThresh        = v;    return *this; }
        Builder& setGroupD_spacerWindowSize         (const NumSequence::size_type v)    { this->groupD_spacerWindowSize        = v;    return *this; }
        Builder& setGroupD_extendedSD               (const string v)                    { this->groupD_extendedSD              = v;    return *this; }
        Builder& setGroupD_minMatchToExtendedSD     (const unsigned v)                  { this->groupD_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupD_allowAGSubstitution      (const bool v)                      { this->groupD_allowAGSubstitution     = v;    return *this; }
        
        Builder& setGroupC_widthPromoter            (const unsigned v)                  { this->groupC_widthPromoter           = v;    return *this; }
        Builder& setGroupC_widthRBS                 (const unsigned v)                  { this->groupC_widthRBS                = v;    return *this; }
        Builder& setGroupC_upstreamLengthPromoter   (const NumSequence::size_type v)    { this->groupC_upstreamLengthPromoter  = v;    return *this; }
        Builder& setGroupC_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupC_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupC_spacerScoreThresh        (const double v)                    { this->groupC_spacerScoreThresh       = v;    return *this; }
        Builder& setGroupC_spacerDistThresh         (const NumSequence::size_type v)    { this->groupC_spacerDistThresh        = v;    return *this; }
        Builder& setGroupC_spacerWindowSize         (const NumSequence::size_type v)    { this->groupC_spacerWindowSize        = v;    return *this; }
        Builder& setGroupC_extendedSD               (const string v)                    { this->groupC_extendedSD              = v;    return *this; }
        Builder& setGroupC_minMatchToExtendedSD     (const unsigned v)                  { this->groupC_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupC_allowAGSubstitution      (const bool v)                      { this->groupC_allowAGSubstitution     = v;    return *this; }
        Builder& setGroupB_widthRBS                 (const unsigned v)                  { this->groupB_widthRBS                = v;    return *this; }
        Builder& setGroupB_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupB_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupB_upstreamRegion3Prime     (const NumSequence::size_type v)    { this->groupB_upstreamRegion3Prime    = v;    return *this; }
        Builder& setGroupB_minMatchRBSPromoter      (const unsigned v)                  { this->groupB_minMatchRBSPromoter     = v;    return *this; }
        Builder& setGroupB_minMatchToExtendedSD     (const unsigned v)                  { this->groupB_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupB_extendedSD               (const string v)                    { this->groupB_extendedSD              = v;    return *this; }
        
        // Group-C-2
        Builder& setGroupB2_widthSDRBS              (const unsigned v)                  { this->groupB2_widthSDRBS             = v;    return *this; }
        Builder& setGroupB2_widthNonSDRBS           (const unsigned v)                  { this->groupB2_widthNonSDRBS          = v;    return *this; }
        Builder& setGroupB2_upstreamLengthSDRBS     (const NumSequence::size_type v)    { this->groupB2_upstreamLengthSDRBS    = v;    return *this; }
        Builder& setGroupB2_upstreamLengthNonSDRBS  (const NumSequence::size_type v)    { this->groupB2_upstreamLengthNonSDRBS = v;    return *this; }
        Builder& setGroupB2_upstreamRegion3Prime    (const NumSequence::size_type v)    { this->groupB2_upstreamRegion3Prime   = v;    return *this; }
        Builder& setGroupB2_minMatchToExtendedSD    (const unsigned v)                  { this->groupB2_minMatchToExtendedSD   = v;    return *this; }
        Builder& setGroupB2_extendedSD              (const string v)                    { this->groupB2_extendedSD             = v;    return *this; }
        
        Builder& setGroupA_widthRBS                 (const unsigned v)                  { this->groupA_widthRBS                = v;    return *this; }
        Builder& setGroupA_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupA_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupA_percentMatchRBS          (const double v)                    { this->groupA_percentMatchRBS         = v;    return *this; }
        Builder& setGroupA_extendedSD               (const string v)                    { this->groupA_extendedSD              = v;    return *this; }
        Builder& setGroupA_minMatchToExtendedSD     (const unsigned v)                  { this->groupA_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupA_allowAGSubstitution      (const bool v)                      { this->groupA_allowAGSubstitution     = v;    return *this; }
        Builder& setGroupX_widthRBS                 (const unsigned v)                  { this->groupX_widthRBS                = v;    return *this; }
        Builder& setGroupX_upstreamLengthRBS        (const NumSequence::size_type v)    { this->groupX_upstreamLengthRBS       = v;    return *this; }
        Builder& setGroupX_lengthUpstreamSignature  (const NumSequence::size_type v)    { this->groupX_lengthUpstreamSignature = v;    return *this; }
        Builder& setGroupX_orderUpstreamSignature   (const unsigned v)                  { this->groupX_orderUpstreamSignature  = v;    return *this; }
        Builder& setGroupX_extendedSD               (const string v)                    { this->groupX_extendedSD              = v;    return *this; }
        Builder& setGroupX_minMatchToExtendedSD     (const unsigned v)                  { this->groupX_minMatchToExtendedSD    = v;    return *this; }
        Builder& setGroupX_allowAGSubstitution      (const bool v)                      { this->groupX_allowAGSubstitution     = v;    return *this; }
        
        
        
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
        Builder& setOptionsMFinder                  (const OptionsMFinder &v)           {  optionsMFinder            = &v;     return *this; }
    };
    
    
    
    
}


#endif /* GMS2Trainer_hpp */











































