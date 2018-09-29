//
//  OptionsGMS2Training.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OptionsGMS2Training_hpp
#define OptionsGMS2Training_hpp

#include <stdio.h>
#include <string>

#include "Options.hpp"
#include "NumSequence.hpp"
#include "GeneticCode.hpp"
#include "OptionsMFinder.hpp"
#include "ProkGeneStartModel.hpp"

using std::string;

namespace gmsuite {
    
    /**
     * @class OptionsGMS2Training
     * @brief A class that deals with parsing command-line options for GMS2Training module
     */
    class OptionsGMS2Training : public Options {
        
        
    public:
        
        OptionsGMS2Training(string mode="gms2-training");
        
        /**
         * Parse the command-line words into arguments
         *
         * @param argc the number of command-line words
         * @param argv a vector of the words
         * @return true if the parse is successful, false otherwise
         */
        bool parse(int argc, const char *argv[]);
        
        
        static void addProcessOptions(OptionsGMS2Training &options, po::options_description &processOptions);
        
        // Below, create a variable for each parameter, to make for easy access
    public:
        
        typedef ProkGeneStartModel::genome_group_t genome_group_t;
        typedef GeneticCode::gcode_t gcode_t;
        
        string fn_sequence;             /**< Input filename containing DNA sequence */
        string fn_labels;               /**< Input filename containing labels */
        string fn_outmod;               /**< Output model file */
        string fn_settings;             /**< Settings to place in output mod file */
        
        // prediction parameters
        double nonProbN;
        double codProbN;
        double nonDurationDecay;
        double codDurationDecay;
        NumSequence::size_type geneMinLengthPrediction;
        
        // GMS2 model parameters
        
        unsigned codingOrder;                       /**< Coding model's Markov order */
        unsigned noncodingOrder;                    /**< Noncoding model's Markov order */
        unsigned startContextOrder;                 /**< Start-context Markov order */
        int startContextMargin;                     /**< Start-context margin (distance from start) */
        NumSequence::size_type upstreamLength;      /**< Upstream length used for motif search */
        NumSequence::size_type startContextLength;  /**< length of start context model */
        NumSequence::size_type MIN_GENE_LEN;        /**< minimum gene length required for using gene in training */ 
        bool trainOnNativeOnly;                     /**< train on genes predicted by native model */
        OptionsMFinder optionsMFinder;              /**< Options for Motif Finder */
        
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
        NumSequence::size_type  groupB2_upstreamLengthSDRBS         ;
        unsigned                groupB2_widthNonSDRBS               ;
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
        
    };
}


#endif /* OptionsGMS2Training_hpp */
