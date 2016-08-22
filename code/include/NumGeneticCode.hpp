//
//  NumGeneticCode.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NumGeneticCode_hpp
#define NumGeneticCode_hpp

#include <stdio.h>

#include "GeneticCode.hpp"
#include "CharNumConverter.hpp"

namespace gmsuite {
    
    /**
     * @class NumGeneticCode
     * @brief A class that represents a genetic code in numeric form
     *
     * A numeric form of a genetic code is basically the semantic equivalent
     * of a standard genetic code representation, but for sequences where the 
     * DNA letter alphabet has been mapped to numbers instead. E.g.
     *
     * A -> 0
     * C -> 1
     * G -> 2
     * T -> 3
     *
     * This provides a convenient method for identifying components in numeric sequences,
     * without having first to convert them back to their letter counterpart.
     */
    class NumGeneticCode {
        
    public:
        
        typedef GeneticCode::gcode_t gcode_t;           /**< Genetic code type @see GeneticCode */
        
        /**
         * Constructor: Create a numeric-based genetic code based off of 
         * a standard genetic code instance.
         *
         * @param gcode a genetic code to be mapped to its "numeric" equivalent
         */
        NumGeneticCode(const GeneticCode &gcode, const CharNumConverter &converter);
        
        
        /******************** Request info on Starts, Stops, and Translating AA ********************/
        
        /**
         * Check if codon is a start codon
         *
         * @param codon candidate codon
         * @return true if codon is a start; false otherwise
         */
        bool isStart(CharNumConverter::seq_t codon) const;
        
        
        /**
         * Check if codon is a stop codon
         *
         * @param codon candidate codon
         * @return true if codon is a stop; false otherwise
         */
        bool isStop(CharNumConverter::seq_t codon) const;
        
        
        /**
         * Get the list of start codons
         *
         * @return a vector of CharNumConverter::seq_t, where each element represents a start
         */
        vector<CharNumConverter::seq_t> getStarts() const;
        
        
        /**
         * Get the list of stop codons
         *
         * @return a vector of CharNumConverter::seq_t, where each element represents a stop
         */
        vector<CharNumConverter::seq_t> getStops() const;
        
        
        /**
         * Get the genetic code value that represents this object.
         *
         * @return the genetic code.
         */
        gcode_t getGCode() const;
        
        
    private:
        
        gcode_t gcode;                                      /**< Genetic code value         */
        vector<CharNumConverter::seq_t> starts;             /**< List of start codons       */
        vector<CharNumConverter::seq_t> stops;              /**< List of stop codons        */
    };
}

#endif /* NumGeneticCode_hpp */
