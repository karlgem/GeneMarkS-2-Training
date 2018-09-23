//
//  CodingMarkov.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 12/1/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef CodingMarkov_hpp
#define CodingMarkov_hpp

#include <stdio.h>
#include "PeriodicMarkov.hpp"
#include "NumGeneticCode.hpp"

namespace gmsuite {
    
    /**
     * @class CodingMarkov
     * @brief A Markov model for coding regions.
     */
    class CodingMarkov : public PeriodicMarkov {
        
    public:
        
        /**
         * Constructor:
         *
         * @param order the model's order
         * @param period the model's period
         * @param alph the alphabet used by the model
         * @param geneticCode the genetic code
         */
        CodingMarkov(unsigned order, size_t period, const NumAlphabetDNA &alph, const NumGeneticCode &geneticCode);
        
        CodingMarkov(const vector<vector<pair<string, double> > >  &keyValue, const NumAlphabetDNA &alph, const CharNumConverter &cnc);
        
    protected:
        
        const NumGeneticCode* geneticCode;      /**< Numeric version of the genetic code */
        
        //@override
        void addPseudocounts(int pcount);
        
    };
}

#endif /* CodingMarkov_hpp */
