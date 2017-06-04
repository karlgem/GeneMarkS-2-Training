//
//  NonCodingMarkov.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 11/3/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef NonCodingMarkov_hpp
#define NonCodingMarkov_hpp

#include <stdio.h>
#include "UniformMarkov.hpp"
#include <string>
#include <vector>

using std::vector;
using std::string;
using std::pair;

namespace gmsuite {
    
    class NonCodingMarkov : public UniformMarkov {
        
    public:
        
        NonCodingMarkov(const vector<pair<string, double> > &keyValue, const NumAlphabetDNA &alph, const CharNumConverter &cnc);
        
        
        
    private:
        
    };
    
}


#endif /* NonCodingMarkov_hpp */
