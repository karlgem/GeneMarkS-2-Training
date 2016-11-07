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

using std::string;

namespace gmsuite {
    
    class NonCodingMarkov : public UniformMarkov {
        
    public:
        
        NonCodingMarkov(string filename, const NumAlphabetDNA &alph);
        
        
        
    private:
        
    };
    
}


#endif /* NonCodingMarkov_hpp */
