//
//  MotifMarkov.hpp
//  Biogem CPP
//
//  Created by Karl Gemayel on 8/20/18.
//  Copyright © 2018 Georgia Institute of Technology. All rights reserved.
//

#ifndef MotifMarkov_hpp
#define MotifMarkov_hpp

#include <stdio.h>
#include "NonUniformMarkov.hpp"
#include <string>
#include <vector>

using std::vector;
using std::string;
using std::pair;

namespace gmsuite {
    
    class MotifMarkov : public NonUniformMarkov {
        
    public:
        
        MotifMarkov(const vector<vector<pair<string, double> > > &keyValue, size_t length, const NumAlphabetDNA &alph, const CharNumConverter &cnc);
    };
}



#endif /* MotifMarkov_hpp */
