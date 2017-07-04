//
//  ProkGeneStartModel.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/13/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ProkGeneStartModel_hpp
#define ProkGeneStartModel_hpp

#include <stdio.h>

namespace gmsuite {
    
    class ProkGeneStartModel {
    
    public:
        
        // class of genome
        typedef enum {C1, C2, C3, C4, C5} genome_class_t;           // C5: non-canonical rbs
    };
}

#endif /* ProkGeneStartModel_hpp */
