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
        
        // Genome groups:
        // A: SD RBS
        // B: non-SD RBS
        // C: Bacteria leaderless
        // D: Archaea leaderless
        // X: Upstream signature + RBS
        // B2: non-SD RBS + SD RBS
        // D2: Archaea leaderless step 2
        typedef enum {A, B, C, D, X, D2, B2} genome_group_t;
    };
}

#endif /* ProkGeneStartModel_hpp */
