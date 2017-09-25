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
        // A: Archaea leaderless
        // B: Bacteria leaderless
        // C: non-SD RBS
        // D: SD RBS
        // E: Upstream signature + RBS
        // A2: Archaea leaderless step 2
        // C2: non-SD RBS + SD RBS
        typedef enum {A, B, C, D, E, A2, C2} genome_group_t;
    };
}

#endif /* ProkGeneStartModel_hpp */
