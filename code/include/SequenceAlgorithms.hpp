//
//  SequenceAlgorithms.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef SequenceAlgorithms_hpp
#define SequenceAlgorithms_hpp

#include <stdio.h>

#include "NumSequence.hpp"


namespace gmsuite {
    
    class SequenceAlgorithms {
        
    public:
        
        static NumSequence longestCommonSubstring(const NumSequence &A, const NumSequence &B);
        
    };
    
    
}

#endif /* SequenceAlgorithms_hpp */
