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
#include <vector>
#include "NumSequence.hpp"

namespace gmsuite {
    
    class SequenceAlgorithms {
        
    public:
        
        static NumSequence longestCommonSubstring(const NumSequence &A, const NumSequence &B,
                                                  const std::vector<std::pair<NumSequence::num_t, NumSequence::num_t> >& subs = std::vector<std::pair<NumSequence::num_t, NumSequence::num_t> > ());
        
        
        static NumSequence longestMatchTo16S(const NumSequence &A, const NumSequence &B,
                                            std::pair<NumSequence::size_type, NumSequence::size_type>& positionsOfMatches,
                                            const std::vector<std::pair<NumSequence::num_t, NumSequence::num_t> >& subs = std::vector<std::pair<NumSequence::num_t, NumSequence::num_t> > ()
                                             );
        
    };
    
    
}

#endif /* SequenceAlgorithms_hpp */
