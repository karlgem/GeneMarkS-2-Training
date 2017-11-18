//
//  SequenceParser.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/13/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef SequenceParser_hpp
#define SequenceParser_hpp

#include <stdio.h>

#include "Label.hpp"
#include "NumSequence.hpp"
#include "CharNumConverter.hpp"


namespace gmsuite {
    
    class SequenceParser {
        
    public:
        
        static NumSequence extractUpstreamSequence(const NumSequence& sequence, const Label &label, const CharNumConverter &cnc, NumSequence::size_type upstrLength);
        
        static void extractUpstreamSequences(const NumSequence& sequence, const vector<Label*> &labels, const CharNumConverter &cnc, NumSequence::size_type upstrLength, vector<NumSequence> &upstreams, bool allowOverlapWithCDS = false, size_t minimumGeneLength=0, const vector<bool> &use = vector<bool>());
        
        static void extractStartContextSequences(const NumSequence& sequence, const vector<Label*> &labels, const CharNumConverter &cnc, long long posRelToStart, NumSequence::size_type length, vector<NumSequence> &contexts, const vector<bool> &use = vector<bool>());
        
        /**
         * Extract the sequence around the provided label's start-codon (based on strand).
         * @param posRelToStart is the number of nucleotides from the begining of the start-context is extracted. Set to 0 means it starts at the gene-start. Negative numbers indicate movement upstream, while positive numbers indicate movement downstream.
         */
        static NumSequence extractStartContextSequence(const NumSequence& sequence, const Label &label, const CharNumConverter &cnc, long long posRelToStart, NumSequence::size_type length);
        
    };
    
    
}
#endif /* SequenceParser_hpp */
