//
//  SequenceParser.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/13/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "SequenceParser.hpp"
#include <stdexcept>

using namespace std;
using namespace gmsuite;


void SequenceParser::extractUpstreamSequences(const NumSequence& sequence, const vector<Label*> &labels, const CharNumConverter &cnc, NumSequence::size_type upstrLength, vector<NumSequence> &upstreamRegions) {
    
    // extract upstream region, and reverse complement when on negative strand
    upstreamRegions.resize(labels.size());
    
    size_t numSkipped = 0;          // number of labels skipped: i.e. won't be used in genome classification (because of no upstream sequence)
    
    for (size_t n = 0; n < labels.size(); n++) {
        try {
            upstreamRegions[n-numSkipped] = extractUpstreamSequence(sequence, *labels[n], cnc, upstrLength);
        }
        catch (out_of_range) {
            numSkipped++;
        }
    }
    
    // resize vector to remove skipped (empty) numsequences
    upstreamRegions.resize(upstreamRegions.size()-numSkipped);
}

// Extract upstream sequence: throws exception when 'not enough' sequence to extract
NumSequence SequenceParser::extractUpstreamSequence(const NumSequence& sequence, const Label &label, const CharNumConverter &cnc, NumSequence::size_type upstrLength) {
    
    if (upstrLength == 0)
        throw out_of_range("Cannot extract upstream sequence of length 0");
    
    if (label.strand == Label::POS) {           // positive strand
        
        if (label.left < upstrLength)
            throw out_of_range("Not enough sequence for position " + to_string(label.left) + " to extract upstream of length " + to_string(upstrLength));
        
        size_t left = label.left - upstrLength;                     // left idx of upstream sequence
        size_t right = label.left-1;                                // right idx of upstream sequence
        size_t length = right - left + 1;                           // length of upstream sequence
        
        NumSequence upstream = sequence.subseq(left, length);       // get upstream subsequence
        return upstream;
    }
    else {                                      // negative strand; reverse complement
        
        if (label.right + upstrLength >= sequence.size())
            throw out_of_range("Not enought sequence for position " + to_string(label.right) + " on negative strand to extract upstream of length " + to_string(upstrLength));
        
        size_t left = label.right + 1;                              // left idx of upstream sequence
        size_t right = label.right + upstrLength;                   // right idx of upstream sequence
        size_t length = right - left + 1;                           // length of upstream sequence
        
        NumSequence upstream = sequence.subseq(left, length);       // get upstream subsequence
        upstream.reverseComplement(cnc);                            // reverse complement
        return upstream;
    }
}