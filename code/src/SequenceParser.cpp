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


void SequenceParser::extractUpstreamSequences(const NumSequence& sequence, const vector<Label*> &labels, const CharNumConverter &cnc, NumSequence::size_type upstrLength, vector<NumSequence> &upstreamRegions, bool allowOverlapWithCDS) {
    
    if (sequence.size() == 0)
        return;
    
    // extract upstream region, and reverse complement when on negative strand
    upstreamRegions.resize(labels.size());
    
    size_t numSkipped = 0;          // number of labels skipped: i.e. won't be used in genome classification (because of no upstream sequence)
    
    for (size_t n = 0; n < labels.size(); n++) {
        try {
            bool skip = false;
            
            if (!allowOverlapWithCDS) {
                // skip if overlapping with CDS
                if (labels[n]->strand == Label::NEG) {
                    // boundaryRight is either end of sequence (if this is the last label), or the left-end of the gene to the "right" of the current gene
                    size_t boundaryRight = (n == labels.size()-1 ? sequence.size()-1 : labels[n+1]->left-1);
                    size_t boundaryLeft = labels[n]->right + 1;
                    
                    size_t length = boundaryRight - boundaryLeft + 1;
                    if (length < upstrLength)
                        skip = true;
                }
                else {
                    // boundaryLeft is either start of sequence (if this is the first label), or the right-end of the gene to the "left" of the current gene
                    size_t boundaryLeft = (n == 0 ? 0 : labels[n-1]->right+1);
                    size_t boundaryRight = labels[n]->left - 1;
                    
                    size_t length = boundaryRight - boundaryLeft + 1;
                    if (length < upstrLength)
                        skip = true;
                }
            }
            
            if (skip)
                numSkipped++;
            else
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
