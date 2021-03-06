//
//  SequenceParser.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/13/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "SequenceParser.hpp"
#include <cstdlib>
#include <stdexcept>
#include <iostream>

using namespace std;
using namespace gmsuite;


void SequenceParser::extractUpstreamSequences(const NumSequence& sequence, const vector<Label*> &labels, const CharNumConverter &cnc, NumSequence::size_type upstrLength, vector<NumSequence> &upstreamRegions, bool allowOverlapWithCDS, size_t minimumGeneLength, const vector<bool> &use) {
    
    if (sequence.size() == 0)
        return;
    
    bool useAll = true;
    if (use.size() == labels.size()) {
        useAll = false;
    }
    
    // extract upstream region, and reverse complement when on negative strand
    upstreamRegions.resize(labels.size());
    
    size_t numSkipped = 0;          // number of labels skipped: i.e. won't be used in genome classification (because of no upstream sequence)
    
    for (size_t n = 0; n < labels.size(); n++) {
        try {
            bool skip = false;
            
            if (labels[n]->right - labels[n]->left + 1 < minimumGeneLength)
                skip = true;
            
            if (!allowOverlapWithCDS) {
                // skip if overlapping with CDS
                if (labels[n]->strand == Label::NEG) {
                    // boundaryRight is either end of sequence (if this is the last label), or the left-end of the gene to the "right" of the current gene
                    size_t boundaryRight = (n == labels.size()-1 ? sequence.size()-1 : labels[n+1]->left-1);
                    size_t boundaryLeft = labels[n]->right + 1;
                    
                    size_t length = (boundaryRight > boundaryLeft ? boundaryRight - boundaryLeft + 1 : 0);
                    if (length < upstrLength)
                        skip = true;
                }
                else {
                    // boundaryLeft is either start of sequence (if this is the first label), or the right-end of the gene to the "left" of the current gene
                    size_t boundaryLeft = (n == 0 ? 0 : labels[n-1]->right+1);
                    size_t boundaryRight = labels[n]->left - 1;
                    
                    size_t length = (boundaryRight > boundaryLeft ? boundaryRight - boundaryLeft + 1 : 0);
                    if (length < upstrLength)
                        skip = true;
                }
            }
            
            if (!useAll && !use[n])
                skip = true;
            
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
            //throw out_of_range("Not enough sequence for position " + to_string(label.left) + " to extract upstream of length " + to_string(upstrLength));
            throw out_of_range("Not enough sequence for upstream.");
        
        size_t left = label.left - upstrLength;                     // left idx of upstream sequence
        size_t right = label.left-1;                                // right idx of upstream sequence
        size_t length = right - left + 1;                           // length of upstream sequence
        
        NumSequence upstream = sequence.subseq(left, length);       // get upstream subsequence
        return upstream;
    }
    else {                                      // negative strand; reverse complement
        
        if (label.right + upstrLength >= sequence.size())
            //throw out_of_range("Not enought sequence for position " + to_string(label.right) + " on negative strand to extract upstream of length " + to_string(upstrLength));
            throw out_of_range("Not enough sequence for position.");
        
        size_t left = label.right + 1;                              // left idx of upstream sequence
        size_t right = label.right + upstrLength;                   // right idx of upstream sequence
        size_t length = right - left + 1;                           // length of upstream sequence
        
        NumSequence upstream = sequence.subseq(left, length);       // get upstream subsequence
        upstream.reverseComplement(cnc);                            // reverse complement
        return upstream;
    }
}







void SequenceParser::extractStartContextSequences(const NumSequence& sequence, const vector<Label*> &labels, const CharNumConverter &cnc, long long posRelToStart, NumSequence::size_type length, vector<NumSequence> &contexts, const vector<bool> &use) {
    
    contexts.clear();
    
    bool useAll = true;
    if (use.size() == labels.size())
        useAll = false;
    
    // for each label
    for (size_t n = 0; n < labels.size(); n++) {
        
        // check to skip
        if (!useAll && !use[n])
            continue;
        
        size_t left = labels[n]->left;
        size_t right = labels[n]->right;
        Label::strand_t strand = labels[n]->strand;
        
        if (right >= sequence.size())
            throw out_of_range("Label 'right' larger than sequence length");
        
        // positive strand
        if (strand == Label::POS) {
            
            // if left of fragment doesn't reach past the beginning or end of "sequence"
            if ((posRelToStart < 0 && abs(posRelToStart) <= left) || (posRelToStart >= 0 && posRelToStart + left < sequence.size())) {
                // left and right of fragment
                size_t fragLeft = left + posRelToStart;
                size_t fragRight = fragLeft + length-1;
                size_t fragLength = fragRight - fragLeft + 1;
                // if the right of fragment doesn't reach past the end of "sequence"
                if (fragRight < sequence.size()) {
                    contexts.push_back(sequence.subseq(fragLeft,fragLength));
                }
            }
        }
        else if (strand == Label::NEG) {
            
            // if left of fragment doesn't reach past the beginning or end of "sequence"
            if ((posRelToStart >= 0 && abs(posRelToStart) <= right) || (posRelToStart < 0 && abs(posRelToStart) + right < sequence.size())) {
                
                // left and right of fragment
                size_t fragLeft = right - posRelToStart - length + 1;
                size_t fragRight = fragLeft + length-1;
                
                if (fragRight < sequence.size()) {
                    NumSequence subseq = sequence.subseq(fragLeft, fragRight - fragLeft + 1);
                    subseq.reverseComplement(cnc);
                    contexts.push_back(subseq);
                }
            }
        }
        else
            contexts.push_back(NumSequence());
    }
}

NumSequence SequenceParser::extractStartContextSequence(const NumSequence& sequence, const Label &labels, const CharNumConverter &cnc, long long posRelToStart, NumSequence::size_type length, const vector<bool> &use) {
    throw std::logic_error("Function has no implementation");
}





















