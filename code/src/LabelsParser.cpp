//
//  LabelsParser.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "LabelsParser.hpp"

using namespace gmsuite;

void LabelsParser::partitionBasedOnOperonStatus(const vector<Label*> &labels, size_t fgioThresh, size_t nfgioThresh,
                                                vector<operon_status_t> &status) {
    
    
    status.resize(labels.size(), AMBIG);            // allocate space, and set all as initially ambiguous
    
    
    // loop over all labels
    for (size_t n = 0; n < labels.size(); n++) {
        
        Label* currLab = labels[n];         // current label
        
        // If current gene on positive strand
        if (currLab->strand == Label::POS) {
            
            // if first gene, set as fgio
            if (n == 0)
                status[n] = FGIO;
            
            // if not first gene in genome
            else {
                Label* prevLab = labels[n-1];
                
                // if previous gene is on negative strand, then this is a fgio
                if (prevLab->strand == Label::NEG)
                    status[n] = FGIO;
                
                // otherwise previous gene is on positive strand: check distance
                else {
                    // if genes overlapping, then set as nfgio
                    if (prevLab->right > currLab->left)
                        status[n] = NFGIO;
                    // if genes not overlapping but too close, set as nfgio
                    else if (currLab->left - prevLab->right < nfgioThresh)
                        status[n] = NFGIO;
                    // if distance is "too far" then fgio
                    else if (currLab->left - prevLab->right > fgioThresh)
                        status[n] = FGIO;
                }
            }
            
        }
        // If current gene on negative strand
        else {
            
            // if "first gene in genome" from right, then set as fgio
            if (n == labels.size() - 1)
                status[n] = FGIO;
            
            // if not first gene in genome
            else {
                Label* prevLab = labels[n+1];
                
                // if previous gene is on positive strand, then this is a fgio
                if (prevLab->strand == Label::POS)
                    status[n] = FGIO;
                
                // otherwise, previous gene is on negative strand; check distance
                else {
                    // if genes overlapping, then set as nfgio
                    if (prevLab->left < currLab->right)
                        status[n] = NFGIO;
                    // if genes not overlapping but too close, set as nfgio
                    else if (prevLab->left - currLab->right < nfgioThresh)
                        status[n] = NFGIO;
                    // if distance is too far, then fgio
                    else if (prevLab->left - currLab->right > fgioThresh)
                        status[n] = FGIO;
                }
            }
        }
    }
}



void LabelsParser::splitBasedOnPartition(const vector<Label*> &labels, const vector<operon_status_t> &status, vector<Label*> &labelsFGIO, vector<Label*> &labelsIGIO, vector<Label*> &labelsAMBIG ) {
    
    size_t numFGIO = 0, numIG = 0, numUNK = 0;
    for (size_t n = 0; n < status.size(); n++) {
        if (status[n] == LabelsParser::FGIO)        numFGIO++;
        else if (status[n] == LabelsParser::NFGIO)  numIG++;
        else
            numUNK++;
    }
    
    // get FGIO and IG upstreams and run motif search
    labelsFGIO.resize(numFGIO);
    labelsIGIO.resize(numIG);
    size_t currFGIO = 0, currIG = 0;        // indices
    
    for (size_t n = 0; n < status.size(); n++) {
        if (status[n] == LabelsParser::FGIO)        labelsFGIO[currFGIO++] = labels[n];
        else if (status[n] == LabelsParser::NFGIO)  labelsIGIO[currIG++] = labels[n];
    }
    
}


