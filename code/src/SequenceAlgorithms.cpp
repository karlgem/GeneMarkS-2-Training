//
//  SequenceAlgorithms.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/4/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "SequenceAlgorithms.hpp"
#include "LabelsParser.hpp"

using namespace gmsuite;

NumSequence SequenceAlgorithms::longestCommonSubstring(const NumSequence &A, const NumSequence &B, const std::vector<std::pair<NumSequence::num_t, NumSequence::num_t> >& subs) {
    
    // allocate space for LCS matrix
    int** LCS = new int*[A.size()+1];
    for (size_t i = 0; i <= A.size(); i++)
        LCS[i] = new int[B.size()+1];
    
    // if A is empty, LCS of A,B=0
    for (size_t i = 0; i <= B.size(); i++)
        LCS[0][i] = 0;
    
    // if B is empty, LCS of A,B = 0
    for (size_t i = 0; i <= A.size(); i++)
        LCS[i][0] = 0;
    
    // fill the rest of the matrix
    for (size_t i = 1; i <= A.size(); i++) {
        for (size_t j = 1; j <= B.size(); j++) {
            // match
            bool match = (A[i-1] == B[j-1]);
            for (size_t n = 0; n < subs.size(); n++) {
                if (A[i-1] == subs[n].first && B[j-1] == subs[n].second)
                    match |= true;
            }
            
            if (match) {
                LCS[i][j] = LCS[i-1][j-1] + 1;
            }
            // mismatch
            else {
                LCS[i][j] = 0;
            }
        }
    }
    
    size_t maxSubstringSize = 0;
    size_t posOfMaxInA = 0;
    size_t posOfMaxInB = 0;
    // find the maximum element in the matrix
    for (size_t i = 1; i <= A.size(); i++) {
        for (size_t j = 1; j <= B.size(); j++) {
            if (maxSubstringSize < LCS[i][j]) {
                maxSubstringSize = LCS[i][j];
                posOfMaxInA = i-1;
                posOfMaxInB = j-1;
            }
        }
    }
    
    // extract substring from either sequence
    vector<NumSequence::num_t> result (maxSubstringSize);
    
    for (size_t i = 0; i < maxSubstringSize; i++) {
        result[i] = A[posOfMaxInA - (maxSubstringSize-1) + i];
    }
    
    // delete LCS allocation
    for (size_t i = 0; i < A.size(); i++)
        delete [] LCS[i];
    delete [] LCS;
    
    
    return NumSequence(result);
}



NumSequence SequenceAlgorithms::longestMatchTo16S(const NumSequence &A, const NumSequence &B,
                                                  std::pair<NumSequence::size_type, NumSequence::size_type>& positionsOfMatches,
                                                  const std::vector<std::pair<NumSequence::num_t, NumSequence::num_t> >& subs) {
    
    // allocate space for LCS matrix
    int** LCS = new int*[A.size()+1];
    for (size_t i = 0; i <= A.size(); i++)
        LCS[i] = new int[B.size()+1];
    
    // if A is empty, LCS of A,B=0
    for (size_t i = 0; i <= B.size(); i++)
        LCS[0][i] = 0;
    
    // if B is empty, LCS of A,B = 0
    for (size_t i = 0; i <= A.size(); i++)
        LCS[i][0] = 0;
    
    // fill the rest of the matrix
    for (size_t i = 1; i <= A.size(); i++) {
        for (size_t j = 1; j <= B.size(); j++) {
            // match
            bool match = (A[i-1] == B[j-1]);
            for (size_t n = 0; n < subs.size(); n++) {
                if (A[i-1] == subs[n].first && B[j-1] == subs[n].second)
                    match |= true;
            }
            
            if (match) {
                LCS[i][j] = LCS[i-1][j-1] + 1;
            }
            // mismatch
            else {
                LCS[i][j] = 0;
            }
        }
    }
    
    size_t maxSubstringSize = 0;
    size_t posOfMaxInA = 0;
    size_t posOfMaxInB = 0;
    // find the maximum element in the matrix
    for (size_t i = 1; i <= A.size(); i++) {
        for (size_t j = 1; j <= B.size(); j++) {
            if (maxSubstringSize < LCS[i][j]) {
                maxSubstringSize = LCS[i][j];
                posOfMaxInA = i-1;
                posOfMaxInB = j-1;
            }
        }
    }
    
    // extract substring from either sequence
    vector<NumSequence::num_t> result (maxSubstringSize);
    
    for (size_t i = 0; i < maxSubstringSize; i++) {
        result[i] = B[posOfMaxInB - (maxSubstringSize-1) + i];
    }
    
    // set matched positions
    positionsOfMatches.first = posOfMaxInA - (maxSubstringSize-1);
    positionsOfMatches.second = posOfMaxInB - (maxSubstringSize-1);
    
    // delete LCS allocation
    for (size_t i = 0; i < A.size(); i++)
        delete [] LCS[i];
    delete [] LCS;
    
    
    return NumSequence(result);
}



// compute GC for entire sequence
double SequenceAlgorithms::computeGC(const Sequence &seq) {
    
    size_t numGC = 0;
    for (size_t n = 0; n < seq.size(); n++) {
        if (seq[n] == 'G' || seq[n] == 'C')
            numGC++;
    }
    
    if (seq.size() == 0)
        return 0;
    
    return 100 * numGC / (double) seq.size();
}

// compute GC per gene
void SequenceAlgorithms::computeGC(const Sequence &seq, const vector<Label*> &labels, vector<double> &gcs) {
    
    gcs.clear();
    
    for (size_t n = 0; n < labels.size(); n++) {
        Label* lab = labels[n];
        
        size_t left = lab->left;
        size_t right = lab->right;
    
        // check for valid coordinates
        if (left <= right && right < seq.size()) {
            
            size_t length = right - left + 1;
            double numGC = 0;
            for (size_t i = left; i <= right; i++) {
                if (seq[i] == 'G' || seq[i] == 'C')
                    numGC++;
            }
            
            double GC = 100 *  numGC / (double) length;
            gcs.push_back(GC);
        }
        else {
            gcs.push_back(0);           // GC=0 for invalid ranges (to maintain equal vector lengths of labels and gcs)
        }
    }
}














