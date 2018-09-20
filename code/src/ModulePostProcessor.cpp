//
//  ModulePostProcessor.cpp
//  Biogem CPP
//
//  Created by Karl Gemayel on 9/18/18.
//  Copyright © 2018 Georgia Institute of Technology. All rights reserved.
//

#include "ModulePostProcessor.hpp"

#include "ModelFile.hpp"
#include "LabelFile.hpp"
#include "SequenceFile.hpp"
#include "NumSequence.hpp"
#include "NumGeneticCode.hpp"
#include "ModelFile.hpp"
#include "NonCodingMarkov.hpp"
#include "CodingMarkov.hpp"
#include <assert.h>
#include <iostream>

using namespace std;
using namespace gmsuite;

double computeConfigurationScore(const NumSequence &numSeq, size_t labelLeft, size_t labelRight, Label::strand_t strand, size_t windowDownstream, size_t windowUpstream, const CharNumConverter &cnc, const NonCodingMarkov &nonCodingMarkov, const CodingMarkov &codingMarkov);

ModulePostProcessor::ModulePostProcessor(const OptionsPostProcessor& opt) : options(opt) {
    
}

void ModulePostProcessor::run() {
    
    AlphabetDNA alph;
    GeneticCode geneticCode (options.gcode);
    CharNumConverter cnc(&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    NumGeneticCode numGeneticCode(geneticCode, cnc);
    
    // read sequence from file
    SequenceFile seqFile(options.fnsequence, SequenceFile::READ);
    Sequence strSeq = seqFile.read();
    NumSequence numSeq (strSeq, cnc);
    
    // read labels from file
    LabelFile labelFile (options.fnlabels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // read model file
    ModelFile modFile(options.fnmod, ModelFile::READ);
    map<string, string> keyValPair;
    modFile.read(keyValPair);
    
    // construct models: non-coding
    vector<pair<string, double> > nonCodingProbs;
    istringstream ssm(keyValPair["NON_MAT"]);
    string line;
    while (std::getline(ssm, line)) {
        
        stringstream lineStream (line);
        string codon;
        double prob;
        
        lineStream >> codon >> prob;
        
        nonCodingProbs.push_back(pair<string, double> (codon, prob));
    }
    
    // build non-coding model
    NonCodingMarkov nonCodingMarkov(nonCodingProbs, numAlph, cnc);
    
    // construct models: coding
    vector<vector<pair<string, double> > > codingProbs(3);
    istringstream ssmCoding(keyValPair["COD_MAT"]);
    line = "";
    while (std::getline(ssmCoding, line)) {
        
        stringstream lineStream (line);
        string codon;
        double prob1;
        double prob2;
        double prob3;
        
        lineStream >> codon >> prob1 >> prob2 >> prob3;
        
        codingProbs[0].push_back(pair<string, double> (codon, prob1));
        codingProbs[1].push_back(pair<string, double> (codon, prob2));
        codingProbs[2].push_back(pair<string, double> (codon, prob3));
    }
    
    // build coding model
    CodingMarkov codingMarkov(codingProbs, numAlph, cnc);
    
    // make sure search boundaries follow in-frame ("round" up)
    size_t neighUpstrInFrame = options.neighborhoodUpstream + (3 - options.neighborhoodUpstream%3);
    size_t neighDownstrInFrame = options.neighborhoodDownstream + (3 - options.neighborhoodDownstream%3);
   

    // for each label
    for (vector<Label*>::const_iterator iter = labels.begin(); iter!= labels.end(); iter++) {
        // get search window around labelled start
        size_t labelLeft = (*iter)->left;
        size_t labelRight = (*iter)->right;
        Label::strand_t strand = (*iter)->strand;
        
        double score = 0;
        
        Label* bestLabel = new Label(*(*iter)); // copy label
        double bestScore = -numeric_limits<double>::infinity();
        
        if (numSeq.size() > 0) {
            
            
            size_t searchLeftBoundary = (abs((int)labelLeft-3)%3);
            size_t searchRightBoundary = numSeq.size()-1 - ((numSeq.size()-1 - labelRight)%3) ;
            
            if (strand == Label::POS) {
                
                if (labelLeft >= neighUpstrInFrame)
                    searchLeftBoundary = labelLeft - neighUpstrInFrame;
                if (labelLeft + 2 + neighDownstrInFrame < numSeq.size())
                    searchRightBoundary = labelLeft + 2 + neighDownstrInFrame;
                
                for (NumSequence::size_type n = searchLeftBoundary; n <= searchRightBoundary-2; n+=3) {
                    
                    NumSequence codonNumSeq = numSeq.subseq(n, 3);
                    CharNumConverter::seq_t codon = CharNumConverter::seq_t (codonNumSeq.begin(), codonNumSeq.end());

                    if (numGeneticCode.isStart(codon)) {
                        
                        size_t newLabelLeft = n;
                        
                        assert(abs((int)n - (int)labelLeft + 1)%3 );
                        

                        score = computeConfigurationScore(numSeq, newLabelLeft, labelRight, strand, options.windowDownstream, options.windowUpstream, cnc, nonCodingMarkov, codingMarkov);
                        
                        if (score > bestScore) {
                            bestLabel->left = newLabelLeft;
                            bestScore = score;
                        }
                    }
                    
                }
                
            }
            else if (strand == Label::NEG) {
                if (labelRight >= neighDownstrInFrame+3)
                    searchLeftBoundary = labelRight - 2 - neighDownstrInFrame;
                if (labelRight + neighUpstrInFrame < numSeq.size())
                    searchRightBoundary = labelRight + neighUpstrInFrame;
                
                for (NumSequence::size_type n = searchLeftBoundary+2; n <= searchRightBoundary; n+=3) {
                    NumSequence codonNumSeqRevComp = numSeq.subseq(n-2, 3);
                    codonNumSeqRevComp.reverseComplement(cnc);
                    CharNumConverter::seq_t codon = CharNumConverter::seq_t (codonNumSeqRevComp.begin(), codonNumSeqRevComp.end());
                    
                    if (numGeneticCode.isStart(codon)) {
                        size_t newLabelRight = n;
                        
                        assert(abs((int)n - (int)labelRight + 1)%3 );
                        
                        score = computeConfigurationScore(numSeq, labelLeft, newLabelRight, strand, options.windowDownstream, options.windowUpstream, cnc, nonCodingMarkov, codingMarkov);
                        
                        if (score > bestScore) {
                            bestLabel->right = newLabelRight;
                            bestScore = score;
                        }
                    }
                }
            }
            else {
                throw logic_error("Invalid strand value.");
            }
            
            
            
            
            
            
        
//            if (strand == Label::POS) {
//
//                size_t windowLeft = 0;
//                size_t windowRight = numSeq.size()-1;
//
//                if (labelLeft >= options.windowUpstream)
//                    windowLeft = labelLeft - options.windowUpstream;
//                if (labelLeft + 2 + options.windowDownstream < numSeq.size())
//                    windowRight = labelLeft + 2 + options.windowDownstream;
//
//
//
//                score += nonCodingMarkov.evaluate(numSeq.begin()+windowLeft, numSeq.begin()+labelLeft, true);
//                score += codingMarkov.evaluate(numSeq.begin()+labelLeft+3, numSeq.begin()+windowRight+1);
//            }
//            // negative strand
//            else if (strand == Label::NEG) {
//                size_t windowLeft = 0;
//                size_t windowRight = numSeq.size()-1;
//
//                if (labelRight >= options.windowDownstream+3)
//                    windowLeft = labelRight - 3 - options.windowDownstream;
//                if (labelLeft + 2 + options.windowUpstream < numSeq.size())
//                    windowRight = labelRight + 1 + options.windowUpstream;
//
//                size_t nonCodingLen = windowRight - labelRight;
//
//                NumSequence revFrag = numSeq.subseq(windowLeft, windowRight-windowLeft+1);
//                revFrag.reverseComplement(cnc);
//
//                score += nonCodingMarkov.evaluate(revFrag.begin(), revFrag.begin() + nonCodingLen);
//                score += codingMarkov.evaluate(revFrag.begin()+nonCodingLen+3, revFrag.end());
//            }
//            // no strand value
//            else {
//                throw logic_error("Invalid strand value.");
//            }
        }
        
        // print score
        cout << bestLabel->toString(true) << endl;
        
    }
    
}




// compute the value of the gene configuration
double computeConfigurationScore(const NumSequence &numSeq, size_t labelLeft, size_t labelRight, Label::strand_t strand, size_t windowDownstream, size_t windowUpstream, const CharNumConverter &cnc, const NonCodingMarkov &nonCodingMarkov, const CodingMarkov &codingMarkov) {
    
    double score = 0;
    
    if (numSeq.size() > 0) {
        
        if (strand == Label::POS) {
            
            size_t windowLeft = 0;
            size_t windowRight = numSeq.size()-1;
            
            if (labelLeft >= windowUpstream)
                windowLeft = labelLeft - windowUpstream;
            if (labelLeft + 2 + windowDownstream < numSeq.size())
                windowRight = labelLeft + 2 + windowDownstream;
            
            
            
            score += nonCodingMarkov.evaluate(numSeq.begin()+windowLeft, numSeq.begin()+labelLeft, true);
            score += codingMarkov.evaluate(numSeq.begin()+labelLeft+3, numSeq.begin()+windowRight+1, true);
        }
        // negative strand
        else if (strand == Label::NEG) {
            size_t windowLeft = 0;
            size_t windowRight = numSeq.size()-1;
            
            if (labelRight >= windowDownstream+3)
                windowLeft = labelRight - 3 - windowDownstream;
            if (labelLeft + 2 + windowUpstream < numSeq.size())
                windowRight = labelRight + 1 + windowUpstream;
            
            size_t nonCodingLen = windowRight - labelRight;
            
            NumSequence revFrag = numSeq.subseq(windowLeft, windowRight-windowLeft+1);
            revFrag.reverseComplement(cnc);
            
            score += nonCodingMarkov.evaluate(revFrag.begin(), revFrag.begin() + nonCodingLen, true);
            score += codingMarkov.evaluate(revFrag.begin()+nonCodingLen+3, revFrag.end(), true);
        }
        // no strand value
        else {
            throw logic_error("Invalid strand value.");
        }
    }
    else {
        return - numeric_limits<double>::infinity();
    }
    
    return score;
}
