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
#include "NonUniformMarkov.hpp"
#include "MotifMarkov.hpp"
#include "CodingMarkov.hpp"
#include "UnivariatePDF.hpp"
#include <assert.h>
#include <iostream>

using namespace std;
using namespace gmsuite;

typedef struct {
    boost::shared_ptr<NonCodingMarkov> nonCodingMarkov;
    boost::shared_ptr<CodingMarkov> codingMarkov;
    boost::shared_ptr<MotifMarkov> rbsMarkov;
    boost::shared_ptr<MotifMarkov> promoterMarkov;
    boost::shared_ptr<UnivariatePDF> rbsSpacer;
    boost::shared_ptr<UnivariatePDF> promSpacer;
} Models;


double computeConfigurationScore(const NumSequence &numSeq, size_t labelLeft, size_t labelRight, Label::strand_t strand, size_t windowDownstream, size_t windowUpstream, const CharNumConverter &cnc, const Models &models, bool printWindow=false);

ModulePostProcessor::ModulePostProcessor(const OptionsPostProcessor& opt) : options(opt) {
    
}


size_t roundUpModulo3(size_t x) {
    size_t mod = 3;
    
    if (x == 0)
        return 0;
    
//    if (x < mod)
//        return (x + x % mod) % mod;
    
    return x + (x + x % mod) % mod;
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
    boost::shared_ptr<NonCodingMarkov> nonCodingMarkov ( new NonCodingMarkov(nonCodingProbs, numAlph, cnc));;
    
    
    
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
    boost::shared_ptr<CodingMarkov> codingMarkov ( new CodingMarkov(codingProbs, numAlph, cnc));
    
    
    // check if models enabled
    bool rbsEnabled = false, promoterEnabled = false;
    
    boost::shared_ptr<MotifMarkov> rbsMarkov, promoterMarkov;
    boost::shared_ptr<UnivariatePDF> rbsSpacerPDF, promoterSpacerPDF;
    
    // RBS enabled?
    if (keyValPair.count("RBS") > 0 && keyValPair["RBS"] == "1")
        rbsEnabled = true;
    
    if (keyValPair.count("PROMOTER") > 0 && keyValPair["PROMOTER"] == "1")
        promoterEnabled = true;
    
    istringstream motifEnabled(keyValPair["RBS"]);
    
    // RBS model
    if (rbsEnabled) {
    
        istringstream ssmMotifWidth(keyValPair["RBS_WIDTH"]);    // FIXME: any motif type
        std::getline(ssmMotifWidth, line);
        int motifWidth = (int) strtol(line.c_str(), NULL, 10);

        
        vector<vector<pair<string, double> > > rbsProbs (motifWidth);
        istringstream ssmRBS(keyValPair["RBS_MAT"]);
        line = "";
        while (std::getline(ssmRBS, line)) {
            stringstream lineStream (line);
            string nt;
            vector<double> prob (motifWidth);
            
            lineStream >> nt;
            
            for (int i = 0; i < motifWidth; i++) {
                lineStream >> prob[i];
                rbsProbs[i].push_back(pair<string, double> (nt, prob[i]));
            }
        }
        

        rbsMarkov = boost::shared_ptr<MotifMarkov>  (new MotifMarkov(rbsProbs, (size_t) motifWidth, numAlph, cnc));
        
        // rbs spacer
        istringstream rbsDurMax(keyValPair["RBS_MAX_DUR"]);        // FIXME: any motif spacer
        std::getline(rbsDurMax, line);
        int maxRBSDur = (int) strtol(line.c_str(), NULL, 10);
        int rbsSpacerLen = maxRBSDur+1;
        
        istringstream rbsDur(keyValPair["RBS_POS_DISTR"]);
        vector<double> rbsPositionProbs (rbsSpacerLen, 0);
        line = "";
        while (std::getline(rbsDur, line)) {
            stringstream lineStream (line);
            int pos;
            double prob;
            
            lineStream >> pos >> prob;
            
            rbsPositionProbs[pos] = prob;
        }
        rbsSpacerPDF = boost::shared_ptr<UnivariatePDF> (new UnivariatePDF(rbsPositionProbs, false));
    }
    
    if (promoterEnabled) {
        
        if (keyValPair.count("PROMOTER_MAX_DUR") > 0) {
            istringstream promDurMax(keyValPair["PROMOTER_MAX_DUR"]);
            std::getline(promDurMax, line);
            int maxPromDur = (int) strtol(line.c_str(), NULL, 10);
            int promSpacerLen = maxPromDur + 1;
            
            istringstream promDur (keyValPair["PROMOTER_POS_DISTR"]);
            vector<double> promPositionsProbs (promSpacerLen, 0);
            line = "";
            while(std::getline(promDur, line)) {
                stringstream lineStream (line);
                int pos;
                double prob;
                
                lineStream >> pos >> prob;
                promPositionsProbs[pos] = prob;
            }
            promoterSpacerPDF = boost::shared_ptr<UnivariatePDF>  (new UnivariatePDF(promPositionsProbs, false));
        }
        
        // Promoter Model
        if (keyValPair.count("PROMOTER_WIDTH") < 0) {
            istringstream ssmPromWidth(keyValPair["PROMOTER_WIDTH"]);
            std::getline(ssmPromWidth, line);
            int promWidth = (int) strtol(line.c_str(), NULL, 10);
            
            vector<vector<pair<string, double> > > promProbs (promWidth);
            istringstream ssmProm(keyValPair["PROMOTER_MAT"]);
            line = "";
            while (std::getline(ssmProm, line)) {
                stringstream lineStream (line);
                string nt;
                vector<double> prob (promWidth);
                
                lineStream >> nt;
                
                for (int i = 0; i < promWidth; i++) {
                    lineStream >> prob[i];
                    promProbs[i].push_back(pair<string, double> (nt, prob[i]));
                }
            }
            
            promoterMarkov = boost::shared_ptr<MotifMarkov>  (new MotifMarkov(promProbs, (size_t) promWidth, numAlph, cnc));
        }
    }
    
    Models models;
    models.nonCodingMarkov = nonCodingMarkov;
    models.codingMarkov = codingMarkov;
    models.rbsMarkov = rbsMarkov;
    models.promoterMarkov = promoterMarkov;
    models.rbsSpacer = rbsSpacerPDF;
    models.promSpacer = promoterSpacerPDF;
    
    // make sure search boundaries follow in-frame ("round" up)
//    size_t neighUpstrInFrame = options.neighborhoodUpstream + (3 - options.neighborhoodUpstream%3);
//    size_t neighDownstrInFrame = options.neighborhoodDownstream + (3 - options.neighborhoodDownstream%3);
    
    size_t neighUpstrInFrame = roundUpModulo3(options.neighborhoodUpstream);
    size_t neighDownstrInFrame = roundUpModulo3(options.neighborhoodDownstream);
   
    if (options.neighborhoodUpstream == 0)
        neighUpstrInFrame = 0;
    if (options.neighborhoodDownstream == 0)
        neighDownstrInFrame = 0;

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
            
            
            size_t searchLeftBoundary = labelLeft % 3;
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
                        
                        
                        score = computeConfigurationScore(numSeq, newLabelLeft, labelRight, strand, options.windowDownstream, options.windowUpstream, cnc, models, options.printWindow);
                        
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
                        
                        score = computeConfigurationScore(numSeq, labelLeft, newLabelRight, strand, options.windowDownstream, options.windowUpstream, cnc, models, options.printWindow);
                        
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
        cout << bestLabel->toString(true) << "\t" << bestScore << endl;
        
    }
    
}


pair<NumSequence::size_type, double> findBestMotifLocation(NumSequence::const_iterator begin, NumSequence::const_iterator end, boost::shared_ptr<const MotifMarkov> motifMarkov, boost::shared_ptr<const UnivariatePDF> motifSpacer) {
    
    NumSequence::size_type bestLoc = 0;
    double bestScore = -numeric_limits<double>::infinity();
    NumSequence::size_type motifWidth = motifMarkov->getLength();
    NumSequence::size_type fragLen = std::distance(begin, end);
    
    // if motif can't fit in fragment, return -infinity score
    if (fragLen < motifWidth)
        return pair<NumSequence::size_type, double> (NumSequence::npos, bestScore);
    
    // loop over all valid positions from motif
    NumSequence::size_type currLoc = 0;
    for (NumSequence::const_iterator currElement = begin; currElement < end-motifWidth; currElement++, currLoc++) {
        
        double score = motifMarkov->evaluate(currElement, currElement+motifWidth, true);
        
        NumSequence::size_type posFromRight = fragLen - motifWidth - currLoc;
        score += log2(motifSpacer->operator[](posFromRight));
        
        if (score > bestScore) {
            bestScore = score;
            bestLoc = currLoc;
        }
    }
    
    return pair<NumSequence::size_type, double> (bestLoc, bestScore);
}


// compute the value of the gene configuration
double computeConfigurationScore(const NumSequence &numSeq, size_t labelLeft, size_t labelRight, Label::strand_t strand, size_t windowDownstream, size_t windowUpstream, const CharNumConverter &cnc, const Models &models, bool printWindow) {
    
    
    boost::shared_ptr<NonCodingMarkov> nonCodingMarkov = models.nonCodingMarkov;
    boost::shared_ptr<CodingMarkov> codingMarkov = models.codingMarkov;
    
    double negInf = -numeric_limits<double>::infinity();
    
    
    double score = 0;
    
    if (numSeq.size() > 0) {
        pair<NumSequence::size_type, double> bestRBS (NumSequence::npos, negInf);
        pair<NumSequence::size_type, double> bestProm (NumSequence::npos, negInf);
        
        if (strand == Label::POS) {
            
            size_t windowLeft = 0;
            size_t windowRight = numSeq.size()-1;
            size_t rbsWindowLeft = 0;
            size_t promWindowLeft = 0;
            
            if (labelLeft >= windowUpstream)
                windowLeft = labelLeft - windowUpstream;
            if (labelLeft + 2 + windowDownstream < numSeq.size())
                windowRight = labelLeft + 2 + windowDownstream;
            if (models.rbsSpacer && labelLeft >= models.rbsSpacer->size())
                rbsWindowLeft = labelLeft - models.rbsSpacer->size();
            if (models.promSpacer && labelLeft >= models.promSpacer->size())
                promWindowLeft = labelLeft - models.promSpacer->size();
            
            if (models.rbsMarkov) {
                bestRBS = findBestMotifLocation(numSeq.begin() + rbsWindowLeft, numSeq.begin()+labelLeft, models.rbsMarkov, models.rbsSpacer);
            }
            if (models.promoterMarkov) {
                bestProm = findBestMotifLocation(numSeq.begin() + promWindowLeft, numSeq.begin()+labelLeft, models.promoterMarkov, models.promSpacer);
            }
            
            if (printWindow) {
                // check if length of sequence is length of windows
                size_t windowLength = windowDownstream + windowUpstream + 3;
                size_t fragLength = windowRight - windowLeft + 1;
                if (fragLength == windowLength) {
            
                    cout << ">" << labelLeft+1 << "\t" << labelRight+1 << "\t" << (strand == Label::POS ? "+" : "-") << endl;
                    NumSequence tmpFrag = numSeq.subseq(windowLeft, fragLength);
                    cout << cnc.convert(numSeq.begin()+windowLeft, numSeq.begin()+windowLeft+windowLength) << endl;
                }
                
            }
            
            pair<NumSequence::size_type, double> bestMotif (NumSequence::npos, -numeric_limits<double>::infinity());
            NumSequence::size_type bestMotifWidth = 0;
            
            if (bestRBS.second != negInf || bestProm.second != negInf){
                if (bestRBS.second > bestProm.second) {
                    bestMotif = bestRBS;
                    bestMotifWidth = models.rbsMarkov->getLength();
                }
                else {
                    bestMotif = bestProm;
                    bestMotifWidth = models.promoterMarkov->getLength();
                }
            }
            
            
            if (nonCodingMarkov) {
                // if motif exists
                if (bestMotif.second > -numeric_limits<double>::infinity()) {
                    score += nonCodingMarkov->evaluate(numSeq.begin()+windowLeft, numSeq.begin()+windowLeft+bestMotif.first, true);
                    score += bestMotif.second;
                    score += nonCodingMarkov->evaluate(numSeq.begin()+windowLeft+bestMotif.first+bestMotifWidth, numSeq.begin() + labelLeft, true);
                    // FIXME: add remaining noncoding?
                }
                // else assume no motif
                else {
                    score += nonCodingMarkov->evaluate(numSeq.begin()+windowLeft, numSeq.begin()+labelLeft, true);
                }
            }
            if (codingMarkov)
                score += codingMarkov->evaluate(numSeq.begin()+labelLeft+3, numSeq.begin()+windowRight+1, true);
        }
        // negative strand
        else if (strand == Label::NEG) {
            size_t windowLeft = 0;
            size_t windowRight = numSeq.size()-1;
            size_t rbsWindowLeft = 0;
            size_t promWindowLeft = 0;
            
            if (labelRight >= windowDownstream+3)
                windowLeft = labelRight - 2 - windowDownstream;
            if (labelRight + 2 + windowUpstream < numSeq.size())
                windowRight = labelRight + windowUpstream;
            
            size_t nonCodingLen = windowRight - labelRight;
            
            NumSequence revFrag = numSeq.subseq(windowLeft, windowRight-windowLeft+1);
            revFrag.reverseComplement(cnc);
            
            if (printWindow) {
                // check if length of sequence is length of windows
                size_t windowLength = windowDownstream + windowUpstream + 3;
                size_t fragLength = windowRight - windowLeft + 1;
                if (fragLength == windowLength) {
                    cout << ">" << labelLeft+1 << "\t" << labelRight+1 << "\t" << (strand == Label::POS ? "+" : "-") << endl;
                    cout << cnc.convert(revFrag.begin(), revFrag.end()) << endl;
                }
                
            }
            
            if (models.rbsSpacer && nonCodingLen >= models.rbsSpacer->size())
                rbsWindowLeft = nonCodingLen - models.rbsSpacer->size();
            if (models.promSpacer && nonCodingLen >= models.promSpacer->size())
                promWindowLeft = nonCodingLen - models.promSpacer->size();

            
            if (models.rbsMarkov) {
                bestRBS = findBestMotifLocation(revFrag.begin()+rbsWindowLeft, revFrag.begin()+nonCodingLen, models.rbsMarkov, models.rbsSpacer);
            }
            if (models.promoterMarkov) {
                bestProm = findBestMotifLocation(revFrag.begin()+promWindowLeft, revFrag.begin()+nonCodingLen, models.promoterMarkov, models.promSpacer);
            }
            
            pair<NumSequence::size_type, double> bestMotif (NumSequence::npos, -numeric_limits<double>::infinity());
            
            NumSequence::size_type bestMotifWidth = 0;
            
            if (bestRBS.second != negInf || bestProm.second != negInf) {
                if (bestRBS.second > bestProm.second) {
                    bestMotif = bestRBS;
                    bestMotifWidth = models.rbsMarkov->getLength();
                }
                else {
                    bestMotif = bestProm;
                    bestMotifWidth = models.promoterMarkov->getLength();
                }
            }
            
            
            if (nonCodingMarkov) {
                // if motif exists
                if (bestMotif.second > -numeric_limits<double>::infinity()) {
                    score += nonCodingMarkov->evaluate(revFrag.begin(), revFrag.begin()+bestMotif.first, true);
                    score += bestMotif.second;
                    score += nonCodingMarkov->evaluate(revFrag.begin()+bestMotif.first+bestMotifWidth, revFrag.end(), true);
                }
                // else assume no motif
                else {
                    score += nonCodingMarkov->evaluate(revFrag.begin(), revFrag.begin() + nonCodingLen, true);
                }
                
            }
            if (codingMarkov)
                score += codingMarkov->evaluate(revFrag.begin()+nonCodingLen+3, revFrag.end(), true);
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
