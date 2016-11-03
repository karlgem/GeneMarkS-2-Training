//
//  ModuleExperiment.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/14/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleExperiment.hpp"

#include "AlphabetDNA.hpp"
#include "NumSequence.hpp"
#include "CharNumConverter.hpp"
#include "SequenceFile.hpp"
#include "LabelFile.hpp"
#include "SequenceParser.hpp"
#include "SequenceAlgorithms.hpp"
#include "MotifFinder.hpp"
#include "GMS2Trainer.hpp"
#include "NonUniformCounts.hpp"
#include "UniformCounts.hpp"
#include "LabelsParser.hpp"
#include "NumGeneticCode.hpp"

#include <algorithm>
#include <iostream>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace gmsuite;
using boost::shared_ptr;


ModuleExperiment::ModuleExperiment(const OptionsExperiment& opt) : options(opt) {
    
}



void ModuleExperiment::run() {
    
    if (options.experiment == OptionsExperiment::MATCH_SEQ_TO_UPSTREAM)
        runMatchSeqToUpstream();
    else if (options.experiment == OptionsExperiment::MATCH_SEQ_TO_NONCODING)
        runMatchSeqToNoncoding();
    else if (options.experiment == OptionsExperiment::BUILD_START_MODELS)
        runBuildStartModels();
    else if (options.experiment == OptionsExperiment::BUILD_START_MODELS2)
        runBuildStartModels2();
    else if (options.experiment == OptionsExperiment::BUILD_START_MODELS3)
        runBuildStartModels3();
    else if (options.experiment == OptionsExperiment::SCORE_STARTS)
        runScoreStarts();
    
}

// match query to upstream regions
void ModuleExperiment::runMatchSeqToUpstream() {
    
    OptionsExperiment::MatchSeqToUpstreamOptions expOptions = options.matchSeqToUpstream;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extractlk upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreams, expOptions.allowOverlaps, expOptions.minGeneLength);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    unsigned matchThresh = 4;            // threshold for nonmatches
    vector<NumSequence> nonMatch;           // keep track of 'non matching'
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        NumSequence sub = upstreams[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, sub);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(upstreams[i]);
    }
    
    // run motif search on non-match
    MotifFinder::Builder b;
    b.setAlign(MFinderModelParams::RIGHT);
    
    MotifFinder mfinder = b.build();
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonMatch, positions);
    
    cout << "The number of remaining sequences is: " << nonMatch.size() << endl;
    for (size_t n = 0; n < nonMatch.size(); n++) {
        cout << cnc.convert(nonMatch[n].begin() + positions[n], nonMatch[n].begin() + positions[n] + 6);
        cout << "\t" << positions[n] + 1 << "\t" << nonMatch[n].size() << endl;
    }
}

// match query to simulated noncoding sequences
void ModuleExperiment::runMatchSeqToNoncoding() {
    
    OptionsExperiment::MatchSeqToNoncodingOptions expOptions = options.matchSeqToNoncoding;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    
    NumSequence numSequence (strSequence, cnc);
    
    OptionsMFinder optionsMFinder;
    
    // run training step
    GMS2Trainer trainer((unsigned) expOptions.pcounts, 0, expOptions.order, 0, 40, 0, ProkGeneStartModel::C1, optionsMFinder, numAlph, 300);
    
    trainer.estimateParameters(numSequence, labels);
    
    // Generate non-coding sequences
    vector<NumSequence> simNonCoding (expOptions.numNoncoding);
    
    for (size_t n = 0; n < simNonCoding.size(); n++)
        simNonCoding[n] = trainer.noncoding->emit(expOptions.length);
    
    // get sequence to match with
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    // for each noncoding sequence, match it
    for (size_t i = 0; i < simNonCoding.size(); i++) {
        NumSequence match = SequenceAlgorithms::longestCommonSubstring(matchSeq, simNonCoding[i]);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
    }
        
}





// match query to upstream regions
void ModuleExperiment::runBuildStartModels() {
    
    OptionsExperiment::BuildStartModelsOptions expOptions = options.buildStartModels;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extractlk upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreams, expOptions.allowOverlaps, expOptions.minGeneLength);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    unsigned matchThresh = expOptions.min16SMatch;          // threshold for nonmatches
    vector<NumSequence> nonMatch;                           // keep track of 'non matching'
    
    vector<pair<NumSequence::size_type, NumSequence::size_type> > positionsOfMatches (upstreams.size());

    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        NumSequence sub = upstreams[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, sub, positionsOfMatches[i], substitutions);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(upstreams[i]);
    }
    
    // run motif search on non-match
    MotifFinder::Builder b;
    
    MotifFinder mfinder = b.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonMatch, positions);
    
    cout << "The number of remaining sequences is: " << nonMatch.size() << endl;
    for (size_t n = 0; n < nonMatch.size(); n++) {
        cout << cnc.convert(nonMatch[n].begin() + positions[n], nonMatch[n].begin() + positions[n] + expOptions.mfinderOptions.width);
        cout << "\t" << positions[n] + 1 << "\t" << nonMatch[n].size() << endl;
    }
}

bool isNull(Label* l) { return l == NULL; }

// match query to upstream regions
void ModuleExperiment::runBuildStartModels2() {
    
    OptionsExperiment::BuildStartModels2Options expOptions = options.buildStartModels2;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    
    if (expOptions.minGeneLength > 0) {
        for (size_t i = 0; i < labels.size(); i++) {
            size_t length = labels[i]->right - labels[i]->left + 1;
            if (length < expOptions.minGeneLength) {
                delete labels[i];
                labels[i] = NULL;
            }
        }
        
        size_t sizeBefore = labels.size();
        
        // erase from vector
        
        vector<Label*>::iterator toRem = remove_if(labels.begin(), labels.end(), isNull);
        labels.erase(toRem, labels.end());
        
        cout << "Num Filtered by Gene Length: " << sizeBefore - labels.size() << endl;
    }
    
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extractlk upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreams, expOptions.allowOverlaps, expOptions.minGeneLength);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    unsigned matchThresh = expOptions.min16SMatch;          // threshold for nonmatches
    vector<NumSequence> nonMatch;                           // keep track of 'non matching'
    vector<NumSequence> matchUpstreams;                     // keep track of matched upstreams
    
    vector<pair<NumSequence::size_type, NumSequence::size_type> > positionsOfMatches (upstreams.size());
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (expOptions.allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    // for each upstream sequence, match it against strMatchSeq
    for (size_t i = 0; i < upstreams.size(); i++) {
        NumSequence sub = upstreams[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, sub, positionsOfMatches[i], substitutions);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(sub.begin(), sub.end()) << "\t" << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << "\t" << positionsOfMatches[i].first << "\t" << positionsOfMatches[i].second << "\t" << 20 << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(upstreams[i]);
        else
            matchUpstreams.push_back(sub);
    }
    
    cout << "TOTAL = MATCHED + UNMATCHED" << endl;
    cout << upstreams.size() << " = " << matchUpstreams.size() << " + " << nonMatch.size() << endl;
    
    
    /************************************ MATCH ************************************/
    
    cout << "EXPERIMENT: MOTIF SEARCH FOR MATCHED UPSTREAMS" << endl;
    
    MotifFinder::Builder bMatch;
    MotifFinder mfinderForMatch = bMatch.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positionsForMatch;
    mfinderForMatch.findMotifs(matchUpstreams, positionsForMatch);
    
    for (size_t n = 0; n < matchUpstreams.size(); n++) {
        cout << cnc.convert(matchUpstreams[n].begin(), matchUpstreams[n].end()) << "\t";
        cout << cnc.convert(matchUpstreams[n].begin() + positionsForMatch[n], matchUpstreams[n].begin() + positionsForMatch[n] + expOptions.mfinderOptions.width);
        cout << "\t" << positionsForMatch[n] + 1 << "\t" << matchUpstreams[n].size() << endl;
    }
    
    
    /************************************ NON MATCH ************************************/
    
    cout << "EXPERIMENT: MOTIF SEARCH FOR UNMATCHED UPSTREAMS" << endl;

    // run motif search on non-match
    MotifFinder::Builder b;
    
    MotifFinder mfinder = b.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonMatch, positions);
    
    for (size_t n = 0; n < nonMatch.size(); n++) {
        cout << cnc.convert(nonMatch[n].begin(), nonMatch[n].end()) << "\t";
        cout << cnc.convert(nonMatch[n].begin() + positions[n], nonMatch[n].begin() + positions[n] + expOptions.mfinderOptions.width);
        cout << "\t" << positions[n] + 1 << "\t" << nonMatch[n].size() << endl;
    }
}




void ModuleExperiment::runBuildStartModels3() {
    
    OptionsExperiment::BuildStartModels3Options expOptions = options.buildStartModels3;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    
    if (expOptions.minGeneLength > 0) {
        for (size_t i = 0; i < labels.size(); i++) {
            size_t length = labels[i]->right - labels[i]->left + 1;
            if (length < expOptions.minGeneLength) {
                delete labels[i];
                labels[i] = NULL;
            }
        }
        
        size_t sizeBefore = labels.size();
        
        // erase from vector
        
        vector<Label*>::iterator toRem = remove_if(labels.begin(), labels.end(), isNull);
        labels.erase(toRem, labels.end());
        
        cout << "Num Filtered by Gene Length: " << sizeBefore - labels.size() << endl;
    }
    
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extractlk upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    vector<NumSequence> upstreamsForProm;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreams, expOptions.allowOverlaps, expOptions.minGeneLength);
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, expOptions.length, upstreamsForProm, false, expOptions.minGeneLength);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);
    
    unsigned matchThresh = expOptions.min16SMatch;          // threshold for nonmatches
    vector<NumSequence> nonMatch;                           // keep track of 'non matching'
    vector<NumSequence> matchUpstreams;                     // keep track of matched upstreams
    
    vector<pair<NumSequence::size_type, NumSequence::size_type> > positionsOfMatches (upstreams.size());
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (expOptions.allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));
    
    size_t upstrLen_RBSSearch = 20;
    
    vector<int> distanceToPreviousGene;
    
    // extract Set of inner-genes in operon
    vector<NumSequence> nonFGIO_upstreams;
    vector<NumSequence> FGIO_upstreams;
    
    for (size_t n = 0; n < upstreams.size(); n++) {
        if (labels[n]->strand == Label::NEG) {
            NumSequence::size_type start = labels[n]->right;
            
            // if last gene
            if (n == upstreams.size()-1) {
                if (numSequence.size() > start + upstrLen_RBSSearch)
                    nonFGIO_upstreams.push_back(upstreams[n]);
            }
            // if not last gene
            else {
                int distToPrev = (int)labels[n+1]->left - (int)start - 1;
                distanceToPreviousGene.push_back(distToPrev);
                
                // if next gene close by, include
                if (start + expOptions.nfgioThresh >= labels[n+1]->left)
                    nonFGIO_upstreams.push_back(upstreams[n]);
                else
                    FGIO_upstreams.push_back(upstreams[n]);
            }
        }
        // positive strand
        else {
            NumSequence::size_type start = labels[n]->left;
            
            // if first gene
            if (n == 0) {
                if (start > upstrLen_RBSSearch)
                    nonFGIO_upstreams.push_back(upstreams[n]);
            }
            // if not first gene
            else {
                int distToPrev = (int)start - (int)labels[n-1]->right - 1;
                distanceToPreviousGene.push_back(distToPrev);
                
                if (labels[n-1]->right + expOptions.nfgioThresh >= start)
                    nonFGIO_upstreams.push_back(upstreams[n]);
                else
                    FGIO_upstreams.push_back(upstreams[n]);
            }
        }
    }
    
    /* Print Distance To Previous Gene */
    for (size_t n = 0; n < distanceToPreviousGene.size(); n++) {
        cout << distanceToPreviousGene[n] << endl;
    }
    
    // for each upstream sequence for non-FGIO, match it against strMatchSeq
    for (size_t i = 0; i < nonFGIO_upstreams.size(); i++) {
        NumSequence sub = nonFGIO_upstreams[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, sub, positionsOfMatches[i], substitutions);
        
        // print match and size
        if (match.size() > 0)
            cout << cnc.convert(sub.begin(), sub.end()) << "\t" << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << "\t" << positionsOfMatches[i].first << "\t" << positionsOfMatches[i].second << "\t" << 20 << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(nonFGIO_upstreams[i]);
        else
            matchUpstreams.push_back(sub);
    }
    
    cout << "TOTAL = MATCHED + UNMATCHED" << endl;
    cout << upstreams.size() << " = " << matchUpstreams.size() << " + " << nonMatch.size() << endl;
    
    
    /************************************ MATCH ************************************/
    
    cout << "EXPERIMENT: MOTIF SEARCH FOR MATCHED UPSTREAMS" << endl;
    
    MotifFinder::Builder bMatch;
    MotifFinder mfinderForMatch = bMatch.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positionsForMatch;
    mfinderForMatch.findMotifs(matchUpstreams, positionsForMatch);
    
    for (size_t n = 0; n < matchUpstreams.size(); n++) {
        cout << cnc.convert(matchUpstreams[n].begin(), matchUpstreams[n].end()) << "\t";
        cout << cnc.convert(matchUpstreams[n].begin() + positionsForMatch[n], matchUpstreams[n].begin() + positionsForMatch[n] + expOptions.mfinderOptions.width);
        cout << "\t" << positionsForMatch[n] + 1 << "\t" << matchUpstreams[n].size() << endl;
    }
    
    
    /************************************ NON MATCH ************************************/
    
    nonMatch.clear();
    positionsOfMatches.clear();
    
    // for each upstream sequence for non-FGIO, match it against strMatchSeq
    for (size_t i = 0; i < upstreamsForProm.size(); i++) {
        NumSequence sub = upstreamsForProm[i].subseq(expOptions.length - 20, 20);
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, sub, positionsOfMatches[i], substitutions);
        
//        // print match and size
//        if (match.size() > 0)
//            cout << cnc.convert(sub.begin(), sub.end()) << "\t" << cnc.convert(match.begin(), match.end()) << "\t" << match.size() << "\t" << positionsOfMatches[i].first << "\t" << positionsOfMatches[i].second << "\t" << 20 << endl;
        
        // keep track of nonmatches
        if (match.size() < matchThresh)
            nonMatch.push_back(upstreamsForProm[i]);
        else
            matchUpstreams.push_back(sub);
    }
    
    cout << "EXPERIMENT: MOTIF SEARCH FOR UNMATCHED UPSTREAMS" << endl;
    
    // run motif search on non-match
    MotifFinder::Builder b;
    
    MotifFinder mfinder = b.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positions;
    mfinder.findMotifs(nonMatch, positions);
    
    for (size_t n = 0; n < nonMatch.size(); n++) {
        cout << cnc.convert(nonMatch[n].begin(), nonMatch[n].end()) << "\t";
        cout << cnc.convert(nonMatch[n].begin() + positions[n], nonMatch[n].begin() + positions[n] + expOptions.mfinderOptions.width);
        cout << "\t" << positions[n] + 1 << "\t" << nonMatch[n].size() << endl;
    }
}


typedef struct MotifModel {
    boost::shared_ptr<NonUniformMarkov> motif;
    boost::shared_ptr<UnivariatePDF> spacer;
    boost::shared_ptr<UniformMarkov> background;
} MotifModel;

MotifModel estimateMotifModel(const vector<NumSequence> &upstreams, const vector<NumSequence::size_type> &positions,
                              unsigned motifOrder, unsigned backOrder,
                              size_t width, double pcounts,
                              const NumAlphabetDNA &alph) {
    
    // get max upstream length
    size_t maxUpstreamLength = 0;
    for (size_t n = 0; n < upstreams.size(); n++)
        if (upstreams[n].size() > maxUpstreamLength)
            maxUpstreamLength = upstreams[n].size();
    
    // get counts
    NonUniformCounts motifCounts(motifOrder, width, alph);
    UniformCounts backCounts(backOrder, alph);
    vector<double> positionCounts(maxUpstreamLength, 0);
    
    for (size_t n = 0; n < upstreams.size(); n++) {
        size_t pos = positions[n];
        
        motifCounts.count(upstreams[n].begin()+pos, upstreams[n].begin()+pos+width);        // motif count
        
        // count background
        if (pos > 0)
            backCounts.count(upstreams[n].begin(), upstreams[n].begin() + pos);
        if (pos < upstreams[n].size() - width)
            backCounts.count(upstreams[n].begin() + pos + width, upstreams[n].end());
        
        // count spacer
        positionCounts[upstreams[n].size() - width - pos]++;
    }
    
    // construct probability from counts
    MotifModel model;
    
    model.motif = boost::shared_ptr<NonUniformMarkov> (new NonUniformMarkov(motifOrder, width, alph));
    model.motif->construct(&motifCounts, pcounts);
    
    model.background = boost::shared_ptr<UniformMarkov> (new UniformMarkov(backOrder, alph));
    model.background->construct(&backCounts, pcounts);
    
    model.spacer = boost::shared_ptr<UnivariatePDF> (new UnivariatePDF(positionCounts, pcounts));

    return model;
}


double scoreMotifAtAllPositions(const MotifModel &model, const NumSequence &sequence, size_t width) {
    
    if (sequence.size() < width)
        throw std::invalid_argument("Sequence length cannot be shorter than motif width.");
    
    Sequence::size_type numPositions = sequence.size() - width + 1;
    
    vector<double> scores (numPositions, 0);
    
    for (Sequence::size_type pos = 0; pos < numPositions; pos++) {
        
        // compute motif score
        double motifScore = model.motif->evaluate(sequence.begin()+pos, sequence.begin()+pos+width);
        
        // compute background scores of motif
        double backScore = model.background->evaluate(sequence.begin()+pos, sequence.begin()+pos+width);;
        
        // compute combined score
        if (backScore == 0)
            continue;
        
        double score = motifScore / backScore;
        score *= (*model.spacer)[sequence.size() - width - pos];
        
        scores[pos] = score;
    }
    
    return *std::max_element(scores.begin(), scores.end());
}


void scoreAllStarts(const NumSequence &sequence, const NumGeneticCode &gc, const MotifModel &mRBS, const MotifModel &mPromoter,
                    size_t upstreamLengthRBS, size_t upstreamLengthPromoter) {
    
    for (size_t n = 0; n < sequence.size(); n++) {
        
        // get upstreams used for RBS and promoter
        NumSequence upstreamRBS, upstreamPromoter;
        
        // check if start on positive strand
        if (n < sequence.size() - 3 && gc.isStart(CharNumConverter::seq_t(sequence.begin()+n, sequence.begin()+n+3))) {
            
            // get RBS upstream
            size_t upstrLeft = 0;
            if (n > upstreamLengthRBS)
                upstrLeft = n - upstreamLengthRBS;
            
             upstreamRBS = sequence.subseq(upstrLeft, upstreamLengthRBS);
            
            // get promoter upstream
            upstrLeft = 0;
            if (n > upstreamLengthPromoter)
                upstrLeft = n - upstreamLengthPromoter;
            
            upstreamPromoter = sequence.subseq(upstrLeft, upstreamLengthPromoter);
            
        }
        // check if start on negative strand
        if (n >= 3 && gc.isStart(CharNumConverter::seq_t(sequence.begin()+n-3, sequence.begin()+n))) {
            
            // get RBS upstream
            size_t upstrLeft = n+1;
            size_t upstrLen = upstreamLengthRBS;
            size_t distToEnd = sequence.size() - n - 1;
            if (distToEnd < upstreamLengthRBS)
                upstrLen = distToEnd;
            
            upstreamRBS = sequence.subseq(upstrLeft, upstrLen);
            
            // get Promoter upstream
            upstrLeft = n+1;
            upstrLen = upstreamLengthPromoter;
            distToEnd = sequence.size() - n - 1;
            if (distToEnd < upstreamLengthPromoter)
                upstrLen = distToEnd;
            
            upstreamPromoter = sequence.subseq(upstrLeft, upstrLen);
        }
        
        double rbsScore = scoreMotifAtAllPositions(mRBS, upstreamRBS, mRBS.motif->getLength());
        double promoterScore = scoreMotifAtAllPositions(mPromoter, upstreamPromoter, mPromoter.motif->getLength());
        
        double maxScore = std::max(rbsScore, promoterScore);
        
        cout << n+1 << "\t" << maxScore << endl;
    }
    
}


void ModuleExperiment::runScoreStarts() {
    
    
    OptionsExperiment::ScoreStarts expOptions = options.scoreStarts;
    
    // read sequence file
    SequenceFile sequenceFile (expOptions.fn_seqeuence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (expOptions.fn_labels, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // filter short genes
    if (expOptions.minGeneLength > 0) {
        for (size_t i = 0; i < labels.size(); i++) {
            size_t length = labels[i]->right - labels[i]->left + 1;
            if (length < expOptions.minGeneLength) {
                delete labels[i];
                labels[i] = NULL;
            }
        }
        
        size_t sizeBefore = labels.size();
        
        // erase from vector
        
        vector<Label*>::iterator toRem = remove_if(labels.begin(), labels.end(), isNull);
        labels.erase(toRem, labels.end());
        
        cout << "Num Filtered by Gene Length: " << sizeBefore - labels.size() << endl;
    }
    
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    NumAlphabetDNA numAlph(alph, cnc);
    NumSequence numSequence (strSequence, cnc);
    
    // get query
    Sequence strMatchSeq (expOptions.matchTo);
    NumSequence matchSeq (strMatchSeq, cnc);

    
    // separate labels of first genes in operon from rest
    vector<LabelsParser::operon_status_t> operonStatus;
    LabelsParser::partitionBasedOnOperonStatus(labels, expOptions.fgioThresh, expOptions.nfgioThresh, operonStatus);
    
    // separate genes to those for promoter VS rbs training
    typedef enum {RBS, PROMOTER, NONE} training_class_t;
    vector<training_class_t> geneTrainingClass (labels.size(), NONE);
    
    vector<pair<NumSequence::size_type, NumSequence::size_type> > positionsOfMatches (labels.size());
    vector<pair<NumSequence::num_t, NumSequence::num_t> > substitutions;
    if (expOptions.allowAGSubstitution)
        substitutions.push_back(pair<NumSequence::num_t, NumSequence::num_t> (cnc.convert('A'), cnc.convert('G')));

    // For each gene, figure out whether it should be used for RBS or promoter building
    for (size_t n = 0; n < labels.size(); n++) {
        
        // extract upstream for match
        NumSequence sub = SequenceParser::extractUpstreamSequence(numSequence, *labels[n], cnc, expOptions.searchUpstrLen);
        // match against 16S tail
        NumSequence match = SequenceAlgorithms::longestMatchTo16S(matchSeq, sub, positionsOfMatches[n], substitutions);
        
        // if gene is FGIO
        if (operonStatus[n] == LabelsParser::FGIO) {
            
            // if match length less than threshold, then this is used for promoter search
            if (match.size() < expOptions.min16SMatch)
                geneTrainingClass[n] = PROMOTER;
        }
        // if gene is NFGIO
        else if (operonStatus[n] == LabelsParser::NFGIO) {
            
            // if match length greater than threshold, then use for RBS search
            if (match.size() >= expOptions.min16SMatch)
                geneTrainingClass[n] = RBS;
        }
        // if gene status is ambiguous
        else {
            
        }
    }
    
    
    vector<NumSequence> upstreamsForRBS;
    vector<NumSequence> upstreamsForPromoter;
    
    for (size_t n = 0; n < labels.size(); n++) {
        // for RBS genes
        if (geneTrainingClass[n] == RBS) {
            // get upstream
            NumSequence upstream = SequenceParser::extractUpstreamSequence(numSequence, *labels[n], cnc, expOptions.upstreamLenRBS);
            upstreamsForRBS.push_back(upstream);
        }
        // for promoter genes
        else if (geneTrainingClass[n] == PROMOTER) {
            
            // get upstream
            NumSequence upstream = SequenceParser::extractUpstreamSequence(numSequence, *labels[n], cnc, expOptions.upstreamLenPromoter);
            upstreamsForPromoter.push_back(upstream);
        }
    }
    
    
    //  "EXPERIMENT: MOTIF SEARCH FOR RBS"
    
    MotifFinder::Builder bRBS;
    MotifFinder mfinderForRBS = bRBS.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positionsForRBS;
    mfinderForRBS.findMotifs(upstreamsForRBS, positionsForRBS);
    
    // build model RBS
    MotifModel rbsModel = estimateMotifModel(upstreamsForRBS, positionsForRBS, expOptions.mfinderOptions.motifOrder, expOptions.mfinderOptions.bkgdOrder, expOptions.mfinderOptions.width, expOptions.mfinderOptions.pcounts, numAlph);

    
    //  "EXPERIMENT: MOTIF SEARCH FOR Promoter"
    
    MotifFinder::Builder bPromoter;
    MotifFinder mfinderForPromoter = bPromoter.build(expOptions.mfinderOptions);
    
    vector<NumSequence::size_type> positionsForPromoter;
    mfinderForPromoter.findMotifs(upstreamsForPromoter, positionsForPromoter);
    
    // build model RBS
    MotifModel promoterModel = estimateMotifModel(upstreamsForPromoter, positionsForPromoter, expOptions.mfinderOptions.motifOrder, expOptions.mfinderOptions.bkgdOrder, expOptions.mfinderOptions.width, expOptions.mfinderOptions.pcounts, numAlph);

    
    GeneticCode gc(GeneticCode::ELEVEN);
    NumGeneticCode numGC(gc, cnc);
    
    scoreAllStarts(numSequence, numGC, rbsModel, promoterModel, expOptions.upstreamLenRBS, expOptions.upstreamLenPromoter);
    
    
    
    
    
}



































