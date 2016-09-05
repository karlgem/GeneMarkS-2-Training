//
//  ProbabilityModelsV1.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ProbabilityModelsV1.hpp"
#include "CountModelsV1.hpp"

using namespace gmsuite;


/**
 * Constructor
 */
ProbabilityModelsV1::ProbabilityModelsV1(const AlphabetDNA &alphabet, NumSequence::size_type width, unsigned motifOrder, unsigned backOrder, double pcounts, MFinderModelParams::align_t align) {
    this->alphabet = &alphabet;
    this->width = width;
    this->motifOrder = motifOrder;
    this->backOrder = backOrder;
    this->pcounts = pcounts;
    this->align = align;
    
    // allocate models
    mMotif = new NonUniformMarkov(motifOrder, width, alphabet);
    mMotifCounts = new NonUniformCounts(motifOrder, width, alphabet);
    mBack = new UniformMarkov(motifOrder, alphabet);
    
    // if alignment is set, create length distribution
    positionDistribution = NULL;
    if (align != MFinderModelParams::NONE) {
        positionDistribution = new UnivariatePDF();      // empty distribution
    }
}

/**
 * Copy constructor
 */
ProbabilityModelsV1::ProbabilityModelsV1(const ProbabilityModelsV1 &obj) {
    this->alphabet = obj.alphabet;
    width = obj.width;
    motifOrder = obj.motifOrder;
    backOrder = obj.backOrder;
    pcounts = obj.pcounts;
    align = obj.align;
    positionCounts = obj.positionCounts;
    
    // deep copy
    mMotif = new NonUniformMarkov(*obj.mMotif);
    mMotifCounts = new NonUniformCounts(*obj.mMotifCounts);
    mBack = new UniformMarkov(*obj.mBack);
    
    if (align != MFinderModelParams::NONE)
        positionDistribution = new UnivariatePDF(*obj.positionDistribution);

}

/**
 * Copy assignment operator
 */
ProbabilityModelsV1& ProbabilityModelsV1::operator=(const ProbabilityModelsV1& other) {
    this->alphabet = other.alphabet;
    width = other.width;
    motifOrder = other.motifOrder;
    backOrder = other.backOrder;
    pcounts = other.pcounts;
    align = other.align;
    positionCounts = other.positionCounts;
    
    // deep copy
    if (mMotif != NULL)
        delete mMotif;
    if (mMotifCounts != NULL)
        delete mMotifCounts;
    if (mBack != NULL)
        delete mBack;
    if (positionDistribution != NULL)
        delete positionDistribution;
    
    mMotif = new NonUniformMarkov(*other.mMotif);
    mMotifCounts = new NonUniformCounts(*other.mMotifCounts);
    mBack = new UniformMarkov(*other.mBack);
    
    if (align != MFinderModelParams::NONE)
        positionDistribution = new UnivariatePDF(*other.positionDistribution);
    
    return *this;
}

/**
 * Destructor
 */
ProbabilityModelsV1::~ProbabilityModelsV1() {
    delete mMotif;
    delete mBack;
    delete mMotifCounts;
    
    if (positionDistribution != NULL)
        delete positionDistribution;
}


/**
 * Construct probabilities from a set of motifs.
 *
 * @param sequences the sequences
 * @param positions the positions of motifs in these sequences
 */
void ProbabilityModelsV1::construct(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions) {
    // create counts
    CountModelsV1 *counts = new CountModelsV1(*alphabet, width, motifOrder, backOrder, align);
    counts->construct(sequences, positions);
    
    construct(counts);
    
    delete counts;
}


/**
 * Construct probabilities from counts
 *
 * @param sequences the sequences
 * @param positions the positions of motifs in these sequences
 */
void ProbabilityModelsV1::construct(const CountModels* counts) {
    const CountModelsV1* countsV1 = dynamic_cast<const CountModelsV1*>(counts);
    
    mMotif->construct(countsV1->mMotif, pcounts);               // create motif probabilities
    (*mMotifCounts) = (*countsV1->mMotif);                      // copy motif counts
    mBack->construct(countsV1->mBack, pcounts);                 // create background probabilities
    
    positionCounts = countsV1->positionCounts;                  // copy position counts
    
    if (align != MFinderModelParams::NONE)
        positionDistribution->construct(countsV1->positionCounts, false, pcounts);         // create position distribution

}


/**
 * Compute conditional log-likelihood
 */
double ProbabilityModelsV1::computeCLL() {
    double score = 0;
    
    // for each position in motif
    for (Sequence::size_type pos = 0; pos < width; pos++) {
        
//        // loop over all keys at current position in motif model
//        for (NonUniformMarkov::const_iterator_keys itKey = mMotif->beginKeys(pos); itKey != mMotif->endKeys(pos); itKey++) {
//            
//            double ratio = 0;
//            if (mBack->valueAtKey(*itKey) != 0)
//                ratio = mMotif->valueAtKey(pos, *itKey) / mBack->valueAtKey(*itKey);
//            
//            if (ratio != 0)
//                score += mMotifCounts->valueAtKey(pos, *itKey) * log(ratio);
//            
//        }
    }
    
    
    if (align != MFinderModelParams::NONE) {
        
//        // for each valid position in model
//        for (vector<double>::size_type p = 0; p < positionCounts.size(); p++) {
//            if ((*positionDistribution)[p] != 0)
//                score += positionCounts[p] * log((*positionDistribution)[p]);
//            
//        }
        
    }
    
    return score;
}


/**
 * Compute the score for a given motif position
 *
 * @param sequence the sequence
 * @param pos the position of the motif in the sequence
 */
double ProbabilityModelsV1::computePositionScore(const NumSequence &sequence, NumSequence::size_type pos) {
    // copute motif score
    double motifScore = mMotif->evaluate(sequence.begin()+pos, sequence.begin()+pos+width);
    
    // compute background scores of motif
    double backScore = mBack->evaluate(sequence.begin()+pos, sequence.begin()+pos+width);;
    
    // compute combined score
    if (backScore == 0)
        return 0;
    
    double score = motifScore / backScore;
    
    if (align != MFinderModelParams::NONE) {
        if (align == MFinderModelParams::LEFT)
            score *= (*positionDistribution)[pos];
        else if (align == MFinderModelParams::RIGHT)
            score *= (*positionDistribution)[sequence.size() - width - pos];
    }
    
    return score;
}


/**
 * Compute the scores for each valid motif position in the sequence
 *
 * @param sequence the sequence
 * @param scores the output scores of all valid positions in the sequence
 */
void ProbabilityModelsV1::computePositionScores(const NumSequence &sequence, vector<double> &scores) {
    if (sequence.size() < width)
        throw std::invalid_argument("Sequence length cannot be shorter than motif width.");
    
    Sequence::size_type numPositions = sequence.size() - width + 1;
    
    scores.resize(numPositions, 0);
    
    for (Sequence::size_type pos = 0; pos < numPositions; pos++)
        scores[pos] = computePositionScore(sequence, pos);

}


/**
 * Sample the position for motif in sequences
 *
 * @param sequence the sequence
 * @param getMax if set, the position with the highest probability is returned; otherwise it is sampled.
 *
 * @return the position of a motif
 */
NumSequence::size_type ProbabilityModelsV1::samplePosition(const NumSequence &sequence, bool getMax) {
    // compute all position scores
    vector<double> scores;
    computePositionScores(sequence, scores);
    
    // if get max is on, simply get position of max score
    if (getMax) {
        return std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));
    }
    
    // otherwise, build a distribution and then sample value
    UnivariatePDF dist (scores);
    return dist.sample();
}

/**
 * Get string representation of counting models
 */
string ProbabilityModelsV1::toString() const {
    string output = "";
    
    output += mMotif->toString() + "\n\n";
    output += mBack->toString();
    
    return output;
}
