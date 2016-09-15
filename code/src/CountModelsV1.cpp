//
//  CountModelsV1.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "CountModelsV1.hpp"

using namespace gmsuite;

// constructor
CountModelsV1::CountModelsV1(const AlphabetDNA &alphabet, Sequence::size_type width, unsigned motifOrder, unsigned backOrder, MFinderModelParams::align_t align) {
    this->alphabet = &alphabet;
    this->width = width;
    this->motifOrder = motifOrder;
    this->backOrder = backOrder;
    this->align = align;
    
    
    // allocate models
    cnc = new CharNumConverter(this->alphabet);
    mMotif = new NonUniformCounts(motifOrder, width, alphabet, *cnc);
    mBack = new UniformCounts(motifOrder, alphabet, *cnc);
}

// copy constructor
CountModelsV1::CountModelsV1(const CountModelsV1 &obj) {
    this->alphabet = obj.alphabet;
    width = obj.width;
    motifOrder = obj.motifOrder;
    backOrder = obj.backOrder;
    align = obj.align;
    positionCounts = obj.positionCounts;
    
    // deep copy
    cnc = new CharNumConverter(*cnc);
    mMotif = new NonUniformCounts(*obj.mMotif);
    mBack = new UniformCounts(*obj.mBack);
}

// copy assignment operator
CountModelsV1& CountModelsV1::operator=(const CountModelsV1& other) {
    this->alphabet = other.alphabet;
    width = other.width;
    motifOrder = other.motifOrder;
    backOrder = other.backOrder;
    align = other.align;
    positionCounts = other.positionCounts;
    
    // deep copy
    if (cnc != NULL)
        delete cnc;
    if (mMotif != NULL)
        delete mMotif;
    if (mBack != NULL)
        delete mBack;
    
    cnc = new CharNumConverter(*other.cnc);
    mMotif = new NonUniformCounts(*other.mMotif);
    mBack = new UniformCounts(*other.mBack);
    return *this;
}

/**
 * Destructor
 */
CountModelsV1::~CountModelsV1() {
    delete cnc;
    delete mMotif;
    delete mBack;
}



/**
 * Construct counts from a set of motifs.
 *
 * @param sequences the sequences
 * @param positions the positions of motifs in these sequences
 */
void CountModelsV1::construct(const vector<NumSequence> &sequences, const vector<NumSequence::size_type> &positions) {
    
    mMotif->resetCounts();
    mBack->resetCounts();
    
    if (align != MFinderModelParams::NONE) {
        // get maximum motif position
        Sequence::size_type maxSequenceSize = 0;
        for (vector<NumSequence>::size_type n = 0; n < sequences.size(); n++) {
            if (sequences[n].size() > maxSequenceSize) {
                maxSequenceSize = sequences[n].size();
            }
        }
        
        // allocate space for position counts
        positionCounts.resize(maxSequenceSize-width+1, 0);
    }
    
    
    // add all sequences to counts
    for (vector<NumSequence>::size_type n = 0; n < sequences.size(); n++) {
        count(sequences[n], positions[n]);
    }
    
}


// Decount a motif
void CountModelsV1::decount(const NumSequence &sequence, NumSequence::size_type pos) {
    
    // decount motif model
    mMotif->decount(sequence.begin() + pos, sequence.begin() + pos + width);
    
    // decount background model
    if (pos > 0)
        mBack->decount(sequence.begin(), sequence.begin() + pos);
    if (pos < sequence.size() - width)
        mBack->decount(sequence.begin()+pos+width, sequence.end());
    
    // decount position
    if (align != MFinderModelParams::NONE) {
        // for left aligned
        if (align == MFinderModelParams::LEFT) {
            if (positionCounts[pos] == 0)
                throw std::invalid_argument("Cannot decount position for sequence.");
            
            positionCounts[pos]--;                          // increment position from left
        }
        // for right aligned
        else if (align == MFinderModelParams::RIGHT) {
            if (positionCounts[sequence.size() - width - pos] == 0)
                throw std::invalid_argument("Cannot decount position for sequence.");
            
            positionCounts[sequence.size() - width - pos]--;  // increment position from right
        }
    }
    
}



/**
 * Count a motif
 *
 * @param sequence the sequence to count
 * @param pos the position of the motif in the sequence
 *
 * @exception invalid_argument if pos is not a valid motif location in the sequence
 */
void CountModelsV1::count(const NumSequence &sequence, NumSequence::size_type pos) {
    // count motif model
    mMotif->count(sequence.begin() + pos, sequence.begin() + pos + width);
    
    // count position model
    if (pos > 0)
        mBack->count(sequence.begin(), sequence.begin() + pos);
    if (pos < sequence.size() - width)
        mBack->count(sequence.begin()+pos+width, sequence.end());
    
    // count position
    if (align != MFinderModelParams::NONE) {
        // for left aligned
        if (align == MFinderModelParams::LEFT) {
            positionCounts[pos]++;                          // increment position from left
        }
        // for right aligned
        else if (align == MFinderModelParams::RIGHT) {
            positionCounts[sequence.size() - width - pos]++;  // increment position from right
        }
    }
}


/**
 * Get string representation of counting models
 */
string CountModelsV1::toString() const {
    string output = "";
    
    output += mMotif->toString() + "\n\n";
    output += mBack->toString();
    
    return output;
}








