//
//  Label.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/24/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Label.hpp"
#include <sstream>

using namespace std;
using namespace gmsuite;

Label::Label(size_t left, size_t right, strand_t strand, string geneClass, string meta) {
    this->left = left;
    this->right = right;
    this->strand = strand;
    this->geneClass = geneClass;
    this->meta = meta;
}


/**
 * Get a string representation of the label. By default, label indeces are 0-indexed.
 * However, for compatibility with existing datafile formats, the indexFromOne option
 * allows the labels to be indexed from 1 (i.e. shifted up by 1).
 */
string Label::toString(bool indexFromOne) const {
    
    // determine amount of shift
    size_t shift = 0;
    if (indexFromOne)
        shift = 1;
    
    stringstream ssm;
    
    // create string
    ssm << left+shift << "\t" << right+shift << "\t" << (strand == POS ? "+" : "-");
    
    return ssm.str();
}

bool Label::compareByLeftAndStrand(Label* a, Label* b) {
    if (a->left < b->left)
        return true;
    if (a->left == b->left)
        return a->strand < b->strand;
    return false;
}


size_t Label::get5Prime() const {
    
    if (strand == Label::NONE)
        invalid_argument("Cannot return 5-prime end. Unknown strand value");
    
    return (strand == Label::POS ? left : right);
}

size_t Label::get3Prime() const {
    
    if (strand == Label::NONE)
        invalid_argument("Cannot return 3-prime end. Unknown strand value");
    
    return (strand == Label::POS ? right : left);
}

string Label::strandToStr() const {
    if (strand == Label::NONE)
        invalid_argument("Cannot return string representation of strand. Unknown strand value");
    
    return (strand == Label::POS ? "+" : "-");
}



