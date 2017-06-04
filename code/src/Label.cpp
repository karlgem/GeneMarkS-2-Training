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
