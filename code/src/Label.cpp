//
//  Label.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/24/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Label.hpp"

using namespace gmsuite;

Label::Label(size_t left, size_t right, strand_t strand) {
    this->left = left;
    this->right = right;
    this->strand = strand;
}