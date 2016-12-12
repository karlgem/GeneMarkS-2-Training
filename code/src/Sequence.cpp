//
//  Sequence.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/21/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Sequence.hpp"
#include <stdexcept>

using namespace gmsuite;

// empty sequence constructor
Sequence::Sequence() {
    
}

// constructor from string
Sequence::Sequence(const string &str, const string &meta, const shape_t shape) {
    this->data = str;
    this->meta = meta;
    this->shape = shape;
}

// Destructor
Sequence::~Sequence() {
    
}



/*************** Basic Operators (length, toString, ...) *******************/


// get sequence size
Sequence::size_type Sequence::size() const {
    return data.size();
}


// get sequence length
Sequence::size_type Sequence::length() const {
    return data.length();
}


// get data
string Sequence::toString(size_type start, size_type length) const {
    if (start != npos) {
        if (start + length > data.length())
            throw std::invalid_argument("Length cannot exceed data size");
        else
            return data.substr(start, length);
    }
    
    return data;
}


// get meta data
string Sequence::getMetaData() const {
    return meta;
}



/*************** Sequence Indexing and Substrings *******************/

// subsequence
Sequence* Sequence::subseq(size_type n, size_type length) const {
    return NULL;
}

// check for char in sequence
bool Sequence::contains(char c, size_type startSearch, size_type searchLength) const {
    return false;
}

// find position of char
Sequence::size_type Sequence::find(char c) const {
    return data.find(c);
}


// get character at idx
const char& Sequence::operator[](size_type idx) const {
    return data[idx];
}


char& Sequence::operator[](size_type idx) {
    return data[idx];
}




/*************** Sequence Iterators *******************/

// Iterators
Sequence::iterator Sequence::begin() {
    return data.begin();
}

Sequence::iterator Sequence::end() {
    return data.end();
}

Sequence::const_iterator Sequence::begin() const {
    return data.begin();
}

Sequence::const_iterator Sequence::end() const {
    return data.end();
}



/*************** Sequence Comparisons *******************/

// overload equality/inequality operators
inline bool Sequence::operator==(const Sequence& rhs) {
    return data == rhs.data;
}

inline bool Sequence::operator!=(const Sequence& rhs) {
    return data != rhs.data;
}

// overload less/greater than operators
inline bool Sequence::operator<(const Sequence& rhs) {
    return data < rhs.data;
}

inline bool Sequence::operator>(const Sequence& rhs) {
    return data > rhs.data;
}

// overload less/greater than or equal operators
inline bool Sequence::operator<=(const Sequence& rhs) {
    return data <= rhs.data;
}

inline bool Sequence::operator>=(const Sequence& rhs) {
    return data >= rhs.data;
}











