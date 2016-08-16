//
//  Genome.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/21/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "Genome.hpp"

using namespace gmsuite;

/*************** Constructors and Destructors *******************/

// default constructor
Genome::Genome(string name) {
    this->name = name;
}


// construct genome with genetic code
Genome::Genome(GeneticCode gcode, string name) {
    this->gcode = gcode;
    this->name = name;
}


// construct genome with sequences and genetic code
Genome::Genome(const vector<Sequence> &sequences, GeneticCode gcode, string name) {
    this->gcode = gcode;
    this->name = name;
    
    // create sequences from strings
    this->sequences = sequences;
}


// destructor
Genome::~Genome() {

}



/*************** Indexing genome (sequences) and adding sequences *******************/

// access i'th sequence
const Sequence&  Genome::operator[](size_type idx) const {
    return this->sequences[idx];
}


// add sequence to genome
void Genome::add(Sequence seq) {
    this->sequences.push_back(seq);
}

// number of sequences
size_t Genome::numOfSeq() const {
    return sequences.size();
}










