//
//  NumGeneticCode.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "NumGeneticCode.hpp"
#include <algorithm>

using namespace gmsuite;


// Constructor: Create a numeric-based genetic code
NumGeneticCode::NumGeneticCode(const GeneticCode &gcode, const CharNumConverter &converter) {
    
    // TODO: 
    
}


/******************** Request info on Starts, Stops, and Translating AA ********************/

// Check if codon is a start codon
bool NumGeneticCode::isStart(int codon) const {
    return std::find(starts.begin(), starts.end(), codon) != starts.end();
}


// Check if codon is a stop codon
bool NumGeneticCode::isStop(int codon) const {
    return std::find(stops.begin(), stops.end(), codon) != stops.end();
}


// Get the list of start codons
vector<int> NumGeneticCode::getStarts() const {
    return starts;
}


// Get the list of stop codons
vector<int> NumGeneticCode::getStops() const {
    return stops;
}


// Translate a codon into an amino acid
int NumGeneticCode::translateCodon(int codon) const {
    return translationTable.at(codon);
}


// Get the genetic code value that represents this object.
NumGeneticCode::gcode_t NumGeneticCode::getGCode() const {
    return this->gcode;
}
