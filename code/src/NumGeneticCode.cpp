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
    
    // copy genetic code
    this->gcode = gcode.getGCode();
    
    // convert starts to numeric
    vector<string> stringStarts = gcode.getStarts();
    for (size_t i = 0; i < stringStarts.size(); i++) {
        CharNumConverter::seq_t newStart;
        converter.convert(stringStarts[i], newStart);           // convert start
        starts.push_back(newStart);
    }
    
    
    // convert stops to numeric
    vector<string> stringStops = gcode.getStops();
    for (size_t i = 0; i < stringStops.size(); i++) {
        CharNumConverter::seq_t newStop;
        converter.convert(stringStops[i], newStop);             // convert stop
        stops.push_back(newStop);
    }
    
    
    
    
}


/******************** Request info on Starts, Stops, and Translating AA ********************/

// Check if codon is a start codon
bool NumGeneticCode::isStart(CharNumConverter::seq_t codon) const {
    return std::find(starts.begin(), starts.end(), codon) != starts.end();
}


// Check if codon is a stop codon
bool NumGeneticCode::isStop(CharNumConverter::seq_t codon) const {
    return std::find(stops.begin(), stops.end(), codon) != stops.end();
}


// Get the list of start codons
vector<CharNumConverter::seq_t> NumGeneticCode::getStarts() const {
    return starts;
}


// Get the list of stop codons
vector<CharNumConverter::seq_t> NumGeneticCode::getStops() const {
    return stops;
}



// Get the genetic code value that represents this object.
NumGeneticCode::gcode_t NumGeneticCode::getGCode() const {
    return this->gcode;
}
