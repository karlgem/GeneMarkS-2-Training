//
//  GeneticCode.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/22/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "GeneticCode.hpp"
#include <algorithm>

using namespace gmsuite;

// empty constructor
GeneticCode::GeneticCode() {
    
}

// constructor with gcode value
GeneticCode::GeneticCode(gcode_t gcode) {
    this->gcode = gcode;
    
    initialize();
}


GeneticCode::GeneticCode(const vector<pair<string, char> > &table, const vector<string> &starts, const vector<string> &stops) {
    gcode = CUSTOM;
    
    // initialize custom genetic code
    
}



/***************************************************************************************\
 *                                                                                     *
 *                  Request info on Starts, Stops, and Translating AA                  *
 *                                                                                     *
 \***************************************************************************************/

/**
 * Check if codon is a start codon
 *
 * @param codon candidate codon
 * @return true if codon is a start { } false otherwise
 */
bool GeneticCode::isStart(string codon) const {
    return std::find(starts.begin(), starts.end(), codon) != starts.end();
}


/**
 * Check if codon is a stop codon
 *
 * @param codon candidate codon
 * @return true if codon is a stop { } false otherwise
 */
bool GeneticCode::isStop(string codon) const {
    return std::find(stops.begin(), stops.end(), codon) != stops.end();
}


/**
 * Get the list of start codons
 *
 * @return a vector of strings, where each string is a 3-letter codon representing a start
 */
vector<string> GeneticCode::getStarts() const {
    return starts;
}


/**
 * Get the list of stop codons
 *
 * @return a vector of strings, where each string is a 3-letter codon representing a stop
 */
vector<string> GeneticCode::getStops() const {
    return stops;
}


/**
 * Translate 3-letter codon to amino acid
 *
 * @param codon the codon to be translated
 */
char GeneticCode::translateCodon(string codon) const {
    return translationTable.at(codon);
}

// get the genetic code value
GeneticCode::gcode_t GeneticCode::getGCode() const {
    return this->gcode;
}



/***************************************************************************************\
 *                                                                                     *
 *                    Initialize Genetic Code and Update Parameters                    *
 *                                                                                     *
\***************************************************************************************/


// genetic code cannot be CUSTOM
void GeneticCode::initialize() {
    if (gcode == FOUR)                      // genetic code 4
        initialize4();
    else if (gcode == ELEVEN)               // genetic code 11
        initialize11();
}




// init genetic code 4
void GeneticCode::initialize4() {
    
}

// init genetic code 11
void GeneticCode::initialize11() {
    
    string AAs      = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    string starts   = "---M-------------------------------M---------------M------------";
    string base1    = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
    string base2    = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
    string base3    = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
    
    fillTranslationTable(AAs, starts, base1, base2, base3);
}







void GeneticCode::fillTranslationTable(string AAs, string starts, string base1, string base2, string base3) {
    
    // make sure strings have the same length
    size_t length = AAs.length();
    
    if (starts.length() != length || base1.length() != length || base2.length() != length || base3.length() != length) {
        // TODO: raise an exception
    }
    
    // reset table and starts
    this->starts.clear();
    this->stops.clear();
    this->translationTable.clear();
    this->translationTableKeys.clear();
    
    string codon = "XXX";       // will hold the codon
    
    for (size_t i = 0; i < length; i++) {
        
        // fill in codon
        codon[0] = base1[i];
        codon[1] = base2[i];
        codon[2] = base3[i];
        
        bool isStart = starts[i] == 'M';
        bool isStop = AAs[i] == '*';
        
        if (isStart) {
            this->starts.push_back(codon);
        }
        else if (isStop) {
            this->stops.push_back(codon);
        }
        
        // push amino acid (including for starts and stops) into translation table
        translationTable[codon] = AAs[i];
        translationTableKeys.push_back(codon);
    }
}






/***************************************************************************************\
 *                                                                                     *
 *                                      Iterators                                      *
 *                                                                                     *
\***************************************************************************************/


// begin iterator
GeneticCode::ttk_const_iterator GeneticCode::begin() const {
    return this->translationTableKeys.begin();
}

// Get end iterator over translation table keys.
GeneticCode::ttk_const_iterator GeneticCode::end() const {
    return this->translationTableKeys.end();
}











