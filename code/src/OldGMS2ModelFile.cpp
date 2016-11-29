//
//  OldGMS2ModelFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 11/7/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OldGMS2ModelFile.hpp"
#include <boost/xpressive/xpressive.hpp>
#include <fstream>
#include <map>

using namespace boost::xpressive;
using namespace std;
using namespace gmsuite;


OldGMS2ModelFile::OldGMS2ModelFile(string filename) {
    this->filename = filename;
}

vector<pair<string, double> > OldGMS2ModelFile::getNoncoding() const {
    vector<pair<string, double> > noncoding;
    
    // open file for reading
    ifstream infile (filename.c_str());
    string line;
    
    sregex eProbWord = sregex::compile("\\$Prob_letter");
    
    // goto probword keyword
    while (getline(infile, line)) {
        smatch match;
        
        if (regex_search(line, match, eProbWord))
            break;
    }
    
    //            Word            COD1            COD2              COD3            NONC
    sregex e = sregex::compile("^\\s*(\\S+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)");

    // read joint probs
    while (getline(infile, line)) {
        smatch match;
        
        if (regex_search(line, match, e) && match.size() > 1) {
            pair<string, double> value (match.str(1), strtod(match.str(5).c_str(), NULL));
            noncoding.push_back(value);
        }
        else
            break;
    }
    
    return noncoding;
}

size_t OldGMS2ModelFile::getNDEC() const {
    // open file for reading
    ifstream infile (filename.c_str());
    string line;
    
    // get model order
    sregex eNDEC = sregex::compile ("NDEC\\s+(\\d+)");
    unsigned NDEC = 0;
    
    
    while (getline(infile, line)) {
        smatch match;
        
        if (regex_search(line, match, eNDEC)) {
            NDEC = (unsigned) strtoul(match.str(1).c_str(), NULL, 10);
            break;
        }
    }
    
    return NDEC;
}



