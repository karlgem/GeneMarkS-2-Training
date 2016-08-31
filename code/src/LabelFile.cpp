//
//  LabelFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/31/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "LabelFile.hpp"

#include <fstream>
#include <boost/xpressive/xpressive.hpp>

using namespace boost::xpressive;
using namespace gmsuite;
using namespace std;


Label* readNextLabelLST(const char*& current, const char* const end);


// constructor
LabelFile::LabelFile(string path, access_t access, format_t format) {
    this->path = path;
    this->access = access;
    this->format = format;
    
    // setup file parameters
    params.path = path;
    if (access == READ)
        params.flags = io::mapped_file::readonly;           // read-only file
    else if (access == WRITE)
        params.flags = io::mapped_file::readwrite;         // read-write file
    
    openFile();                         // open file
    
    // if format not specified (i.e. AUTO), try to guess what it is
    if (this->format == AUTO) {
        this->format = detectFormat();      // detect format
    }
}


// free memory (if any) and close file
LabelFile::~LabelFile() {
    closeFile();
}


// Read labels from file. This behaves differently for separate file formats.
void LabelFile::read(vector<Label*> &output) const {
    
    // read based on format set
    if (this->format == LST)
        read_lst(output);
}


// Write labels to file.
void LabelFile::write(const vector<Label*> &labels) const {
    
}



/**
 * Read labels from LST file.
 *
 * @param output a vector label pointers that have been read from the file.
 */
void LabelFile::read_lst(vector<Label*> &output) const {
    
    
    output.clear();         // clear output vector (sanity check)
    Label* currlabel;       // stores single label information
    
    // point to start of data
    const char* current = begin_read;
    
    // loop over all the file
    while (current != end_read) {
        
        // read next label
        currlabel = readNextLabelLST(current, end_read);
        
        if (currlabel != NULL)
            output.push_back(currlabel);
        else
            throw runtime_error("Could not read label from file.");
    }    
}







Label* readNextLabelLST(const char*& current, const char* const end) {
    string output = "";
    
    // start of line
    const char* startOfLine = current;
    
    // read the remainder of the line
    while (current != end && *current != '\n' && *current != '\r')
        current++;
    
    const char* endOfLine = current;
    
    cmatch match;
    cregex expr = cregex::compile("^\\s*(\\d+)\\s+([+,-])\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+([1,2])\\s*$");
    
    Label* label = NULL;
    
    // parse 'line' for label information
    if (regex_search(startOfLine, endOfLine, match, expr) && match.size() > 1) {
        size_t left;
        size_t right;
        char strandChar;
        
        left = (size_t) strtol(match.str(3).c_str(), NULL, 10);
        right = (size_t) strtol(match.str(4).c_str(), NULL, 10);
        strandChar = match.str(2)[0];
        
        Label::strand_t strand = Label::POS;
        if (strandChar == '-')
            strand = Label::NEG;

        label = new Label(left, right, strand);
    }
    
    
    return label;
}















