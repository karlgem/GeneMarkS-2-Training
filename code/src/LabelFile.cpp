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
Label* readNextLabelGFF(const char*& current, const char* const end);
void gotoKey(const char*& current, const char* const end, string key);


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
    if (access == READ && this->format == AUTO) {
        this->format = detectFormat();      // detect format
    }
    else if (access == WRITE && this->format == AUTO)
        this->format = LST;
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
    if (this->format == GFF)
        read_gff(output);
}


// Write labels to file.
void LabelFile::write(const vector<Label*> &labels) const {
    if (access != WRITE)
        throw logic_error("File not opened for writing.");
    
    if (format == LST)
        write_lst(labels);
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
    
    // ignore all lines until we reach "Predicted Genes"
    gotoKey(current, end_read, "SequenceID");
    
    
    bool foundFirstLabel = false;
    
    // loop over all the file
    while (current != end_read) {
        
        // skip all white spaces
        while (current != end_read && isspace(*current))
            current++;
        
        // skip all commented lines
        while (current != end_read && *current == '#') {
            // skip all white spaces
            while (current != end_read && *current != '\n' && *current != '\r')
                current++;
            
            // skip all white spaces
            while (current != end_read && isspace(*current))
                current++;
        }
        
        if (current == end_read)
            break;
        
        // read next label
        currlabel = readNextLabelLST(current, end_read);
        
        if (currlabel != NULL) {
            output.push_back(currlabel);
            foundFirstLabel = true;
        }
        //        else if (foundFirstLabel)
        //            throw runtime_error("Could not read label from file.");
    }
}



/**
 * Read labels from GFF file.
 *
 * @param output a vector label pointers that have been read from the file.
 */
void LabelFile::read_gff(vector<Label*> &output) const {
    
    
    output.clear();         // clear output vector (sanity check)
    Label* currlabel;       // stores single label information
    
    // point to start of data
    const char* current = begin_read;
    
    bool foundFirstLabel = false;
    
    // loop over all the file
    while (current != end_read) {
        
        // skip all white spaces
        while (current != end_read && isspace(*current))
            current++;
        
        // skip all commented lines
        while (current != end_read && *current == '#') {
            // skip all white spaces
            while (current != end_read && *current != '\n' && *current != '\r')
                current++;
            
            // skip all white spaces
            while (current != end_read && isspace(*current))
                current++;
        }
        
        if (current == end_read)
            break;
        
        // read next label
        currlabel = readNextLabelGFF(current, end_read);
        
        if (currlabel != NULL) {
            output.push_back(currlabel);
            foundFirstLabel = true;
        }
        //        else if (foundFirstLabel)
        //            throw runtime_error("Could not read label from file.");
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
    //                                  gene #     strand      left     right    length     class           meta
    cregex expr = cregex::compile("^\\s*(\\d+)\\s+([+,-])\\s+(<?\\d+)\\s+(>?\\d+)\\s+(\\d+)\\s+(\\S+)(?:\\s+([ACGT]+))?");
    
    Label* label = NULL;
    
    string t (startOfLine, endOfLine);
    
    // parse 'line' for label information
    if (regex_search(startOfLine, endOfLine, match, expr) && match.size() > 1) {
        size_t left;
        size_t right;
        char strandChar;
        string geneClass;
        string meta;
        
        // skip incomplete genes
        if (match.str(3).find("<") != string::npos || match.str(4).find(">") != string::npos)
            return NULL;
        
        left = (size_t) strtol(match.str(3).c_str(), NULL, 10)-1;
        right = (size_t) strtol(match.str(4).c_str(), NULL, 10)-1;
        strandChar = match.str(2)[0];
        geneClass = match.str(6);
        
        // check if meta exists
        if (match[7].matched)
            meta = match.str(7);
        
        Label::strand_t strand;
        
        if (strandChar == '+')
            strand = Label::POS;
        else if (strandChar == '-')
            strand = Label::NEG;
        else
            throw invalid_argument("Invalid gene strand: " + match.str(2));
        
        label = new Label(left, right, strand, geneClass, meta);
    }
    
    
    return label;
}


Label* readNextLabelGFF(const char*& current, const char* const end) {
    string output = "";
    
    // start of line
    const char* startOfLine = current;
    
    // read the remainder of the line
    while (current != end && *current != '\n' && *current != '\r')
        current++;
    
    const char* endOfLine = current;
    
    cmatch match;
    //                                  gene #     strand      left     right    length     class           meta
    cregex expr = cregex::compile("^(\\S+)\t(.+)\t(\\S+)\t(\\d+)\t(\\d+)\t([-+]?\\d+(?:\\.\\d+)?|\\.)\t([+-])\t(\\d+|\\.)\t(.+)$");
    
    Label* label = NULL;
    
    string t (startOfLine, endOfLine);
    
    // parse 'line' for label information
    if (regex_search(startOfLine, endOfLine, match, expr) && match.size() > 1) {
        string seqname;          // chromosome/scaffold name
        string source;           // data source (e.g. project name)
        string feature;          // feature type name (e.g. gene)
        size_t left;             // start position (starting with 1)
        size_t right;            // end position (starting with 1)
        double score;            // a float
        char strandChar;         // + or -
        int frame;               // 0, 1, or 2
        string attribute;        // semicolon-separate list of tag-value pairs
        
        
        seqname         = match.str(1);
        source          = match.str(2);
        feature         = match.str(3);
        left            = (size_t) strtol(match.str(4).c_str(), NULL, 10)-1;
        right           = (size_t) strtol(match.str(5).c_str(), NULL, 10)-1;
        score           = strtod(match.str(6).c_str(), NULL);
        strandChar      = match.str(7)[0];
        frame           = (int) strtol(match.str(8).c_str(), NULL, 10)-1;
        attribute       = match.str(9);
        
        if (source == "RefSeq")
            return NULL;
        if (feature != "CDS")
            return NULL;
        
        Label::strand_t strand;
        
        if (strandChar == '+')
            strand = Label::POS;
        else if (strandChar == '-')
            strand = Label::NEG;
        else
            throw invalid_argument("Invalid gene strand: " + match.str(2));
        
        label = new Label(left, right, strand, "", seqname);
    }
    
    
    return label;
}



void gotoKey(const char*& current, const char* const end, string key) {
    
    cmatch match;
    cregex expr = cregex::compile(key);
    
    while (current != end) {
        
        // skip all new-lines
        while (current != end && (*current == '\n' || *current == '\r'))
            current++;
            
            
        // start of line
        const char* startOfLine = current;
        
        // read the remainder of the line
        while (current != end && *current != '\n' && *current != '\r')
            current++;
        
        const char* endOfLine = current;
        
        // check if line matches key
        if (regex_search(startOfLine, endOfLine, match, expr)) {
            return;
        }
    }
}


/**
 * Detect the file's format. This assumes an already opened file
 *
 * @return the file's format.
 */
LabelFile::format_t LabelFile::detectFormat() const {
    if (detectLST(begin_read, end_read))
        return LST;
    
    if (detectGFF(begin_read, end_read))
        return GFF;
    
    // if no format detected, throw exception
    throw logic_error("File format could not be detected.");
}


bool LabelFile::detectGFF(const char* const begin, const char* const end) const {
    const char* current = begin;
    
    while (current != end) {
        
        // skip whitespaces
        while (current != end && isspace(*current))
            current++;
        
        // start of line
        const char* startOfLine = current;
        
        // read the remainder of the line
        while (current != end && *current != '\n' && *current != '\r')
            current++;
        
        const char* endOfLine = current;
        
        cmatch match;
        cregex expr = cregex::compile("^\\S+\t.+\t\\S+\t\\d+\t\\d+\t([-+]?\\d+(?:\\.\\d+)?|\\.)\t[+-]\t(\\d+|\\.)\t.+$");
        // check if key found
        if (regex_search(startOfLine, endOfLine, match, expr)) {
            return true;
        }
    }
    
    return false;
}

bool LabelFile::detectLST(const char* const begin, const char* const end) const {
    
    const char* current = begin;
    
    while (current != end) {
        
        // skip whitespaces
        while (current != end && isspace(*current))
            current++;
        
        // start of line
        const char* startOfLine = current;
        
        // read the remainder of the line
        while (current != end && *current != '\n' && *current != '\r')
            current++;
        
        const char* endOfLine = current;
        
        cmatch match;
        cregex expr = cregex::compile("SequenceID");
        
        // check if key found
        if (regex_search(startOfLine, endOfLine, match, expr)) {
            return true;
        }
    }
    
    return false;

}





/**
 * Open the file and set start/end pointers to the data.
 */
void LabelFile::openFile() {
    
    // if file is open, close it
    if (mfile.is_open())
        mfile.close();
    
    // (re)set pointers
    if (access == READ) {
        // open file
        mfile.open(params);
        begin_read = mfile.const_data();
        end_read = begin_read + mfile.size();
    }
    else if (access == WRITE) {
//        begin_write = mfile.data();
//        end_write = begin_write + mfile.size();
    }
}

/**
 * Close the file
 */
void LabelFile::closeFile() {
    if (mfile.is_open())
        mfile.close();
}







void LabelFile::write_lst(const vector<Label *> &labels) const {
    
    ofstream out;
    out.open(params.path.c_str());
    
    // FIXME: add correct header information
    out << "SequenceID: 1";
    out << endl;
    
    // FIXME: need to incorporate incomplete genes
    // loop over all labels
    for (size_t n = 0; n < labels.size(); n++) {
        out << n+1 << "\t";                                                 // gene
        out << (labels[n]->strand == Label::POS ? "+" : "-") << "\t";       // strand
        out << labels[n]->left+1 << "\t";                                   // left
        out << labels[n]->right+1 << "\t";                                  // right
        out << labels[n]->right - labels[n]->left + 1 << "\t";              // length
        out << labels[n]->geneClass << "\t";                                // gene class
        out << labels[n]->meta << "\t";                                     // meta
        out << 1 << endl;
    }
    
    out.close();
}



