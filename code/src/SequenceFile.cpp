//
//  SequenceFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/25/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "SequenceFile.hpp"
#include <assert.h>
using namespace gmsuite;


// function prototypes

bool containsFastaDefLine(const char* const begin, const char* const end);      // check for fasta definition at start of file

// read file per format
string readNextFastaDefLine(const char*& current, const char* const end);
string readNextFastaSequence(const char*& current, const char* const end);

// write file per format
void writeNextFastaDefLine(const char*& current, const char* const end);
void writeNextFastaSequence(const char*& current, const char* const end);


SequenceFile::SequenceFile(string path, access_t access, format_t format) {
    this->path = path;
    this->access = access;
    this->format = format;
    
    // if format not specified (i.e. AUTO), try to guess what it is
    if (this->format == AUTO)
        this->format = detectFormat();
    
    // setup file parameters
    params.path = path;
    if (access == READ)
        params.flags = io::mapped_file::readonly;           // read-only file
    else if (access == WRITE)
        params.flags = io::mapped_file::readwrite;          // read-write file
    
    // open new file
    openfile();
    
    
    
}



SequenceFile::format_t SequenceFile::mapToFormat_t(string format) const {
    // fasta format
    if (format == "fasta")
        return FASTA;
    
    // if unknown format is set (or empty string), then set as AUTO
    else {
        return AUTO;
    }
}


SequenceFile::format_t SequenceFile::getFormat() const {
    return this->format;
}




void SequenceFile::write(const vector<Sequence> &sequences) const {
    
    assert(false);                      // don't use: still not complete
    if (this->format == FASTA)
        write_fasta(sequences);
}

// read sequences from file
void SequenceFile::read(vector<Sequence> &output) const {
    
    // read based of format set
    if (this->format == FASTA)
        read_fasta(output);
    
}



SequenceFile::format_t SequenceFile::detectFormat() const {
    
    format_t result = PLAIN;            // assume it's PLAIN format
    
    if (containsFastaDefLine(begin_read, end_read))
        result = FASTA;
    
    return result;
}




void SequenceFile::openfile() {
    
    // if file is already open, close it
    if (mfile.is_open()) {
        mfile.close();
    }
    
    // open file
    mfile.open(params);
    
    // reset pointers
    if (access == READ) {
        begin_read = mfile.const_data();
        end_read = begin_read + mfile.size();
    }
    else if (access == WRITE) {
        begin_write = mfile.data();
        end_write = begin_write + mfile.size();
    }
    
    
}




void SequenceFile::write_fasta(const vector<Sequence> &sequences) const {
    
    
    for (size_t n = 0; n < sequences.size(); n++) {
        
        // write fasta definition
        
    }
}






// check for fasta definition at start of file
bool containsFastaDefLine(const char* const begin, const char* const end) {
    
    // point to start of data
    const char *  current = begin;
    
    // skip starting all space characters
    while (current != end && isspace(*current))
        current++;
    
    // if we've reached the end of the file or the first non-space character is NOT a '>'
    // then this is not a FASTA file
    if (current == end || *current != '>')
        return false;
    
    return true;
}


// read sequences from fasta file. This assumes that the file is indeed in FASTA format (i.e. check first)
void SequenceFile::read_fasta(vector<Sequence> &output) const{
    
    string fastaDef;        // stores single fasta definition
    string fastaSeq;        // stores single fasta sequence
    
    output.clear();         // clear output vector (sanity check)
    
    // point to start of data
    const char *  current = begin_read;
    
    // loop over all file
    while (current != end_read) {
        
        // read next fasta definition line
        fastaDef = readNextFastaDefLine(current, end_read);
        
        // read (entire) fasta sequence. This may constitute multiple
        // lines. The function will read until the next fasta header definition
        // is found
        fastaSeq = readNextFastaSequence(current, end_read);
        
        // create new sequence and add it to the output vector
        output.push_back(Sequence(fastaSeq, fastaDef));
    }
    
}











string readNextFastaDefLine(const char*& current, const char* const end) {
    string output = "";
    
    // skip starting all space characters
    while (current != end && isspace(*current))
        current++;
    
    // if we've reached the end of the file or the first non-space character is NOT a '>'
    // then this is not a FASTA file
    if (current == end || *current != '>')
        return output;
    
    // the current character is now '>'. Skip it
    current++;
    const char* startOfDef = current;                   // beginning of fasta definition line
    
    // read the remaining definition line until newline character is reached
    while (current != end && *current != '\n' && *current != '\r')
        current++;
    
    const char* endOfDef = current;
    
    // create string from fasta definition
    output = string(startOfDef, endOfDef);
    
    return output;
}

string readNextFastaSequence(const char*& current, const char* const end) {
    string output = "";
    
    // skip starting all space characters
    while (current != end && isspace(*current))
        current++;
    
    // if we've reached the end of the file or the first non-space character IS a '>'
    // then no sequence exists
    if (current == end || *current == '>')
        return output;
    
    // the current character is now the first character of the sequence
    const char* startOfLine = current;
    
    // read the remaining sequence lines until the next '>' character is reached, or the end of file is reached
    while (current != end && *current != '>') {
        
        // if the character is an whitespace (or newline, ...) then extract the current fragment
        if (isspace(*current)) {
            
            output += string(startOfLine, current);     // extract fragment and append it to existing fragments
            
            current++;
            // skip all whitespaces
            while (current != end && isspace(*current))
                current++;
            
            startOfLine = current;          // new line
        }
        // otherwise, character is element of sequence
        else {
            current++;
        }
    }
    
    // append the last batch of characters (if any)
    output += string(startOfLine, current);
    
    return output;
}









































