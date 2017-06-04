//
//  SequenceFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/25/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "SequenceFile.hpp"
#include <assert.h>

#include <fstream>

using namespace std;
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
    
    // setup file parameters
    params.path = path;
    if (access == READ)
        params.flags = io::mapped_file::readonly;           // read-only file
    else if (access == WRITE)
        params.flags = io::mapped_file::readwrite;          // read-write file
    
    // open new file
    openfile();
    
    if (this->format == AUTO && access == WRITE)
        this->format = PLAIN;
    
    // if format not specified (i.e. AUTO), try to guess what it is
    if (this->format == AUTO && access == READ)
        this->format = detectFormat();
    
    
    
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
    
    if (this->format == FASTA)
        write_fasta(sequences);
    else if (this->format == PLAIN)
        write_plain(sequences);
}

// read sequences from file
void SequenceFile::read(vector<Sequence> &output) const {
    
    // read based of format set
    if (this->format == FASTA)
        read_fasta(output);
    else if (this->format == PLAIN)
        read_plain(output);
    
}


Sequence SequenceFile::read() const {
    
    // read based on format set
    if (this->format == FASTA) {
        
        // FIXME: Make this step more efficient and stable. No need to read all fasta defs
        vector<Sequence> output;
        read_fasta(output);
        if (output.size() > 0)
            return output[0];
    }
    else if (this->format == PLAIN) {
        vector<Sequence> output;
        read_plain(output);
        if (output.size() > 0)
            return output[0];
    }
    
    return Sequence();
}


// file should already be opened
SequenceFile::format_t SequenceFile::detectFormat() const {
    
    format_t result = PLAIN;            // assume it's PLAIN format
    
    if (containsFastaDefLine(begin_read, end_read))                 // check if fasta
        result = FASTA;
    
    return result;
}




void SequenceFile::openfile() {
    
    // if file is already open, close it
    if (mfile.is_open()) {
        mfile.close();
    }
    
    // reset pointers
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




void SequenceFile::write_fasta(const vector<Sequence> &sequences) const {
    
    // open a file using fstream
    ofstream out;
    out.open(params.path.c_str());
    
    for (size_t n = 0; n < sequences.size(); n++) {
        out << ">" << sequences[n].getMetaData() << endl;       // write fasta definition
        out << sequences[n].toString() << endl;                 // write sequence
    }
    
    out.close();
}


void SequenceFile::write_plain(const vector<Sequence> &sequences) const {
    // open a file using fstream
    ofstream out;
    out.open(params.path.c_str());
    
    for (size_t n = 0; n < sequences.size(); n++) {
        out << sequences[n].toString() << endl;                 // write sequence
    }
    
    out.close();
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



void SequenceFile::read_plain(vector<Sequence> &output) const {
    
    string seq;        // stores single sequence
    
    output.clear();         // clear output vector (sanity check)
    
    // point to start of data
    const char *  current = begin_read;
    
    // loop over all file
    while (current != end_read) {
        
        // skip all newline characters
        while (current != end_read && isspace(*current))
            current++;
        
        // if reached end of file, break
        if (current == end_read)
            break;
        
        // otherwise, read line
        const char* startOfLine = current;
        while (current != end_read && *current != '\n' && *current != '\r')
            current++;
        
        const char* endOfLine = current;
        
        // create new sequence and add it to the output vector
        output.push_back(Sequence(string(startOfLine, endOfLine)));
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









































