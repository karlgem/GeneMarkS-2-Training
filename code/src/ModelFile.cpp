//
//
//  ModelFile.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/7/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModelFile.hpp"

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/stream.hpp>

using namespace std;
using namespace gmsuite;


string readNextKey(const char*& current, const char* const end);
string readNextValue(const char*& current, const char* const end);

// constructor
ModelFile::ModelFile(string path, access_t access) {
    this->path = path;
    this->access = access;
    
    // setup file parameters
    params.path = path;
    if (access == READ)
        params.flags = io::mapped_file::readonly;           // read-only file
    else if (access == WRITE) {
        params.flags = io::mapped_file::readwrite;          // read-write file
    }
    
    openFile();                                             // open file
}

// destructor
ModelFile::~ModelFile() {
    closeFile();
}


void ModelFile::read(map<string, string> &output) const {
    output.clear();           // clear output map (sanity check)
    
    // point to start of data
    const char* current = begin_read;
    
    string key = "";
    string value = "";
    
    // loop over file
    while (current != end_read) {
        
        // read next key
        key = readNextKey(current, end_read);
        boost::trim(key);
        
        // no key found (i.e. should have read the end
        if (key.empty())
            continue;
        
        // read next value, which consists of all characters up
        // until the next key
        value = readNextValue(current, end_read);
        boost::trim(value);
        
        // add key-value pair to map (note, value can be empty, but key should not be)
        output[key] = value;
    }
}


void ModelFile::read(map<string, string> &output, const vector<string> &keys) const {
    
    const char* current = begin_read;
    
    string key = "";
    string value = "";
    
    // loop over file
    while (current != end_read) {
        
        // read next key
        key = readNextKey(current, end_read);
        boost::trim(key);
        
        // no key found (i.e. should have read the end
        if (key.empty())
            continue;
        
        // check if key is NOT in the query
        if (find(keys.begin(), keys.end(), key) == keys.end())
            continue;
        
        // read next value, which consists of all characters up
        // until the next key
        value = readNextValue(current, end_read);
        boost::trim(value);
        
        // add key-value pair to map (note, value can be empty, but key should not be)
        output[key] = value;
    }
}



string readNextKey(const char*& current, const char* const end) {
    
    string key = "";
    
    // read until next $ sign
    while (current != end && *current != '$')
        current++;
    
    // if reached end of file, return empty key
    if (current == end)
        return key;
    
    // otherwise, current is now at $, so read next keyword
    
    current++;      // skip $
    
    const char* startOfKey = current;
    while (current != end && !isspace(*current))        // loop until next whitespace character
        current++;
    
    const char *endOfKey = current;
    
    key = string(startOfKey, endOfKey);
    return key;
}


string readNextValue(const char*& current, const char* const end) {
    string value = "";
    
    // skip all whitespace characters
    while (current != end && isspace(*current))
        current++;
    
    // if reached end of file, or reached next key, return empty value
    if (current == end || *current == '$')
        return value;
    
    // otherwise, current is now at start of value, so read up until next key
    
    const char* startOfValue = current;
    while (current != end && *current != '$')        // loop until next key
        current++;
    
    const char *endOfValue = current;
    
    value = string(startOfValue, endOfValue);
    return value;

}





/**
 * Read single value from file, for a given key
 *
 * @param key the key
 * @exception invalid_ar
 */
string ModelFile::readValueForKey(string key) const {
    
    
    // point to start of data
    const char* current = begin_read;
    
    string currKey = "";
    string value = "";
    
    // loop over file
    while (current != end_read) {
        
        // read next key
        currKey = readNextKey(current, end_read);
        boost::trim(currKey);
        
        // no key found (i.e. should have read the end
        if (currKey.empty())
            continue;
        
        // if key found, read value
        if (currKey == key) {
        
            // read next value, which consists of all characters up
            // until the next key
            value = readNextValue(current, end_read);
            boost::trim(value);
            return value;
        }
    }
    
    throw invalid_argument("Key not found: " + key);
}


/**
 * Check if the key exists.
 */
bool ModelFile::keyExists(string key) const {
    
    boost::trim(key);
    
    // point to start of data
    const char* current = begin_read;
    
    string currKey = "";
    
    // loop over file
    while (current != end_read) {
        
        // read next key
        currKey = readNextKey(current, end_read);
        boost::trim(currKey);
        
        if (key == currKey)
            return true;
    }
    
    return false;
}


/**
 * Write model parameters to file in key-value pair format.
 */
void ModelFile::write(const vector<pair<string, string> > &keyValue, string sectionHeader) {

    // open filestream
    ofstream out;
    out.open(params.path.c_str());
    
    if (!sectionHeader.empty()) {
        out << "__" << sectionHeader << endl;
    }
    
    for (vector<pair<string, string> >::const_iterator iter = keyValue.begin(); iter != keyValue.end(); iter++) {
        out << "$" << iter->first;
        
        string second = iter->second;
        boost::trim(second);
        
        // if value doesn't contain new line, print it in the same line as the key
        if (second.find('\n') == string::npos)
            out << " " << second << endl;
        // otherwise, start the value on a new line
        else
            out << endl << second << endl << endl;
        
    }
    
    // close filestream
    out.close();
}




/**
 * (Re)open file and reset parameters
 */
void ModelFile::openFile() {
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
        begin_write = mfile.data();
        end_write = begin_write + mfile.size();
    }
}

/**
 * Close file.
 */
void ModelFile::closeFile() {
    if (mfile.is_open())
        mfile.close();
}






