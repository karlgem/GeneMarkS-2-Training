//
//  ModelFile.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/7/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef ModelFile_hpp
#define ModelFile_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>

#include <boost/iostreams/device/mapped_file.hpp>
namespace io = boost::iostreams;

using std::map;
using std::string;
using std::vector;

namespace gmsuite {
    
    /**
     * @class Model File
     * @brief Represents a file containing model parameters
     *
     * This class is used to handle model files, in which model parameters are stored.
     * The parameters are in a key-value pair format, where keys are signified by
     * the $ sign.
     */
    class ModelFile {
        
    public:
        
        typedef enum {READ, WRITE} access_t;                    /**< Read or write access to file */
        
        /**
         * Constructor: Create a model file instance at a given path, with read/write access.
         *
         * @param path the path to the file
         * @param access indicates READ or WRITE access permission
         */
        ModelFile(string path, access_t access);
        
        /**
         * Destructor: close file
         */
        ~ModelFile();
        
        /**
         * Read all key-value pairs from the file
         *
         * @param output a map of key-value pair strings where the output is placed
         */
        void read(map<string, string> &output) const;
        
        /**
         * Read key-value pairs from the file, for a specific list of keys
         *
         * @param output a map of key-value pair strings where the output is placed
         */
        void read(map<string, string> &output, const vector<string> &keys) const;
        
        /**
         * Read single value from file, for a given key
         *
         * @param key the key
         * @exception invalid_ar
         */
        string readValueForKey(string key) const;
        
        
        /**
         * Check if the key exists.
         */
        bool keyExists(string key) const;
        
        
        /**
         * Write model parameters to file in key-value pair format.
         */
        void write(const map<string, string> &keyValue) const;
        
        
        
    protected:
        
        string path;                /**< full path to file */
        access_t access;            /**< Whether file has READ or WRITE access */
        
        io::mapped_file_params params;      /**< Parameters for mapped file */
        
        io::mapped_file mfile;       /**< Mapped file */
        
        
        const char* begin_read;                 /**< Start of read-only mapped file */
        const char* end_read;                   /**< End of read-only mapped file */
        
        char* begin_write;                /**< Start of readwrite mapped file */
        char* end_write;                  /**< End of readwrite mapped file */
        

        /**
         * (Re)open file and reset parameters
         */
        void openFile();
        
        /**
         * Close file.
         */
        void closeFile();
        
        
    };
}


#endif /* ModelFile_hpp */
