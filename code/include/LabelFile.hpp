//
//  LabelFile.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/31/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef LabelFile_hpp
#define LabelFile_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <boost/iostreams/device/mapped_file.hpp>

#include "Label.hpp"

using std::string;
using std::vector;
namespace io = boost::iostreams;

namespace gmsuite {
    
    /**
     * @class LabelFile
     * @brief Represents a file containing labels of genetic information
     *
     * This class is used to represent a file containing sequence labels; i.e.
     * left/right coordinates of fragments in a larger sequence, along with other 
     * appropriate properties
     */
    class LabelFile {
        
    public:
        
        typedef enum {READ, WRITE} access_t;                /**< Read or write access to file */
        typedef enum {LST, AUTO} format_t;                  /**< File format (defines how labels are read/written) */
        
        /**
         * Constructor: Create a LabelFile instance for a file at a given path, with
         * read/write access, and using a specific format.
         *
         * @param path the path to the file
         * @param access indicates READ or WRITE access permission
         * @param format defines the file's format
         */
        LabelFile(string path, access_t access, format_t format = AUTO);
        
        /**
         * Destructor: free memory and close the file
         */
        ~LabelFile();
        
        
        /**
         * Read labels from file. This behaves differently for separate file formats.
         *
         * @param output a vector label pointers that have been read from the file.
         */
        void read(vector<Label*> &output) const;
        
        
        /**
         * Write labels to file. 
         *
         * @param labels the vector of labels to be written to file
         */
        void write(const vector<Label*> &labels) const;
        
        
    private:
        
        /**
         * Detect the file's format. This assumes an already opened file
         *
         * @return the file's format.
         */
        format_t detectFormat() const;
        
        /**
         * Check if the file is in LST format. This assumes an already opened file
         *
         * @return true if the file is in LST format; false otherwise
         */
        bool detectLST(const char* const begin, const char* const end) const;
        
        /**
         * Open the file and set start/end pointers to the data.
         */
        void openFile();
        
        /**
         * Close the file
         */
        void closeFile();
        
        
        /**
         * Read labels from LST file.
         *
         * @param output a vector label pointers that have been read from the file.
         */
        void read_lst(vector<Label*> &output) const;
        
        /**
         * Write labels to LST file.
         *
         * @param labels the vector of labels to be written to the file
         */
        void write_lst(const vector<Label*> &labels) const;
        
        
        
        
        
        
        string path;                        /**< full path to file */
        access_t access;                    /**< Whether file has READ or WRITE access */
        format_t format;                    /**< File'  (sequence) format */
        
        io::mapped_file mfile;              /**< Mapped file */
        io::mapped_file_params params;      /**< Parameters for mapped file */
        
        
        const char* begin_read;             /**< Start of read-only mapped file */
        const char* end_read;               /**< End of read-only mapped file */
        
        char* begin_write;                  /**< Start of readwrite mapped file */
        char* end_write;                    /**< End of readwrite mapped file */
        
    };
}

#endif /* LabelFile_hpp */
