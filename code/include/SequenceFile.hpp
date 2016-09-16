//
//  SequenceFile.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/25/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef SequenceFile_hpp
#define SequenceFile_hpp

#include <stdio.h>
#include <string>
#include <vector>

#include "Sequence.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
namespace io = boost::iostreams;

using std::string;
using std::vector;

namespace gmsuite {
    
    /**
     * @class SequenceFile
     * @brief Represents a file containing one or more sequences
     *
     * This class is used to represent sequence files.
     */
    class SequenceFile {
        
    public:
        
        typedef enum {READ, WRITE} access_t;                    /**< Read or write access to file */
        typedef enum {FASTA, PLAIN, AUTO} format_t;             /**< File format (defines how sequences are read/written) */
        
        /**
         * Constructor: Create a sequence file at a given path, with read/write access,
         * and using a specific format.
         *
         * @param path the path to the file
         * @param access indicates READ or WRITE access permission
         * @param format defines the file format
         */
        SequenceFile(string path, access_t access, format_t format = AUTO);
        
        /**
         * Get format_t representation of the string "format". If the string is not
         * recognized, NONE is returned. Allowed formats include: fasta, plain
         *
         * @param format a string containing a format type.
         * @return a format_t value representing the string
         */
        format_t mapToFormat_t(string format) const;
        
        /**
         * Get the file's format.
         *
         * @return file's format as format_t type.
         */
        format_t getFormat() const;
        
        /**
         * Read sequences from file. This behaves differently for separate
         * file formats:
         *
         * - FASTA: For each fasta header definition (e.g. >...), all subsequent
         *          lines not starting with ">" will be part of the same sequence.
         *          So if the file contains N lines starting with ">", then this method
         *          will return N sequences.
         */
        void read(vector<Sequence> &output) const;
        
        /**
         * Read single sequence from file. This behaves differently for separate 
         * file formats:
         *
         * - FASTA: For the first fasta header definition (e.g. >...), all subsequent
         *          lines not starting with ">" will be constitute a single sequence.
         *          All remaining fasta definitions/sequences will be ignored, since
         *          this function is meant to read a single sequence from a file.
         */
        Sequence read() const;
        
        
        /**
         * Write sequences to file. This behaves differently for separate
         * file formats:
         *
         * - FASTA: before a sequence is written, a FASTA definition line (i.e. >...)
         *          is written to file. The content of the definition line is
         *          derived fromt he Sequence meta data.
         */
        void write(const vector<Sequence> &sequences) const;
        
        
        
    protected:
        
        string path;                /**< full path to file */
        access_t access;            /**< Whether file has READ or WRITE access */
        format_t format;            /**< File'  (sequence) format */
        
        io::mapped_file_params params;      /**< Parameters for mapped file */
        
        io::mapped_file mfile;       /**< Mapped file */
        
        
        const char* begin_read;                 /**< Start of read-only mapped file */
        const char* end_read;                   /**< End of read-only mapped file */
        
        char* begin_write;                /**< Start of readwrite mapped file */
        char* end_write;                  /**< End of readwrite mapped file */
        
        
        
        /**
         * Detect the file's format
         */
        format_t detectFormat() const;
        
        // read file in fasta format
        void read_fasta(vector<Sequence> &output) const;
        
        // write to file in fasta format
        void write_fasta(const vector<Sequence> &sequences) const;
        
        // write to file in plain format (i.e. just sequences each on a line)
        void write_plain(const vector<Sequence> &sequences) const;
        
        
        /**
         * Reopen file and reset parameters
         */
        void openfile();
        
        
    };
}

#endif /* SequenceFile_hpp */
