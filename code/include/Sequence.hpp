//
//  Sequence.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/21/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Sequence_hpp
#define Sequence_hpp

#include <stdio.h>
#include <string>

using std::string;

namespace gmsuite {
    
    /**
     * @class Sequence
     * @brief This class handles sequences of DNA
     *
     * A Sequence is a container for a real sequence of characters. This interface
     * defines a set of common instructions that are applied to sequences, such as creating
     * a subsequence, getting a string representation, etc...
     *
     * A Sequence is usually made up of an alphabet, the sequence data, and the metadata,
     * where
     *  - alphabet: use for performance reasons
     *  - data: a the sequence's actual data (i.e. representative of its characters)
     *  - metadata: extra information about the sequence (e.g. fasta definitions, ...)
     */
    class Sequence {
        
    public:
        
        typedef string::size_type size_type;                /**< Type of sequence length */
        static const size_type npos = string::npos;         /**< Returned to indicate no matches */
        typedef enum {LINEAR, CIRCLE, NONE} shape_t;        /**< Sequence shape; e.g. of chromosomes */
        
        
        /*************** Constructors and Destructors *******************/
        
        /**
         * Default Constructor: create an empty sequence.
         */
        Sequence();
        
        /**
         * Constructor: create a sequence object from a string
         *
         * @param str the string representation of the sequence
         */
        Sequence(const string &str, const string &meta = "", const shape_t shape = NONE);
        
        /**
         * Destructor
         */
        virtual ~Sequence();
        
        
        
        /*************** Basic Operators (length, toString, ...) *******************/
        
        
        /**
         * Get the size of the sequence (equivalent to sequence length()).
         *
         * @return the size of the sequence.
         */
        virtual size_type size() const;
        
        
        /**
         * Get the length of the sequence (number of characters).
         *
         * @return the length of the sequence.
         */
        virtual size_type length() const;
        
        
        /**
         * Get a string representation of the sequence. E.g. ACCGAT...
         *
         * @return a string representation of the sequence.
         */
        virtual string toString() const;
        
        
        /**
         * Get the meta data (e.g. FASTA definition).
         *
         * @return meta data in string format.
         */
        virtual string getMetaData() const;
        
        
        
        /*************** Sequence Indexing and Substrings *******************/
        
        /**
         * Get a subsequence of the current sequence. This creates a copy
         * of the selected fragment, and a new sequence object is created.
         *
         * Note: The return value is a pointer to an object created by "new".
         * Therefore, it is the caller's responsibility to delete this object
         * When they're done using it, in order to avoid memory leaks.
         *
         * @param n the location of the first character in the subsequence
         * @param length the length of the subsequence
         *
         * @exception std::invalid_argument thrown if n is larger than the sequence's
         * length, or if n + length is larger than the sequence's length.
         *
         * @return a pointer to the subsequence
         */
        virtual Sequence* subseq(size_type n, size_type length) const;
        
        /**
         * Check if string contains a character (can search within specified location)
         *
         * @param c the character
         * @param startSearch where to start searching
         * @param searchLength number of characters to search, starting from startSearch
         *
         * @exception out_of_range if search goes out of range
         *
         * @return true if the sequence contains c
         */
        virtual bool contains(char c, size_type startSearch = 0, size_type searchLength = npos) const;
        
        
        /**
         * Find the first occurrence of a character in the sequence and return its
         * position. If the character is not found, return npos.
         *
         * @param c the character
         * @return the position of the first occurrence of c; npos otherwise.
         */
        virtual size_type find(char c) const;
        
        
        /**
         * Access an element from a const sequence (i.e. cannot be modified)
         * 
         * @param idx the index of the element
         */
        virtual const char&  operator[](size_type idx) const;
        
        /**
         * Access an element from a sequence
         *
         * @param idx the index of the element
         */
        virtual char& operator[](size_type idx);
        
        
        
        
        /*************** Sequence Iterators *******************/
        
        // Iterators
        typedef string::iterator iterator;                  /**< Sequence iterator */
        
        virtual iterator begin();                           /**< Start of iterator */
        virtual iterator end();                             /**< End of iterator   */
        
        
        // Const Iterators
        typedef string::const_iterator const_iterator;      /**< Const sequence iterator */
        
        virtual const_iterator begin() const;               /**< Start of const iterator */
        virtual const_iterator end() const;                 /**< End of const iterator   */
        
        
        
        /*************** Sequence Comparisons *******************/
        
        // overload equality/inequality operators
        virtual inline bool operator==(const Sequence& rhs);    /**< Check if sequence is equal to rhs     */
        virtual inline bool operator!=(const Sequence& rhs);    /**< Check if sequence is not equal to rhs */
        
        // overload less/greater than operators
        virtual inline bool operator<(const Sequence& rhs);     /**< Check if sequence precedes rhs lexicographically */
        virtual inline bool operator>(const Sequence& rhs);     /**< Check if sequence follows rhs lexicographically  */
        
        // overload less/greater than or equal operators
        virtual inline bool operator<=(const Sequence& rhs);    /**< Check if sequence precedes or equals rhs lexicographically */
        virtual inline bool operator>=(const Sequence& rhs);    /**< Check if sequence follows or equals rhs lexicographically  */
        
        
        
    protected:
        
        string data;                /**< The sequence's data */
        string meta;                /**< The sequence's meta information (e.g. fasta definition) */
        shape_t shape;              /**< The sequence's shape */
        string name;                /**< The sequence's name */
        unsigned id;                /**< Automatically generated unique ID for each sequence instance */
        
        
    };
    
}



#endif /* Sequence_hpp */
