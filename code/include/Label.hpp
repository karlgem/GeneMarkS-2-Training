//
//  Label.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 8/24/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Label_hpp
#define Label_hpp

#include <stdio.h>
#include <string>

using std::string;

namespace gmsuite {
    
    /**
     * @class Label
     * @brief A class for handling sequence labels (e.g. gene locations in DNA)
     *
     * A Label holds information on the location of fragments (e.g. genes) in DNA. It
     * is usually made up of the left, right, and strand information.
     * where
     *  - left: index of left-end of fragment in DNA (inclusive)
     *  - right: index of right-end of fragment in DNA (inclusive)
     *  - strand: strand of fragment
     *
     * Note: left and right indices start from 0.
     */
    class Label {
        
    public:
        
        typedef enum {POS, NEG, NONE} strand_t;       /**< strand type - positive, negative, none */
        
        /**
         * Constructor: Create a label by defining the left, right, and strand factors
         *
         * @param left the index of the left-end of the fragment (inclusive)
         * @param right the index of the right-end of the fragment (inclusive)
         * @param strand the strand (+,-)
         * @param geneClass the class of the gene
         */
        Label(size_t left, size_t right, strand_t strand, string geneClass = "", string meta = "");
        
        
        /**
         * Get a string representation of the label. By default, label indeces are 0-indexed.
         * However, for compatibility with existing datafile formats, the indexFromOne option
         * allows the labels to be indexed from 1 (i.e. shifted up by 1).
         *
         * @param indexFromOne if set to true, output indeces start from 1 instead of 0
         * @return a string representation of the label
         */
        string toString(bool indexFromOne) const;
        
        
        size_t left;            /**< left-end of the fragment (inclusive) */
        size_t right;           /**< right-end of the fragment (inclusive) */
        strand_t strand;        /**< Strand of the fragment */
        string geneClass;       /**< The gene's class */
        string meta;            /**< Meta information about label */
        
    };
}

#endif /* Label_hpp */
