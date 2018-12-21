//
//  GeneticCode.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/22/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef GeneticCode_hpp
#define GeneticCode_hpp

#include <stdio.h>
#include <vector>
#include <utility>
#include <string>
#include <map>

using std::vector;
using std::string;
using std::pair;
using std::map;

namespace gmsuite {
    
    /**
     * @class GeneticCode
     * @brief A class accessing genetic code based functions
     *
     * A GeneticCode describes the start/stop codons, as well as codon to amino acid
     * translation table, based on a genetic code description.
     */
    class GeneticCode {
        
    public:
        
        typedef enum {FOUR, ELEVEN, CUSTOM} gcode_t;                        /**< Genetic code type */
        
        typedef vector<string>::const_iterator ttk_const_iterator;          /**< Iterator over translation table keys (e.g. codons) */
        
        
        /**
         * Constructor: initialize genetic code 11
         */
        GeneticCode();
        
        /**
         * Constructor: initialize genetic code based on the gcode value. This method cannot be used
         * for custom genetic codes. Please use the relevant constructor for custom configurations.
         *
         * @param gcode the genetic code (cannot be CUSTOM)
         */
        GeneticCode(gcode_t gcode);
        
        
        /**
         * Constructor: initialize custom genetic code with amino acid translation
         * table, as well as the start and stop codons.
         *
         * @param table a vector of string-char pairs, where string is the 3 letter codon, and char is the amino acid
         * @param starts vector of strings, where each string is a 3-letter start codon
         * @param stops vector of strings, where each string is a 3-letter stop codon
         */
        GeneticCode(const vector<pair<string, char> > &table, const vector<string> &starts, const vector<string> &stops);
        
        
        /******************** Request info on Starts, Stops, and Translating AA ********************/
        
        /**
         * Check if codon is a start codon
         *
         * @param codon candidate codon
         * @return true if codon is a start; false otherwise
         */
        bool isStart(string codon) const;
        
        
        /**
         * Check if codon is a stop codon
         *
         * @param codon candidate codon
         * @return true if codon is a stop; false otherwise
         */
        bool isStop(string codon) const;
        
        
        /**
         * Get the list of start codons
         *
         * @return a vector of strings, where each string is a 3-letter codon representing a start
         */
        vector<string> getStarts() const;
        
        
        /**
         * Get the list of stop codons
         *
         * @return a vector of strings, where each string is a 3-letter codon representing a stop
         */
        vector<string> getStops() const;
        
        
        /**
         * Get the translation table for the current genetic code
         *
         * @return a map of string-char pairs, showing relationship between codon and amino acid
         */
        map<string, char> getTranslationTable() const;
        
        /**
         * Translate 3-letter codon to amino acid
         *
         * @param codon the codon to be translated
         * @return amino acid character
         */
        char translateCodon(string codon) const;
        
        
        /**
         * Get the genetic code value that represents this object.
         *
         * @return the genetic code.
         */
        gcode_t getGCode() const;
        
        
        /**
         * Get first iterator over translation table keys. This iterates over
         * the codons (e.g. ATT, AGT, TGA, ...). You can get the corresponding
         * amino acid by using the conversion method. @see translateCodon
         *
         * @return first iterator over translation table keys.
         */
        ttk_const_iterator begin() const;
        
        /**
         * Get end iterator over translation table keys. 
         *
         * @return last (end) iterator over translation table keys.
         */
        ttk_const_iterator end() const;
        
        /**
         * Get name of genetic code
         */
        string getName() const;
        
        
    protected:
        
        gcode_t gcode;                                      /**< Genetic code value     */
        string name;                                        /**< Genetic code name      */
        vector<string> starts;                              /**< List of start codons   */
        vector<string> stops;                               /**< List of stop codons    */
        map<string, char> translationTable;                 /**< Translation table of codons to amino-acids */
        vector<string>  translationTableKeys;               /**< Stores the keys (codons) of the translation table in vector format (for iterator) */
        
        
        /**
         * Initialize translation table and start/stop codons based on genetic code
         */
        void initialize();
        
        /**
         * Initialize translation table and start/stop codons for genetic code 11
         */
        void initialize11();
        
        /**
         * Initialize translation table and start/stop codons for genetic code 4
         */
        void initialize4();
        
        /**
         * Fill the translation table and identify start and stop codons. To do this, 5 (parallel) strings
         * are provided: (1) amino acids, (2) starts, (3) base1, (4) base2, and (5) base3. All strings should
         * have the same length. An example input of the strings is:
         *
         * (1) FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
         * (2) ---M-------------------------------M---------------M------------
         * (3) TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
         * (4) TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
         * (5) TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
         *
         * Codons are identified by the base strings (i.e. strings (3), (4), and (5)). For example, codon TTT
         * has amino acid F. Codon TAC has amino acid C, etc.... Stop codons are identified by * in string (1).
         * Finally, starts are identified in string (2), where M indicates that this codon is a start codon.
         *
         * @param AAs string of amino acids. i'th location indicates i'th amino acid. The "*" character
         *            indicates a stop codon.
         * @param starts if i'th location is a 'M' character, then it indicates a start.
         *
         */
        void fillTranslationTable(string AAs, string starts, string base1, string base2, string base3);
    };
    
}

#endif /* GeneticCode_hpp */
