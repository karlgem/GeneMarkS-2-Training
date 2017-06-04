//
//  Genome.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 7/21/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Genome_hpp
#define Genome_hpp

#include <stdio.h>
#include <string>
#include <vector>

#include "Sequence.hpp"
#include "GeneticCode.hpp"

using std::string;
using std::vector;

namespace gmsuite {
    
    /**
     * @class Genome
     * @brief A class for storing full genomes
     *
     * A Genome is a container that stores the sequences representing the genome. Note that a genome
     * can be composed of multiple, separate sequences; each of these is represented
     * as a Sequence object.
     */
    class Genome {
        
    public:
        
        typedef string::size_type size_type;            /**< Type of genome length */
        static const size_type npos = string::npos;     /**< Returned to indicate no matches */
        
        
        /*************** Constructors and Destructors *******************/
        
        /**
         * Default constructor: create an empty genome
         *
         * @param name the genome's name (default: empty)
         */
        Genome(string name = "");
        
        
        /**
         * Default constructor: create an empty genome with genetic code
         *
         * @param gcode the genome's genetic code
         * @param name the genome's name (default: empty)
         */
        Genome(GeneticCode gcode, string name = "");
        
        
        /**
         * Default constructor: create a genome with sequences and genetic code
         *
         * @param sequences the sequences making up the genome
         * @param gcode the genome's genetic code
         * @param name the genome's name (default: empty)
         */
        Genome(const vector<Sequence> &sequences, GeneticCode gcode, string name = "");
        
        
        /**
         * Destructor
         */
        ~Genome();
        
        
        
        /*************** Indexing genome (sequences) and adding sequences *******************/
        
        /**
         * Access an a const Sequence from the genome (i.e. cannot be modified)
         *
         * @param idx the index of the element
         */
        const Sequence&  operator[](size_type idx) const;
        
        
        /**
         * Add a sequence to the genome.
         *
         * @param seq a new sequence
         */
        void add(Sequence seq);
        
        
        /**
         * Get the number of sequences in the genome
         *
         * @return the number of sequences.
         */
        size_t numOfSeq() const;
        
        
        /*************** Constructors and Destructors *******************/
        


        
    protected:
        
        unsigned id;                        /**< Automatically generated unique ID for each Genome instance */
        string name;                        /**< The genome's name */
        GeneticCode gcode;                  /**< The genome's genetic code */
        
        vector<Sequence> sequences;         /**< Vector of sequences constituting the genome */
        
    };
    
}

#endif /* Genome_hpp */
