//
//  Labels.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef Labels_hpp
#define Labels_hpp

#include <vector>
#include <stdio.h>

#include <boost/shared_ptr.hpp>

#include "Label.hpp"

using std::vector;

namespace gmsuite {
    
    /**
     * @class Labels
     * @brief Holds and manages a list of labels
     */
    class Labels {
        
    private:
        
        typedef vector<boost::shared_ptr<Label> > labels_t;     /**< type of labels container */
        labels_t labels;                                        /**< holds the labels */
        
        
    public:
        
        typedef labels_t::size_type size_type;                  /**< Type of labels list size */
        
        /**
         * Create an empty labels list
         */
        Labels();
        
        /**
         * Constructor: create a labels instance from a set of labels
         */
        Labels(const vector<boost::shared_ptr<Label> >& labels);
        
        
        /******************************\
         *          Iterators         *
        \******************************/
         
        
        typedef labels_t::iterator iterator;                    /**< Type of iterator */
        typedef labels_t::const_iterator const_iterator;        /**< Type of const iterator */
        
        iterator begin();                                       /**< Iterator begin */
        iterator end();                                         /**< Iterator end */
        
        const_iterator begin() const;                           /**< Const iterator begin */
        const_iterator end() const;                             /**< Const iterator end */
        
        Label*&  operator[](size_type idx);                     /**< Access an element in the list */
        const Label*  operator[](size_type idx) const;          /**< Access an element in the list */
        
    };
}

#endif /* Labels_hpp */
