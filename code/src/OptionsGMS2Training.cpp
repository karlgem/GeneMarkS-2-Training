//
//  OptionsGMS2Training.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsGMS2Training.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>


using namespace boost::xpressive;

using namespace std;
using namespace gmsuite;
namespace po = boost::program_options;

OptionsGMS2Training::OptionsGMS2Training(string mode) : Options(mode), optionsMFinder(mode) {
    
}

namespace gmsuite {
    std::istream& operator>>(std::istream& in, ProkGeneStartModel::genome_class_t& unit)
    {
        std::string token;
        in >> token;
        if (token == "1")
            unit = ProkGeneStartModel::C1;      // class D
        else if (token == "2")
            unit = ProkGeneStartModel::C2;      // class A
        else if (token == "3")
            unit = ProkGeneStartModel::C3;      // class B
        else if (token == "4")
            unit = ProkGeneStartModel::C4;      // class E
        else if (token == "5")
            unit = ProkGeneStartModel::C5;      // class C
        else if (token == "6")
            unit = ProkGeneStartModel::C6;      // arhcaea step 2
    //    else
    //        throw boost::program_options::validation_error("Invalid genome class");
        
        return in;
    }
    
    std::istream& operator>>(std::istream& in, GeneticCode::gcode_t& unit)
    {
        std::string token;
        in >> token;
        if (token == "11")
            unit = GeneticCode::ELEVEN;
        else if (token == "4")
            unit = GeneticCode::FOUR;
        //    else
        //        throw boost::program_options::validation_error("Invalid genome class");
        
        return in;
    }
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              ProkGeneStartModel::genome_class_t * target_type, int)
{
    static sregex expr = sregex::compile("^\\s*(\\d)*");
    using namespace boost::program_options;
    
    // Make sure no previous assignment to 'a' was made.
    validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const string& s = validators::get_single_string(values);
    
    // Do regex match and convert the interesting part to
    // int.
    smatch match;
    if (regex_search(s, match, expr)) {
//        v = boost::any(ProkGeneStartModel::genome_class_t(boost::lexical_cast<int>(match[1])));
    } else {
        throw validation_error(validation_error::invalid_option_value);
    }
}



// parse CMD options
bool OptionsGMS2Training::parse(int argc, const char *argv[]) {
    
    
    try {
        vector<string> config_fnames;                           // holds names of all configuration files (specified by user)
        
        // Declare a group of options that will be allowed
        // only on the command line (CML)
        po::options_description generic("General Options");
        generic.add_options()
        ("version", "Print version string")
        ("help,h", "Display help message")
        ("config", po::value(&config_fnames), "Config file where options may be specified (can be specified more than once)")
        ;
        
        // Declare a group of options that will be allowed
        // on both CML and in the config files
        po::options_description config("Configuration");
        config.add_options()
        ("verbose,v", po::value<int>(&verbose)->default_value(0), "Verbose level")
        ("fn-sequence,s", po::value<string>(&fn_sequence)->required(), "Name of sequence file")
        ("fn-labels,l", po::value<string>(&fn_labels)->required(), "Name of labels file")
        ("fn-mod,m", po::value<string>(&fn_outmod)->required(), "Name of output model file")
        
        
        // GMS2 model parameters
        ("genome-class", po::value<ProkGeneStartModel::genome_class_t>(&genomeClass)->required(), "The genome's class: 1,2,3,4,5,6")
        ("pcounts", po::value<double>(&pcounts)->default_value(1), "Pseudocounts for gms2 models")
        ("coding-order", po::value<unsigned>(&codingOrder)->default_value(4), "Order for coding Markov model")
        ("noncoding-order", po::value<unsigned>(&noncodingOrder)->default_value(2), "Order for noncoding Markov model")
        ("sc-order", po::value<unsigned>(&startContextOrder)->default_value(0), "Order for start-context model")
        ("sc-length", po::value<NumSequence::size_type>(&startContextLength)->default_value(18), "Length of start-context model")
        ("upstream-length", po::value<NumSequence::size_type>(&upstreamLength)->default_value(40), "Length of upstream region for motif search")
        ("MIN_GENE_LEN", po::value<NumSequence::size_type>(&MIN_GENE_LEN)->default_value(300), "Minimum gene length allowed in training")
        ("sc-margin", po::value<int>(&startContextMargin)->default_value(-15), "Margin for start context matrix")
        ("genetic-code", po::value<gcode_t>(&geneticCode)->default_value(GeneticCode::ELEVEN), "Genetic code")
        ("train-on-native-only", po::value<bool>(&trainOnNativeOnly)->default_value(false), "Train on native genes only")
        ("fn-settings", po::value<string>(&fn_settings), "Settings to put in output model file.")
        ("run-motif-search", po::value<bool>(&runMotifSearch)->default_value(true), "Enable/disable motif search.")
        ("upstream-length-fgio", po::value<size_t>(&upstrLenFGIO)->default_value(40), "Upstream length of first-gene-in-operon.")
        ("width-archaea-promoter", po::value<unsigned>(&widthArchaeaPromoter)->default_value(12), "Width for promoters in Archaea.")
        ("match-to", po::value<string>(&matchTo)->default_value("TAAGGAGGTGA"), "16S tail")
        ("allow-ag-substitution", po::bool_switch(&allowAGSubstitution)->default_value(true), "Allow AG substitution.")
        ("match-thresh", po::value<unsigned>(&matchThresh)->default_value(4), "Match threshold for 16S tail.")
        ("upstream-sig-length", po::value<NumSequence::size_type>(&upstreamSignatureLength)->default_value(35), "Length of full upstream signature model")
        ("upstream-sig-order", po::value<unsigned> (&upstreamSignatureOrder)->default_value(2), "Order of upstream signature model")
        ("train-noncoding-on-genome", po::bool_switch(&trainNonCodingOnFullGenome)->default_value(false), "If set, non-coding is trained on full genome, rather than only on intergenic regionns")
        ("fgio-dist-thresh", po::value<unsigned>(&fgioDistThresh)->default_value(25), "Distance threshold to separate first-genes-in-operon")
//        // MFinder options
//        ("pcounts-mfinder", po::value<double>(&optionsMFinder.pcounts)->default_value(1), "Pseudocounts for mfinder models")
//        ("width", po::value<unsigned>(&optionsMFinder.width)->default_value(6), "Width of motif in MFinder")
//        ("motif-order", po::value<unsigned>(&optionsMFinder.motifOrder)->default_value(0), "Order for motif Markov model")
//        ("bkgd-order", po::value<unsigned>(&optionsMFinder.bkgdOrder)->default_value(0), "Order for background model in MFinder")
//        ("align", po::value<string>(&optionsMFinder.align)->default_value("none"), "Set to left or right to allow use of spacer distribution in mfinder")
//        ("tries", po::value<unsigned>(&optionsMFinder.tries)->default_value(10), "Number of tries in mfinder")
        ;
        
        po::options_description prediction("Prediction Parameters");
        prediction.add_options()
        ("NON_P_N", po::value<double>(&nonProbN)->default_value(0.6), "Probability for N value in non-coding region")
        ("COD_P_N", po::value<double>(&codProbN)->default_value(0.4), "Probability for N value in coding region")
        ("NON_DURATION_DECAY", po::value<double>(&nonDurationDecay)->default_value(150), "Duration decay for non-coding region")
        ("COD_DURATION_DECAY", po::value<double>(&codDurationDecay)->default_value(300), "Duration decay for coding region")
        ("GENE_MIN_LENGTH", po::value<NumSequence::size_type>(&geneMinLengthPrediction)->default_value(89), "Minimum length for genes in prediction step")
        ;
        
        config.add(prediction);
        
        // mfinder options
        po::options_description mfinder("Motif Finder");
        OptionsMFinder::addProcessOptions(optionsMFinder, mfinder);
        config.add(mfinder);
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
        ("mode", po::value<string>(&mode)->required(), "Program Mode")
        ;
        
        // Congregate options into further groups
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);
        
        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);
        
        po::options_description visible("Allowed options");
        visible.add(generic).add(config);
        
        // Specify positional arguments
        po::positional_options_description pos;
        pos.add("mode",-1);
        
        
        
        // create storage component for storing names and values of arguments
        po::variables_map vm;
        
        // store command-line options
        po::store(po::command_line_parser(argc, argv).          // pass in input
                  options(cmdline_options).                     // specify options list
                  positional(pos).                              // specify which are positional
                  run(),                                        // parse options
                  vm);                                          // specify storage container
        if (vm.count("config") > 0) {
            config_fnames = vm["config"].as<std::vector<std::string> >();

        
            for (size_t n = 0; n < config_fnames.size(); n++) {
                std::ifstream file(config_fnames[n].c_str());
                
                if(file.fail())
                {
                    std::cerr << "Error opening config file: " << config_fnames[n] << std::endl;
                    return false;
                }
                
                po::store(po::parse_config_file(file, config_file_options), vm);
                file.close();
            }
        }

        
        // if help specified, print usage message and quit
        if (vm.count("help")) {
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // try parsing arguments.
        po::notify(vm);
        
//        // get genome sequence
//        switch (vm["genome-class"].as<int>()) {
//            case 1:
//                genomeClass = ProkGeneStartModel::C1;
//                break;
//            case 2:
//                genomeClass = ProkGeneStartModel::C2;
//            case 3:
//                genomeClass = ProkGeneStartModel::C3;
//                
//            default:
//                throw invalid_argument("Unknown genome class");
//        }
        
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
    
    
}



void OptionsGMS2Training::addProcessOptions(OptionsGMS2Training &options, po::options_description &processOptions) {
 
    processOptions.add_options()
    ("genome-class", po::value<ProkGeneStartModel::genome_class_t>(&options.genomeClass)->required(), "The genome's class: 1,2,3")
    ("pcounts", po::value<double>(&options.pcounts)->default_value(1), "Pseudocounts for gms2 models")
    ("coding-order", po::value<unsigned>(&options.codingOrder)->default_value(4), "Order for coding Markov model")
    ("noncoding-order", po::value<unsigned>(&options.noncodingOrder)->default_value(2), "Order for noncoding Markov model")
    ("sc-order", po::value<unsigned>(&options.startContextOrder)->default_value(0), "Order for start-context model")
    ("sc-length", po::value<NumSequence::size_type>(&options.startContextLength)->default_value(18), "Length of start-context model")
    ("upstream-length", po::value<NumSequence::size_type>(&options.upstreamLength)->default_value(40), "Length of upstream region for motif search")
    ("MIN_GENE_LEN", po::value<NumSequence::size_type>(&options.MIN_GENE_LEN)->default_value(300), "Minimum gene length allowed in training")
    ("sc-margin", po::value<int>(&options.startContextMargin)->default_value(-15), "Margin for start context matrix")
    ("genetic-code", po::value<gcode_t>(&options.geneticCode)->default_value(GeneticCode::ELEVEN), "Genetic code")
    ("train-on-native-only", po::value<bool>(&options.trainOnNativeOnly)->default_value(false), "Train on native genes only")
    ("fn-settings", po::value<string>(&options.fn_settings), "Settings to put in output model file.")
    ("run-motif-search", po::value<bool>(&options.runMotifSearch)->default_value(true), "Enable/disable motif search.")
    ("upstream-length-fgio", po::value<size_t>(&options.upstrLenFGIO)->default_value(40), "Upstream length of first-gene-in-operon.")
    ("width-archaea-promoter", po::value<unsigned>(&options.widthArchaeaPromoter)->default_value(12), "Width for promoters in Archaea.")
    ("upstream-sig-length", po::value<NumSequence::size_type>(&options.upstreamSignatureLength)->default_value(35), "Length of full upstream signature model")
    ("upstream-sig-order", po::value<unsigned> (&options.upstreamSignatureOrder)->default_value(2), "Order of upstream signature model")
    ("train-noncoding-on-genome", po::bool_switch(&options.trainNonCodingOnFullGenome)->default_value(false), "If set, non-coding is trained on full genome, rather than only on intergenic regionns")
    ("fgio-dist-thresh", po::value<unsigned>(&options.fgioDistThresh)->default_value(25), "Distance threshold to separate first-genes-in-operon")
    ;
    
    
    po::options_description prediction("Prediction Parameters");
    prediction.add_options()
    ("NON_P_N", po::value<double>(&options.nonProbN)->default_value(0.6), "Probability for N value in non-coding region")
    ("COD_P_N", po::value<double>(&options.codProbN)->default_value(0.4), "Probability for N value in coding region")
    ("NON_DURATION_DECAY", po::value<double>(&options.nonDurationDecay)->default_value(150), "Duration decay for non-coding region")
    ("COD_DURATION_DECAY", po::value<double>(&options.codDurationDecay)->default_value(300), "Duration decay for coding region")
    ("GENE_MIN_LENGTH", po::value<NumSequence::size_type>(&options.geneMinLengthPrediction)->default_value(89), "Minimum length for genes in prediction step")
    ;
    
    processOptions.add(prediction);
    
    // mfinder options
    po::options_description mfinder("Motif Finder");
    OptionsMFinder::addProcessOptions(options.optionsMFinder, mfinder);
    processOptions.add(mfinder);
    
}

