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
    std::istream& operator>>(std::istream& in, ProkGeneStartModel::genome_group_t& unit)
    {
        
        using namespace boost::program_options;
        
        std::string token;
        in >> token;
        if (token == "D")
            unit = ProkGeneStartModel::D;      // class D
        else if (token == "A")
            unit = ProkGeneStartModel::A;      // class A
        else if (token == "B")
            unit = ProkGeneStartModel::B;      // class B
        else if (token == "E")
            unit = ProkGeneStartModel::E;      // class E
        else if (token == "C")
            unit = ProkGeneStartModel::C;      // class C
        else if (token == "A2")
            unit = ProkGeneStartModel::D2;      // arhcaea step 2
        else if (token == "C2")
            unit = ProkGeneStartModel::B2;      // class C2
        else
            throw validation_error(validation_error::invalid_option_value);
        
        return in;
    }
    
//    std::istream& operator>>(std::istream& in, GeneticCode::gcode_t& unit)
//    {
//        std::string token;
//        in >> token;
//        if (token == "11")
//            unit = GeneticCode::ELEVEN;
//        else if (token == "4")
//            unit = GeneticCode::FOUR;
//        //    else
//        //        throw boost::program_options::validation_error("Invalid genome class");
//        
//        return in;
//    }
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              ProkGeneStartModel::genome_group_t * target_type, int)
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
//        v = boost::any(ProkGeneStartModel::genome_group_t(boost::lexical_cast<int>(match[1])));
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
        ;
        
        addProcessOptions(*this, config);
        
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
        
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
    
    
}



void OptionsGMS2Training::addProcessOptions(OptionsGMS2Training &options, po::options_description &processOptions) {
 
    typedef NumSequence::size_type numseqsize;
    
    processOptions.add_options()
    // Coding and NonCoding Models
    ("order-coding",        po::value<unsigned>     (&options.orderCoding                        )->default_value(5),    "Order of coding model")
    ("order-noncoding",     po::value<unsigned>     (&options.orderNonCoding                     )->default_value(2),    "Order of non-coding model")
    ("order-start-context", po::value<unsigned>     (&options.orderStartContext                  )->default_value(2),    "Order of start-context model")
    ("len-start-context",   po::value<numseqsize>   (&options.lengthStartContext                 )->default_value(18),   "Length of start-context model")
    ("margin-start-context",po::value<int>          (&options.marginStartContext                 )->default_value(-15),   "3' Position of start-context model relative to start (negative numbers move downstream into the gene")
    // Misc Variables
    ("fgio-dist-thr",       po::value<numseqsize>   (&options.fgioDistanceThresh                 )->default_value(25),   "Minimum distance between FGIO and upstream gene on same strand")
    ("igio-dist-thr",       po::value<numseqsize>   (&options.igioDistanceThresh                 )->default_value(22),   "Maximum distance between IGIO and upstream gene on same strand")
    ("pcounts",             po::value<unsigned>     (&options.pcounts                            )->default_value(1),    "Pseudocounts")
    ("genome-group",        po::value<genome_group_t> (&options.genomeGroup                      )->required(),          "The genome's group: A,B,C,D,E,A2")
    ("genetic-code",        po::value<gcode_t>      (&options.gcode                              )->default_value(GeneticCode::ELEVEN), "Genetic code")
    ("min-gene-len",        po::value<numseqsize>   (&options.minimumGeneLengthTraining          )->default_value(300),  "Minimym gene length used in training parameters")
    ("only-train-on-native",po::value<bool>         (&options.onlyTrainOnNativeGenes             )->default_value(false),"Only train on native genes")
    ("run-motif-search",    po::value<bool>         (&options.runMotifSearch                     )->default_value(true), "Run motif search")
    // Group-D
    ("gd-width-prom",       po::value<unsigned>     (&options.groupD_widthPromoter               )->default_value(12),   "Group D: promoter width")
    ("gd-width-rbs",        po::value<unsigned>     (&options.groupD_widthRBS                    )->default_value(6),    "Group D: rbs width")
    ("gd-upstr-len-prom",   po::value<numseqsize>   (&options.groupD_upstreamLengthPromoter      )->default_value(40),   "Group D: upstream length for promoter training")
    ("gd-upstr-len-rbs",    po::value<numseqsize>   (&options.groupD_upstreamLengthRBS           )->default_value(20),   "Group D: upstream length for rbs training")
    ("gd-spacer-score-thr", po::value<double>       (&options.groupD_spacerScoreThresh           )->default_value(0.1),  "Group D: minimum peak for valid motif position distribution")
    ("gd-spacer-dist-thr",  po::value<numseqsize>   (&options.groupD_spacerDistThresh            )->default_value(14),   "Group D: minimum distance from start for valid motif")
    ("gd-spacer-window",    po::value<numseqsize>   (&options.groupD_spacerWindowSize            )->default_value(1),    "Group D: size of window to determine peak value")
    ("gd-extended-sd",      po::value<string>       (&options.groupD_extendedSD                  )->default_value("TAAGGAGGTGA"), "Group D: extended SD sequence")
    ("gd-min-match-to-sd",  po::value<unsigned>     (&options.groupD_minMatchToExtendedSD        )->default_value(4),    "Group D: minimum number of consecutive matches to SD")
    ("gd-allow-ag-sub",     po::value<bool>         (&options.groupD_allowAGSubstitution         )->default_value(true), "Group D: allow A-G substitution for valid match to SD")
    // Group C
    ("gc-width-prom",       po::value<unsigned>     (&options.groupC_widthPromoter               )->default_value(6),    "Group C: promoter width")
    ("gc-width-rbs",        po::value<unsigned>     (&options.groupC_widthRBS                    )->default_value(6),    "Group C: rbs width")
    ("gc-upstr-len-prom",   po::value<numseqsize>   (&options.groupC_upstreamLengthPromoter      )->default_value(20),   "Group C: upstream length for promoter training")
    ("gc-upstr-len-rbs",    po::value<numseqsize>   (&options.groupC_upstreamLengthRBS           )->default_value(20),   "Group C: upstream length for rbs training")
    ("gc-spacer-score-thr", po::value<double>       (&options.groupC_spacerScoreThresh           )->default_value(0.25), "Group C: minimum peak for valid motif position distribution")
    ("gc-spacer-dist-thr",  po::value<numseqsize>   (&options.groupC_spacerDistThresh            )->default_value(14),   "Group C: minimum distance from start for valid motif")
    ("gc-spacer-window",    po::value<numseqsize>   (&options.groupC_spacerWindowSize            )->default_value(1),    "Group C: size of window to determine peak value")
    ("gc-extended-sd",      po::value<string>       (&options.groupC_extendedSD                  )->default_value("TAAGGAGGTGA"), "Group C: extended SD sequence")
    ("gc-min-match-to-sd",  po::value<unsigned>     (&options.groupC_minMatchToExtendedSD        )->default_value(4),    "Group C: minimum number of consecutive matches to SD")
    ("gc-allow-ag-sub",     po::value<bool>         (&options.groupC_allowAGSubstitution         )->default_value(true), "Group C: allow A-G substitution for valid match to SD")
    // Group B
    ("gb-width-rbs",        po::value<unsigned>     (&options.groupB_widthRBS                    )->default_value(6),    "Group B: rbs width")
    ("gb-upstr-len-rbs",    po::value<numseqsize>   (&options.groupB_upstreamLengthRBS           )->default_value(20),   "Group B: upstream length for rbs training")
    ("gb-upstr-reg-3-prime",po::value<numseqsize>   (&options.groupB_upstreamRegion3Prime        )->default_value(0),    "Group B: 3'prime end of upstream region used for rbs training")
    ("gb-min-match-rbs-prom",  po::value<unsigned>  (&options.groupB_minMatchRBSPromoter         )->default_value(3),    "Group B: minimum number of consecutive matches between rbs and promoter")
    ("gb-min-match-to-sd",  po::value<unsigned>     (&options.groupB_minMatchToExtendedSD        )->default_value(4),    "Group B: minimum number of consecutive matches to SD")
    ("gb-extended-sd",      po::value<string>       (&options.groupB_extendedSD                  )->default_value("TAAGGAGGTGA"), "Group B: extended SD sequence")
    // Group B2
    ("gb2-width-sd-rbs",        po::value<unsigned>     (&options.groupB2_widthSDRBS             )->default_value(6),    "Group B2: sd rbs width")
    ("gb2-width-non-sd-rbs",    po::value<unsigned>     (&options.groupB2_widthNonSDRBS          )->default_value(6),    "Group B2: non-sd rbs width")
    ("gb2-upstr-len-sd-rbs",    po::value<numseqsize>   (&options.groupB2_upstreamLengthSDRBS    )->default_value(20),   "Group B2: upstream length for SD RBS training")
    ("gb2-upstr-len-non-sd-rbs",po::value<numseqsize>   (&options.groupB2_upstreamLengthNonSDRBS )->default_value(20),   "Group B2: upstream length for non-SD RBS training")
    ("gb2-upstr-reg-3-prime",po::value<numseqsize>   (&options.groupB2_upstreamRegion3Prime      )->default_value(0),    "Group B2: 3'prime end of upstream region used for rbs training")
    ("gb2-min-match-to-sd",  po::value<unsigned>     (&options.groupB2_minMatchToExtendedSD      )->default_value(4),    "Group B2: minimum number of consecutive matches to SD")
    ("gb2-extended-sd",      po::value<string>       (&options.groupB2_extendedSD                )->default_value("TAAGGAGGTGA"), "Group B2: extended SD sequence")
    // Group A
    ("ga-width-rbs",        po::value<unsigned>     (&options.groupA_widthRBS                    )->default_value(6),    "Group A: rbs width")
    ("ga-upstr-len-rbs",    po::value<numseqsize>   (&options.groupA_upstreamLengthRBS           )->default_value(20),   "Group A: upstream length for rbs training")
    ("ga-percent-match-rbs",  po::value<double>     (&options.groupA_percentMatchRBS             )->default_value(0.5),  "Group A: minimum percentage of predicted motifs that match SD")
    ("ga-extended-sd",      po::value<string>       (&options.groupA_extendedSD                  )->default_value("TAAGGAGGTGA"), "Group A: extended SD sequence")
    ("ga-min-match-to-sd",  po::value<unsigned>     (&options.groupA_minMatchToExtendedSD        )->default_value(4),    "Group A: minimum number of consecutive matches to SD")
    ("ga-allow-ag-sub",     po::value<bool>         (&options.groupA_allowAGSubstitution         )->default_value(true), "Group A: allow A-G substitution for valid match to SD")
    // Group X
    ("gx-width-rbs",        po::value<unsigned>     (&options.groupX_widthRBS                    )->default_value(6),    "Group X: rbs width")
    ("gx-upstr-len-rbs",    po::value<numseqsize>   (&options.groupX_upstreamLengthRBS           )->default_value(20),   "Group X: upstream length for rbs training")
    ("gx-len-upstr-sig",    po::value<numseqsize>   (&options.groupX_lengthUpstreamSignature     )->default_value(35),   "Group X: length of upstream signature model")
    ("gx-order-upstr-sig",  po::value<unsigned>     (&options.groupX_orderUpstreamSignature      )->default_value(2),    "Group X: order of upstream signature model")
    ("gx-extended-sd",      po::value<string>       (&options.groupX_extendedSD                  )->default_value("TAAGGAGGTGA"), "Group X: extended SD sequence")
    ("gx-min-match-to-sd",  po::value<unsigned>     (&options.groupX_minMatchToExtendedSD        )->default_value(4),    "Group X: minimum number of consecutive matches to SD")
    ("gx-allow-ag-sub",     po::value<bool>         (&options.groupX_allowAGSubstitution         )->default_value(true), "Group X: allow A-G substitution for valid match to SD")
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

