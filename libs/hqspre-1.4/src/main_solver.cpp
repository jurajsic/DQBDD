/*
 * This file is part of HQSpre.
 *
 * Copyright 2016/17 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
 * Albert-Ludwigs-Universitaet Freiburg, Freiburg im Breisgau, Germany
 *
 * HQSpre is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <easylogging++.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>

#include "formula.hpp"
#include "literal.hpp"
#include "settings.hpp"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char** argv)
{
    hqspre::Timer timer;
    timer.start();

    std::cout << "This is HQSpre_solver built on " << __DATE__ << " at " << __TIME__ << "\n";
#ifndef NDEBUG
    std::cout << "This is the debug version.\n";
#endif

    hqspre::Formula formula;
    hqspre::Settings& settings = formula.settings();
    std::string in_name;
    unsigned int num_exp = 1;
    unsigned int num_pre_iter = 2;
    int subst_increase = 0;
    int resol_increase = 0;
    el::base::type::VerboseLevel verbosity = 0;


    boost::program_options::options_description public_options("Options");
    public_options.add_options()
      ("help,h", "Show available options")
      ("verbosity,v", boost::program_options::value<el::base::type::VerboseLevel>(&verbosity)->default_value(verbosity), "Set verbosity of the preprocessor")
#ifndef NDEBUG
      ("cons_check",  boost::program_options::value<bool>(&settings.consistency_check)->default_value(settings.consistency_check),    "Enable/disable consistency checks in preprocessor")
#endif
      ("num_exp",    boost::program_options::value<unsigned int>(&num_exp)->default_value(num_exp), "Number of universal expansions per iteration")
      ("num_pre_iter", boost::program_options::value<unsigned int>(&num_pre_iter)->default_value(num_pre_iter), "Number of preprocessing iteration between expansions")
      ("max_resolution_cost", boost::program_options::value<int>(&settings.max_resolution_cost)->default_value(settings.max_resolution_cost),
            "Maximal costs of resolutions")
      ("resolution_cost_increase", boost::program_options::value<int>(&resol_increase)->default_value(resol_increase),
            "Increase of resolution costs per iteration")
      ("max_substitution_cost", boost::program_options::value<int>(&settings.max_substitution_cost)->default_value(settings.max_substitution_cost),
            "Maximal costs of substitutions")
      ("substitution_cost_increase", boost::program_options::value<int>(&subst_increase)->default_value(subst_increase),
            "Increase of substitution costs per iteration")
      ("random_patterns", boost::program_options::value<std::size_t>(&settings.num_random_patterns)->default_value(100), "Number of random patterns for incomplete UNSAT checks")
    ;

    boost::program_options::options_description hidden_options("Options");
    hidden_options.add_options()
      ("infile,i", boost::program_options::value<std::string>(&in_name), "Input file name")
    ;

    boost::program_options::options_description all_options("All options");
    all_options.add(public_options).add(hidden_options);

    boost::program_options::positional_options_description p;
    p.add("infile", -1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cout << "Call with\n  " << argv[0] << " <options> <input file>\n\n" << public_options << "\n";
        std::exit(0);
    }

    if (!vm.count("infile")) {
       LOG(ERROR) << "Invalid syntax. Call with\n\t" << argv[0] << " <options> <input file>\n\n" << public_options;
       std::exit(0);
    }

    try {
        const std::string ext = in_name.substr(in_name.length() - 3, 3);

        if (ext == ".gz" || ext == ".GZ") {
        std::ifstream in(in_name, std::ios_base::in | std::ios_base::binary);
         if (!in) {
                LOG(ERROR) << "Could not open input file '" << in_name << "'!";
                std::exit(-1);
            }
            boost::iostreams::filtering_istream in_filter;
            in_filter.push(boost::iostreams::gzip_decompressor());
            in_filter.push(in);
            in_filter >> formula;
            in.close();
        } else {
            std::ifstream in(in_name);
            if (!in) {
                LOG(ERROR) << "Could not open input file '" << in_name << "'!";
                std::exit(-1);
            }
            in >> formula;
            in.close();
        }

        // Start the solving routine
        settings.max_loops = 10;
        settings.rewrite = true;
        settings.pure_sat_timeout = 0;
        settings.univ_expand = 2;

        do {
            formula.preprocess();

            for (unsigned int i = 0; i < num_exp; ++i) {
                hqspre::Variable next = 0;
                long int costs = std::numeric_limits<long int>::max();
                for (hqspre::Variable var = formula.minVarIndex(); var <= formula.maxVarIndex(); ++var) {
                    if (!formula.isUniversal(var)) continue;

                    std::set<hqspre::Variable> pseudo_deps;
                    formula.sstdQuadDep(var, true, true, &pseudo_deps);
                    const long int new_cost = formula.computeExpansionCosts2(hqspre::var2lit(var), pseudo_deps);
                    pseudo_deps.clear();
                    if (new_cost < costs) {
                        costs = new_cost;
                        next = var;
                    }
                }
                if (next > 0) {
//                    VLOG(2) << formula.numUVars() << " universal variables left (formula size = " << formula.numLiterals() << " literals). Expanding variable " << next;
                    formula.elimUVar(hqspre::var2lit(next, false));
                }
            }

            settings.max_loops = num_pre_iter;
            settings.max_resolution_cost += resol_increase;
            settings.max_substitution_cost += subst_increase;
        } while (formula.numUVars() > 0);

        formula.checkSAT();
    } catch (hqspre::SATException&) {
        std::cout << "SAT\n";
        timer.stop();
        VLOG(1) << "Time: " << timer.read();
        std::exit(10);
    } catch (hqspre::UNSATException&) {
        std::cout << "UNSAT\n";
        timer.stop();
        VLOG(1) << "Time: " << timer.read();
        std::exit(20);
    }

    std::cout << "UNKNOWN\n";
    timer.stop();
    VLOG(1) << "Time: " << timer.read();
    return 0;
}
