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

#define ELPP_STL_LOGGING

#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>

#include <easylogging++.hpp>
#include "formula.hpp"
#include "literal.hpp"
#include "settings.hpp"

INITIALIZE_EASYLOGGINGPP

static hqspre::Formula formula;

void
signalHandler(int /* signum*/)
{
    formula.setInterrupt(true);

    // Remove the signal handler.
    std::signal(SIGINT, nullptr);
    std::signal(SIGUSR1, nullptr);
}

int
main(int argc, char** argv)
{
    hqspre::Timer timer;
    timer.start();

    const std::string sat   = "p cnf 0 0\n";
    const std::string unsat = "p cnf 0 1\n0\n";

    // Configure logging
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Enabled, "true");
    defaultConf.setGlobally(el::ConfigurationType::Format, "[%level] %msg");
    defaultConf.set(el::Level::Verbose, el::ConfigurationType::Format, "[%level-%vlevel] %msg");
    defaultConf.setGlobally(el::ConfigurationType::ToFile, "false");
    defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
#ifdef ELPP_DISABLE_VERBOSE_LOGS
    defaultConf.set(el::Level::Verbose, el::ConfigurationType::ToStandardOutput, "false");
#endif

    hqspre::Settings&            settings = formula.settings();
    std::string                  in_name;
    std::string                  log_name;
    std::string                  out_name;
    std::string                  elim_method("no");
    bool                         enforce_dqbf = false;
    el::base::type::VerboseLevel verbosity    = 1;

    // Preprocessing options
    bool compress_output = true;
#ifdef linux
    unsigned int timeout = 0;
#endif
    bool use_limits = true;

    // Declaration of parameters:
    boost::program_options::options_description public_options("Options");
    public_options.add_options()("help,h", "Show available options")(
        "dqbf", boost::program_options::value<bool>(&enforce_dqbf), "Treat the formula always as a DQBF")(
        "verbosity,v",
        boost::program_options::value<el::base::type::VerboseLevel>(&verbosity)->default_value(verbosity),
        "Set verbosity of the preprocessor")
#ifdef linux
        ("timeout", boost::program_options::value<unsigned int>(&timeout)->default_value(timeout),
         "Set a time limit (0 = none)")
#endif
            ("outfile,o", boost::program_options::value<std::string>(&out_name), "Write the result to the given file")(
                "logfile", boost::program_options::value<std::string>(&log_name), "Write the output to a log file")(
                "compress", boost::program_options::value<bool>(&compress_output)->default_value(true),
                "Rename variables in the output file")(
                "preserve_gates",
                boost::program_options::value<bool>(&settings.preserve_gates)->default_value(settings.preserve_gates),
                "Enable/disable gates protection")("pipe",
                                                   "Suppress status messages and print the formula to the "
                                                   "standard output");

    boost::program_options::options_description pre_options("Preprocessor options");
    pre_options.add_options()(
        "max_loops", boost::program_options::value<std::size_t>(&settings.max_loops)->default_value(settings.max_loops),
        "Set maximal number of preprocessor loops")(
        "univ_reduct",
        boost::program_options::value<bool>(&settings.univ_reduction)->default_value(settings.univ_reduction),
        "Use universal reduction")(
        "bce", boost::program_options::value<unsigned int>(&settings.bce)->default_value(settings.bce),
        "Use blocked clause elimination, 0: no, 1: old, 2: new")(
        "hidden", boost::program_options::value<unsigned int>(&settings.hidden)->default_value(settings.hidden),
        "Add hidden literals before BCE, 0: no, 1: incomplete, 2: complete")(
        "covered", boost::program_options::value<bool>(&settings.covered)->default_value(settings.covered),
        "Add covered literals before BCE")(
        "ble", boost::program_options::value<bool>(&settings.ble)->default_value(settings.ble),
        "Use blocked literal elimination for universal literals")(
        "bla", boost::program_options::value<bool>(&settings.bla)->default_value(settings.bla),
        "Use blocked literal addition for universal literals")(
        "bia", boost::program_options::value<bool>(&settings.bia)->default_value(settings.bia),
        "Use blocked implication addition")(
        "max_clause_size",
        boost::program_options::value<std::size_t>(&settings.max_clause_size)->default_value(settings.max_clause_size),
        "Maximal size of clause using hidden and covered literals")(
        "hse", boost::program_options::value<bool>(&settings.hse)->default_value(settings.hse),
        "Use hidden subsumption elimination")(
        "hec", boost::program_options::value<bool>(&settings.hec)->default_value(settings.hec),
        "Find hidden equivalences and contradictions")(
        "ic", boost::program_options::value<std::uint32_t>(&settings.impl_chains)->default_value(settings.impl_chains),
        "Eliminate implication chains (0=no, 1=strong, 2=weak")(
        "contra", boost::program_options::value<bool>(&settings.contradictions)->default_value(settings.contradictions),
        "Find contradictions")(
        "substitute", boost::program_options::value<bool>(&settings.substitution)->default_value(settings.substitution),
        "Use gate substitution")(
        "equiv_gates", boost::program_options::value<bool>(&settings.equiv_gates)->default_value(settings.equiv_gates),
        "Detect structurally equivalent gates")(
        "sem_gates",
        boost::program_options::value<bool>(&settings.semantic_gates)->default_value(settings.semantic_gates),
        "Use SAT-based semantic gate detection")(
        "extend_gates",
        boost::program_options::value<bool>(&settings.extend_gates)->default_value(settings.extend_gates),
        "Add blocked or transitive implications to find more AND gates")(
        "max_substitution_cost",
        boost::program_options::value<int>(&settings.max_substitution_cost)
            ->default_value(settings.max_substitution_cost),
        "Maximal substitution costs")("max_substitution_loops",
                                      boost::program_options::value<std::size_t>(&settings.max_substitution_loops)
                                          ->default_value(settings.max_substitution_loops),
                                      "Maximal substitution loops")(
        "rewrite", boost::program_options::value<bool>(&settings.rewrite)->default_value(settings.rewrite),
        "Use gate rewriting")(
        "self_sub",
        boost::program_options::value<bool>(&settings.self_subsumption)->default_value(settings.self_subsumption),
        "Eliminate self-subsuming literals")(
        "subsumption", boost::program_options::value<bool>(&settings.subsumption)->default_value(settings.subsumption),
        "Eliminate subsumed clauses")(
        "resolution", boost::program_options::value<bool>(&settings.resolution)->default_value(settings.resolution),
        "Eliminate variables by resolution")(
        "max_resolution_cost",
        boost::program_options::value<int>(&settings.max_resolution_cost)->default_value(settings.max_resolution_cost),
        "Maximal resolution costs")(
        "sat_const", boost::program_options::value<bool>(&settings.sat_const)->default_value(settings.sat_const),
        "Detect constants with SAT-based techniques")(
        "sat_impl", boost::program_options::value<bool>(&settings.sat_impl)->default_value(settings.sat_impl),
        "Detect implications with SAT-based techniques")(
        "pure_sat_timeout", boost::program_options::value<unsigned int>(&settings.pure_sat_timeout)->default_value(20),
        "Timeout for solving pure SAT instances")(
        "univ_exp",
        boost::program_options::value<unsigned int>(&settings.univ_expand)->default_value(settings.univ_expand),
        "Apply universal expansion")(
        "incomplete",
        boost::program_options::value<bool>(&settings.sat_incomplete)->default_value(settings.sat_incomplete),
        "Apply incomplete decision procedures")(
        "random_patterns",
        boost::program_options::value<std::size_t>(&settings.num_random_patterns)
            ->default_value(settings.num_random_patterns),
        "Number of random patterns for incomplete UNSAT checks")(
        "cons_check",
        boost::program_options::value<bool>(&settings.consistency_check)->default_value(settings.consistency_check),
        "Enable/disable consistency checks in preprocessor")(
        "use_limits", boost::program_options::value<bool>(&use_limits)->default_value(use_limits),
        "Enable/disable process limits")(
        "upla", boost::program_options::value<bool>(&settings.upla)->default_value(settings.upla),
        "Use unit propagation lookahead (UPLA)")(
        "vivify", boost::program_options::value<bool>(&settings.vivify)->default_value(settings.vivify),
        "Use vivification to shorten clauses")(
        "vivify_delete",
        boost::program_options::value<bool>(&settings.vivify_delete)->default_value(settings.vivify_delete),
        "Delete clauses using vivification if possible")(
        "vivify_fp", boost::program_options::value<bool>(&settings.vivify_fp)->default_value(settings.vivify_fp),
        "Apply vivification until a fixed point is reached")(
        "vivify_min_size",
        boost::program_options::value<std::size_t>(&settings.vivify_min_size)->default_value(settings.vivify_min_size),
        "Minimal size of clauses to be vivified")(
        "to_qbf", boost::program_options::value<std::string>(&elim_method)->default_value("no"),
        "Try to convert the DQBF into a QBF (no|var|dep)");

    boost::program_options::options_description hidden_options("Options");
    hidden_options.add_options()("infile,i", boost::program_options::value<std::string>(&in_name), "Input file name");

    boost::program_options::options_description all_options("All options");
    all_options.add(public_options).add(hidden_options).add(pre_options);

    boost::program_options::positional_options_description p;
    p.add("infile", -1);

    boost::program_options::variables_map vm;
    boost::program_options::store(
        boost::program_options::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("pipe") > 0) {
        verbosity = 0;
    }

    if (verbosity > 0) {
        std::cout << "This is HQSpre built on " << __DATE__ << " at " << __TIME__ << "\n";
#ifndef NDEBUG
        std::cout << "This is the debug version.\n";
#endif
    }

    // Evaluation of parameters
    if (vm.count("help") > 0) {
        std::cout << "Call with\n  " << argv[0] << "<options> <input file>\n\n"
                  << public_options << "\n"
                  << pre_options << "\n";
        std::exit(0);
    }

    if (vm.count("logfile") > 0) {
        defaultConf.setGlobally(el::ConfigurationType::ToFile, "true");
        defaultConf.setGlobally(el::ConfigurationType::Filename, log_name);
    }

    if (vm.count("pipe") > 0) {
        defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "false");
    }

    el::Loggers::reconfigureAllLoggers(defaultConf);
    el::Loggers::setVerboseLevel(verbosity);

    if (in_name.empty()) {
        LOG(ERROR) << "No input file given.";
        std::exit(-1);
    }

    if (settings.vivify_min_size <= 2) {
        LOG(ERROR) << "The value of vivify_min_size must be at least 3!";
        std::exit(-1);
    }

    if (elim_method == "var") {
        settings.convert_to_qbf = hqspre::DQBFtoQBFmethod::VAR_ELIM;
    } else if (elim_method == "dep") {
        settings.convert_to_qbf = hqspre::DQBFtoQBFmethod::DEP_ELIM;
    } else {
        settings.convert_to_qbf = hqspre::DQBFtoQBFmethod::NONE;
    }

    formula.enforceDQBF(enforce_dqbf);

    std::signal(SIGINT, signalHandler);

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
            try {
                in_filter >> formula;
            } catch (hqspre::FileFormatException& e) {
                LOG(ERROR) << "Bad file format in '" << in_name << "' (" << e.what() << ").";
                in.close();
                return -1;
            }
            in.close();
        } else {
            std::ifstream in(in_name);
            if (!in) {
                LOG(ERROR) << "Could not open input file '" << in_name << "'!";
                std::exit(-1);
            }
            try {
                in >> formula;
            } catch (hqspre::FileFormatException& e) {
                LOG(ERROR) << "Bad file format in '" << in_name << "' (" << e.what() << ").";
                in.close();
                return -1;
            }
            in.close();
        }

        // Start the solving routine
        const auto exist_vars   = formula.numEVars();
        const auto univ_vars    = formula.numUVars();
        const auto clauses      = formula.numClauses();
        const auto dependencies = formula.numDependencies();
        const auto literals     = formula.numLiterals();
        const auto before_qbf   = formula.isQBF();

        formula.preprocess();
#ifdef linux
        if (timeout > 0) {
            hqspre::createTimeout(timeout, signalHandler);
        }
#endif

        if (!formula.isQBF() && settings.convert_to_qbf == hqspre::DQBFtoQBFmethod::DEP_ELIM) {
            const auto elim_set = formula.computeDepElimSet();
            VLOG(2) << "We need to eliminate the following dependencies to obtain a QBF:";
            for (auto y = formula.minVarIndex(); y <= formula.maxVarIndex(); ++y) {
                if (!formula.isExistential(y)) {
                    continue;
                }
                if (y >= elim_set.size()) {
                    break;
                }
                if (elim_set[y].empty()) {
                    continue;
                }
                VLOG(2) << " * " << y << " <-- " << elim_set[y];

                std::vector<hqspre::Variable> current_vars(1, y);
                std::vector<hqspre::Variable> next_vars;

                for (auto x : elim_set[y]) {
                    if (formula.varDeleted(x)) {
                        continue;
                    }
                    while (!current_vars.empty()) {
                        const hqspre::Variable var = current_vars.back();
                        current_vars.pop_back();
                        const auto result = formula.elimDependency(x, var);
                        if (result.first > 0) {
                            next_vars.push_back(result.first);
                        }
                        if (result.second > 0) {
                            next_vars.push_back(result.second);
                        }
                    }
                    std::swap(current_vars, next_vars);
                    formula.applySubstitution();
                }
            }
            if (formula.isQBF()) {
                formula.convertToQBF();
                formula.preprocess();
            } else {
                LOG(ERROR) << "Formula is still not a QBF!";
                std::exit(-1);
            }
        }

        // Output some statistics.
        if (verbosity > 0) {
            formula.printStatistics();
            std::cout << "\nc preprocess -> after / before"
                      << "\nc exist_vars    = " << formula.numEVars() << " / " << exist_vars
                      << "\nc univ_vars     = " << formula.numUVars() << " / " << univ_vars
                      << "\nc literals      = " << formula.numLiterals() << " / " << literals
                      << "\nc clauses       = " << formula.numClauses() << " / " << clauses
                      << "\nc dependies     = " << formula.numDependencies() << " / " << dependencies
                      << "\nc maxVarIndex() = " << formula.maxVarIndex() << "\nc is_qbf        = " << formula.isQBF()
                      << " / " << before_qbf << '\n';
        }

    } catch (hqspre::SATException& e) {
        if (verbosity > 0) {
            formula.printStatistics();
            LOG(INFO) << "Preprocessor determined satisfiability.";
            VLOG(1) << e.what();
        }
        if (verbosity == 0 && vm.count("pipe") == 0) {
            std::cout << "SAT\n";
        }

        if (vm.count("pipe") > 0) {
            std::cout << sat;
            return 0;
        }

        std::ofstream out;
        if (!out_name.empty()) {
            const std::string                   ext = out_name.substr(out_name.length() - 3, 3);
            boost::iostreams::filtering_ostream out_filter;
            if (ext == ".gz" || ext == ".GZ") {
                out.open(out_name, std::ios_base::out | std::ios_base::binary);
                if (!out) {
                    LOG(ERROR) << "Could not open output file '" << out_name << "'!";
                    std::exit(-1);
                }
                out_filter.push(boost::iostreams::gzip_compressor());
            } else {
                out.open(out_name);
                if (!out) {
                    LOG(ERROR) << "Could not open output file '" << out_name << "'!";
                    std::exit(-1);
                }
            }
            out_filter.push(out);
            out_filter << sat;
        }
        return 10;

    } catch (hqspre::UNSATException& e) {
        formula.printStatistics();
        LOG(INFO) << "Preprocessor determined unsatisfiability.";
        VLOG(1) << e.what();
        if (verbosity > 0 && vm.count("pipe") == 0) {
            std::cout << "UNSAT\n";
        }

        if (vm.count("pipe") > 0) {
            std::cout << unsat;
            return 0;
        }

        std::ofstream out;
        if (!out_name.empty()) {
            const std::string                   out_ext = out_name.substr(out_name.length() - 3, 3);
            boost::iostreams::filtering_ostream out_filter;
            if (out_ext == ".gz" || out_ext == ".GZ") {
                out.open(out_name, std::ios_base::out | std::ios_base::binary);
                if (!out) {
                    LOG(ERROR) << "Could not open output file '" << out_name << "'!";
                    std::exit(-1);
                }
                out_filter.push(boost::iostreams::gzip_compressor());
            } else {
                out.open(out_name);
                if (!out) {
                    LOG(ERROR) << "Could not open output file '" << out_name << "'!";
                    std::exit(-1);
                }
            }
            out_filter.push(out);
            out_filter << unsat;
        }
        return 20;
    }

    if (verbosity > 0) {
        LOG(INFO) << "Preprocessor could not solve the formula.";
    }
    if (verbosity == 0 && vm.count("pipe") == 0) {
        std::cout << "UNKNOWN\n";
    }

    if (vm.count("pipe") > 0) {
        formula.write(std::cout, compress_output);
        return 0;
    }

    if (!out_name.empty()) {
        const std::string ext = out_name.substr(out_name.length() - 3, 3);
        if (ext == ".gz" || ext == ".GZ") {
            std::ofstream out(out_name, std::ios_base::out | std::ios_base::binary);
            if (!out) {
                LOG(ERROR) << "Could not open output file '" << out_name << "'!";
                std::exit(-1);
            }
            boost::iostreams::filtering_ostream out_filter;
            out_filter.push(boost::iostreams::gzip_compressor());
            out_filter.push(out);
            formula.write(out_filter, compress_output);
        } else {
            std::ofstream out(out_name);
            if (!out) {
                LOG(ERROR) << "Could not open output file '" << out_name << "'!";
                std::exit(-1);
            }
            formula.write(out, compress_output);
        }
    }

    return 0;
}
