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

#ifndef HQSPRE_SETTINGS_HPP_
#define HQSPRE_SETTINGS_HPP_

#include <iosfwd>
#include <cstddef>

/**
 * \file settings.hpp
 * \brief Settings for DQBF preprocessor.
 * \author Sven Reimer
 * \date 2016
 */

namespace hqspre {

enum class DQBFtoQBFmethod
{
    NONE,
    VAR_ELIM,
    DEP_ELIM
};

/**
 * \brief Represents the settings of the preprocessor
 */
struct Settings
{
    //@{
    /**\name Construction, destruction, and assigment
     */

    /**
     * \brief Constructs an default settings.
     */
    Settings() = default;

    /**
     * \brief Copy constructor
     */
    Settings(const Settings&) = default;

    /**
     * \brief Move constructor
     */
    Settings(Settings&&) = default;

    /**
     * \brief Frees the memory occupied by the settings.
     */
    ~Settings() noexcept = default;

    /**
     * \brief Assignment of settings
     */
    Settings& operator=(const Settings&) = default;

    /**
     * \brief Move assignment
     */
    Settings& operator=(Settings&&) = default;
    //@}

    //@{
    /**
     * \name Variables for storing the current preprocessor configuration
     */
    int          max_num_vars    = 1000000;      ///< Maximal number of variables in the input file
    int          max_num_clauses = 10000000;     ///< Maximal number of clauses in the input file
    std::size_t  max_loops       = 20;           ///< Maximal number of preprocessing iterations
    std::size_t  max_fast_loops  = 10;           ///< Maximal number of iterations in fastPreprocess
    bool         univ_reduction  = true;         ///< Use universal reduction
    unsigned int bce             = 2;            ///< Use blocked clause elimination (BCE)
    std::size_t  min_bce_size    = 3;            ///< Minimum size of a clause for which BCE is performed
    unsigned int hidden          = 1;            ///< Improve BCE by hidden literal addition (HLA)
    bool         covered         = true;         ///< Improve BCE by covered literal addition (CLA)
    bool         ble             = true;         ///< Use blocked literal elimination (BLE) for universal literals
    bool         bla             = true;         ///< Use blocked literal addition (BLA) for universal literals
    bool         bia             = false;        ///< Add blocked implications
    std::size_t  max_clause_size = 256;          ///< Maximal size for an extended clause by
                                                 ///< hidden or covered extension
    bool hse                   = true;           ///< Use hidden subsumption elimination (HSE)
    bool hec                   = false;          ///< Find hidden equivalences and contradiction
    bool impl_chains           = true;           ///< Eliminate implication chains
    bool contradictions        = true;           ///< Detect contradicting implication chains
    bool substitution          = true;           ///< Eliminate Tseitin variables by substitution
    int  max_substitution_cost = 100;            ///< Maximal increase of the number of literals in the formula by
                                                 ///< substitution of gates
    std::size_t max_substitution_loops = 2;      ///< Maximal number of substitution iterations
    bool        rewrite                = true;   ///< Replace Tseitin variables by double Plaisted encoding
    bool        self_subsumption       = true;   ///< Use self-subsuming resolution
    bool        subsumption            = true;   ///< Use subsumption checks
    bool        resolution             = true;   ///< Eliminate variables using resolution
    int         max_resolution_cost    = 1;      ///< Maximal increase of the number of literals
                                                 ///< in the formula by resolution
    bool sat_const      = true;                  ///< Use constant detection with SAT
    bool sat_impl       = false;                 ///< Find implications using SAT
    bool sat_incomplete = true;                  ///< Use incomplete SAT-based checks to determine
                                                 ///< (un)satisfiability of the (D)QBF
    bool            sat_decide          = true;  ///< If the formula is a SAT problem, call a SAT solver to decide it
    std::size_t     num_random_patterns = 100;   ///< Number of random patters for incomplete UNSAT checks
    unsigned int    sat_timeout         = 3;     ///< Timeout for SAT checks if the formula is not a pure SAT problem
    unsigned int    pure_sat_timeout    = 20;    ///< Timeout for SAT check if the formula is a pure SAT problem
    unsigned int    univ_expand         = 2;     ///< Apply universal expansion if reasonable
    DQBFtoQBFmethod convert_to_qbf      = DQBFtoQBFmethod::NONE;  ///< Which method to use to turn a DQBF into an
                                                                  ///< equisatisfiable QBF
    std::size_t max_expansion_size = 0;   ///< Maximal expansion size. This value is set and used during
                                          ///< universal expansion
    bool        vivify          = true;   ///< Apply vivification to shorten clauses
    bool        vivify_delete   = false;  ///< Delete clauses during vivification if possible
    bool        vivify_fp       = false;  ///< Apply vivification until a fixed point is reached
    std::size_t vivify_min_size = 4;      ///< Minimal size of clauses to be vivified
    bool        upla            = false;  ///< Apply unit propagation lookahead

    // Switches related to gate detection
    bool equiv_gates    = true;   ///< Detect structurally equivalent gates
    bool semantic_gates = false;  ///< Use SAT-based semantic gate detection
    bool preserve_gates = false;  ///< Apply techniques like BCE only to clauses
                                  ///< which do not encode gates
    bool extend_gates = false;    ///< Add blocked implications and transitive
                                  ///< implications to get more gates

    bool consistency_check = true;  ///< Check consistence of database after each prepro-function
#ifdef SKOLEM
    bool skolem = false;  ///< Enable the computation of Skolem functions
#endif
    //@}
};

std::ostream& operator<<(std::ostream& stream, const Settings& s);

}  // end namespace hqspre

#endif
