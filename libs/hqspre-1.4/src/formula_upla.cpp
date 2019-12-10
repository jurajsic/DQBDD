#include "formula.hpp"

#include <algorithm>
#include <vector>
#include <easylogging++.hpp>

#include "propagator.hpp"

/**
 * \file formula_upla.cpp
 * \author Ralf Wimmer
 * \date 02/2018
 * \brief Implementation of unit propagation lookahead
 */

namespace hqspre {


/**
 * \brief Performs unit propagation lookahead to find units, equivalences and subsuming binary clauses
 *
 * \return true if the formula was changed by UPLA
 */
bool Formula::upla()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer upla(getTimer(WhichTimer::UPLA));

    std::vector<Literal> units;
    std::vector<std::vector<Literal>> to_add;

    SatPropagator propagate_pos(_prefix, _clauses, _implications);
    SatPropagator propagate_neg = propagate_pos;

    std::vector<Literal> bin_clause(2, 0);

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (_interrupt) break;
        if (varDeleted(var)) continue;

        const Literal pos_lit = var2lit(var, false);
        const Literal neg_lit = var2lit(var, true);
        bool conflict = false;

        propagate_pos.clear();
        try {
            propagate_pos.bcp(pos_lit);
        } catch (ConflictException&) {
            // an exception has occurred -> the negation of the current literal,
            // i.e., the negative literal, is unit.
            units.push_back(neg_lit);
            conflict = true;
        }
        if (conflict) continue;

        std::vector<Literal>& pos_impl = propagate_pos.getUnits();
        std::sort(pos_impl.begin(), pos_impl.end());

        propagate_neg.clear();
        try {
            propagate_neg.bcp(neg_lit);
        } catch (ConflictException&) {
            // an exception has occurred -> the negation of the current literal,
            // i.e., the positive literal, is unit.
            units.push_back(pos_lit);
            conflict = true;
        }
        if (conflict) continue;

        std::vector<Literal>& neg_impl = propagate_neg.getUnits();
        std::sort(neg_impl.begin(), neg_impl.end());

        // check for subsuming binary clauses for the positive implications
        bin_clause[0] = neg_lit;
        for (const Literal implied: pos_impl) {
            if (implied == pos_lit) continue;
            bin_clause[1] = implied;
            std::uint64_t sig = 0;
            addSignatureLit(sig, bin_clause[0]);
            addSignatureLit(sig, bin_clause[1]);
            if (!hasImplication(pos_lit, implied) && isBackwardSubsuming(bin_clause, sig, -1, false)) {
                to_add.push_back(bin_clause);
            }
        }

        // check for subsuming binary clauses for the negative implications
        bin_clause[0] = pos_lit;
        for (const Literal implied: neg_impl) {
            if (implied == pos_lit) continue;

            bin_clause[1] = implied;
            std::uint64_t sig = 0;
            addSignatureLit(sig, bin_clause[0]);
            addSignatureLit(sig, bin_clause[1]);
            if (!hasImplication(neg_lit, implied) && isBackwardSubsuming(bin_clause, sig, -1, false)) {
                to_add.push_back(bin_clause);
            }
        }


        // Check for units and equivalences
        auto pos_ptr = pos_impl.cbegin();
        auto neg_ptr = neg_impl.cbegin();
        while (pos_ptr != pos_impl.cend() && neg_ptr != neg_impl.cend()) {
            if (lit2var(*pos_ptr) == lit2var(*neg_ptr)) {
                if (lit2var(*pos_ptr) != var) {
                    if (*pos_ptr == *neg_ptr) {
                        units.push_back(*pos_ptr);
                    }  else {
                        bin_clause[0] = neg_lit;
                        bin_clause[1] = *pos_ptr;
                        to_add.push_back(bin_clause);

                        bin_clause[0] = pos_lit;
                        bin_clause[1] = *neg_ptr;
                        to_add.push_back(bin_clause);
                    }
                }
                ++pos_ptr;
                ++neg_ptr;
            } else {
                if (*pos_ptr < *neg_ptr) ++pos_ptr;
                else ++neg_ptr;
            }
        }


    } // end for (Variable var ... )

    VLOG(2) << __FUNCTION__ << " found " << units.size() << " unit literals and added " << to_add.size() << " binary clauses.";

    for (auto& clause: to_add) addClause(std::move(clause));
    for (const Literal lit: units) pushUnit(lit, PureStatus::UNIT);
    fastPreprocess();

    return !units.empty() && !to_add.empty();
}

} // end namespace hqspre
