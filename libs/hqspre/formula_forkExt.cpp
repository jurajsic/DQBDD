#include <queue>

#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"

#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#define TESTFORK


namespace hqspre
{

/**
 * \brief Perform Forksplitting
 * \return true if the formula has been splitted and became a 3 QBF
 */
bool Formula::forkExtension()
{
	// d = distinct
	// e = equal
	// dec = for all pairs of existential variables y1, y2
	// it holds that D1 = D2, D1 \cup D2 = \emptyset or D1 = UnivVars or D2 = UnivVars
	// will not split, if the formula is in dec


	if (_prefix->getRMB().size() > 0) {
		return false;
	}

	if (!inDEPlus()) {
		return false;
	}


	// if we get here formula is in de or de+
	// -> perform forksplitting until all clauses are splitted properly
	
	// store original existential vars to put them into the rmb in the end
	const std::vector<Variable> origExistVars = _prefix->getExistVars();

	// for each clause in the formula, we get a vector of new clauses
	std::vector<std::vector<std::vector<Literal>>> newClauses(_clauses.size(),std::vector<std::vector<Literal>>());
	std::vector<Literal> tmpVec;
	std::vector<ClauseID> clausesToDelete;
	size_t numberOfNewVars = 0;

	for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) { // split each clause
		if (_clauses[c_nr].getStatus() == ClauseStatus::DELETED) continue;

		Clause& currentClause = _clauses[c_nr];


		for (size_t pos = 0; pos < currentClause.size(); ++pos) {
			Literal currentLit = currentClause[pos];
			bool inserted = false;

			// if there is already a clause where the literal fits in, we add the literal to it
			// otherwise we start a new clause
			for (auto& splitClause : newClauses[c_nr]) {
				if (belongIntoOneClause[lit2var(currentLit)][lit2var(splitClause[0])] || 
					belongIntoOneClause[lit2var(splitClause[0])][lit2var(currentLit)]) {
					splitClause.push_back(currentLit);
					inserted = true;
					break;
				}
			}

			if (!inserted) {
				tmpVec.push_back(currentLit);
				newClauses[c_nr].push_back(tmpVec);
				tmpVec.clear();
			}
		}

		// check if the clause has been splitted and if so, mark for deletion
		if (newClauses[c_nr].size() > 1) {
			clausesToDelete.insert(clausesToDelete.end(),c_nr);
		}
		else if (newClauses[c_nr].size() == 1) {
			newClauses[c_nr].clear();
		}
		else {
			std::cout << "   didn't get new clause, something wrong" << std::endl;
			exit(0);
		}
	} // end iterating all clauses


	belongIntoOneClause.clear();
	belongIntoOneClause.shrink_to_fit();


	// delete all splitted clauses
	for (const auto c_nr : clausesToDelete) {
		removeClause(c_nr);
	}


	Variable newVar;

	std::set<Literal> emptyDepSet;

	// add new vars to new Clauses
	for (auto& clauses : newClauses) {
		if (clauses.size() == 0) continue;
		if (clauses.size() == 1) {
			std::cout << "unit left, something wrong" << std::endl;
			exit(0);
		}

		for (size_t c1 = 0; c1 < clauses.size()-1; ++c1) {
			for (size_t c2 = c1+1; c2 < clauses.size(); ++c2) {
				newVar = addEVar(emptyDepSet);
				clauses[c1].insert(clauses[c1].end(),var2lit(newVar,false));
				clauses[c2].insert(clauses[c2].end(),var2lit(newVar,true));
			}
		}

		for (const auto& clause : clauses) addClause(clause);
	}


	// moving all original ex vars into rmb
	for (const auto& ex : origExistVars) {
		_prefix->moveToRMB(ex);
	}

	return true;
	
} // end forkExtension



/** \brief d = distinct
 *         e = equal
 *  	   de+ = for all pairs of existential variables that occur together in a clause it holds that 
 *  	   their dependency sets are either equal or distinct
 *  	   determine whether the formula is in de+
 *  \return true if so
 * requires that the formula is universally reduced and that the clauses are sorted
 */
bool Formula::inDEPlus()
{
	belongIntoOneClause = std::vector<std::vector<bool>>(maxVarIndex()+1,std::vector<bool>(maxVarIndex()+1,false));
	
	// check for each pair of variables
	for (Variable var1 = 1; var1 <= maxVarIndex(); ++var1) {
		if (varDeleted(var1)) continue;
		
		const auto depSet1 = getDependencies(var1); 

		for (Variable var2 = var1+1; var2 <= maxVarIndex(); ++var2) {
			if (varDeleted(var2)) continue;

			if (isExistential(var1)) { // var1 is existential
				if (isExistential(var2)) { // var2 is existential

					const auto depSet2 = getDependencies(var2);
					
					if (depSet1.size() == 0) {
						if (depSet2.size() == 0) {
							belongIntoOneClause[var1][var2] = true;
						}
						continue;
					}
					else if (depSet2.size() == 0) {
						continue;
					}
					else {
 						if (depSet1 == depSet2) {
							belongIntoOneClause[var1][var2] = true;
						}
						// if the dependency sets are not the same or distinct, it is still ok
						// if the variables do not occur in one clause
						else {
							std::set<Variable> intersection;
							std::set_intersection(depSet1.begin(), depSet1.end(), depSet2.begin(), depSet2.end(), std::inserter(intersection, intersection.begin()));

							if (intersection.size() == 0) continue;

							if (!inOneClause(var1,var2)) continue;

							belongIntoOneClause.clear();
							return false;
						}
					}
				}
				else { // var2 is universal
					if (depSet1.find(var2) != depSet1.end()) {
						belongIntoOneClause[var1][var2] = true;
					}
				
				}
			} // end if var1 is existential
			else { // if var1 is universal
				if (isExistential(var2)) { // if var2 is existential
					if (depSet1.find(var2) != depSet1.end()) {
						belongIntoOneClause[var1][var2] = true;
					}
					else {
						continue;
					}
				}
				else { // if var2 is universal
					const auto depSet2 = getDependencies(var2);

					if (depSet1 == depSet2) {
						belongIntoOneClause[var1][var2] = true;
						continue;
					}
				}
			} // end if var1 is universal
		}
	}

	return true;

} // end inDEPlus



/**
 * \return true if there is a clause where both variables occur
 */
bool 
Formula::inOneClause(Variable var1, Variable var2) {
	std::set<ClauseID> occ1;

	// collect clauses for both polarities
	for (const auto c : _occ_list[var2lit(var1,false)]) {
		if (_clauses[c].getStatus() == ClauseStatus::DELETED) continue;

		occ1.insert(c);
	}
	for (const auto c : _occ_list[var2lit(var1,true)]) {
		if (_clauses[c].getStatus() == ClauseStatus::DELETED) continue;

		occ1.insert(c);
	}

	if (occ1.size() == 0) return false;

	std::set<ClauseID> occ2;

	for (const auto c : _occ_list[var2lit(var2,false)]) {
		if (_clauses[c].getStatus() == ClauseStatus::DELETED) continue;

		occ2.insert(c);
	}
	for (const auto c : _occ_list[var2lit(var2,true)]) {
		if (_clauses[c].getStatus() == ClauseStatus::DELETED) continue;

		occ2.insert(c);
	}

	if (occ2.size() == 0) return false;


	// compare occurrence lists
	std::set<Variable> intersection;
	std::set_intersection(occ1.begin(), occ1.end(), occ2.begin(), occ2.end(), std::inserter(intersection, intersection.begin()));

	if (intersection.size() > 0) return true;

	return false;
} // end inOneClause





/** 
 * \brief create copy of formula and try to resolve away all literals in rmb. 
 * If this succeeded then take new formula, else stick with old formula
 */
bool Formula::resolveAfterForkSplit()
{
	if (numClauses() > 80000) return false;

	Formula formula(*this);

	bool emptyRmb = formula.resolveRmb();

	if (emptyRmb) {
		*this = formula;
		return true;
	}

	return false;
} // end resolveAfterForkSplit




/** 
 * \return true if rmb is empty in the end 
 * can only be used after forksplitting, because it uses another limit than applyresolution
 */
bool Formula::resolveRmb()
{
	max_resolveRmb_cost = static_cast<int>(numClauses()) * 5;

	const std::set<Variable>& rmbSet = _prefix->getRMB();

	std::vector<Variable> rmb;
	for (const auto r : rmbSet) {
		rmb.push_back(r);
	}

	// eliminate variables
	for (const auto var : rmb) {
		if (varDeleted(var)) continue;

		if (!isResolvable(var)) return false;

		const Literal lit_pos = var2lit(var, false);
		const Literal lit_neg = var2lit(var, true);
	
		// Subtract the cost of those clauses that will be deleted
		max_resolveRmb_cost += _occ_list[lit_pos].size();
		max_resolveRmb_cost += _occ_list[lit_neg].size();

		std::vector<Clause> to_add;
		to_add.reserve(_occ_list[lit_pos].size() * _occ_list[lit_neg].size());

		for (const ClauseID clause_pos: _occ_list[lit_pos]) {
			// Ignore the optional clauses
			if (clauseDeleted(clause_pos)) continue;
			if (clauseOptional(clause_pos)) continue;

			for (const ClauseID clause_neg: _occ_list[lit_neg]) {
				// Ignore the optional clauses
				//if (clauseDeleted(clause_pos)) break;
				if (clauseDeleted(clause_neg)) continue;
				if (clauseOptional(clause_neg)) continue;

				auto resolvent = resolve(_clauses[clause_pos], _clauses[clause_neg], var);

				if (!resolvent.isTautology()) {
					universalReduction(resolvent, -1);
					// add cost of added resolvent
					max_resolveRmb_cost -= 1;
					if (max_resolveRmb_cost <= 0 || numClauses() > 80000) {
						return false;
					}
					to_add.push_back(std::move(resolvent));
				}
			}
		}

	
		// Delete all clauses (incl. the optional ones) that contain the eliminated variable
		while (!_occ_list[lit_pos].empty()) {
			const Clause& clause = _clauses[_occ_list[lit_pos].front()];
			removeClause(_occ_list[lit_pos].front());
		}

		while (!_occ_list[lit_neg].empty()) {
			const Clause& clause = _clauses[_occ_list[lit_neg].front()];
			removeClause(_occ_list[lit_neg].front());
		}

		removeVar(var);

		// add the new clauses
		for (Clause& resolvent: to_add) {
			addClause(std::move(resolvent));
		}
	
		// Propagate possible new units
		if (!_unit_stack.empty()) unitPropagation();
	}

	return true;
} // end resolvermb



} // end namespace
