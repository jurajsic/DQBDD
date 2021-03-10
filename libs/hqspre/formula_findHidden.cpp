#include <queue>
#include <map>
#include <iostream>
#include <fstream>

#include "aux.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"

namespace hqspre
{

/**
 * assume clause set C
 * if we can find a complete undirected bipartite graph G where
 * nodes represent literals and
 * an edge (l1,l2) means that we have the clause {l1,l2} in C
 * assume all literals in G are l1,...,ln and we partition G into 
 * the nodes Vl(eft) (= ll1,...,llj) and Vr(ight) = (= lr1,...,lrk)
 * and C/C' (C' contains all binary clauses which are involved in G) does not contain l1,...,ln
 * then we can pick one representative, assume ll1, and 
 * replace ll2,...,llj by ll1 and lr1,...,lrk by ~ll1.
 */
void 
Formula::findHiddenEquivalences()
{
	// go through all clauses and collect binary clauses, occurrences
	// find node with most adjacent nodes to check for equivalences first
	nodes = std::vector<Node*>((maxVarIndex()+1)*2,NULL);

	// go find all binary clauses and build the graphs
	for (size_t c = 0; c < _clauses.size(); ++c) {
		if (_clauses[c].getStatus() == ClauseStatus::DELETED) continue;
		if (_clauses[c].size() != 2) continue;

		Literal l1 = _clauses[c].front();
		Literal l2 = _clauses[c].back();

		if (nodes[l1] == NULL) {
			Node* newNode = new Node();
			newNode->lit = l1;
			nodes[l1] = newNode;
		}
		if (nodes[l2] == NULL) {
			Node* newNode = new Node();
			newNode->lit = l2;
			nodes[l2] = newNode;
		}

		nodes[l1]->adjacents.push_back(nodes[l2]);
		nodes[l2]->adjacents.push_back(nodes[l1]);
	}

	// sort nodes by number of their adjacents
	std::multimap<size_t, Node*> sortByAdjacents;
	for (size_t i = 2; i < nodes.size(); ++i) {
		if (nodes[i] == NULL) continue;
		
		sortByAdjacents.insert(std::pair<size_t,Node*>(nodes[i]->adjacents.size(),nodes[i])); 
	}

	// find equivalent literals beginning with the node with most adjacents
	std::vector<Node*> equivNodesLeft;
	std::vector<Node*> equivNodesRight;
	for (std::multimap<size_t,Node*>::reverse_iterator n = sortByAdjacents.rbegin(); n != sortByAdjacents.rend(); ++n) {
		if (n->second->seen) continue;

		Node* node = n->second;
		if (checkForEquivalence(node, equivNodesLeft, equivNodesRight)) {
			replaceEquivalentLiterals(equivNodesLeft, equivNodesRight);
		}
		equivNodesLeft.clear();
		equivNodesRight.clear();
	}

	// eliminate all variables that only occur once positive and once negative by resolution
	for (Variable i = 1; i <= maxVarIndex(); ++i) {
		if (!isExistential(i)) continue;
		if (_occ_list[i*2].size() == 1 && _occ_list[i*2+1].size() == 1) {
			elimEVarLimit(i, std::numeric_limits<long int>::max(), NULL);
		}
	}

	unitPropagation();

} // end findHiddenEquivalences



/**
 * \brief for a node, follow its adjacents in the graph to find out whether it is bipartite and complete
 * if so, fill equivNodesLeft and Right with the nodes such that we can replace the literals
 */
bool 
Formula::checkForEquivalence(Node* node, std::vector<Node*>& equivNodesLeft, std::vector<Node*>& equivNodesRight)
{
	val_assert(node->side == Node::Side::unknown);
	
	node->side = Node::Side::left;
		
	//bool bipartite = true;
	unsigned int leftSideNodes = 1;
	unsigned int rightSideNodes = 0;
	unsigned int edges = 0;
	std::vector<Variable> seenVars(maxVarIndex()+1,false); // check if every variable only occur once in component
	std::vector<Variable> varsTwiceInCC;				   // because we might have v and -v in the graph

	equivNodesLeft.push_back(node);

	std::stack<Node*> pending;
	pending.push(node);

	// for a node turn seen to true at the first time it is taking out from pending
	while (!pending.empty()) {
		if (pending.top()->seen) {
			pending.pop();
			continue;
		}

		Node* currentNode = pending.top();
		pending.pop();

		currentNode->seen = true;

		// store var if occur for the second time
		if (seenVars[lit2var(currentNode->lit)])
			varsTwiceInCC.push_back(lit2var(currentNode->lit));
		else
			seenVars[lit2var(currentNode->lit)] = true;

		// go through all adjacent nodes
		for (auto& a : currentNode->adjacents) {
			if (a->seen) continue;
			
			pending.push(a);

			if (a->side == currentNode->side) {
				return false;
			}

			++edges;

			if (a->side != Node::Side:: unknown)
				continue;

			if (currentNode->side == Node::Side::left) {
				a->side = Node::Side::right;

				equivNodesRight.push_back(a);
				++rightSideNodes;
			}
			else if (currentNode->side == Node::Side::right) {
				a->side = Node::Side::left;

				equivNodesLeft.push_back(a);
				++leftSideNodes;
			}

		} // end for all adjacents 

		// check for purity
		if (_occ_list[currentNode->lit].size() != currentNode->adjacents.size()) {
			return false;
		}

	} // end while

	// check for completeness
	if (edges != (rightSideNodes * leftSideNodes)) {
		return false;
	}


	if (varsTwiceInCC.size() == 0)
		return true;
	else {
		// we assume that the literals of one variable are on one side,
		// because we do not have tautological clauses.


		// TODO: Do not stop here, but find test cases
		std::cout << "found var in different polarities in one graph -> untested case" << std::endl;
		return false;

		// TODO: find test cases
		Node::Side firstSide = nodes[var2lit(varsTwiceInCC[0],false)]->side;

		// if we can find two pairs of literals on opposite sides, then the formula is unsatisfiable.
		for (const auto var : varsTwiceInCC) {
			if (nodes[var2lit(var,false)]->side != firstSide) {
				for (auto n : nodes) delete n;
				throw UNSATException("found unsatisfiable bipartite graph");
			}
		}

		// if all literals are on the same side, then we can assign the literals 
		// on the other side of the graph as units.

		if (firstSide == Node::Side::left) {
			for (const auto r : equivNodesRight) pushUnit(r->lit, PureStatus::UNIT);
			for (const auto l : equivNodesLeft) pushUnit(negate(l->lit), PureStatus::UNIT);

		}
		else {
			for (const auto l : equivNodesLeft) pushUnit(l->lit, PureStatus::UNIT);
			for (const auto r : equivNodesRight) pushUnit(negate(r->lit), PureStatus::UNIT);
		}		

		return false; // because these literals are declared unit already

		// TODO: find cases where we can add tautological clauses to complete the graph
	}
} // end checkForEquivalence




void 
Formula::replaceEquivalentLiterals(std::vector<Node*>& equivNodesLeft, std::vector<Node*>& equivNodesRight)
{
	const Literal representative = equivNodesLeft[0]->lit;

	// delete all clauses involved in CC
	for (const auto l : equivNodesLeft) {
		for (auto o : _occ_list[l->lit]) {
			if (_clauses[o].getStatus() == ClauseStatus::DELETED) continue;

			removeClause(o);
		}

		if (l->lit == representative) continue;
		replaceLiteral(l->lit, representative);
	}

	for (const auto r : equivNodesRight) {
		for (auto o : _occ_list[r->lit]) {
			if (_clauses[o].getStatus() == ClauseStatus::DELETED) continue;

			removeClause(o);
		}
		replaceLiteral(negate(r->lit), representative);
	}

	// replace in rest all literals by one representative
	
} // end replaceEquivalentNodes




/** 
 * \brief calls MPhaseSAT64 SAT solver to solve forksplitted instances
 */
void
Formula::solveSAT()
{
	std::cout << "calling external SAT solver MPhaseSAT64-static" << std::endl;
	val_assert(numUVars() == 0);

	std::string dimacsname = std::string(std::tmpnam(NULL))+".dimacs";
	std::ofstream out(dimacsname);
	std::cout << "file = " << dimacsname << std::endl;

	updateVars();
	val_assert(checkConsistency());

	std::vector<Variable> translation_table(maxVarIndex()+1, 0);
	Variable current = 0;
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var)) continue;
        translation_table[var] = var;
    }

	out << "p cnf " << maxVarIndex() << " " << numClauses() << std::endl;
	writeClauses(out);
	out.close();

	const std::string solverPath = "./MPhaseSAT64-static " + dimacsname;

	const int result = std::system(solverPath.c_str());

	std::remove(dimacsname.c_str());

	if ( result == 2560 ) throw SATException("Solved by MPhase");
	else if ( result == 5120 ) throw UNSATException("Solved by Mphase");
}	

Formula::Node::Node()
	:
	side(unknown),
	seen(false)
{}


Formula::Node::~Node()
{}


}

