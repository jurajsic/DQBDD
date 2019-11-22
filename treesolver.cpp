#include <string>
#include <sstream>
#include "treesolver.hpp"

TreeSolver::TreeSolver(const Cudd &mgr) : Solver(mgr) {}

TreeSolver::~TreeSolver() {
    delete root;
}

void TreeSolver::readFile(std::ifstream& file) {
    std::string line;

    QuantifierTree *root = new QuantifierTree(true, std::list<QuantifierTreeNode*>{}, qvMgr);

    while(std::getline(file, line)) {
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        if (token == "p") {
            continue;
            // TODO maybe initialize manager here based on the size??
            /*
            streamline >> token; // ignore "cnf"
            streamline >> token; // number of variables
            // TODO decide initial number of variables
            int numOfVariables = std::stoi(token) + 1;
            streamline >> token; // number of CNF conjuncts
            // TODO decide iniatiliazon of BDD based on the size of formula
            bddProcessor.initialize(100000,10000);
            bddProcessor.setNumOfVars(numOfVariables);
            */
        } else if (token == "a") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token), mgr);
                root->addUnivVar(univVar);
            }
        } else if (token == "e") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable existVar(std::stoi(token), mgr);
                root->addExistVar(existVar, root->getUnivVars());
            }    
        } else if (token == "d") {
            streamline >> token;
            Variable existVar(std::stoi(token), mgr);
            root->addExistVar(existVar);
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token), mgr);
                root->addDependency(existVar, univVar);
            }
        } else { // parse clause (disjunction of literals)
            // n -> returns BDD variable with index n
            // -n -> returns negated BDD variable with index n
            auto getVarFromStr = [&](std::string tok) {
                int i = std::stoi(tok);
                QuantifierTreeFormula *qtf = new QuantifierTreeFormula(mgr, qvMgr);
                if (i < 0) {
                    qtf->setMatrix(!Variable(-i, mgr));
                } else {
                    qtf->setMatrix(Variable(i, mgr));
                }
                return qtf;
            };
            QuantifierTree *clause = new QuantifierTree(false, {}, qvMgr);
            clause->addChild(getVarFromStr(token));
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                clause->addChild(getVarFromStr(token));
            }
            root->addChild(clause);
        }
    }
    this->root = root;
}

void TreeSolver::setTest1Formula() {
    delete root;
    QuantifierTreeFormula *f1 = new QuantifierTreeFormula(mgr, qvMgr);
    f1->setMatrix(Variable(0,mgr));
    QuantifierTreeFormula *f2 = new QuantifierTreeFormula(mgr, qvMgr);
    f2->setMatrix(!Variable(1,mgr));

    root = new QuantifierTree(true, std::list<QuantifierTreeNode*>{ f1, f2 }, qvMgr);
    root->addExistVar(Variable(0,mgr));
    root->addExistVar(Variable(1,mgr));
}

bool TreeSolver::solve() {
    root->localise();
    QuantifierTreeFormula *f = root->getFormula(mgr);
    root = f;
    return (f->getMatrix().IsOne());
}