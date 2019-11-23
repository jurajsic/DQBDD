#include <string>
#include <sstream>
#include <iostream>
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

void TreeSolver::setTest2Formula() {
    delete root;
    QuantifierTreeFormula *f1 = new QuantifierTreeFormula(mgr, qvMgr);
    f1->setMatrix(Variable(0,mgr));
    QuantifierTreeFormula *f2 = new QuantifierTreeFormula(mgr, qvMgr);
    f2->setMatrix(!Variable(1,mgr));
    QuantifierTreeFormula *f3 = new QuantifierTreeFormula(mgr, qvMgr);
    f3->setMatrix(!Variable(2,mgr));
    QuantifierTreeFormula *f4 = new QuantifierTreeFormula(mgr, qvMgr);
    f4->setMatrix(Variable(3,mgr));

    auto qf1 = new QuantifierTree(false, std::list<QuantifierTreeNode*>{ f1, f2 }, qvMgr);
    auto qf2 = new QuantifierTree(false, std::list<QuantifierTreeNode*>{ f3, f4 }, qvMgr);
    root = new QuantifierTree(true, std::list<QuantifierTreeNode*>{ qf1, qf2 }, qvMgr);
    root->addUnivVar(Variable(1,mgr));
    root->addUnivVar(Variable(3,mgr));
    root->addExistVar(Variable(0,mgr), VariableSet{ Variable(1,mgr) });
    root->addExistVar(Variable(2,mgr), VariableSet{ Variable(3,mgr) });
}


void TreeSolver::setTest3Formula() {
    delete root;
    QuantifierTreeFormula *f1 = new QuantifierTreeFormula(mgr, qvMgr);
    f1->setMatrix(Variable(0,mgr));
    QuantifierTreeFormula *f2 = new QuantifierTreeFormula(mgr, qvMgr);
    f2->setMatrix(!Variable(1,mgr));

    root = new QuantifierTree(false, std::list<QuantifierTreeNode*>{ f1, f2 }, qvMgr);
    root->addUnivVar(Variable(0,mgr));
    root->addUnivVar(Variable(1,mgr));
}

void TreeSolver::runTests() {
    int numOfCorrect = 0;
    int numOfAll = 0;
    ++numOfAll;
    std::cout << "Test number :" << numOfAll << std::endl;
    setTest1Formula();
    if (solve() == true) {
        ++numOfCorrect;
    }
    ++numOfAll;
    std::cout << "Test number :" << numOfAll << std::endl;
    setTest3Formula();
    if (solve() == false) {
        ++numOfCorrect;
    }
    ++numOfAll;
    std::cout << "Test number :" << numOfAll << std::endl;
    setTest2Formula();
    if (solve() == true) {
        ++numOfCorrect;
    }
    std::cout << numOfCorrect << "/" << numOfAll << std::endl;
}

bool TreeSolver::solve() {
    std::cout << *root << std::endl;
    root->localise();
    QuantifierTreeFormula *f = root->getFormula(mgr);
    root = f;
    std::cout << *f << std::endl;
    return (f->getMatrix().IsOne());
}