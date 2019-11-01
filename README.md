# DQBF-BDD

DQBF (Dependency quantifier boolean formulas) solver using BDDs (binary decision diagrams). This implementation uses cudd BDD library (https://github.com/ivmai/cudd). Put cudd in the same folder where you put folder with this implementation (that is cudd-release/ and DQBF-BDD/ are in the same folder). Compile cudd using the manual on the github, then go to DQBF-BDD/ folder and run make. This should result in program solve which can be called as

./solve
runs test example, should result in UNSAT

./solve "name of file in .dqdimacs format"
reads DQBF from the file and checks its satisfiability