# DQBF-BDD

DQBF (Dependency quantifier boolean formulas) solver using BDDs (binary decision diagrams). This implementation uses BuDDy BDD library (https://github.com/jgcoded/BuDDy). Put Buddy in the same folder where you put folder with this implementation (that is BuDDy/ and DQBF-BDD/ are in the same folder). Compile BuDDy using the manual on the github, then go to DQBF-BDD/ folder and run make. This should result in program solve which can be called as

./solve "name of file in .dqdimacs format"