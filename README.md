# DQBF-BDD

DQBF (Dependency quantifier boolean formulas) solver using BDDs (binary decision diagrams). This implementation uses cudd BDD library (https://github.com/ivmai/cudd). Just run make to compile.

## Usage

```
./solve 0
```
runs test examples using simple solver

```
./solve 1
```
runs test examples using solver that pushes quantifiers inside formula

```
./solve 0 "name of file in .dqdimacs format"
```
reads DQBF from the file and checks its satisfiability using simple solver

```
./solve 1 "name of file in .dqdimacs format"
```
reads DQBF from the file and checks its satisfiability using solver that pushes quantifiers inside formula