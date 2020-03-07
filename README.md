# DQBF-BDD

Copyright 2020 Juraj Síč  
This program is released under the version 3 of the
GNU Lesser General Public License  
(LGPL v3, see https://www.gnu.org/licenses/lgpl-3.0.en.html)  



TODO: rewrite this, this is old info

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