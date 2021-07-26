# DQBDD

DQBDD is a dependency quantified Boolean formula (DQBF) solver that uses binary decision diagrams (BDDs) as an underlying representation of formulas. It is written in C++ and it reads DQBFs encoded in DQDIMACS or prenex cleansed DQCIR format for which it checks their satisfiability using quantifier elimination. For an explanation of the techniques used in DQBDD, see [this paper](https://doi.org/10.1007/978-3-030-80223-3_36) and my [master's thesis](https://is.muni.cz/th/prexv/).

## Installation

You can find binaries in [tagged release versions](https://github.com/jurajsic/DQBDD/releases). If you want to compile it yourself, you need C++ compiler supporting C++14 standard (GCC 7.4.0 on Ubuntu 18.04.1 and Clang 11.0.3 on macOS 11.1 are known to work) and [CMake](https://cmake.org/) (3.5 or later). Execute 
```
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
to build DQBDD which will be located in `Release/src/`. However, using `master` branch is not recommended as it is generally a work in progress.

## Dependencies
There is no need to install any dependency, all of them are in `libs/` and are compiled with DQBDD. They are these:
- [antom](https://projects.informatik.uni-freiburg.de/projects/antom) - SAT solver used in HQSpre
- [CUDD v3.0.0](https://github.com/ivmai/cudd) - BDD library
- [cxxopts v2.2.0](https://github.com/jarro2783/cxxopts) - argument parser
- [Easylogging++ v9.96.7](https://github.com/zuhd-org/easyloggingpp) - C++ logger used in HQSpre
- [HQSpre](https://abs.informatik.uni-freiburg.de/src/projects_view.php?projectID=21) - DQBF preprocessor
- [PicoSAT](http://fmv.jku.at/picosat/) - SAT solver used in HQSpre

## Input Formats

DQBDD accepts two input formats: [DQDIMACS](https://doi.org/10.29007/1s5k) and prenex cleansed DQCIR.

### DQDIMACS

DQDIMACS is used to represent DQBFs in prenex CNF form. It is an extension of [QDIMACS](http://www.qbflib.org/qdimacs.html) format where not only `a` and `e` quantifiers can be used (representing universal and existential variables respectively), but also `d` quantifier is introduced (representing existential variable with its dependencies). [Appendix C](https://is.muni.cz/th/prexv/thesis.pdf#appendix.C) of my [thesis](https://is.muni.cz/th/prexv/) contains a full definition of DQDIMACS format.

Example of DQDIMACS format:
```
p cnf 6 2
a 1 2 0
e 3 0
a 4 0
e 5 0
d 6 1 4 0
-1 6 -3 5 0
-2 -4 5 0
```
This represents DQBF `∀x1 ∀x2 ∀x4 ∃x3(x1, x2) ∃x5(x1, x2, x4) ∃x6(x1, x4).(¬x1 ∨ x6 ∨ ¬x3 ∨ x5) ∧ (¬x2 ∨ ¬x4 ∨ x5)`

### Prenex Cleansed DQCIR

Prenex cleansed DQCIR is an extension of [prenex cleansed QCIR](http://www.qbflib.org/qcir.pdf) format. This format is used for DQBFs in prenex form, however, it does not need to be in CNF. It adds quantifiers of type `depend(v, v1, ..., vn)` which represents existential variable `v` with the dependency set `v1, ..., vn`. It is assumed that `v1, ..., vn` were already defined as `forall` variables and that `v` was not yet defined. The quantifier `exists` can still be used and it represents existential variables with the dependency set equal to the set of universal variables which were already defined with `forall` quantifier.

Example of prenex cleansed DQCIR format:
```
#QCIR-G14 10
forall(1, 2)
exists(3)
forall(4)
exists(5)
depend(6, 1, 4)
output(7)
8 = and(-6, 3)
9 = or(-1, -8, 5)
10 = or(-2, -4, 5)
7 = and(9, 10)
```
This represents DQBF `∀x1 ∀x2 ∀x4 ∃x3(x1, x2) ∃x5(x1, x2, x4) ∃x6(x1, x4).(¬x1 ∨ ¬(¬x6 ∧ x3) ∨ x5) ∧ (¬x2 ∨ ¬x4 ∨ x5)`

## Usage

    DQBDD [OPTION...] <input file>

`<input file>` should be a path to formula in DQDIMACS or DQCIR format. For the list of options run `DQBDD --help`. If the formula is satisfiable, the return value is 10, otherwise it is 20.

## Examples

```
DQBDD file.dqdimacs
```
solves the formula in `file.dqdimacs` with the default settings.

```
DQBDD --preprocess 0 --dyn-reordering 0 file.dqdimacs
```
solves the formula in `file.dqdimacs` without running the preprocessor HQSpre first (the default for DQDIMACS files is `preprocess 1`, i.e. preprocessing is enabled) and without using dynamic reordering of variables in BDDs as implemented in CUDD (dynamic reordering is enabled by default, DQBDD uses [Rudell's sifting algorithm](https://ieeexplore.ieee.org/document/580029)).

```
DQBDD --localise 0 --uvar-choice 1 file.dqdimacs
```
solves the formula in `file.dqdimacs` without localising quantifiers (or creating quantifier tree, the default behaviour is `--localise 1`) where the next universal variable for universal expansion is always the one that has the minimal number of dependent existential variables (`--uvar-choice 1`). The other options of `--uvar-choice` are:
- 0 - the order of universal variables to expand is set at beginning from the smallest to the largest number of dependencies (this is the default),
- 2 - the next variable is chosen by the number of variables in BDDs representing the two conjucts of universal expansion.

```
DQBDD --localise 1 --elimination-choice 2 file.dqdimacs
```
solves the formula in `file.dqdimacs` with localising quantifiers (this is the default behaviour) where it eliminates all universal and possible existential variables while creating the final BDD from the quantifier tree (`--elimination-choice 2`). The other options of `--elimination-choice` are:
- 0 - does not eliminate any quantifiers,
- 1 - eliminates only universal variables which do not have any dependencies and all possible existential variables (this is the default).

```
DQBDD file.dqcir
```
solves the formula in DQCIR format with the default settings. The default settings are the same as for DQDIMACS format except for preprocessing. As HQSpre can only accept files in DQDIMACS format, there is no preprocessing for DQCIR files.

```
DQBDD --hqspre-dqcir-output file.dqdimacs
```
does not solve the formula, but runs HQSpre preprocessor and saves the result of HQSpre with gate extraction to `file.dqcir`. This is basically the same output that DQBDD takes when it processes DQDIMACS file with preprocessing.

## Licence

- **[LGPL v3](https://www.gnu.org/licenses/lgpl-3.0.en.html)**
- Copyright 2020, 2021 Juraj Síč
