# DQBDD

DQBDD is a dependency quantified Boolean formula (DQBF) solver that uses binary decision diagrams (BDDs) as an underlying representation of formulas. It is written in C++ and it reads DQBFs encoded in [DQDIMACS](https://doi.org/10.29007/1s5k) format for which it checks their satisfiability using quantifier elimination. For an explanation of the techniques used in DQBDD, see my [master's thesis](https://is.muni.cz/th/prexv/).

## Installation

You can find binaries in [tagged release versions](https://github.com/jurajsic/DQBDD/releases). If you want to compile it yourself, you need C++ compiler supporting C++14 standard and [CMake](https://cmake.org/). Execute 
```
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
to build DQBDD which will be located in `Release/src/`. However, do not use `master` branch as it is generally a work in progress.

## Usage

    DQBDD [OPTION...] <input file>

`<input file>` should be formula to solve in DQDIMACS format, for the list of options run `DQBDD --help`. If the formula is satisfiable, the return value is 10, otherwise it is 20.

## Examples

```
DQBDD file.dqdimacs
```
solves the formula in `file.dqdimacs` with the default settings.

```
DQBDD --preprocess 0 --dyn-reordering 0 file.dqdimacs
```
solves the formula in `file.dqdimacs` without running the preprocessor HQSpre first and without using dynamic reordering of variables in BDDs as implemented in CUDD (the default behaviour is to use both preprocessing and dynamic reordering).

```
DQBDD --localise 0 --uvar-choice 1 file.dqdimacs
```
solves the formula in `file.dqdimacs` without localising quantifiers (or creating quantifier tree) where the next universal variable for universal expansion is always the one that has the minimal number of dependent existential variables. The other options of `--uvar-choice` are:
- 0 - the order of universal variables to expand is set at beginning from the smallest to the largest number of dependencies (this is the default),
- 2 - the next variable is chosen by the number of variables in BDDs representing the two conjucts of universal expansion.

```
DQBDD --localise 1 --elimination-choice 2 file.dqdimacs
```
solves the formula in `file.dqdimacs` with localising quantifiers (this is the default behaviour) where it eliminates all universal and possible existential variables while creating the final BDD from the quantifier tree. The other options are:
- 0 - does not eliminate any quantifiers,
- 1 - eliminates only universal variables which do not have any dependencies and all possible existential variables (this is the default).

## Dependencies
There is no need to install any dependency, all of them are in `libs/` and are compiled with DQBDD. They are these:
- [antom](https://projects.informatik.uni-freiburg.de/projects/antom) - SAT solver used in HQSpre
- [CUDD v3.0.0](https://github.com/ivmai/cudd) - BDD library
- [cxxopts v2.2.0](https://github.com/jarro2783/cxxopts) - argument parser
- [Easylogging++ v9.96.7](https://github.com/zuhd-org/easyloggingpp) - C++ logger used in HQSpre
- [HQSpre](https://abs.informatik.uni-freiburg.de/src/projects_view.php?projectID=21) - DQBF preprocessor
- [PicoSAT](http://fmv.jku.at/picosat/) - SAT solver used in HQSpre

## Licence

- **[LGPL v3](https://www.gnu.org/licenses/lgpl-3.0.en.html)**
- Copyright 2020 Juraj Síč
