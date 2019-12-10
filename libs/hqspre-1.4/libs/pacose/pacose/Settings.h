#ifndef SETTINGS_H
#define SETTINGS_H

namespace Glucose{

  enum EncodingType {
      warn,               // adder warners 1996, [Warners 1998?]
      bail,               // totalizer [Bailleux & Boufkhad 2003]
      ASIN,               // [Asin et. al 2011] Robert Asin, Robert Nieuwenhuis, Albert Oliveras, Enric Rodriguez-Carbonell "Cardinality Networks: a theoretical and empirical study"
      ogaw,               // modulo Totalizer [Ogawa et. al 2013]
      bailw2,             // BailleuxW2
      wmto,               // Weighted MaxSAT Totalizer
      mrwto,              // Mixed Radix Weighted Totalizer
      mrwto2,             // Mixed Radix Weighted Totalizer
      heuristicQMaxSAT,   // auto heuristic [Koshi, 2014]
                          // selects between "warn, bail, ogaw"
      heuristicDGPW18     // heuristic described in DGPW SAT paper [Paxian & Reimer 2018]
                          // selects between "warn, DGPW"
  };

  enum SolverType {
      SOLVER,
      SIMPSOLVER,
      PARALLELSOLVER,
      MAXSATSOLVER
  };


struct Settings
{
  public:
    Settings():
          _encoding(heuristicQMaxSAT),
          _solverType(SOLVER)
    {}

    EncodingType _encoding;
    SolverType _solverType;

    /**
     * @brief ResetCore
     *          sets all settings to standard values.
     */
    void ResetCore()
    {
        _encoding = heuristicQMaxSAT;
        _solverType = SOLVER;
    }
};

}

#endif // SETTINGS_H
