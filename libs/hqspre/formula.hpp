/*
 * This file is part of HQSpre.
 *
 * Copyright 2016/17 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
 * Albert-Ludwigs-Universitaet Freiburg, Freiburg im Breisgau, Germany
 *
 * HQSpre is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HQSPRE_FORMULA_HPP_
#define HQSPRE_FORMULA_HPP_

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <set>
#include <stack>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "aux.hpp"
#include "clause.hpp"
#include "exceptions.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "process_limits.hpp"
#include "settings.hpp"
#include "timer.hpp"
#include "varheap.hpp"

#ifdef SKOLEM
#    include "skolem.hpp"
#endif

/**
 * \file formula.hpp
 * \brief Main header file for DQBF formulas.
 * \author Ralf Wimmer
 * \date 2015-2016
 * \details Contains all functions for manipulating literals and DQBFs.
 */

namespace hqspre {

/**
 * \brief Type of the dependency scheme that is to be applied.
 */
enum class DependencyScheme : int
{
    TRIVIAL = 0,         ///< Trivial dependency scheme as given in the prefix
    STANDARD,            ///< Standard dependency scheme
    STRICT_STANDARD,     ///< Strict standard dependency scheme
    REF_TRIANGLE,        ///< Reflexive triangle dependency scheme
    REF_QUADRANGLE,      ///< Reflexive quadrangle dependency scheme
    RP_STANDARD,         ///< Standard resolution path dependency scheme
    RP_STRICT_STANDARD,  ///< Strict standard resolution path dependency scheme
    RP_REF_TRIANGLE,     ///< Reflexive triangle dependency scheme
    RP_REF_QUADRANGLE,   ///< Reflexive quadrangle dependency scheme
    GATE                 ///< Gate detection
};

/**
 * \brief How to modify the dependencies
 */
enum class DependencyOperation : int
{
    ADD,        ///< Add dependencies using the selected dependency scheme
    REMOVE,     ///< Remove dependencies using the selected dependency scheme
    DO_NOTHING  ///< Do not modify the dependencies.
};

/// Timers for different preprocessing techniques
enum class WhichTimer : unsigned int
{
    UNIT_SAT,                 ///< Finding unit literals (a.k.a. backbones) by SAT calls
    IMPLICATIONS_SAT,         ///< Finding implications (binary clauses) by SAT calls
    INCOMPLETE_SAT,           ///< Incomplete (un)satisfiability checks using a SAT solver
    IMPLICATION_CHAINS,       ///< Implications chains (special case of resolution)
    DEPENDENCY_SCHEMES,       ///< Application of dependency schemes
    DEP_ELIM_SET,             ///< Time for computing a dependency elimination set
    EQUIV_LITS,               ///< Finding equivalent literals using the binary implication
                              ///< graph
    EQUIV_GATES,              ///< Finding equivalent gate definitions (a.k.a. structural
                              ///< hashing)
    CONTRADICTIONS,           ///< Finding contradicting implications (a =>* ~a means that
                              ///< ~a is unit)
    HIDDEN_EQUIV_CONTRA,      ///< Finding hidden equivalences and contradictions
    UNIV_REDUCTION,           ///< Universal reduction
    UNIV_EXPANSION,           ///< Universal expansion
    BCE,                      ///< Blocked clause elimination
    HSE,                      ///< Hidden subsumption elimination
    GATE_DETECTION,           ///< Gate detection (accumulated time for the different types
                              ///< of gates)
    GATE_MUX_DETECTION,       ///< Detection of MUX gates
    GATE_AND_DETECTION,       ///< Detection of AND gates
    GATE_XOR_DETECTION,       ///< Detection of binary XOR gates
    GATE_SEMANTIC_DETECTION,  ///< Semantic gate detection using a SAT solver
    SUBSTITUTION,             ///< Gate substitution
    REWRITING,                ///< Gate rewriting
    RESOLUTION,               ///< Variable elimination by resolution
    SUBSUMPTION,              ///< Subsumption checks
    SELF_SUBSUMPTION,         ///< Self-subsuming resolution
    UPLA,                     ///< Unit propagation lookahead
    VIVIFICATION,             ///< Vivification
    TOTAL_TIME,               ///< Total preprocessing time
    _last_                    ///< Dummy element. Must always be the last one!
};

enum class Statistics : unsigned int
{
    UNIT,
    PURE,
    UNIV_REDUCTION,
    UNIV_EXPANSION,
    BCE,
    BIA,
    HTE,
    HSE,
    BLE,
    BLA,
    IMPLICATION_CHAINS,
    CONTRADICTIONS,
    GATES_AND,
    GATES_XOR,
    GATES_MUX,
    GATES_SEMANTIC,
    SUBSUMPTION,
    RESOLUTION,
    SELF_SUBSUMPTION,
    EQUIV_LITS,
    EQUIV_GATES,
    HIDDEN_EQUIV_LITS,
    HIDDEN_UNIT,
    SUBSTITUTION,
    REWRITING,
    UNIT_SAT,
    IMPL_SAT,
    SAT_CALLS,
    VIVIFIED_LITERALS,  ///< Number of literals removed by vivification
    VIVIFIED_CLAUSES,   ///< Number of clauses deleted by vivification
    ADD_DEPENDENCY_SCHEMES,
    REM_DEPENDENCY_SCHEMES,
    PREPRO_LOOPS,
    _last_
};

// Forward declaration
class SatPropagator;

/**
 * \brief Represents a DQBF with a number of operations on it.
 */
class Formula
{
   public:
    using ImplicationSet = std::unordered_set<BinaryClause>;

    //@{
    /**
     * \name Construction, destruction, and assigment
     */

    /**
     * \brief Constructs an empty formula of a given type
     */
    Formula();

    /**
     * \brief Creates a copy of an existing formula.
     */
    Formula(const Formula& other);

    /**
     * \brief Moves an existing formula into a new one.
     */
    Formula(Formula&& other) noexcept;

    /**
     * \brief Frees the memory occupied by a formula.
     */
    ~Formula() noexcept;

    /**
     * \brief Assigns a new formula to the current object.
     */
    Formula& operator=(Formula& other); // = delete;

    /**
     * \brief Assigns a new formula to the current object (rvalue version)
     */
    Formula& operator=(Formula&& other) = delete;

    //@}

    //@{
    /**\name Formula input and output
     */
    void read(std::istream& stream);
    void readSAT(std::istream& stream);
    void write(std::ostream& stream, bool compact = false) const;
    //@}

    //@{
    /**\name Variable handling
     */
    Variable addUVar();
    Variable addEVar();
    template <typename Container>
    Variable addEVar(const Container& deps);
    Variable addEVar(std::set<Variable>&& deps);

    std::size_t               numVars() const noexcept;
    std::size_t               numUVars() const noexcept;
    std::size_t               numEVars() const noexcept;
    std::size_t               numLiterals() const noexcept;
    bool                      isUniversal(Variable var) const noexcept;
    bool                      isExistential(Variable var) const noexcept;
    void                      setMaxVarIndex(Variable index);
    constexpr static Variable minVarIndex() noexcept;
    Variable                  maxVarIndex() const noexcept;
    constexpr static Literal  minLitIndex() noexcept;
    Literal                   maxLitIndex() const noexcept;
    void                      removeVar(Variable var);
    void                      updateVars();
    bool                      varDeleted(Variable var) const noexcept;
    void                      setDontTouch(Variable var, bool status) noexcept;
    bool                      isDontTouch(Variable var) const noexcept;
    bool                      isGateOutput(Variable var) const noexcept;

    std::size_t getLevel(Variable var) const noexcept;
    //@}

    void reset();

    //@{
    /**\name Clause handling
     */
    std::size_t numClauses() const noexcept;
    ClauseID    maxClauseIndex() const noexcept;

    template <typename Container>
    int addClause(const Container& clause, bool needs_sorting = true, bool check_subsumption = true,
                  ClauseStatus status = ClauseStatus::MANDATORY);

    int addClause(Clause::ClauseData&& clause, bool needs_sorting = true, bool check_subsumption = true,
                  ClauseStatus status = ClauseStatus::MANDATORY);
    int addClause(const Clause& clause);
    int addClause(Clause&& clause);
    template <typename Container>
    int           findClause(const Container& clause);
    const Clause& getClause(ClauseID c_nr) const noexcept;
    void          removeClause(ClauseID c_nr);
    bool          removeOptionalClauses();
    bool          clauseDeleted(ClauseID c_nr) const noexcept;
    bool          clauseOptional(ClauseID c_nr) const noexcept;
    bool          clauseMandatory(ClauseID c_nr) const noexcept;
    bool          isGateClause(ClauseID c_nr) const noexcept;
    void          writeClauses(std::ostream& stream, std::vector<Variable>* translation_table = nullptr) const;
    //@}

    //@{
    /** \name Gate handling
     */
    Literal addAndGate(Literal input1, Literal input2);
    Literal addNandGate(Literal input1, Literal input2);
    Literal addOrGate(Literal input1, Literal input2);
    Literal addNorGate(Literal input1, Literal input2);
    Literal addXorGate(Literal input1, Literal input2);
    Literal addMux(Literal x0, Literal x1, Literal select);
    Literal addAndGate(const std::vector<Literal>& inputs);
    Literal addNandGate(const std::vector<Literal>& inputs);
    Literal addOrGate(const std::vector<Literal>& inputs);
    Literal addNorGate(const std::vector<Literal>& inputs);
    Literal addXorGate(const std::vector<Literal>& inputs);
    //@}

    //@{
    /**\name Dependency handling
     */
    std::size_t               numDependencies() const noexcept;
    std::size_t               numDependencies(Variable var) const noexcept;
    bool                      depends(Variable var1, Variable var2) const;
    void                      removeDependency(Variable var1, Variable var2);
    void                      addDependency(Variable var1, Variable var2);
    const std::set<Variable>& getDependencies(Variable var) const;
    bool                      dependenciesSubset(Variable var1, Variable var2) const;
    //@}

    //@{
    /** \name SAT/QBF-based test for satisfiability
     */
    void checkSAT();
    void checkUNSAT(std::size_t num_random_patterns = 0);
    //@}

    //@{
    /** \name Access to timers and statistics
     */
    double getTime(WhichTimer type) const noexcept { return _timers[static_cast<unsigned int>(type)].read(); }

    std::size_t getStatistics(Statistics which) const noexcept { return _statistics[static_cast<unsigned int>(which)]; }
    //@}

    //@{
    /**
     * \name Changing the preprocessor settings
     */
    Settings&       settings() noexcept { return _settings; }
    const Settings& settings() const noexcept { return _settings; }

    void enforceDQBF(bool val) noexcept;
    void setInterrupt(bool val) noexcept;

    //@}

    //@{
    /**\name Preprocessing techniques
     * \todo bool hyperBinaryResolution()
     * \todo bool removeRedundantClausesBySAT()
     */

    void                     preprocess();
    std::vector<Gate>        preprocessGates();
    const std::vector<Gate>& getGates() const noexcept { return _gates.getGates(); }

    // The preprocessor methods
    void                               fastPreprocess(bool until_fixedpoint = true);
    bool                               unitPropagation();
    bool                               findPure();
    bool                               findContradictions();
    bool                               findEquivalences();
    bool                               findHiddenEquivAndContraDefinitions();
    bool                               removeBlockedAndFriends();
    Literal                            checkImplicationChain(Literal lit);
    bool                               findImplicationChains();
    bool                               removeSubsumedClauses();
    bool                               universalReduction();
    bool                               applyResolution();
    bool                               selfSubsumingResolution();
    bool                               applySubstitution();
    bool                               applyDependencyScheme(DependencyScheme scheme, DependencyOperation operation);
    std::vector<std::vector<Variable>> identifyDontCares();
    bool                               applyUniversalExpansion();
    bool                               applyUniversalExpansion2();
    bool                               applyUniversalExpansionDQBF();
    bool                               addBlockedImplications();
    bool                               findConstantsBySAT();
    bool                               findImplicationsBySAT();
    bool                               findEquivalentGates();
    bool                               vivifyClause(ClauseID c_nr, SatPropagator& prop);
    bool                               applyVivification();
    bool                               upla();
    //@}

    //@{
    /**\name Gate detection
     */
    bool determineGates(bool and_gates = true, bool xor_gates = true, bool mux_gates = true,
                        bool semantic_gates = false);
    //@}

    //@{
    /**\name Elimination routines
     */
    std::vector<std::vector<Variable>> computeDepElimSet(bool unit_costs = false);
    std::vector<Variable>              computeVarElimSet();
    void                               elimEVar(Variable var, std::unordered_set<Variable>* recalc_vars = nullptr);
    bool elimEVarLimit(Variable var, std::int64_t max_cost, std::unordered_set<Variable>* recalc_vars = nullptr);
    void elimUVar(Literal lit);
    std::pair<Variable, Variable> elimDependency(Variable univ, Variable exist);
    //@}

    //{@
    /**\name Miscellaneous
     */
    bool              isQBF() const noexcept;
    void              convertToQBF();
    void              printStatistics() const;
    void              printSettings() const;
    const Prefix*     getPrefix() const noexcept { return _prefix; }
    const QBFPrefix*  getQBFPrefix() const noexcept { return _qbf_prefix; }
    const DQBFPrefix* getDQBFPrefix() const noexcept { return _dqbf_prefix; }
    //@}

    //@{
    /**\name Debugging functions
     */
    bool checkConsistency() const;
    //@}

    //@{
    /**
     \name Dependency schemes
    */
    bool stdTriDep(Variable u_var, bool resolution_paths, bool triangle, std::set<Variable>* pseudo_deps = nullptr);
    bool stdTriDep(bool resolution_paths, bool triangle);
    bool sstdQuadDep(Variable u_var, bool resolution_paths, bool quadrangle, std::set<Variable>* pseudo_deps = nullptr);
    bool sstdQuadDep(bool resolution_paths, bool quadrangle);
    bool invStdTriDep(Variable u_var, bool resolution_paths, bool triangle, std::set<Variable>* pseudo_deps = nullptr);
    bool invStdTriDep(bool resolution_paths, bool triangle);
    bool invSstdQuadDep(Variable u_var, bool resolution_paths, bool quadrangle,
                        std::set<Variable>* pseudo_deps = nullptr);
    bool invSstdQuadDep(bool resolution_paths, bool quadrangle);
    bool gateDependencies(DependencyOperation operation);
    //@}

   private:
    //@{
    /**
    \name Internal variable handling

    These functions are used internally to create variables with a predefined index.
    Normally, a newly created variable simply gets an arbitrary index that is not yet
    in use. After using these function, the list of available variable indicies
    Formula::_deleted_var_numbers needs to be re-created.
    */
    Variable setEVar(Variable index);
    template <typename Container>
    Variable setEVar(Variable index, const Container& deps);
    Variable setEVar(Variable index, std::set<Variable>&& deps);
    Variable setUVar(Variable index);
    //@}

    //@{
    /**
     \name Dependency schemes
    */
    template <typename Function>
    void searchPath(const std::vector<Literal>& start_lits, Function forbidden, std::vector<unsigned char>& seen) const;

    template <typename Function>
    void searchResolutionPath(const std::vector<Literal>& start_lits, Function forbidden,
                              std::vector<unsigned char>& seen) const;

    //@}

    //@{
    /// \name Gate detection
    bool        checkGatePrecondition(Variable output_var, const Clause& clause) const;
    std::size_t findMUXGates(std::vector<std::vector<Variable>>& gateDep);
    std::size_t findANDGates(std::vector<std::vector<Variable>>& gateDep, bool extend = true);
    std::size_t findXORGates(std::vector<std::vector<Variable>>& gateDep);
    std::size_t findSemanticGates(std::vector<std::vector<Variable>>& gateDep);
    //@}

    //@{
    /// \name Unit propagation and pure literals
    enum class PureStatus
    {
        UNIT,
        PURE
    };
    void pushUnit(Literal lit, PureStatus status = PureStatus::UNIT);
    bool checkPure(Literal lit) const;
    void checkPure(const Clause& clause, Literal except_lit);
    //@}

    //@{
    /// \name Internal functions for blocked clause elimination
    template <typename Container>
    Literal clauseBlocked(const Container& current_clause) const;
    bool    clauseBlockedByLit(Literal blocking_lit) const;
    template <typename Container>
    bool checkResolventTautology(const Container& clause, Variable pivot_var) const;
    bool addHiddenLiterals(int c_nr, Clause::ClauseData& clause, std::uint64_t& sign) const;
    bool addHiddenLiteralsBinary(int c_nr, Clause::ClauseData& clause, std::uint64_t& sign) const;
    bool addCoveredLiterals(Clause::ClauseData& clause, std::uint64_t& sign) const;
    template <typename Container>
    bool addBlockingLiterals(const Container& clause);
    void removeClauseAndUpdateCandidates(ClauseID c_nr);

    // Experimental function:
    bool pairBlockedClauses();
    //@}

    // Universal expansion
   public:
    std::pair<int, int> computeExpansionCosts(Variable uvar) const;
    std::int64_t        computeExpansionCosts2(Literal ulit, const std::set<Variable>& pseudo_deps);

   private:
    void markTransitiveUnits(std::stack<Literal>& units, std::vector<bool>& marked) const;

    // Subsumption checks
    template <typename Container>
    std::size_t isBackwardSubsuming(const Container& short_clause, std::uint64_t signature, int c_nr = -1,
                                    bool delete_subsumed = true);

    std::size_t isBackwardSubsuming(const Clause& short_clause, int c_nr = -1, bool delete_subsumed = true);

    template <typename Container>
    bool isForwardSubsumed(const Container& clause, std::uint64_t sign, int except = -1);

    bool isForwardSubsumed(const Clause& clause, int except = -1);

    template <typename Container>
    bool isForwardSubsumedByBinary(const Container& clause, int except = -1);
    template <typename Container>
    Literal getMinOccLit(const Container& clause) const;

    // Universal reduction
    bool universalReduction(Clause& clause, int c_nr = -1);

    // Resolution
    int                   computeResolutionCosts(Variable var) const;
    std::vector<Variable> getResolvableVariables() const;
    bool                  isResolvable(Variable var) const;

    // Gate substitution
    bool substituteGate(const Gate& g);
    int  computeSubstitutionCosts(const Gate& g) const;

    // Gate rewriting, done (according to sQueezeBF's way) when substitution would
    // be too costly
    void rewriteGate(Gate& g);

    // Equivalence Reduction
    std::size_t replaceSCC(const std::vector<Literal>& scc);

    // Clause modification
    int                   addClauseToLists(ClauseID c_nr, bool check_subsumption = true);
    bool                  removeLiteral(ClauseID c_nr, Literal lit);
    void                  replaceLiteral(Literal lit, Literal replacement);
    void                  replaceLiteralMono(Literal lit, Literal replacement);
    void                  addImplications(Literal lit1, Literal lit2,
                                          ClauseID c_nr);  // adds both implications, check for duplicates
    void                  addImplication(Literal lit1, Literal lit2, ClauseID c_nr);
    std::vector<bool>     impliedLiterals(Literal start_lit) const;
    std::pair<bool, bool> hasImplicationTransitive(Literal start_lit, Variable target_var) const;
    int                   hasImplication(Literal lit1, Literal lit2) const noexcept;
    int                   getImplicationClause(Literal lit1, Literal lit2) const noexcept;
    void                  removeImplication(Literal lit1, Literal lit2);
    void                  removeFromOccList(Literal lit, ClauseID c_nr);

    // Variable functions
    Variable nextVar();
    Variable copyVar(Variable var);

    /**
     * \brief Get access to the timer for a certain operation.
     */
    Timer& getTimer(WhichTimer type) { return _timers[static_cast<unsigned int>(type)]; }

    std::size_t& stat(Statistics which);

    // Helper functions
    void    initCandidateLists();
    void    countOccurences();
    void    countOccurencesVar();
    void    countImplications();
    int32_t varAssignment(Variable var) const noexcept;

    //@{
    /**\name Debugging functions
     */
    bool    checkSeen() const;
    Literal getAssignment(Literal lit) const;
    void    printFullOccurenceList(Variable var, std::ostream& stream) const;
    void    printOccurenceList(Literal lit, std::ostream& stream) const;
    void    printImplications(Literal lit, std::ostream& stream) const;
    void    printAllClauses(std::ostream& stream, bool print_implications = false) const;
    //@}

    // DATA

    /**
     * \brief The quantifier prefix.
     *
     * There are three versions for the quantifier prefix: Formula::_prefix
     * is the prefix representation shared by both QBFs and DQBFs. It provides
     * access to all operations that are common for both formula types.
     * In contrast, Formula::_qbf_prefix is available only for QBFs and
     * Formula::_dqbf_prefix only for DQBFs.
     */
    Prefix* _prefix{nullptr};

    /**
     * \brief The quantifier prefix for DQBFs only.
     *
     * It provides access to all functions to manipulate a DQBF prefix.
     * For all operations that can be applied to QBFs as well, one should
     * use Formula::_prefix. Formula::_dqbf_prefix is different from nullptr
     * if and only if the formula is a DQBF.
     */
    DQBFPrefix* _dqbf_prefix{nullptr};

    /**
     * \brief The quantifier prefix for QBFs only.
     *
     * It provides access to all functions to manipulate a QBF prefix.
     * For all operations that can be applied to DQBFs as well, one should
     * use Formula::_prefix. Formula::_qbf_prefix is different from nullptr
     * if and only if the formula is a QBF.
     */
    QBFPrefix* _qbf_prefix{nullptr};

    /**
     * \brief The quantifier prefix for SAT only.
     *
     * It provides access to all functions to manipulate a SAT prefix.
     * For all operations that can be applied to (D)QBFs as well, one should
     * use Formula::_prefix. Formula::_sat_prefix is different from nullptr
     * if and only if the formula is a SAT problem.
     */
    SATPrefix* _sat_prefix{nullptr};

    /**
     * \brief If true, every formula is treated as a DQBF, even it is actually a
     * QBF. \sa Formula::enforceDQBF(bool)
     */
    bool _enforce_dqbf = false;

    /**
     * \brief The list of clauses of the formula.
     *
     * Formula::_clauses contains all clauses of the formula
     * with the exception of unit clauses, which are stored in
     * Formula::unit_stack until they are propagated. The binary
     * clauses are not only stored in Formula::clauses, but also
     * redundantly in Formula::_implications. Note that some of
     * the entries of Formula::_clauses might not correspond to
     * valid clauses because they might have been deleted. This
     * can be checked using Formula::clauseDeleted(unsigned int).
     */
    std::vector<Clause> _clauses;

    /**
     * \brief Information about the detected gates in the formula.
     * \sa Formula::determineGates()
     */
    GateInfo _gates;

    /**
     * \brief Occurrence lists for all literals.
     *
     * For a literal \f$\ell\f$, Formula::_occ_list[\f$\ell\f$] contains the
     * IDs of all clauses that contain \f$\ell\f$.
     */
    std::vector<std::vector<ClauseID>> _occ_list;

    /**
     * \brief Stores the implication graph of the binary clauses.
     *
     * Each binary clause \f$(a,b)\f$ corresponds to the two implications
     * \f$\neg a\to b\f$ and \f$\neg b\to a\f$. These implications are stored
     * redundantly in Formula::_implications. This is useful not only for
     * equivalence reasoning, but speeds up subsumption checking, hidden
     * literal addition etc.
     */
    std::vector<ImplicationSet> _implications;

    /**
     * \brief Have new implications been added since the last call to
     * findEquivalences()?
     */
    bool _implications_added = false;

    /**
     * \brief Contains the IDs of those clauses which are deleted and therefore
     * available
     *
     * The IDs of clauses that got deleted are recycled when new clauses are
     * added. Formula::_deleted_clause_numbers is a list of these IDs.
     */
    std::vector<ClauseID> _deleted_clause_numbers;

    /**
     * \brief Contains the IDs of those variables which are deleted and therefore
     * available
     *
     * The IDs of variables that got deleted are recycled when new variables are
     * added. Formula::_deleted_var_numbers is a list of these IDs.
     */
    std::vector<Variable> _deleted_var_numbers;

    /**
     * \brief Stack with unit literals that need to be propagated
     *
     * This stack is used by Formula::pushUnit() and Formula::unitPropagation()
     * for propagating unit literals.
     */
    std::vector<Literal> _unit_stack;

    /**
     * \brief Assignment of unit literals; used by Formula::pushUnit() to detect
     * conflicts
     *
     * For each unit literal that is propagated, we store its assignment.
     * We can then recognize conflicting unit clauses (which indicate that the
     * formula is unsatisfiable).
     */
    std::vector<TruthValue> _assignment;

    /**
     * \brief States whether a variable may not be touched.
     *
     * This enables us to perform preprocessing for incremental (D)QBF solving.
     * Variables which will appear in clauses that are added later may not be
     * touched (i.e. eliminated, serve as blocking literal of clauses etc.)
     * \note Currently this is not used.
     */
    std::vector<bool> _dont_touch;

    /**
     * \brief Vector with costs for each variable, literal, or clause
     *
     * This vector is used together with Formula::candidates to implement
     * a priority queue (min heap). It is used, e.g., to eliminate variables
     * with low costs first.
     */
    std::vector<int> _variable_score;

    /**
     * \brief Vector with the length of each clause.
     *
     * This is required for usage in a priority queue.
     */
    std::vector<std::size_t> _clause_sizes;

    /**
     * \brief Priority queue (min heap) for variables, literals, or clauses
     *
     * This priority queue is used together with Formula::_variable_score.
     * It is used, e.g., to eliminate variables with low costs first.
     */
    VarHeap<AscendingOrder, Variable, int> _candidates;

    /**
     * \brief Priority queue (max heap) for usage during "removeBlockedAndFriends"
     */
    VarHeap<DescendingOrder, ClauseID, std::size_t> _blocked_candidates;

    /**
     * \brief Used to temporarily store removed literals during universal
     * reduction
     *
     * This vector is used during universal reduction to store the removed
     * literals. This avoid creating the vector each time a clause is added
     * or modified.
     */
    std::vector<Literal> _removed_lits;

    /**
     * \brief A helper vector to store seen literals
     *
     * This vector is used by various methods to memorize which
     * literals have already be seen.
     * Two data structures are necessary since there are some methods where two
     * independent literals sets have to remembered. \note Due to performance
     * reasons the content of the vector is kept within every hidden/covered
     * literal addition, block clause and friends methods as well as subsumption
     * methods => You have to call "setSeen" and "cleanSeen" manually outside
     * these methods. \note All methods assume that the vector is cleared at the
     * beginning of its call (apart from mentioned methods above). \note seen is
     * currently used for all blocked clause and friends methods (+hidden,
     * covered), dependency schemes, hidden equivalence and contradiction checks,
     * and constant SAT checks \note seen2 is currently used in subsumption checks
     * and intersection for covered literals \note The user is responsible for
     * clearing the vector, otherwise bad things can happen.
     */
    mutable std::vector<unsigned char> _seen;
    mutable std::vector<unsigned char> _seen2;

    void clearSeen() const;
    template <typename Container>
    void clearSeen(const Container& container) const;
    template <typename Container>
    void setSeen(const Container& container) const;
    template <typename Container>
    void clearSeen2(const Container& container) const;
    template <typename Container>
    void setSeen2(const Container& container) const;

#ifdef SKOLEM
    std::vector<std::unique_ptr<SkolemEntry>> _skolem_data;
#endif

    /**
     * \brief Preprocessor settings
     *
     * Specifies which techniques should be applied to optimize
     * the formula.
     */
    Settings _settings;

    /**
     * \brief Used to abort certain preprocessor methods if they take too long
     *
     * For each method, Formula::_process_limit contains a counter which is
     * decremented each time a certain operation is performed (like resolving
     * a clause). If the counter becomes zero, the preprocessing method is
     * aborted.
     */
    mutable ProcessLimit _process_limit;

    /**
     * \brief If set to true, the preprocessor stops as soon as possible.
     */
    bool _interrupt = false;

    mutable std::array<std::size_t, static_cast<unsigned int>(Statistics::_last_)>
        _statistics;  ///< Statistics for the different operations
    mutable std::array<Timer, static_cast<unsigned int>(WhichTimer::_last_)>
        _timers;  ///< Timers for the different operations

    // fork splitting
    bool forkExtension();
    bool inDEPlus();
    bool inOneClause(Variable var1, Variable var2);

 	std::vector<std::vector<bool>> belongIntoOneClause;

	bool resolveAfterForkSplit();
	bool resolveRmb();

	long max_resolveRmb_cost = 0; // cost used while performing resolveRmb 
	
	void findHiddenEquivalences();

	
	class Node
	{
	public:

		enum Side {unknown, left, right};

		Literal lit = 0;
		Side side = unknown;
		std::vector<Node*> adjacents;
		bool seen = false;
	};

	std::vector<Node*> nodes;
	//std::vector<Node*> equivNodesLeft;
	//std::vector<Node*> equivNodesRight;

	// TODO: find test cases for case that var occur on one side but in different polarities
	bool checkForEquivalence(Node* node, std::vector<Node*>& equivNodesLeft, std::vector<Node*>& equivNodesRight);

	void replaceEquivalentLiterals(const std::vector<Node*>& equivNodesLeft, const std::vector<Node*>& equivNodesRight);

	void solveSAT();
};

/**
 * \brief Output operator for DQBFs.
 *
 * Writes a DQBF in DQDIMACS format to the given stream.
 * \param stream the output stream the formula is written to
 * \param formula the formula that is to be written
 * \return the stream
 * \relates Formula
 */
std::ostream& operator<<(std::ostream& stream, const Formula& formula);

/**
 * \brief Input operator for (D)QBFs (in (D)QDIMACS format)
 *
 * Reads a formula in (D)QDIMACS format from an input stream.
 * \param stream the input stream the formula is read from
 * \param formula the data structure the read formula is stored in
 * \return the stream
 * \relates Formula
 */
std::istream& operator>>(std::istream& stream, Formula& formula);

}  // end namespace hqspre

#include "formula.ipp"
#include "formula_blocked_et_al.ipp"
#include "formula_subsumption.ipp"

#endif
