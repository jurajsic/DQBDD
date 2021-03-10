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

#ifndef HQSPRE_PREFIX_HPP_
#define HQSPRE_PREFIX_HPP_

#include <iosfwd>
#include <set>
#include <vector>

#include "aux.hpp"
#include "literal.hpp"

/**
 * \file prefix.hpp
 * \brief Classes for representing a formula's prefix
 * \author Ralf Wimmer
 * \date 06/2016
 */

namespace hqspre {

/**
 * \todo Don't store the variables in the RMB explicitly in the universal
 * variables' dependency sets.
 */

/**
 * \brief Status of a variable (quantifier/deleted)
 */
enum class VariableStatus : short
{
    EXISTENTIAL,  ///< Variable is existential
    UNIVERSAL,    ///< Variable is universal
    DELETED       ///< Variable has been deleted
};

/**
 * \brief Type of the prefix (QBF/DQBF).
 */
enum class PrefixType : short
{
    SAT,  ///< Formula is purely existential
    QBF,  ///< Formula is a standard QBF
    DQBF  ///< Formula is a DQBF
};

class Prefix;
class SATPrefix;
class QBFPrefix;
class DQBFPrefix;

/**
 * \brief Representation of a quantifed formula's prefix.
 */
class Prefix
{
   public:
    Prefix() : _var_status(1, VariableStatus::DELETED) {}
    Prefix(const Prefix& other) = default;
    Prefix(Prefix&& other)      = default;
    virtual ~Prefix() noexcept  = default;
    Prefix& operator=(const Prefix& other) = default;
    Prefix& operator=(Prefix&& other) = default;

    /**
     * \brief Returns the type of the prefix (QBF vs. DQBF).
     */
    virtual PrefixType type() const noexcept = 0;

    virtual void clear();

    //@{
    /**
     * \name Creating and deleting variables.
     */
    virtual void          setMaxVarIndex(Variable var);
    virtual void          addEVar(Variable var);
    virtual void          addUVar(Variable var);
    virtual void          removeVar(Variable var);
    std::vector<Variable> getExistVars() const noexcept;
    std::vector<Variable> getUnivVars() const noexcept;
    virtual void          makeCopy(Variable var, Variable original);
    virtual void          updateVars();

    virtual bool       isSAT() const noexcept { return _num_u_vars == 0; }
    virtual SATPrefix* convertToSAT() const;
    //@}

    //@{
    /**
     * \name Querying the variable status
     */
    bool           isUniversal(Variable var) const noexcept;
    bool           isExistential(Variable var) const noexcept;
    bool           varDeleted(Variable var) const noexcept;
    VariableStatus getStatus(Variable var) const noexcept;
    //@}

    //@{
    /**
     * \name Number of variables and minimal/maximal indices
     */

    /// Returns the minimal possible index of a variable.
    constexpr static Variable minVarIndex() noexcept { return 1; }

    /// Returns the maximal possible index of a variable.
    Variable maxVarIndex() const noexcept { return static_cast<Variable>(_var_status.size() - 1); }

    /// Returns the minimal possible index of a literal.
    constexpr static Literal minLitIndex() noexcept { return 2u; }

    /// Returns the maximal possible index of a literal.
    Literal maxLitIndex() const noexcept { return static_cast<Literal>(2 * _var_status.size() - 1); }

    /// Returns the number of (non-deleted) variables.
    std::size_t numVars() const noexcept { return _num_e_vars + _num_u_vars; }

    /// Returns the number of existential variables.
    std::size_t numEVars() const noexcept { return _num_e_vars; }

    /// Returns the number of universal variables.
    std::size_t numUVars() const noexcept { return _num_u_vars; }
    //@}

    //@{
    /**
     * \name Access to dependencies
     */
    /// Returns the number of dependencies of variable `var`
    virtual std::size_t numDependencies(Variable var) const noexcept = 0;

    // Returns the accumulated number of dependencies of all existential variables
    virtual std::size_t numDependencies() const noexcept = 0;

    /// Returns whether `var1` depends on `var2`
    virtual bool depends(Variable var1, Variable var2) const = 0;

    /// Returns whether `var1`'s dependencies are a subset of `var2`'s
    /// dependencies
    virtual bool dependenciesSubset(Variable var1, Variable var2) const = 0;

    /// Replaces `var1`'s dependencies by the intersection of `var1` and `var2`'s
    /// dependencies
    virtual void intersectDependencies(Variable var1, Variable var2) = 0;

    /// Replaces `var1`'s dependencies by the union of `var1` and `var2`'s
    /// dependencies
    virtual void uniteDependencies(Variable var1, Variable var2) = 0;

    /// Moves an existential variable to the right-most block (i.e., makes it
    /// depend on all universals)
    virtual void moveToRMB(Variable index) = 0;

    /// Checks if a variable is contained in the right-most block
    virtual bool inRMB(Variable var) const noexcept = 0;

    /// Checks if an existential variable is in the left-most block (i.e. has no dependencies)
    virtual bool inLMB(Variable var) const noexcept = 0;

    /// Returns the variables in the right-most block
    virtual const std::set<Variable>& getRMB() const noexcept = 0;
    //@}

    //@{
    /**
     * \name Input/output
     */

    /*
     * \brief Writes the prefix to an output stream.
     * \param stream the stream to write to
     * \param translation_table can be used to re-number the variable ID such that
     *        there are no holes due to deleted variables. The table needs to be
     * an injective partial mapping that assigns to each non-deleted variable
     * index a new unique index.
     */
    virtual void write(std::ostream& stream, std::vector<Variable>* translation_table = nullptr) const = 0;
    //@}

    //@{
    /**
     * \name Miscellaneous
     */
    virtual bool checkConsistency() const;
    //@}

   private:
    std::vector<VariableStatus> _var_status;      ///< Status of all variables (existential/universal/deleted)
    std::size_t                 _num_u_vars = 0;  ///< Number of universal variables
    std::size_t                 _num_e_vars = 0;  ///< Number of existential variables
};

/**
 * \brief Representation of a SAT problem (purely existential prefix).
 */
class SATPrefix : public Prefix
{
   public:
    SATPrefix()                       = default;
    SATPrefix(const SATPrefix& other) = default;
    SATPrefix(SATPrefix&& other)      = default;
    virtual ~SATPrefix() noexcept     = default;
    SATPrefix& operator=(const SATPrefix& other) = default;
    SATPrefix& operator=(SATPrefix&& other) = default;

    PrefixType type() const noexcept override;
    void       clear() override;

    //@{
    /**
     * \name Creating and deleting variables.
     */
    void addEVar(Variable var) override;
    void addUVar(Variable) override;
    void makeCopy(Variable var, Variable original) override;
    void updateVars() override;
    void removeVar(Variable var) override;
    //@}

    //@{
    /**
     * \name Access to dependencies
     */
    std::size_t               numDependencies(Variable) const noexcept override;
    std::size_t               numDependencies() const noexcept override;
    bool                      depends(Variable, Variable) const override;
    bool                      dependenciesSubset(Variable, Variable) const override;
    void                      intersectDependencies(Variable, Variable) override;
    void                      uniteDependencies(Variable, Variable) override;
    void                      moveToRMB(Variable) override;
    bool                      inRMB(Variable) const noexcept override { return true; }
    bool                      inLMB(Variable) const noexcept override { return true; }
    const std::set<Variable>& getRMB() const noexcept override;
    //@}

    //@{
    /**
     * \name Input/output
     */

    void write(std::ostream& stream, std::vector<Variable>* translation_table = nullptr) const override;
    //@}

    bool checkConsistency() const override;

   private:
    std::set<Variable> _vars;
};

/**
 * \brief Representation of a QBF's prefix.
 */
class QBFPrefix : public Prefix
{
   public:
    QBFPrefix() : Prefix(), _depth(1, 0), _blocks(1, std::set<Variable>()) {}
    QBFPrefix(const QBFPrefix& other) = default;
    QBFPrefix(QBFPrefix&& other)      = default;
    ~QBFPrefix() noexcept override    = default;
    QBFPrefix& operator=(const QBFPrefix& other) = default;
    QBFPrefix& operator=(QBFPrefix&& other) = default;

    PrefixType type() const noexcept override;
    void       clear() override;
    bool       checkConsistency() const override;

    //@{
    /**
     * \name Creating and deleting variables
     */
    void addEVar(Variable var) override;
    void addUVar(Variable var) override;
    void makeCopy(Variable var, Variable original) override;
    void setMaxVarIndex(Variable var) override;
    void removeVar(Variable var) override;
    void updateVars() override;
    //@}

    //@{
    /**
     * \name Access to dependencies
     */
    bool        depends(Variable var1, Variable var2) const override;
    bool        dependenciesSubset(Variable var1, Variable var2) const override;
    void        write(std::ostream& stream, std::vector<Variable>* translation_table = nullptr) const override;
    void        intersectDependencies(Variable var1, Variable var2) override;
    void        uniteDependencies(Variable var1, Variable var2) override;
    std::size_t numDependencies(Variable var) const noexcept override;
    std::size_t numDependencies() const noexcept override;
    void        moveToRMB(Variable var) override;
    const std::set<Variable>& getRMB() const noexcept override;
    bool inRMB(Variable var) const noexcept override { return _depth[var] == (_blocks.size() - 1); }
    bool inLMB(Variable var) const noexcept override { return _depth[var] == 0 && isExistential(var); }
    //@}

    //@{
    /**
     * \name QBF-specific functions
     */
    DQBFPrefix*               convertToDQBF() const;
    std::size_t               getLevel(Variable var) const noexcept;
    std::size_t               getMaxLevel() const noexcept { return _blocks.size() - 1; }
    VariableStatus            getLevelQuantifier(std::size_t level) const noexcept;
    void                      setLevel(Variable var, std::size_t level);
    const std::set<Variable>& getVarBlock(std::size_t level) const noexcept;
    //@}

   private:
    std::vector<std::size_t> _depth;          ///< Assigns to each variable its quantifier
                                              ///< level (starting from 0)
    std::vector<std::set<Variable>> _blocks;  ///< Contains the variables in each quantifier block
};

/**
 * \brief Representation of a DQBF's prefix.
 * \todo Don't store the variables in the RMB explicitly in the universal
 * variables' dependency sets.
 */
class DQBFPrefix : public Prefix
{
   public:
    DQBFPrefix() : _rmb(), _dependencies(1), _in_rmb(1, false), _univ_vars() {}
    DQBFPrefix(const DQBFPrefix& other) = default;
    DQBFPrefix(DQBFPrefix&& other)      = default;
    ~DQBFPrefix() noexcept override     = default;
    DQBFPrefix& operator=(const DQBFPrefix& other) = default;
    DQBFPrefix& operator=(DQBFPrefix&& other) = default;

    PrefixType type() const noexcept override;
    void       clear() override;

    //@{
    /**
     * \name Variable handling
     */
    void addEVar(Variable var) override;
    template <typename Container>
    void addEVar(Variable var, const Container& dependencies);
    void addEVar(Variable var, std::set<Variable>&& dependencies);
    void addUVar(Variable var) override;
    void removeVar(Variable var) override;
    void makeCopy(Variable var, Variable original) override;
    void setMaxVarIndex(Variable var) override;
    //@}

    //@{
    /**
     * \name Dependency handling
     */
    bool                      depends(Variable var1, Variable var2) const override;
    bool                      dependenciesSubset(Variable var1, Variable var2) const override;
    void                      intersectDependencies(Variable var1, Variable var2) override;
    void                      uniteDependencies(Variable var1, Variable var2) override;
    std::size_t               numDependencies() const noexcept override;
    std::size_t               numDependencies(Variable var) const noexcept override;
    const std::set<Variable>& getRMB() const noexcept override;
    bool                      inRMB(Variable var) const noexcept override { return _in_rmb[var]; }
    bool                      inLMB(Variable var) const noexcept override { return isExistential(var) && _dependencies[var].empty(); }
    void                      moveToRMB(Variable var) override;

    template <typename Container>
    void                      setDependencies(Variable var, const Container& deps);
    void                      setDependencies(Variable var, std::set<Variable>&& deps);
    void                      removeDependency(Variable var1, Variable var2);
    void                      addDependency(Variable var1, Variable var2);
    const std::set<Variable>& getDependencies(Variable var) const noexcept;
    void                      clearDependencies(Variable var);
    //@}

    bool checkConsistency() const override;
    void write(std::ostream& stream, std::vector<Variable>* translation_table = nullptr) const override;

    bool       isQBF() const;
    QBFPrefix* convertToQBF() const;

   private:
    /**
     * \brief The variables in the right-most block
     *
     * The variables which depend on all universal variables (like
     * the Tseitin variables) are stored separately in order to save
     * memory. `rmb` (for ``right-most block'') is the set of all these
     * variables. Their entry in `dependencies` is empty; however
     * they appear in the dependency lists of all universal variables.
     * \sa DQBFPrefix::in_rmb
     * \sa DQBFPrefix::dependencies
     */
    std::set<Variable> _rmb;

    /**
     * \brief The dependency sets of all variables except those in the right-most
     * block
     *
     * Apart from the variables in the right-most block, each variable has an
     * entry in dependencies. For an existential variable `evar`,
     * `dependencies[evar]` contains all universal variables on which `evar`
     * depends. For a universal variable `uvar`, `dependencies[uvar]` contains all
     * existential variables which depend on `uvar'. \sa DQBFPrefix::rmb \sa
     * DQBFPrefix::in_rmb
     */
    std::vector<std::set<Variable>> _dependencies;

    /**
     * \brief For a variable \f$v\f$, in_rmb[\f$v\f$] is true iff \f$v\f$ depends
     * on all universal variables.
     *
     * Only existential variables can be in the right-most block.
     * \sa DQBFPrefix::rmb
     * \sa DQBFPrefix::dependencies
     */
    std::vector<bool> _in_rmb;

    /**
     * \brief The set of all universal variables.
     */
    std::set<Variable> _univ_vars;
};

}  // end namespace hqspre

#include "prefix.ipp"

#endif
