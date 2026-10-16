/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSets.h"
#include "oops/util/Logger.h"

namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
}

namespace oops {
  class FieldSet3D;
  class FieldSet4D;
  template <class MODEL> class Geometry;
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// @brief Inputs prepared by the caller and handed to a block chain at
/// construction, for cases where the caller owns information the block chain
/// cannot obtain on its own.
///
/// The motivating case is Parallel Hybrid: only the parent communicator can see
/// the whole ensemble, so only it can form perturbations about a single common
/// mean. Each component is then handed its own members rather than reading them.
///
/// Move-only. Ownership of the supplied ensemble transfers to the block chain,
/// and an empty context means "source your own inputs, as before".
class BlockChainContext {
 public:
  BlockChainContext() = default;
  explicit BlockChainContext(std::unique_ptr<oops::FieldSets> suppliedEnsemble)
    : suppliedEnsemble_(std::move(suppliedEnsemble)) {}

  BlockChainContext(BlockChainContext &&) = default;
  BlockChainContext & operator=(BlockChainContext &&) = default;
  BlockChainContext(const BlockChainContext &) = delete;
  BlockChainContext & operator=(const BlockChainContext &) = delete;

  /// @brief Whether the context carries anything at all.
  bool empty() const {return suppliedEnsemble_ == nullptr;}

  /// @brief Ensemble perturbations already prepared on the geometry passed to
  /// the block chain, unnormalized. Null when the chain should read its own.
  std::unique_ptr<oops::FieldSets> & suppliedEnsemble() {return suppliedEnsemble_;}
  const std::unique_ptr<oops::FieldSets> & suppliedEnsemble() const {return suppliedEnsemble_;}

 private:
  std::unique_ptr<oops::FieldSets> suppliedEnsemble_;
};

// -----------------------------------------------------------------------------
/// Base class for SABER block chains that have a self-adjoint central block
/// (ensemble and non-ensemble).
class BlockChainBase {
 public:
  BlockChainBase() = default;
  virtual ~BlockChainBase() = default;

  virtual void randomize(oops::FieldSet4D &) const = 0;
  virtual void multiply(oops::FieldSet4D &) const = 0;
  virtual size_t ctlVecSize() const = 0;
  virtual void randomCtlVec(atlas::Field &, const size_t &) const = 0;
  virtual void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const = 0;
  virtual void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const = 0;
  virtual const atlas::FunctionSpace & outerFunctionSpace() const = 0;
  virtual const oops::Variables & outerVariables() const = 0;

  /// @brief Diagonal variance of this block chain on its outer function space.
  virtual oops::FieldSet3D variance() const = 0;
};

template<typename MODEL>
class BlockChainFactory {
 public:
  typedef oops::Geometry<MODEL> Geometry_;

  static std::unique_ptr<BlockChainBase> create(const Geometry_ &,
                                                const oops::Variables &,
                                                oops::FieldSet4D &,
                                                oops::FieldSet4D &,
                                                const eckit::Configuration &,
                                                BlockChainContext && = BlockChainContext());

  virtual ~BlockChainFactory() = default;

 protected:
  explicit BlockChainFactory(const std::string &);

 private:
  virtual std::unique_ptr<BlockChainBase> make(const Geometry_ &,
                                               const oops::Variables &,
                                               oops::FieldSet4D &,
                                               oops::FieldSet4D &,
                                               const eckit::Configuration &,
                                               BlockChainContext &&) = 0;

  static std::map <std::string, BlockChainFactory<MODEL> *> & getMakers() {
    static std::map <std::string, BlockChainFactory<MODEL> *> makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class BlockChainMaker : public BlockChainFactory<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;

  /// @brief Forward the construction context only to block chains that accept
  /// one. Chains that do not keep their existing constructor untouched, and
  /// being handed a non-empty context is an error rather than a silent drop.
  std::unique_ptr<BlockChainBase> make(const Geometry_ & geom,
                                       const oops::Variables & outerVars,
                                       oops::FieldSet4D & fset4dXb,
                                       oops::FieldSet4D & fset4dFg,
                                       const eckit::Configuration & conf,
                                       BlockChainContext && context) override {
    if constexpr (std::is_constructible_v<T, const Geometry_ &, const oops::Variables &,
                                          oops::FieldSet4D &, oops::FieldSet4D &,
                                          const eckit::Configuration &,
                                          BlockChainContext &&>) {
      return std::make_unique<T>(geom, outerVars, fset4dXb, fset4dFg, conf, std::move(context));
    } else {
      if (!context.empty()) {
        throw eckit::UserError("Block chain '" + name_ + "' does not accept a supplied "
                               "ensemble, but one was provided.", Here());
      }
      return std::make_unique<T>(geom, outerVars, fset4dXb, fset4dFg, conf);
    }
  }

  /// @brief Registered name, kept for diagnostics.
  std::string name_;

 public:
  explicit BlockChainMaker(const std::string & name)
    : BlockChainFactory<MODEL>(name), name_(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
BlockChainFactory<MODEL>::BlockChainFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(name + " already registered in saber::BlockChainFactory.",
                              Here());
  getMakers()[name] = this;
}

template <typename MODEL>
std::unique_ptr<BlockChainBase>
BlockChainFactory<MODEL>::create(const Geometry_ & geom,
                                 const oops::Variables & outerVars,
                                 oops::FieldSet4D & fset4dXb,
                                 oops::FieldSet4D & fset4dFg,
                                 const eckit::Configuration & conf,
                                 BlockChainContext && context) {
  oops::Log::trace() << "BlockChainFactory<MODEL>::create starting" << std::endl;
  std::string name = "parametric";
  if (conf.has("covariance type")) {
    name = conf.getString("covariance type");
  }
  typename std::map<std::string, BlockChainFactory<MODEL>*>::iterator jbc =
    getMakers().find(name);
  if (jbc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in saber::BlockChainFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<BlockChainBase> ptr =
    jbc->second->make(geom, outerVars, fset4dXb, fset4dFg, conf, std::move(context));
  oops::Log::trace() << "BlockChainFactory<MODEL>::create done" << std::endl;
  return ptr;
}

}  // namespace saber
