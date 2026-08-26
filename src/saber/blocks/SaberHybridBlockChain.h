/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetSubCommunicators.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

class CovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters, oops::Parameters)
 public:
  oops::ConfigurationParameter saberBlockChainParams{this};
};

// -----------------------------------------------------------------------------

class WeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WeightParameters, oops::Parameters)
 public:
  // Scalar weight
  oops::Parameter<double> value{"value", 1.0, this};

  // File-base weight
  oops::OptionalParameter<eckit::LocalConfiguration> file{"file", this};
};

// -----------------------------------------------------------------------------

class ComponentParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ComponentParameters, oops::Parameters)
 public:
  // Covariance
  oops::RequiredParameter<CovarianceParameters> covariance{"covariance", this};
  // Weight
  oops::RequiredParameter<WeightParameters> weight{"weight", this};
};

// -----------------------------------------------------------------------------

/// @brief An ensemble read once, on the parent communicator, and shared by all
/// components of a parallel hybrid.
///
/// Only the parent communicator can see every member, so only it can form
/// perturbations about a single common mean. Each component is then handed the
/// members assigned to it. The perturbations are passed on unnormalized: the
/// 1/(N-1) factor is carried by the component weights, which are checked
/// against `normalization denominator`.
class SharedEnsembleParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SharedEnsembleParameters, oops::Parameters)
 public:
  /// Ensemble of states. This is the only source form currently supported,
  /// because it is the only one for which readEnsemble subtracts a mean.
  oops::RequiredParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};

  /// Denominator of the ensemble covariance normalization, conventionally N-1.
  oops::RequiredParameter<double> denominator{"normalization denominator", this};

  /// Relative tolerance when checking component weights against 1/denominator.
  oops::Parameter<double> weightTolerance{"weight tolerance", 1.0e-10, this};
};

// -----------------------------------------------------------------------------

class SaberHybridBlockChainParameters: public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(SaberHybridBlockChainParameters,
                           ErrorCovarianceParametersBase)
 public:
  // Optional outer blocks
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocksParams{"saber outer blocks", this};
  // Vector of components
  oops::RequiredParameter<std::vector<ComponentParameters>> components{"components", this};
  // Geometry [optional]
  oops::OptionalParameter<eckit::LocalConfiguration> hybridGeometry{"geometry", this};
  // Switch to run components in parallel
  oops::Parameter<bool> runInParallel{"run in parallel", false, this};
  // Switch to run components recursively (for diagnostics)
  oops::Parameter<bool> runComponentsRecursively{"run components recursively", false, this};
  // Resource weighting for each hybrid component.
  oops::OptionalParameter<std::vector<double>> parallelCovarRelativeCPUWeight{
      "parallel covariance relative cpu weight", this};

  // Which method to use in the atlas field redistribution between subcommunicator
  // and parent communicator.
  oops::Parameter<std::string> commRedistributionMethod{"comm redistribution method",
      "straight", this};

  // Ensemble read once on the parent communicator and shared across components.
  // Requires `run in parallel`.
  oops::OptionalParameter<SharedEnsembleParameters> sharedEnsemble{"shared ensemble", this};
};

// -----------------------------------------------------------------------------

/// @brief Collect and validate the member-to-component mapping for a shared
/// ensemble.
///
/// Every rank holds the full component configuration, so the whole mapping can
/// be checked before a single file is opened. Returns one ascending member list
/// per component, using the 1-based numbering of the shared ensemble.
inline std::vector<std::vector<size_t>> sharedEnsembleMemberMap(
    const std::vector<ComponentParameters> & components,
    const SharedEnsembleParameters & sharedParams) {
  std::vector<std::vector<size_t>> memberMap;
  memberMap.reserve(components.size());

  for (size_t component = 0; component < components.size(); ++component) {
    const eckit::LocalConfiguration cmpConf =
      components[component].covariance.value().toConfiguration();
    const std::string label = "Component " + std::to_string(component + 1);

    if (!cmpConf.has("shared ensemble members")) {
      throw eckit::UserError(label + " of a hybrid block declaring `shared ensemble` does not "
                             "configure `shared ensemble members`.", Here());
    }
    const std::vector<int> raw = cmpConf.getIntVector("shared ensemble members");
    if (raw.empty()) {
      throw eckit::UserError(label + " has an empty `shared ensemble members` list.", Here());
    }
    std::vector<size_t> members;
    members.reserve(raw.size());
    for (const int member : raw) {
      if (member < 1) {
        throw eckit::UserError(label + " lists member " + std::to_string(member) + " in "
                               "`shared ensemble members`; members are numbered from 1.",
                               Here());
      }
      members.push_back(static_cast<size_t>(member));
    }
    // Ascending order keeps the local member index monotonic in the global one.
    std::sort(members.begin(), members.end());
    memberMap.push_back(std::move(members));
  }

  // Every member exactly once, numbered 1..N with no gaps and no duplicates.
  std::vector<size_t> all;
  for (const auto & members : memberMap) {
    all.insert(all.end(), members.begin(), members.end());
  }
  std::sort(all.begin(), all.end());
  for (size_t jj = 0; jj < all.size(); ++jj) {
    if (all[jj] != jj + 1) {
      throw eckit::UserError("`shared ensemble members` across all components must list every "
                             "member from 1 to " + std::to_string(all.size()) + " exactly once. "
                             "Found " + std::to_string(all[jj]) + " where " +
                             std::to_string(jj + 1) + " was expected.", Here());
    }
  }

  // Every component weight must carry the global normalization, since a supplied
  // ensemble bypasses the scaling applied at read time.
  const double expectedWeight = 1.0/sharedParams.denominator.value();
  const double tolerance = sharedParams.weightTolerance.value();
  for (size_t component = 0; component < components.size(); ++component) {
    const double weight = components[component].weight.value().value.value();
    if (std::abs(weight - expectedWeight) > tolerance*std::abs(expectedWeight)) {
      throw eckit::UserError("Component " + std::to_string(component + 1) + " has weight " +
                             std::to_string(weight) + ", but a shared ensemble with "
                             "`normalization denominator` " +
                             std::to_string(sharedParams.denominator.value()) +
                             " requires " + std::to_string(expectedWeight) + ".", Here());
    }
  }

  return memberMap;
}

/// Hybrid covariance block chain implementation
template<typename MODEL>
class SaberHybridBlockChain : public SaberBlockChainBase {
 public:
  SaberHybridBlockChain(const oops::Geometry<MODEL> & geom,
                        const oops::Variables & outerVars,
                        oops::FieldSet4D & fset4dXb,
                        oops::FieldSet4D & fset4dFg,
                        const eckit::Configuration & conf);
  ~SaberHybridBlockChain() = default;

  /// @brief Randomize the increment according to this hybrid B matrix.
  void randomize(oops::FieldSet4D &) const override;
  /// @brief Multiply the increment by this hybrid B matrix.
  void multiply(oops::FieldSet4D &) const override;

  /// @brief Control vector size
  size_t ctlVecSize() const override
    {throw eckit::NotImplemented("ctlVecSize not implemented yet", Here());}
  /// @brief Generate a random control vector.
  void randomCtlVec(atlas::Field &, const size_t &) const override
    {throw eckit::NotImplemented("randomCtlVec not implemented yet", Here());}
  /// @brief Square-root multiplication
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const override
    {throw eckit::NotImplemented("multiplySqrt not implemented yet", Here());}
  /// @brief Adjoint of square-root multiplication
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const override
    {throw eckit::NotImplemented("multiplySqrtAD not implemented yet", Here());}

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const override {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::Variables & outerVariables() const override {return outerVariables_;}

 private:
  /// Function space
  const atlas::FunctionSpace & outerFunctionSpace_;
  /// Variables
  const oops::Variables outerVariables_;

  /// Chain of outer blocks applied to all components of hybrid covariances.
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// Vector of hybrid B components.
  std::vector<std::unique_ptr<SaberBlockChainBase>> hybridBlockChain_;
  /// Vector of scalar weights for hybrid B components.
  std::vector<double> hybridScalarWeightSqrt_;
  /// Vector of field weights for hybrid B components.
  std::vector<oops::FieldSet3D> hybridFieldWeightSqrt_;

  /// Whether to run Hybrid in parallel
  bool parallelHybrid_;
  /// Index of component if running Hybrid in parallel
  size_t myComponent_;
  /// local geometry for parallel Hybrid block (type-erased)
  std::shared_ptr<atlas::FunctionSpace> localHybridFs_, globalHybridFs_;

  /// Which method to use for comm redistribution, if any.
  std::string redistributionMethod_;

  std::unique_ptr<oops::Geometry<MODEL>> localHybridGeom_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberHybridBlockChain<MODEL>::SaberHybridBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.functionSpace()), outerVariables_(outerVars),
    parallelHybrid_(false), myComponent_(0),
    redistributionMethod_{""}
{
  oops::Log::trace() << "SaberHybridBlockChain ctor starting" << std::endl;

  // Deserialize parameters and fill configuration with missing values
  SaberHybridBlockChainParameters params;
  params.deserialize(conf);
  eckit::LocalConfiguration fullConf;
  params.serialize(fullConf);

  // Extract ErrorCovarianceParametersBase from fullConf
  ErrorCovarianceParametersBase paramsBase;
  paramsBase.deserialize(fullConf);

  // Initialize current outer variables
  oops::Variables currentOuterVars(outerVars);

  // Build common (for all hybrid components) outer blocks if they exist
  if (params.saberOuterBlocksParams.value()) {
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVariables_,
                          fset4dXb, fset4dFg, fullConf,
                          *params.saberOuterBlocksParams.value());
    currentOuterVars = outerBlockChain_->innerVars();
  }

  // Hybrid central block
  parallelHybrid_ = params.runInParallel;
  redistributionMethod_ = params.commRedistributionMethod;

  if (params.sharedEnsemble.value() && !parallelHybrid_) {
    throw eckit::UserError("`shared ensemble` requires `run in parallel: true`. Run serially "
                           "and each component reads its own ensemble, as before.", Here());
  }

  const eckit::mpi::Comm & defaultSpaceComm = geom.getComm();
  const size_t ntasks = defaultSpaceComm.size();
  const size_t nComponents = params.components.value().size();
  globalHybridFs_.reset(new atlas::FunctionSpace(geom.functionSpace()));

  std::vector<double> parallelCovRelativeCpuWeight;
  if (params.parallelCovarRelativeCPUWeight.value()) {
    parallelCovRelativeCpuWeight = *params.parallelCovarRelativeCPUWeight.value();
  } else {
    parallelCovRelativeCpuWeight = std::vector<double>(nComponents,
                      1.0 / static_cast<double>(nComponents));
  }

  std::vector<size_t> ntasksPerComponent(nComponents, 0);
  std::vector<size_t> globalTaskOffsetPerComponent(nComponents+1, 0);
  if (parallelHybrid_) {
    oops::Log::info() << "Info     : Creating Hybrid block in parallel" << std::endl;
    // checks
    ASSERT(nComponents == parallelCovRelativeCpuWeight.size());

    // need to check to ensure that the total sum of PEs over components is consistent
    // with the MPI size on the default communicator and that each component
    // has a minimum MPI size of 1.
    for (size_t component = 0; component < nComponents; ++component) {
      ntasksPerComponent[component] =
        std::round(parallelCovRelativeCpuWeight[component] * ntasks);
      ASSERT(ntasksPerComponent[component] > 0);
    }
    int discrepencyPE =
      std::accumulate(ntasksPerComponent.begin(), ntasksPerComponent.end(), 0) - ntasks;

    for (size_t component = 0; component < nComponents && discrepencyPE != 0; ++component) {
      if (discrepencyPE > 0 && ntasksPerComponent[component] >= 2) {
        ntasksPerComponent[component] -= 1;
        discrepencyPE -= 1;
      } else if (discrepencyPE < 0) {
        ntasksPerComponent[component] += 1;
        discrepencyPE += 1;
      }
    }

    ASSERT(std::accumulate(ntasksPerComponent.begin(),
                           ntasksPerComponent.end(), 0) - ntasks == 0);

    for (size_t component = 1; component < nComponents; ++component) {
      globalTaskOffsetPerComponent[component] =
        globalTaskOffsetPerComponent[component-1] + ntasksPerComponent[component-1];
    }
    globalTaskOffsetPerComponent[nComponents] = ntasks;

    // Report the full task split. A component whose task count does not match
    // the decomposition its central block expects (MGBF requires nxm*nym) will
    // fail later inside that block, so make the numbers visible up front.
    oops::Log::info() << "Info     : Hybrid task split over " << nComponents
                      << " components:";
    for (size_t component = 0; component < nComponents; ++component) {
      oops::Log::info() << " " << ntasksPerComponent[component];
    }
    oops::Log::info() << " (total " << ntasks << ")" << std::endl;

    // Shared ensemble: read once here, on the parent communicator, while every
    // rank can still see every member. This must happen before the communicator
    // split below, and in particular before setCommDefault: readEnsemble's
    // other-geometry branch captures eckit::mpi::comm() as its geometry
    // communicator, which the split is about to change underneath it.
    std::unique_ptr<oops::FieldSets> sharedPerts;
    std::vector<std::vector<size_t>> sharedMemberMap;
    if (params.sharedEnsemble.value()) {
      const auto & sharedParams = *params.sharedEnsemble.value();

      // Validate the whole mapping before any I/O. Every rank holds the full
      // component configuration, so this costs nothing and turns a silent
      // mis-assignment into a named error.
      sharedMemberMap = sharedEnsembleMemberMap(params.components.value(), sharedParams);

      eckit::LocalConfiguration sharedConf;
      sharedParams.serialize(sharedConf);

      oops::Log::info() << "Info     : Reading shared ensemble on the parent communicator"
                        << std::endl;
      sharedPerts = std::make_unique<oops::FieldSets>(
        readEnsemble(geom, currentOuterVars,
                     fset4dXb.times(), fset4dXb.commTime(), fset4dXb.commEns(),
                     sharedConf));

      size_t nmembers = 0;
      for (const auto & members : sharedMemberMap) nmembers += members.size();
      if (sharedPerts->ens_size() != nmembers) {
        throw eckit::UserError("The shared ensemble holds " +
                               std::to_string(sharedPerts->ens_size()) + " members but "
                               "`shared ensemble members` accounts for " +
                               std::to_string(nmembers) + ".", Here());
      }
      oops::Log::info() << "Info     : Shared ensemble: " << nmembers << " members over "
                        << nComponents << " components, common mean subtracted, "
                        << "normalization 1/" << sharedParams.denominator.value()
                        << " carried by the component weights" << std::endl;
    }

    const eckit::mpi::Comm & initialDefaultComm = eckit::mpi::comm();
    ASSERT(initialDefaultComm.name() == defaultSpaceComm.name());

    // We split the space communicators only, the time parallelization is untouched
    const size_t myTask = defaultSpaceComm.rank();

    // Set myComponent_  tasksPerComponent
    size_t tasksPerComponent = 0;
    for (size_t component = 0; component < nComponents; ++component) {
      if ((myTask >= globalTaskOffsetPerComponent[component]) &&
          (myTask < globalTaskOffsetPerComponent[component+1])) {
        myComponent_ = component;
        tasksPerComponent = ntasksPerComponent[component];
      }
    }

    oops::Log::info() << "Info     : Creating component " << myComponent_ + 1
                      << "/" << nComponents
                      << " of Hybrid block using " << tasksPerComponent
                      << " MPI tasks." << std::endl;

    // Create communicators for same component, for communications in space
    const auto spaceCommName = ("comm_space_" + std::to_string(myComponent_));
    if (eckit::mpi::hasComm(spaceCommName.c_str())) {
      eckit::mpi::deleteComm(spaceCommName.c_str());
    }
    const auto & localSpaceComm = defaultSpaceComm.split(myComponent_, spaceCommName.c_str());

    // Set up default MPI communicator for atlas
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Create block geometry (needed for ensemble reading and local geometries)
    if (params.hybridGeometry.value() == boost::none) {
      throw eckit::UserError("Parallel hybrid block requires geometry key", Here());
    }
    const auto geomConf = *params.hybridGeometry.value();
    // The hybrid Geometry is stored as a class member to ensure it doesn't go
    // out of scope after construction, as it is directly used (not copied) by
    // the hybrid Block Chains.
    localHybridGeom_ =
        std::make_unique<oops::Geometry<MODEL>>(geomConf, localSpaceComm, geom.timeComm());
    localHybridFs_.reset(new atlas::FunctionSpace(localHybridGeom_->functionSpace()));
    // Copy and redistribute the background and first guess
    oops::FieldSet4D localFset4dXb(fset4dXb.times(), fset4dXb.commTime(),
                                   localSpaceComm, *localHybridFs_, fset4dXb.variables());
    oops::FieldSet4D localFset4dFg(fset4dFg.times(), fset4dFg.commTime(),
                                   localSpaceComm, *localHybridFs_, fset4dFg.variables());

    for (size_t jtime = 0; jtime < fset4dXb.size(); jtime++) {
      util::redistributeToSubcommunicator(redistributionMethod_,
                                          fset4dXb[jtime].fieldSet(),
                                          localFset4dXb[jtime].fieldSet(),
                                          *localHybridFs_);
      util::redistributeToSubcommunicator(redistributionMethod_,
                                          fset4dFg[jtime].fieldSet(),
                                          localFset4dFg[jtime].fieldSet(),
                                          *localHybridFs_);
    }
    defaultSpaceComm.barrier();

    // Slice the shared ensemble into the members belonging to this component.
    //
    // The loop runs over every global member on every parent rank, so the
    // sequence of collectives is identical everywhere however the members are
    // distributed; only the decision to keep a member is local. That is what
    // allows components to hold different numbers of members.
    SaberBlockChainContext context;
    if (sharedPerts) {
      const std::vector<size_t> & mine = sharedMemberMap[myComponent_];

      // Original member identifiers are preserved so that diagnostics stay
      // traceable to the ensemble numbering the user configured.
      std::vector<int> memberIds;
      memberIds.reserve(mine.size());
      for (const size_t member : mine) memberIds.push_back(static_cast<int>(member));

      auto myPerts = std::make_unique<oops::FieldSets>(
        fset4dXb.times(), fset4dXb.commTime(), memberIds, fset4dXb.commEns());

      for (size_t member = 0; member < sharedPerts->ens_size(); ++member) {
        const auto itMine = std::find(mine.begin(), mine.end(), member + 1);
        const bool keep = (itMine != mine.end());
        const size_t localIndex = keep
          ? static_cast<size_t>(std::distance(mine.begin(), itMine)) : 0;

        // Allocated per member: the block chain keeps a reference to these
        // fields, so reusing one buffer would alias every member together.
        oops::FieldSet4D localMember(fset4dXb.times(), fset4dXb.commTime(),
                                     localSpaceComm, *localHybridFs_, currentOuterVars);

        for (size_t jtime = 0; jtime < fset4dXb.size(); jtime++) {
          util::redistributeToSubcommunicator(redistributionMethod_,
                                              (*sharedPerts)(jtime, member).fieldSet(),
                                              localMember[jtime].fieldSet(),
                                              *localHybridFs_);
          if (keep) myPerts->emplace_back(jtime, localIndex, localMember[jtime]);
        }
      }

      // The global container is no longer needed; release it before the
      // component block chain is built.
      sharedPerts.reset();
      defaultSpaceComm.barrier();

      oops::Log::info() << "Info     : Component " << myComponent_ + 1 << " received "
                        << mine.size() << " shared ensemble members" << std::endl;
      context = SaberBlockChainContext(std::move(myPerts));
    }

    const auto & cmpParams = params.components.value()[myComponent_];

    // Initialize component outer variables
    const oops::Variables cmpOuterVars(currentOuterVars);

    // Set weight
    const auto & weightParams = cmpParams.weight.value();
    // Scalar weight
    hybridScalarWeightSqrt_.push_back(std::sqrt(weightParams.value.value()));
    // File-base weight
    oops::FieldSet3D fsetWeight(localFset4dXb[0].validTime(), localSpaceComm);
    if (weightParams.file.value()) {
      // File-base weight
      readHybridWeight(*localHybridGeom_,
                       cmpOuterVars,
                       localFset4dXb[0].validTime(),
                       *weightParams.file.value(),
                       fsetWeight);
      fsetWeight.sqrt();
    }
    hybridFieldWeightSqrt_.push_back(fsetWeight);

    // Set covariance parameters
    const auto & cmpCovParams = cmpParams.covariance.value();

    // Merge component configuration with full configuration base (order of arguments matters!)
    const eckit::LocalConfiguration cmpMergedConf =
      util::mergeConfigs(cmpCovParams.toConfiguration(), paramsBase.toConfiguration());

    // Add block chain
    hybridBlockChain_.push_back(
        SaberBlockChainFactory<MODEL>::create
         (*localHybridGeom_,
          cmpOuterVars,
          localFset4dXb,
          localFset4dFg,
          cmpMergedConf,
          std::move(context)));

    ASSERT(hybridBlockChain_.size() > 0);

    // Restore previous default MPI communicator for atlas
    eckit::mpi::setCommDefault(defaultSpaceComm.name().c_str());
  } else {
    oops::Log::info() << "Info     : Creating Hybrid block serially" << std::endl;
    // Create block geometry
    const oops::Geometry<MODEL> * hybridGeom = &geom;
    if (params.hybridGeometry.value()) {
      hybridGeom = new oops::Geometry<MODEL>(
        *params.hybridGeometry.value(),
        geom.getComm());
    }
    for (const auto & cmpParams : params.components.value()) {
      // Initialize component outer variables
      const oops::Variables cmpOuterVars(currentOuterVars);

      // Set weight
      const auto & weightParams = cmpParams.weight.value();
      // Scalar weight
      hybridScalarWeightSqrt_.push_back(std::sqrt(weightParams.value));
      // File-base weight
      oops::FieldSet3D fsetWeight(fset4dXb[0].validTime(), geom.getComm());
      if (weightParams.file.value()) {
        // File-base weight
        readHybridWeight(*hybridGeom,
                         cmpOuterVars,
                         fset4dXb[0].validTime(),
                         *weightParams.file.value(),
                         fsetWeight);
        fsetWeight.sqrt();
      }
      hybridFieldWeightSqrt_.push_back(fsetWeight);

      // Set covariance parameters
      const auto & cmpCovParams = cmpParams.covariance.value();

      // Merge component configuration with full configuration base (order of arguments matters!)
      const eckit::LocalConfiguration cmpMergedConf =
        util::mergeConfigs(cmpCovParams.toConfiguration(), paramsBase.toConfiguration());

      // Add block chain
      hybridBlockChain_.push_back
          (SaberBlockChainFactory<MODEL>::create
           (*hybridGeom,
            cmpOuterVars,
            fset4dXb,
            fset4dFg,
            cmpMergedConf));
    }
    ASSERT(hybridBlockChain_.size() > 0);
  }

  oops::Log::trace() << "SaberHybridBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberHybridBlockChain<MODEL>::randomize(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "SaberHybridBlockChain::randomize starting" << std::endl;
  util::Timer timer("SaberHybridBlockChain", "randomize");

  // Initialize FieldSet4D
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    if (outerBlockChain_) {
      fset4d[jtime].init(outerBlockChain_->innerGeometryData().functionSpace(),
        outerBlockChain_->innerVars());
    } else {
      fset4d[jtime].init(outerFunctionSpace_, outerVariables_);
    }
  }
  fset4d.zero();

  if (parallelHybrid_) {
    // Run components of the central block in parallel
    oops::Log::debug() << "Parallel execution of doRandomize in Hybrid" << std::endl;
    oops::Log::debug() << "Running Hybrid component " << myComponent_ + 1 << std::endl;
    ASSERT(hybridBlockChain_.size() == 1);

    // global communicator and functionSpace
    const auto & defaultSpaceComm = fset4d[0].commGeom();

    // check global communicator is the default one for atlas MPI
    ASSERT(eckit::mpi::comm().name() == defaultSpaceComm.name());

    // subcommunicator within this component
    const auto spaceCommName = "comm_space_" + std::to_string(myComponent_);
    const auto & localSpaceComm = eckit::mpi::comm(spaceCommName.c_str());

    // Set up atlas MPI
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Create temporary FieldSet on subcommunicator
    oops::FieldSet4D fset4dCmp(fset4d.times(), fset4d.commTime(), localSpaceComm);

    hybridBlockChain_[0]->randomize(fset4dCmp);

    // Weight square-root multiplication
    if (hybridScalarWeightSqrt_[0] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[0];
    }
    if (!hybridFieldWeightSqrt_[0].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[0];
    }

    // Add components
    defaultSpaceComm.barrier();

    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      // Redistribute to global communicator and sum
      util::gatherAndSumFromSubcommunicator(redistributionMethod_,
                                            fset4dCmp[jtime].fieldSet(),
                                            fset4d[jtime].fieldSet(),
                                            *localHybridFs_,
                                            *globalHybridFs_);
    }

    // Restore atlas MPI to previous
    eckit::mpi::setCommDefault(defaultSpaceComm.name().c_str());

    fset4d += fset4dCmp;
  } else {
    // Loop over components for the central block
    for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
      // Randomize covariance
      oops::FieldSet4D fset4dCmp(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());
      hybridBlockChain_[jj]->randomize(fset4dCmp);

      // Weight square-root multiplication
      if (hybridScalarWeightSqrt_[jj] != 1.0) {
        // Scalar weight
        fset4dCmp *= hybridScalarWeightSqrt_[jj];
      }
      if (!hybridFieldWeightSqrt_[jj].empty()) {
        // File-based weight
        fset4dCmp *= hybridFieldWeightSqrt_[jj];
      }

      // Add component
      fset4d += fset4dCmp;
    }
  }

  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4d);

  oops::Log::trace() << "SaberHybridBlockChain::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberHybridBlockChain<MODEL>::multiply(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "SaberHybridBlockChain::multiply starting" << std::endl;
  util::Timer timer("SaberHybridBlockChain", "multiply");

  // Apply outer blocks adjoint
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocksAD(fset4d);

  // Initialize sum to zero
  oops::FieldSet4D fset4dSum = oops::copyFieldSet4D(fset4d);
  fset4dSum.zero();

  // Loop over B components
  if (parallelHybrid_) {
    oops::Log::debug() << "Parallel execution of Hybrid::multiply, component "
                       << myComponent_ + 1 << std::endl;
    ASSERT(hybridBlockChain_.size() == 1);
    ASSERT(hybridScalarWeightSqrt_.size() == 1);
    ASSERT(hybridFieldWeightSqrt_.size() == 1);

    // Global communicator
    const auto & defaultSpaceComm = fset4d[0].commGeom();
    ASSERT(defaultSpaceComm.name() == eckit::mpi::comm().name());

    // Subcommunicator within component
    const std::string spaceCommName = "comm_space_" + std::to_string(myComponent_);
    const auto & localSpaceComm = eckit::mpi::comm(spaceCommName.c_str());

    // Create temporary FieldSet copy on communicator of this component
    oops::FieldSet4D fset4dCmp(fset4d.times(), fset4d.commTime(), localSpaceComm);
    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      util::redistributeToSubcommunicator(redistributionMethod_,
                                          fset4d[jtime].fieldSet(),
                                          fset4dCmp[jtime].fieldSet(),
                                          *localHybridFs_);
    }

    // Set up atlas MPI
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Apply weight
    if (hybridScalarWeightSqrt_[0] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[0];
    }
    if (!hybridFieldWeightSqrt_[0].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[0];
    }

    // Apply covariance
    hybridBlockChain_[0]->multiply(fset4dCmp);

    // Apply weight
    if (hybridScalarWeightSqrt_[0] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[0];
    }
    if (!hybridFieldWeightSqrt_[0].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[0];
    }

    // Wait for all components to have finished multiplying
    defaultSpaceComm.barrier();

    // Gather and sum data across components
    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      util::gatherAndSumFromSubcommunicator(redistributionMethod_,
                                            fset4dCmp[jtime].fieldSet(),
                                            fset4dSum[jtime].fieldSet(),
                                            *localHybridFs_,
                                            *globalHybridFs_);
    }

    // Set back default MPI communicator
    eckit::mpi::setCommDefault(defaultSpaceComm.name().c_str());

  } else {
    if (hybridBlockChain_.size() > 1) {
        oops::Log::debug() << "Serial execution of Hybrid::multiply" << std::endl;
    }
    for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
      // Create temporary FieldSet
      oops::FieldSet4D fset4dCmp = oops::copyFieldSet4D(fset4d);

      // Apply weight
      if (hybridScalarWeightSqrt_[jj] != 1.0) {
        // Scalar weight
        fset4dCmp *= hybridScalarWeightSqrt_[jj];
      }
      if (!hybridFieldWeightSqrt_[jj].empty()) {
        // File-based weight
        fset4dCmp *= hybridFieldWeightSqrt_[jj];
      }

      // Apply covariance
      hybridBlockChain_[jj]->multiply(fset4dCmp);

      // Apply weight
      if (hybridScalarWeightSqrt_[jj] != 1.0) {
        // Scalar weight
        fset4dCmp *= hybridScalarWeightSqrt_[jj];
      }
      if (!hybridFieldWeightSqrt_[jj].empty()) {
        // File-based weight
        fset4dCmp *= hybridFieldWeightSqrt_[jj];
      }

      // Add component
      fset4dSum += fset4dCmp;
    }
  }

  // Apply outer blocks forward
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4dSum);

  fset4d.deepCopy(fset4dSum);

  oops::Log::trace() << "SaberHybridBlockChain::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
