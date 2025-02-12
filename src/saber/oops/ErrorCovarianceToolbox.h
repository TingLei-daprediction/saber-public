/*
 * (C) Copyright 2021-2023 UCAR
 * (C) Copyright 2023 Meteorologisk Institutt
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <math.h>
#include <netcdf.h>
#include <omp.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "atlas/functionspace.h"
#include "atlas/util/Earth.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/StateWriter.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/ConfigHelpers.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"
#include "saber/util/HorizontalProfiles.h"

#include "oops/base/PostProcessor.h"
#include "oops/base/StructuredGridPostProcessor.h"
#include  "oops/base/StructuredGridWriter.h"

namespace saber {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ErrorCovarianceToolbox application.
template <typename MODEL> class ErrorCovarianceToolboxParameters :
  public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceToolboxParameters, oops::ApplicationParameters)

 public:
  /// Geometry parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};

  /// Background parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> background{"background", this};

  /// Background error covariance model.
  oops::RequiredParameter<eckit::LocalConfiguration> backgroundError{"background error", this};

  /// Geometry parameters.
  oops::Parameter<bool> parallel{"parallel subwindows", true, this};

  /// Outer variables parameters
  oops::OptionalParameter<oops::Variables> incrementVars{"increment variables", this};

  /// Dirac location/variables parameters.
  oops::OptionalParameter<eckit::LocalConfiguration> dirac{"dirac", this};

  /// Diagnostic location/variables parameters.
  oops::OptionalParameter<eckit::LocalConfiguration> diagnostic{"diagnostic points", this};

  /// Where to write the output(s) of Dirac tests
  oops::OptionalParameter<eckit::LocalConfiguration> outputDirac{"output dirac", this};

  /// Whether and where to write the perturbations generated by the randomization.
  oops::OptionalParameter<eckit::LocalConfiguration> outputPerturbations{"output perturbations",
    this};

  /// Whether and where to write the states generated by the randomization.
  oops::OptionalParameter<eckit::LocalConfiguration> outputStates{"output states", this};

  /// Where to write the output of randomized variance.
  oops::OptionalParameter<eckit::LocalConfiguration> outputVariance{"output variance", this};

  /// Whether and how to compute unidimensional covariance profiles for isotropic cases
  oops::OptionalParameter<eckit::LocalConfiguration> covarianceProfile{
                                    "covariance profile", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ErrorCovarianceToolbox : public oops::Application {
  typedef oops::ModelSpaceCovarianceBase<MODEL>           CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                  CovarianceFactory_;
  typedef ModelSpaceCovarianceParametersBase<MODEL>       CovarianceParametersBase_;
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef oops::Increment4D<MODEL>                        Increment4D_;
  typedef oops::State<MODEL>                              State_;
  typedef oops::State4D<MODEL>                            State4D_;
  typedef oops::Localization<MODEL>                       Localization_;
  typedef ErrorCovarianceToolboxParameters<MODEL>         ErrorCovarianceToolboxParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit ErrorCovarianceToolbox(const eckit::mpi::Comm & comm = eckit::mpi::comm()) :
    Application(comm) {
    oops::instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~ErrorCovarianceToolbox() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
    // Deserialize parameters
    ErrorCovarianceToolboxParameters_ params;
    params.deserialize(fullConfig);

    // Define number of subwindows
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "background");
    size_t nsubwin = 1;
    if (backgroundConfig.has("states")) {
      std::vector<eckit::LocalConfiguration> confs;
      backgroundConfig.get("states", confs);
      nsubwin = confs.size();
    }

    // Define space and time communicators
    const eckit::mpi::Comm * commSpace = &this->getComm();
    const eckit::mpi::Comm * commTime = &oops::mpi::myself();
    if (nsubwin > 1) {
      // Define sub-windows
      const size_t ntasks = this->getComm().size();
      size_t mysubwin = 0;
      size_t nsublocal = nsubwin;
      if (params.parallel && (ntasks % nsubwin == 0)) {
        nsublocal = 1;
        mysubwin = this->getComm().rank() / (ntasks / nsubwin);
        ASSERT(mysubwin < nsubwin);
      } else if (params.parallel) {
        oops::Log::warning() << "Parallel time subwindows specified in yaml "
                             << "but number of tasks is not divisible by "
                             << "the number of subwindows, ignoring." << std::endl;
      }

      // Create a communicator for same sub-window, to be used for communications in space
      const std::string sgeom = "comm_geom_" + std::to_string(mysubwin);
      char const *geomName = sgeom.c_str();
      commSpace = &this->getComm().split(mysubwin, geomName);

      // Create a communicator for same local area, to be used for communications in time
      const size_t myarea = commSpace->rank();
      const std::string stime = "comm_time_" + std::to_string(myarea);
      char const *timeName = stime.c_str();
      commTime = &this->getComm().split(myarea, timeName);
      ASSERT(commTime->size() == (nsubwin / nsublocal));
    }

    // Get number of MPI tasks and OpenMP threads
    size_t ntasks = commSpace->size();
    size_t nthreads = 1;
#ifdef _OPENMP
    # pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }
#endif

    // Replace patterns in full configuration and deserialize parameters
    eckit::LocalConfiguration fullConfigUpdated(fullConfig);
    util::seekAndReplace(fullConfigUpdated, "_MPI_", std::to_string(ntasks));
    util::seekAndReplace(fullConfigUpdated, "_OMP_", std::to_string(nthreads));
    params.deserialize(fullConfigUpdated);

    // Setup geometry
    const Geometry_ geom(params.geometry, *commSpace, *commTime);

    // Setup background
    const State4D_ xx(geom, params.background.value(), *commTime);

    // Setup variables
    oops::Variables tmpVars = xx.variables();
    if (params.incrementVars.value() != boost::none) {
      tmpVars = params.incrementVars.value().value();
    }
    const oops::Variables vars = tmpVars;

    // Setup time
    util::DateTime time = xx[0].validTime();

    const eckit::LocalConfiguration covarConf(fullConfigUpdated, "background error");

    // Dirac test
    const auto & diracParams = params.dirac.value();
    if (diracParams != boost::none) {
      // Setup Dirac field
      Increment4D_ dxi(geom, vars, xx.times(), *commTime);
      dxi.dirac(*diracParams);
      oops::Log::test() << "Input Dirac increment:" << dxi << std::endl;

      // Test configuration
      eckit::LocalConfiguration testConf;

      // Add dirac configuration
      testConf.set("dirac", *diracParams);

      // Add diagnostic print configuration
      const auto & diagnostic = params.diagnostic.value();
      if (diagnostic != boost::none) {
        testConf.set("diagnostic points", *diagnostic);
      }

      // Add output Dirac configuration
      eckit::LocalConfiguration outputDiracUpdated = params.outputDirac.value().value();
      setMPI(outputDiracUpdated, ntasks);
      testConf.set("output dirac", outputDiracUpdated);

      // Add covariance profile configuration
      const auto & profileConfig = params.covarianceProfile.value();
      if (profileConfig != boost::none) {
        testConf.set("covariance profile", *profileConfig);
      }

      // Apply B matrix components recursively
      std::string id;
      dirac(covarConf, testConf, id, geom, vars, xx, dxi);
    }

    // Background error covariance parameters
    CovarianceParametersBase_ covarParams;
    covarParams.deserialize(covarConf);
    const auto & randomizationSize = covarParams.randomizationSize.value();
    if ((diracParams == boost::none) || (randomizationSize != boost::none)) {
      // Background error covariance training
      std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                            geom, vars, covarConf, xx, xx));

      // Randomization
      randomization(params, geom, vars, xx, Bmat, ntasks);
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ErrorCovarianceToolbox<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
// The passed geometry should be consistent with the passed increment
  void print_value_at_positions(const eckit::LocalConfiguration & diagConf,
                                const Geometry_ & geom,
                                const Increment4D_ & data) const {
    oops::Log::trace() << appname() << "::print_value_at_position starting" << std::endl;

    // Create diagnostic field
    Increment4D_ diagPoints(data);
    diagPoints.dirac(diagConf);

    // Get diagnostic values
    for (size_t jj = 0; jj < data.size(); ++jj) {
      util::printDiagValues(diagPoints.commTime(),
                            geom.getComm(),
                            geom.functionSpace(),
                            data[jj].fieldSet().fieldSet(),
                            diagPoints[jj].fieldSet().fieldSet());
    }

    oops::Log::trace() << appname() << "::print_value_at_position done" << std::endl;
  }
// -----------------------------------------------------------------------------
// The passed geometry should be consistent with the passed increment
  void extract_1d_covariances(const eckit::LocalConfiguration & diracConf,
                              const eckit::LocalConfiguration & profileConf,
                              const Geometry_ & geom,
                              const Increment4D_ & data) const {
  oops::Log::trace() << appname() << "::extract_1d_covariances starting" << std::endl;

  if (data.size() > 1) {
    throw eckit::NotImplemented("Not implemented for 4D covariances", Here());
  }

  // Create dirac field
  Increment4D_ diracPoints(data);
  diracPoints.dirac(diracConf);

  // Define maximum length of horizontal profile
  double maxLength = std::numeric_limits<double>::infinity();
  if (profileConf.has("maximum distance")) {
    oops::Log::info() << "Info     : Reading user-provided maximum distance:" << std::endl;
    maxLength = profileConf.getDouble("maximum distance");
  } else {
    // If no maximum length input is given, try to compute a default value from the grid
    const auto & fspace = geom.functionSpace();
    if (fspace.type() == "StructuredColumns") {
      const atlas::functionspace::StructuredColumns fs(fspace);
      if (fs.grid().name().compare(0, 1, "F") == 0) {
        oops::Log::info() << "Info     : maximum distance computed as twice the cell "
                             "length at the Equator:" << std::endl;
        const int n = std::stoi(fs.grid().name().substr(1, std::string::npos));
        maxLength = atlas::util::Earth().radius() * M_PI / n;
      }
    } else if (fspace.type() == "NodeColumns") {
      const atlas::functionspace::NodeColumns fs(fspace);
      if (fs.mesh().grid().name().compare(0, 7, "CS-LFR-") == 0) {
        oops::Log::info() << "Info     : maximum distance computed as twice the cell "
                             "length at the Equator:" << std::endl;
        const int n = std::stoi(fs.mesh().grid().name().substr(7, std::string::npos));
        maxLength = atlas::util::Earth().radius() * M_PI / n;
      }
    }
  }
  oops::Log::info() << "Info     : maximum distance is " << maxLength << " m." << std::endl;

  // Boolean flag to remove duplicate points (is isotropy for instance)
  const bool removeDuplicates = profileConf.getBool("remove duplicate distances",
                                                    false);

  // Get values as a function of separation distance
  auto[distances, covariances, lons, lats, levs, fieldIndexes] =
      util::sortBySeparationDistance(geom.getComm(),
                                     geom.functionSpace(),
                                     data[0].fieldSet().fieldSet(),
                                     diracPoints[0].fieldSet().fieldSet(),
                                     maxLength,
                                     removeDuplicates);

  // Write to file or to test Log
  const auto & names = data[0].fieldSet().fieldSet().field_names();
  if (profileConf.has("output filepath")) {
    // Write to file
    const auto outputPath = profileConf.getString("output filepath");
    util::write_1d_covariances(geom.getComm(),
                               distances,
                               covariances,
                               lons,
                               lats,
                               levs,
                               fieldIndexes,
                               names,
                               outputPath);
  }
  // Output first 10 values of first 10  profiles to test log
  const int numProfiles = std::min(10, static_cast<int>(distances.size()));
  for (int i=0; i < numProfiles; i++) {
    oops::Log::test() << "Covariance profile for variable " << names[fieldIndexes[i]]
                      << ", at level " << levs[i] << ": " << std::endl;
    const int numValues = std::min(10, static_cast<int>(distances[i].size()));
    const std::vector<double> subDist(distances[i].begin(), distances[i].begin() + numValues);
    const std::vector<double> subCovs(covariances[i].begin(), covariances[i].begin() + numValues);
    oops::Log::test() << "Separation distance: " << subDist << std::endl;
    oops::Log::test() << "Covariance: " << subCovs << std::endl;
  }

  oops::Log::trace() << appname() << "::extract_1d_covariances done" << std::endl;
}
// -----------------------------------------------------------------------------
// The passed geometry/variables should be consistent with the passed increment
  void dirac(const eckit::LocalConfiguration & covarConf,
             const eckit::LocalConfiguration & testConf,
             std::string & id,
             const Geometry_ & geom,
             const oops::Variables & vars,
             const State4D_ & xx,
             const Increment4D_ & dxi) const {
    // Define output increment
    // tothinkdo
    oops::Log::trace() <<  "dirac starting" << std::endl;
    Increment4D_ dxo(dxi, false);

    // Covariance
    oops::Log::trace() <<  "dirac Bmat being created" << std::endl;
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                          geom, vars, covarConf, xx, xx));

    // Multiply
    oops::Log::trace() <<  "dirac Bmat 's multiply to be invoked" << std::endl;
    Bmat->multiply(dxi, dxo);

    // Update ID
    if (id != "") id.append("_");
    id.append(Bmat->covarianceModel());

    if (testConf.has("diagnostic points")) {
      oops::Log::test() << "Covariance(" << id << ") diagnostics:" << std::endl;

      // Print variances
      oops::Log::test() << "- Variances at Dirac points:" << std::endl;
      print_value_at_positions(testConf.getSubConfiguration("dirac"), geom, dxo);

      // Print covariances
      oops::Log::test() << "- Covariances at diagnostic points:" << std::endl;
      print_value_at_positions(testConf.getSubConfiguration("diagnostic points"), geom, dxo);
    }

    if (testConf.has("covariance profile")) {
      oops::Log::test() << "Extracting covariances as a function of separation distance"
                        << std::endl;
      eckit::LocalConfiguration covProfileConf(testConf.getSubConfiguration("covariance profile"));

      // Seek and replace %id% with id, recursively
      util::seekAndReplace(covProfileConf, "%id%", id);

      extract_1d_covariances(testConf.getSubConfiguration("dirac"),
                             covProfileConf,
                             geom,
                             dxo);
    }

    // Copy configuration
    eckit::LocalConfiguration outputBConf(testConf.getSubConfiguration("output dirac"));

    // Seek and replace %id% with id, recursively
    util::seekAndReplace(outputBConf, "%id%", id);
     for (auto& fld:dxo[0].fieldSet() ) {
    oops::Log::trace()<<"thinkdeb dxo 0 input fld is "<<fld<<std::endl;
//clt    fld.dump(std::cout);
   };

    // Write output increment
    dxo[0].write(outputBConf);

     for (auto& fld:dxo[0].fieldSet() ) {
    oops::Log::trace()<<"thinkdeb dxo 1 input fld is "<<fld<<std::endl;
//clt    fld.dump(std::cout);
   };
//cltb1 for output on latlon grids
    if (outputBConf.has("analysis to structured grid")) {
       oops::PostProcessor<State_> post;
       const eckit::LocalConfiguration anLatlonConf(outputBConf, "analysis to structured grid");
       oops::Log::trace() << "thinkdeb anaLatlonConf: " << anLatlonConf << std::endl;
       oops::StructuredGridWriter<MODEL> latlonwriter(anLatlonConf,dxo.geometry());
       oops::Log::trace() << "thinkdeb before latlonwrite.interpolatAndWriter " <<  std::endl;
       latlonwriter.interpolateAndWrite(dxo[0]); //clttodo
       
//       post.enrollProcessor(new oops::StructuredGridPostProcessor<MODEL, State_>(
 //           anLatlonConf, dxo.geometry() ));
        };

       
 
//end of b1
    oops::Log::test() << "Covariance(" << id << ") * Increment:" << dxo << std::endl;

    // Look for hybrid or ensemble covariance models
    const std::string covarianceModel(covarConf.getString("covariance model"));
    if (covarianceModel == "hybrid") {
      std::vector<eckit::LocalConfiguration> confs;
      covarConf.get("components", confs);
      size_t componentIndex(1);
      for (const auto & conf : confs) {
        std::string idC(id + std::to_string(componentIndex));
        const eckit::LocalConfiguration componentConfig(conf, "covariance");
        dirac(componentConfig, testConf, idC, geom, vars, xx, dxi);
        ++componentIndex;
      }
    }
    if (covarianceModel == "SABER") {
      const std::string saberCentralBlockName =
        covarConf.getString("saber central block.saber block name");
      if (saberCentralBlockName == "Hybrid") {
        // Check for outer blocks (can't pass the correct geometry/variables in that case)
        if (!covarConf.has("saber outer blocks")) {
          std::vector<eckit::LocalConfiguration> confs;
          covarConf.get("saber central block.components", confs);
          size_t componentIndex(1);
          for (const auto & conf : confs) {
            std::string idC(id + std::to_string(componentIndex));
            eckit::LocalConfiguration componentConfig(conf, "covariance");
            componentConfig.set("covariance model", "SABER");
            if (covarConf.has("adjoint test")) {
              componentConfig.set("adjoint test", covarConf.getBool("adjoint test"));
            }
            if (covarConf.has("inverse test")) {
              componentConfig.set("inverse test", covarConf.getBool("inverse test"));
            }
            if (covarConf.has("square-root test")) {
              componentConfig.set("square-root test", covarConf.getBool("square-root test"));
            }
            dirac(componentConfig, testConf, idC, geom, vars, xx, dxi);
            ++componentIndex;
          }
        }
      }
    }
    if (covarianceModel == "ensemble" && covarConf.has("localization")) {
      // Localization configuration
      eckit::LocalConfiguration locConfig(covarConf.getSubConfiguration("localization"));
      locConfig.set("date", xx[0].validTime().toString());

      // Define output increment
      Increment4D_ dxo(dxi);

      // Setup localization
      Localization_ Lmat(geom, vars, locConfig);

      // Apply localization
      Lmat.multiply(dxo);

      // Update ID
      std::string idL(id);
      idL.append("_localization");

      // Print localization
      oops::Log::test() << "Localization(" << idL << ") diagnostics:" << std::endl;
      oops::Log::test() << "- Localization at zero separation:" << std::endl;

      print_value_at_positions(testConf.getSubConfiguration("dirac"), geom, dxo);
      if (testConf.has("diagnostic points")) {
        oops::Log::test() << "- Localization at diagnostic points:" << std::endl;
        print_value_at_positions(testConf.getSubConfiguration("diagnostic points"), geom, dxo);
      }

      // Copy configuration
      eckit::LocalConfiguration outputLConf(testConf.getSubConfiguration("output dirac"));

      // Seek and replace %id% with id, recursively
      util::seekAndReplace(outputLConf, "%id%", idL);

      // Write output increment
      dxo[0].write(outputLConf);
      oops::Log::test() << "Localization(" << id << ") * Increment:" << dxo << std::endl;
    }
  }
// -----------------------------------------------------------------------------
  void randomization(const ErrorCovarianceToolboxParameters_ & params,
                     const Geometry_ & geom,
                     const oops::Variables & vars,
                     const State4D_ & xx,
                     const std::unique_ptr<CovarianceBase_> & Bmat,
                     const size_t & ntasks) const {
    if (Bmat->randomizationSize() > 0) {
      oops::Log::info() << "Info     : " << std::endl;
      oops::Log::info() << "Info     : Generate perturbations:" << std::endl;
      oops::Log::info() << "Info     : -----------------------" << std::endl;

      // Create increments
      Increment4D_ dx(geom, vars, xx.times(), xx.commTime());
      Increment4D_ dxsq(geom, vars, xx.times(), xx.commTime());
      Increment4D_ variance(geom, vars, xx.times(), xx.commTime());

      // Initialize variance
      variance.zero();

      // Create empty ensemble
      std::vector<Increment_> ens;

      // Output options
      const auto & outputPerturbations = params.outputPerturbations.value();
      const auto & outputStates = params.outputStates.value();
      const auto & outputVariance = params.outputVariance.value();

      for (size_t jm = 0; jm < Bmat->randomizationSize(); ++jm) {
        // Generate member
        oops::Log::info() << "Info     : Member " << jm << std::endl;
        Bmat->randomize(dx);

        if ((outputPerturbations != boost::none) || (outputStates != boost::none)) {
          // Save member
          ens.push_back(dx[0]);
        }

        // Square perturbation
        dxsq = dx;
        dxsq.schur_product_with(dx);

        // Update variance
        variance += dxsq;
      }
      oops::Log::info() << "Info     : " << std::endl;

      if ((outputPerturbations != boost::none) || (outputStates != boost::none)) {
        oops::Log::info() << "Info     : Write states and/or perturbations:" << std::endl;
        oops::Log::info() << "Info     : ----------------------------------" << std::endl;
        oops::Log::info() << "Info     : " << std::endl;
        for (size_t jm = 0; jm < Bmat->randomizationSize(); ++jm) {
          oops::Log::test() << "Member " << jm << ": " << ens[jm] << std::endl;

          if (outputPerturbations != boost::none) {
            // Update config
            auto outputPerturbationsUpdated = *outputPerturbations;
            util::setMember(outputPerturbationsUpdated, jm+1);
            setMPI(outputPerturbationsUpdated, ntasks);

            // Write perturbation
            ens[jm].write(outputPerturbationsUpdated);
          }

          if (outputStates != boost::none) {
            // Update config
            auto outputStatesUpdated = *outputStates;
            util::setMember(outputStatesUpdated, jm+1);
            setMPI(outputStatesUpdated, ntasks);

            // Add background state to perturbation
            State_ xp(xx[0]);
            xp += ens[jm];

            // Write state
            xp.write(outputStatesUpdated);
          }

          oops::Log::info() << "Info     : " << std::endl;
        }
      }

      if (outputVariance != boost::none) {
        oops::Log::info() << "Info     : Write randomized variance:" << std::endl;
        oops::Log::info() << "Info     : --------------------------" << std::endl;
        oops::Log::info() << "Info     : " << std::endl;
        if (Bmat->randomizationSize() > 1) {
          // Normalize variance
          double rk_norm = 1.0/static_cast<double>(Bmat->randomizationSize());
          variance *= rk_norm;
        }

        // Update config
        auto outputVarianceUpdated = *outputVariance;
        setMPI(outputVarianceUpdated, ntasks);

        // Write variance
        variance[0].write(outputVarianceUpdated);
        oops::Log::test() << "Randomized variance: " << variance << std::endl;
      }
    }
  }
// -----------------------------------------------------------------------------
};

}  // namespace saber
