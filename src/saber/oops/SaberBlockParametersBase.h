/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberBlockParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SaberBlockParametersBase, Parameters)
 public:
  // Parameters
  oops::RequiredParameter<std::string> saberBlockName{"saber block name", this};
  oops::OptionalParameter<oops::Variables> activeVars{"active variables", this};
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputFieldConfs{"input fields",
    this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};
  oops::OptionalParameter<double> adjointTolerance{"adjoint tolerance", this};

  // Mandatory active variables
  virtual oops::Variables mandatoryActiveVars() const = 0;
};

// -----------------------------------------------------------------------------

}  // namespace saber
