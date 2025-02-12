/*
 * (C) Copyright 2021-2024 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace quench {

// -------------------------------------------------------------------------------------------------
/// LinearVariableChange parameters class

class LinearVariableChangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, Parameters)
 public:
  oops::OptionalParameter<oops::Variables> inputVariables{"input variables", this};
  oops::OptionalParameter<oops::Variables> outputVariables{"output variables", this};
  // ATLAS file (multiplicative factor)
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFile{"atlas file", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
