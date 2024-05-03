/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "oops/base/GeometryData.h"

namespace saber {
namespace mgbf {

// -------------------------------------------------------------------------------------------------
class mgbf_InterpolationParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, SaberBlockParametersBase)

 public:
  // File containing grid and coefficients
  oops::RequiredParameter<std::string> mgbfFile{"mgbf error covariance file", this};
  oops::RequiredParameter<std::string> mgbfNML{"mgbf berror namelist file", this};
  oops::RequiredParameter<std::string> mgbfGRD{"gsi akbk", this};

  // Handle vertical top-2-bottom and vice-verse wrt to GSI
  oops::Parameter<bool> vflip{"flip vertical grid", true, this};

  // Processor layout

  // Debugging mode
  oops::Parameter<bool> debugMode{"debugging mode", false, this};
  oops::Parameter<bool> bypassGSI{"debugging bypass gsi", false, this};

  // Mandatory active variables
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -------------------------------------------------------------------------------------------------

template <typename T_interpolator>
class mgbfInterpolation : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::mgbf::Interpolation";}

  typedef mgbfInterpolationParameters Parameters_;
  typedef T_Interpolator  Interpolator_;

  mgbfInterpolation( const oops::Geometry & targetGeometry;
                const oops::Variables &,
                const eckit::Configuration &,
                const Parameters_ &,
                const oops::FieldSet3D &,
                const oops::FieldSet3D &);
  virtual ~Interpolation();
// source stuff are corresponding to stuff with the  innner block
// target stuff are corresponding to stuff with the outer block
  const oops::GeometryData & innerGeometryData() const override {return *innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}
  const atlas::functionspace outerFunctionspace ,

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::GeometryData> innerGeometryData_;
  oops::Variables innerVars_;
  const oops::functionspace targetFunctionspace_ ,

  // Interpolation object
  std::unique_ptr<oops::GlobalInterpolator> interpolator_;

  // Inverse interpolation object (need adjoint)
};

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
