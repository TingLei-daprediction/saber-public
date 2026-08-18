/*
 * (C) Copyright 2024 DOC/NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"


using atlas::option::levels;
using atlas::option::name;

namespace oops {
  class Variables;
}

namespace saber {
namespace mgbf {

typedef int MGBF_CovarianceKey;

// -------------------------------------------------------------------------------------------------

class MGBF_CovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MGBF_CovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<std::string> SDL_MGBFNML{"mgbf sdl and vdl init namelist file", this};
  oops::OptionalParameter<std::string> MGBFNML{"mgbf namelist file", this};
  oops::OptionalParameter<bool> debugPrint{"debug print", this};
  oops::OptionalParameter<std::string> timerOutputFile{"mgbf timer output file", this};

  // Mandatory active variables
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -------------------------------------------------------------------------------------------------

class MGBF_Covariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::mgbf::Covariance";}
  typedef MGBF_CovarianceParameters Parameters_;

  MGBF_Covariance(const oops::GeometryData & geometryData,
                  const oops::Variables & centralVars,
                  const eckit::Configuration & covarConf,
                  const Parameters_ & params,
                  const oops::FieldSet3D & xb,
                  const oops::FieldSet3D & fg);

  virtual ~MGBF_Covariance();
  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  // Fortran LinkedList key
  MGBF_CovarianceKey keySelf_;
  // Parameter
  Parameters_ params_;
  // Variables
  std::vector<std::string> variables_;
  // Function space
  atlas::FunctionSpace mgbfGridFuncSpace_;
  const eckit::mpi::Comm * comm_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mgbf
}  // namespace saber
