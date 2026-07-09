/*
 * (C) Copyright 2024 DOC/NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"

#include "saber/mgbf/covariance/MGBF_Covariance.interface.h"
#include "saber/oops/Utilities.h"


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

  void directCalibration(const oops::FieldSets &) override {};

  void iterativeCalibrationInit() override {};
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override{};
  void iterativeCalibrationFinal() override{};

  void write() const override {};
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
    override {};

  void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const override {};
  void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const override {};


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


MGBF_Covariance::MGBF_Covariance(const oops::GeometryData & geometryData,
                                 const oops::Variables & centralVars,
                                 const eckit::Configuration & covarConf,
                                 const Parameters_ & params,
                                 const oops::FieldSet3D & xb,
                                 const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime(), geometryData, centralVars),
     params_(params), variables_(params.activeVars.value().get_value_or(centralVars).variables()),
     mgbfGridFuncSpace_(geometryData.functionSpace()), comm_(&geometryData.comm())
{
  oops::Log::trace() << classname() << "MGBF::Covariance starting" << std::endl;

  util::Timer timer(classname(), "Covariance");
  eckit::LocalConfiguration mgbf_config = params.toConfiguration();
  if (params.doCalibration()) {
    throw eckit::UserError("doCalibration=.true. is not implemented ", Here());
  }

  // Assert that there is no variable change in this block

  // Function space


  // Create covariance module
  mgbf_covariance_create_f90(keySelf_, *comm_, mgbf_config,
                            mgbfGridFuncSpace_.get(), xb.get(), fg.get());

  oops::Log::trace() << classname() << "::Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

MGBF_Covariance::~MGBF_Covariance() {
  oops::Log::trace() << classname() << "::~Covariance starting" << std::endl;
  util::Timer timer(classname(), "~Covariance");
  mgbf_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void MGBF_Covariance::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");


  for (auto sabField : fset) {
    // Get the name
    const auto fieldName = name(sabField.name());

    // Ensure that the field name is in the input/output list
    const std::string fieldNameStr = fieldName.getString("name");
    if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
      ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
    }

  }


  mgbf_covariance_randomize_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void MGBF_Covariance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  int index_member;
  if (fset.fieldSet().metadata().has("ensemble member index")) {
    index_member = fset.fieldSet().metadata().get<int>("ensemble member index");
  } else {
    index_member = 9999;
  }

  mgbf_covariance_multiply_f90(keySelf_, fset.get(), index_member);
  // Mark all fields as having dirty halos after modification
  for (const auto & fieldname : fset.field_names()) {
    atlas::Field field = fset[fieldname];
    field.set_dirty();  // Mark field as having dirty halos that need to be synchronized
  }
  // Perform the actual halo exchange
  fset.fieldSet().haloExchange();
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//
//
// -------------------------------------------------------------------------------------------------


void MGBF_Covariance::print(std::ostream & os) const {
  os << classname();
}



}  // namespace mgbf
}  // namespace saber
