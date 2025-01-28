/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/GpToHpm1.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/eval_geostrophic_to_hydrostatic_pressure_levels_minus_one.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

namespace {

oops::Variables removeOuterOnlyVar(const oops::Variables & vars) {
  oops::Variables innerVars(vars);
  innerVars -= innerVars["hydrostatic_pressure_levels_minus_one"];
  return innerVars;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GpToHpm1>
  makerGpToHpm1_("mo_hydrostatic_pressure_levels_minus_one_from_geostrophic_pressure");

// -----------------------------------------------------------------------------

GpToHpm1::GpToHpm1(const oops::GeometryData & outerGeometryData,
                   const oops::Variables & outerVars,
                   const eckit::Configuration & covarConf,
                   const Parameters_ & params,
                   const oops::FieldSet3D & xb,
                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(removeOuterOnlyVar(getUnionOfInnerActiveAndOuterVars(params, outerVars))),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    params_(params),
    covFieldSet_(),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::GpToHpm1 starting" << std::endl;
  const oops::Variables stateVariables = params.mandatoryStateVars();
  augmentedStateFieldSet_.clear();
  for (const auto & s : stateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s.name()]);
  }

  oops::Log::trace() << classname() << "::GpToHpm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

GpToHpm1::~GpToHpm1() {
  oops::Log::trace() << classname() << "::~GpToHpm1 starting" << std::endl;
  util::Timer timer(classname(), "~GpToHpm1");
  oops::Log::trace() << classname() << "::~GpToHpm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_hydrostatic_pressure_levels_minus_one_tl(fset.fieldSet(),
                                                    augmentedStateFieldSet_);

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_hydrostatic_pressure_levels_minus_one_ad(fset.fieldSet(),
                                                    augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done"  << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  if (!fset.has("geostrophic_pressure_levels_minus_one")) {
    oops::Log::error() << "The inverse of block "
          << classname()
          << " is not correctly defined if geostrophic_pressure_levels_minus_one"
          << " is not provided as an input." << std::endl;
    throw eckit::UserError("Please only use leftInverseMultiply of this block "
                           "within the mo_hydrostatic_pressure_levels_minus_one "
                           "block.", Here());
  }
  //   Allocate inner-only variables except air temperature
  oops::Variables innerOnlyVarsForInversion(innerOnlyVars_);
  innerOnlyVarsForInversion -=
    innerOnlyVarsForInversion["geostrophic_pressure_levels_minus_one"];
  checkFieldsAreNotAllocated(fset, innerOnlyVarsForInversion);
  allocateMissingFields(fset, innerOnlyVarsForInversion, innerOnlyVarsForInversion,
                        innerGeometryData_.functionSpace());

  // Retrieve unbalanced pressure from hydrostatic pressure and geostrophic pressure.
  mo::eval_hydrostatic_pressure_levels_minus_one_tl_inv(fset.fieldSet(),
                                                        augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

void GpToHpm1::read() {
  oops::Log::trace() << classname() << "::read start " << params_ <<  std::endl;
  const auto & readParams = params_.readParams.value();
  if (readParams != boost::none) {
    eckit::LocalConfiguration lconf;
    readParams.value().serialize(lconf);
    const eckit::Configuration & conf = lconf;
    // Covariance FieldSet
    covFieldSet_ = createGpRegressionStats(innerGeometryData_.functionSpace(),
                                           innerVars_,
                                           conf);
    // also copy variables from covariance fieldset if required
    if (covFieldSet_.has("interpolation_weights")) {
      augmentedStateFieldSet_.add(covFieldSet_["vertical_regression_matrices"]);
      augmentedStateFieldSet_.add(covFieldSet_["interpolation_weights"]);
    }
  }
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

void GpToHpm1::directCalibration(const oops::FieldSets & fset) {
  oops::Log::info() << classname() << "::directCalibration start" << std::endl;
  const auto & calibrationReadParams = params_.calibrationReadParams.value();
  if (calibrationReadParams != boost::none) {
    eckit::LocalConfiguration lconf;
    calibrationReadParams.value().serialize(lconf);
    const eckit::Configuration & conf = lconf;
    // Covariance FieldSet
    covFieldSet_ = createGpRegressionStats(innerGeometryData_.functionSpace(),
                                           innerVars_,
                                           conf);

    // also copy variables from covariance fieldset if required
    if (covFieldSet_.has("interpolation_weights")) {
      augmentedStateFieldSet_.add(covFieldSet_["vertical_regression_matrices"]);
      augmentedStateFieldSet_.add(covFieldSet_["interpolation_weights"]);
    }
  }
  oops::Log::info() << classname() << "::directCalibration end" << std::endl;
}

void GpToHpm1::write() const {
  oops::Log::trace() << classname() << "::write start" << std::endl;
  // write regression matrix to file.
  oops::Log::trace() << classname() << "::write end" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
