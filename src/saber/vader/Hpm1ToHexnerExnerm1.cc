/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/Hpm1ToHexnerExnerm1.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"

#include "mo/constants.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<Hpm1ToHexnerExnerm1>
  makerHpm1ToHexnerExnerm1_("mo_hp_lvlsm1_to_hexner_lvls");

// -----------------------------------------------------------------------------

Hpm1ToHexnerExnerm1::Hpm1ToHexnerExnerm1(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & outerVars,
                                         const eckit::Configuration & covarConf,
                                         const Parameters_ & params,
                                         const oops::FieldSet3D & xb,
                                         const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(getUnionOfInnerActiveAndOuterVars(params, outerVars)),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    xb_(xb.validTime(), xb.commGeom())
{
  xb_.shallowCopy(xb);
  oops::Log::trace() << classname() << "::Hpm1ToHexnerExnerm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

Hpm1ToHexnerExnerm1::~Hpm1ToHexnerExnerm1() {
  oops::Log::trace() << classname() << "::~Hpm1ToHexnerExnerm1 starting" << std::endl;
  util::Timer timer(classname(), "~Hpm1ToHexnerExnerm1");
  oops::Log::trace() << classname() << "::~Hpm1ToHexnerExnerm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

void Hpm1ToHexnerExnerm1::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // State fields
  auto exnerView =
    make_view<const double, 2>(xb_.fieldSet()["dimensionless_exner_function_levels_minus_one"]);
  auto hexnerView = make_view<const double, 2>(xb_.fieldSet()["hydrostatic_exner_levels"]);
  auto hpView = make_view<const double, 2>(xb_.fieldSet()["hydrostatic_pressure_levels"]);
  auto pView = make_view<const double, 2>(xb_.fieldSet()["air_pressure_levels"]);

  // Input increments
  const auto hpIncView = atlas::array::make_view<const double, 2>(
    fset["hydrostatic_pressure_levels_minus_one"]);

  // Populate output fields.
  auto hExnerIncView = atlas::array::make_view<double, 2>(
    fset["hydrostatic_exner_levels"]);
  auto exnerIncView = atlas::array::make_view<double, 2>(
    fset["dimensionless_exner_function_levels_minus_one"]);

  const atlas::idx_t sizeOwned =
    util::getSizeOwned(fset["hydrostatic_pressure_levels_minus_one"].functionspace());
  const atlas::idx_t lvls =
    fset["hydrostatic_pressure_levels_minus_one"].shape(1);

  for (atlas::idx_t jn = 0; jn < sizeOwned; ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      hExnerIncView(jn, jl) = hpIncView(jn, jl) *
        (mo::constants::rd_over_cp * hexnerView(jn, jl)) / hpView(jn, jl);
      exnerIncView(jn, jl) = hpIncView(jn, jl) *
        (mo::constants::rd_over_cp * exnerView(jn, jl)) / pView(jn, jl);
    }
    hExnerIncView(jn, lvls) = hpIncView(jn, lvls-1) *
      std::pow(pView(jn, lvls-1)/pView(jn, lvls), mo::constants::rd_over_cp - 1.0) *
      (mo::constants::rd_over_cp * hexnerView(jn, lvls)) / hpView(jn, lvls);
  }

  fset["dimensionless_exner_function_levels_minus_one"].set_dirty();
  fset["hydrostatic_exner_levels"].set_dirty();

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Hpm1ToHexnerExnerm1::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  // State fields
  auto exnerView =
    make_view<const double, 2>(xb_.fieldSet()["dimensionless_exner_function_levels_minus_one"]);
  auto hexnerView = make_view<const double, 2>(xb_.fieldSet()["hydrostatic_exner_levels"]);
  auto hpView = make_view<const double, 2>(xb_.fieldSet()["hydrostatic_pressure_levels"]);
  auto pView = make_view<const double, 2>(xb_.fieldSet()["air_pressure_levels"]);

  // Input increments
  auto hpIncView = atlas::array::make_view<double, 2>(
    fset["hydrostatic_pressure_levels_minus_one"]);

  // Populate output fields.
  auto hExnerIncView = atlas::array::make_view<double, 2>(
    fset["hydrostatic_exner_levels"]);
  auto exnerIncView = atlas::array::make_view<double, 2>(
    fset["dimensionless_exner_function_levels_minus_one"]);

  const atlas::idx_t sizeOwned =
    util::getSizeOwned(fset["hydrostatic_pressure_levels_minus_one"].functionspace());
  const atlas::idx_t lvls =
    fset["hydrostatic_pressure_levels_minus_one"].shape(1);

  for (atlas::idx_t jn = 0; jn < sizeOwned; ++jn) {
    hpIncView(jn, lvls-1) += hExnerIncView(jn, lvls) *
      std::pow(pView(jn, lvls-1)/pView(jn, lvls), mo::constants::rd_over_cp - 1.0) *
      (mo::constants::rd_over_cp * hexnerView(jn, lvls)) / hpView(jn, lvls);
    hExnerIncView(jn, lvls) = 0.0;

    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      hpIncView(jn, jl) += hExnerIncView(jn, jl) *
        (mo::constants::rd_over_cp * hexnerView(jn, jl)) / hpView(jn, jl);
      hpIncView(jn, jl) += exnerIncView(jn, jl) *
        (mo::constants::rd_over_cp * exnerView(jn, jl)) / pView(jn, jl);
      hExnerIncView(jn, jl) = 0.0;
      exnerIncView(jn, jl) = 0.0;
    }
  }

  fset["dimensionless_exner_function_levels_minus_one"].set_dirty();
  fset["hydrostatic_exner_levels"].set_dirty();
  fset["hydrostatic_pressure_levels_minus_one"].set_dirty();

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void Hpm1ToHexnerExnerm1::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  // State fields
  auto hexnerView = make_view<const double, 2>(xb_.fieldSet()["hydrostatic_exner_levels"]);
  auto hpView = make_view<const double, 2>(xb_.fieldSet()["hydrostatic_pressure_levels"]);

  // Input increments
  auto hExnerIncView = atlas::array::make_view<const double, 2>(
    fset["hydrostatic_exner_levels"]);

  // Populate output fields.
  auto hpIncView = atlas::array::make_view<double, 2>(
    fset["hydrostatic_pressure_levels_minus_one"]);

  const atlas::idx_t sizeOwned =
    util::getSizeOwned(fset["hydrostatic_pressure_levels_minus_one"].functionspace());
  const atlas::idx_t lvls =
    fset["hydrostatic_pressure_levels_minus_one"].shape(1);

  for (atlas::idx_t jn = 0; jn < sizeOwned; ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      hpIncView(jn, jl) = hExnerIncView(jn, jl) / (
        (mo::constants::rd_over_cp * hexnerView(jn, jl)) / hpView(jn, jl));
    }
  }
  fset["hydrostatic_pressure_levels_minus_one"].set_dirty();

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Hpm1ToHexnerExnerm1::directCalibration(const oops::FieldSets & fset) {
  oops::Log::trace() << classname() << "::directCalibration start" << std::endl;
  oops::Log::info() << classname() << "::directCalibration (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::directCalibration end" << std::endl;
}

// -----------------------------------------------------------------------------

void Hpm1ToHexnerExnerm1::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
