/*
 * (C) Crown Copyright 2024-2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HydroBalm1.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "mo/constants.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace vader {

namespace {

using atlas::array::make_view;
using atlas::idx_t;

void eval_hydrobalm1_virtual_potential_temperature_tl(atlas::FieldSet & incFlds,
                                                      const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() <<
    "[eval_hydrobalm1_virtual_potential_temperature_tl()] starting ..." << std::endl;
  auto hlView = make_view<const double, 2>(augStateFlds["height_above_mean_sea_level_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto hexnerIncView = make_view<const double, 2>(incFlds["hydrostatic_exner_levels"]);
  auto thetavIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  const idx_t n_levels = incFlds["virtual_potential_temperature"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["virtual_potential_temperature"].functionspace());
  for (idx_t ih = 0; ih < sizeOwned; ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      thetavIncView(ih, ilev) =
        (hexnerIncView(ih, ilev+1) - hexnerIncView(ih, ilev)) *
        (::mo::constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)) /
        (::mo::constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)));
    }
  }
  incFlds["virtual_potential_temperature"].set_dirty();
  oops::Log::trace()
    << "[eval_hydrobalm1_virtual_potential_temperature_tl()] ... done" << std::endl;
}

void eval_hydrobalm1_virtual_potential_temperature_ad(atlas::FieldSet & hatFlds,
                                                      const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() <<
    "[eval_hydrobalm1_virtual_potential_temperature_ad] starting ..." << std::endl;

  auto hlView = make_view<const double, 2>(augStateFlds["height_above_mean_sea_level_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto hexnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  const idx_t tlev = hatFlds["virtual_potential_temperature"].shape(1);
  const idx_t sizeOwned =
        util::getSizeOwned(hatFlds["virtual_potential_temperature"].functionspace());

  atlas_omp_parallel_for(idx_t ih = 0; ih < sizeOwned; ++ih) {
    for (idx_t ilev = tlev-1; ilev > -1; --ilev) {
      hexnerHatView(ih, ilev+1) += thetavHatView(ih, ilev) *
        (::mo::constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)) /
        (::mo::constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)) );
      hexnerHatView(ih, ilev) -= thetavHatView(ih, ilev) *
        (::mo::constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)) /
        (::mo::constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)));
      thetavHatView(ih, ilev) = 0.0;
    }
  }
  hatFlds["hydrostatic_exner_levels"].set_dirty();
  hatFlds["virtual_potential_temperature"].set_dirty();

  oops::Log::trace() << "[eval_hydrobalm1_virtual_potential_temperature_ad] ... done" << std::endl;
}

void eval_hydrobalm1_hydrostatic_exner_levels_tl(atlas::FieldSet & incFlds,
                                                 const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_hydrobalm1_hydrostatic_exner_levels_tl()] ... starting" << std::endl;

  auto hlView = make_view<const double, 2>(augStateFlds["height_above_mean_sea_level_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto thetavIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto exnerIncView = make_view<const double, 2>(
    incFlds["dimensionless_exner_function_levels_minus_one"]);
  auto hexnerIncView = make_view<double, 2>(incFlds["hydrostatic_exner_levels"]);

  const idx_t n_levels = incFlds["dimensionless_exner_function_levels_minus_one"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["dimensionless_exner_function_levels_minus_one"].functionspace());

  for (idx_t ih = 0; ih < sizeOwned; ++ih) {
    hexnerIncView(ih, 0) = exnerIncView(ih, 0);
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      hexnerIncView(ih, ilev+1) = hexnerIncView(ih, ilev) +
        ((::mo::constants::grav * thetavIncView(ih, ilev) *
          (hlView(ih, ilev+1) - hlView(ih, ilev))) /
         (::mo::constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)));
    }
  }
  incFlds["hydrostatic_exner_levels"].set_dirty();
  oops::Log::trace() << "[eval_hydrobalm1_hydrostatic_exner_levels_tl()] ... done" << std::endl;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HydroBalm1> makerHydroBalm1_("mo_hydro_bal2");

// -----------------------------------------------------------------------------

HydroBalm1::HydroBalm1(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::HydroBalm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

HydroBalm1::~HydroBalm1() {
  oops::Log::trace() << classname() << "::~HydroBalm1 starting" << std::endl;
  util::Timer timer(classname(), "~HydroBalm1");
  oops::Log::trace() << classname() << "::~HydroBalm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalm1::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  eval_hydrobalm1_virtual_potential_temperature_tl(fset.fieldSet(), xb_.fieldSet());

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalm1::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  eval_hydrobalm1_virtual_potential_temperature_ad(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalm1::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  eval_hydrobalm1_hydrostatic_exner_levels_tl(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalm1::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  oops::Log::info() << classname() << "::directCalibration (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalm1::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
