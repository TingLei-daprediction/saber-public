/*
 * (C) Crown Copyright 2024-2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryAirDensityFromExnerm1.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "mo/constants.h"
#include "mo/eval_dry_air_density.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace saber {
namespace vader {

namespace {

using atlas::array::make_view;
using atlas::idx_t;

const char specific_humidity_mo[] = "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water";

/// \details Calculate the dry air density from potential temperature,
///          specific humidity and exner pressure (using exner_pressure_levels_minus_one)
void eval_dry_air_density_from_exner_levels_minus_one_tl(atlas::FieldSet & incFlds,
                                                         const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_tl()] starting ..."
                     << std::endl;
  // State Fields
  auto dryrhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);
  auto exnerView = make_view<const double, 2>(
    stateFlds["dimensionless_exner_function_levels_minus_one"]);
  auto hlView = make_view<const double, 2>(stateFlds["height_above_mean_sea_level_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height_above_mean_sea_level"]);
  auto ptView = make_view<const double, 2>(stateFlds["air_potential_temperature"]);
  auto qView = make_view<const double, 2>(stateFlds[specific_humidity_mo]);
  auto qclView = make_view<const double, 2>(
    stateFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto qcfView = make_view<const double, 2>(
    stateFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"]);

  // Increment Fields
  auto exnerIncView = make_view<const double, 2>(
    incFlds["dimensionless_exner_function_levels_minus_one"]);
  auto ptIncView = make_view<const double, 2>(incFlds["air_potential_temperature"]);
  auto qIncView = make_view<const double, 2>(incFlds[specific_humidity_mo]);
  auto qclIncView = make_view<const double, 2>(
    incFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto qcfIncView = make_view<const double, 2>(
    incFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto dryrhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = incFlds["dry_air_density_levels_minus_one"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["dry_air_density_levels_minus_one"].functionspace());

  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydensInc;
  double vptdrydens_jlm1;
  double vptdrydensInc_jlm1;
  double vptdrydens_intp_times_h_minus_hm1;
  double vptdrydensInc_intp_times_h_minus_hm1;

  for (idx_t jn = 0; jn < sizeOwned; ++jn) {
    vptdrydens = ptView(jn, 0) * (1.0 + ::mo::constants::c_virtual * qView(jn, 0)
                 - qclView(jn, 0) - qcfView(jn, 0)) / (1.0 - qView(jn, 0)
                 - qclView(jn, 0) - qcfView(jn, 0));
    vptdrydensInc = ((1.0 + ::mo::constants::c_virtual * qView(jn, 0)
                      - qclView(jn, 0) - qcfView(jn, 0)) * ptIncView(jn, 0)
                     + ((1.0 + ::mo::constants::c_virtual) * qIncView(jn, 0)
                        * (1.0 - qclView(jn, 0) - qcfView(jn, 0)) * ptView(jn, 0)
                        / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)))
                     + (qclIncView(jn, 0) *  qView(jn, 0)
                       * (1.0 + ::mo::constants::c_virtual) * ptView(jn, 0))
                       / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                     + (qcfIncView(jn, 0) *  qView(jn, 0)
                       * (1.0 + ::mo::constants::c_virtual) * ptView(jn, 0))
                       / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)))
                    / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0));

    dryrhoIncView(jn, 0) = dryrhoView(jn, 0) * (
      (1.0 - ::mo::constants::rd_over_cp)  * exnerIncView(jn, 0) /
      (exnerView(jn, 0) * ::mo::constants::rd_over_cp) -
      vptdrydensInc / vptdrydens);

    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      vptdrydens_jlm1 = vptdrydens;
      vptdrydens = ptView(jn, jl) * (1.0 + ::mo::constants::c_virtual * qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl)) / (1.0 - qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl));
      vptdrydensInc_jlm1 = vptdrydensInc;
      vptdrydensInc = ((1.0 + ::mo::constants::c_virtual * qView(jn, jl)
                        - qclView(jn, jl) - qcfView(jn, jl)) * ptIncView(jn, jl)
                       + ((1.0 + ::mo::constants::c_virtual) * qIncView(jn, jl)
                          * (1.0 - qclView(jn, jl) - qcfView(jn, jl)) * ptView(jn, jl)
                          / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                       + (qclIncView(jn, jl) *  qView(jn, jl)
                          * (1.0 + ::mo::constants::c_virtual) * ptView(jn, jl))
                          / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                       + (qcfIncView(jn, jl) *  qView(jn, jl)
                          * (1.0 + ::mo::constants::c_virtual) * ptView(jn, jl))
                          / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                      / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl));
      vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                           hl_minus_hm1 * vptdrydens;
      vptdrydensInc_intp_times_h_minus_hm1 = hl_minus_hm1 * vptdrydensInc +
                                             h_minus_hl * vptdrydensInc_jlm1;

      dryrhoIncView(jn, jl) = dryrhoView(jn, jl) * (
        (1.0 - ::mo::constants::rd_over_cp) *
        (exnerIncView(jn, jl) /
        (exnerView(jn, jl) * ::mo::constants::rd_over_cp)) -
        vptdrydensInc_intp_times_h_minus_hm1 / vptdrydens_intp_times_h_minus_hm1);
    }
  }
  incFlds["dry_air_density_levels_minus_one"].set_dirty();
  oops::Log::trace() << "[eval_dry_air_density_from_exner_levels_minus_one_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_exner_levels_minus_one_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_exner_levels_minus_one_ad()] starting ..."
                     << std::endl;
  // State fields
  auto dryrhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);
  auto exnerView = make_view<const double, 2>(
    stateFlds["dimensionless_exner_function_levels_minus_one"]);
  auto hlView = make_view<const double, 2>(stateFlds["height_above_mean_sea_level_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height_above_mean_sea_level"]);
  auto ptView = make_view<const double, 2>(stateFlds["air_potential_temperature"]);
  auto qView = make_view<const double, 2>(stateFlds[specific_humidity_mo]);
  auto qclView = make_view<const double, 2>(
    stateFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto qcfView = make_view<const double, 2>(
    stateFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"]);

  // Increment (adjoint) fields
  auto exnerHatView = make_view<double, 2>(
    hatFlds["dimensionless_exner_function_levels_minus_one"]);
  auto ptHatView = make_view<double, 2>(
    hatFlds["air_potential_temperature"]);
  auto qHatView = make_view<double, 2>(hatFlds[specific_humidity_mo]);
  auto qclHatView = make_view<double, 2>(
    hatFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto qcfHatView = make_view<double, 2>(
    hatFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto dryrhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydens_intp_times_h_minus_hm1;
  double vptdrydensHat;
  double vptdrydens_jlm1;
  double vptdrydensHat_jlm1;

  const idx_t numLevels = hatFlds["dry_air_density_levels_minus_one"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(hatFlds["dry_air_density_levels_minus_one"].functionspace());

  for (idx_t jn = 0; jn < sizeOwned; ++jn) {
    for (idx_t jl = numLevels-1; jl >= 1; --jl) {
      // Passive fields.
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      vptdrydens = ptView(jn, jl) * (1.0 + ::mo::constants::c_virtual * qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl)) / (1.0 - qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl));

      exnerHatView(jn, jl) += dryrhoView(jn, jl) *
        (1.0 - ::mo::constants::rd_over_cp)  * dryrhoHatView(jn, jl) /
        (exnerView(jn, jl) * ::mo::constants::rd_over_cp);

      vptdrydens_jlm1 = ptView(jn, jl-1) * (1.0 + ::mo::constants::c_virtual * qView(jn, jl-1)
                        - qclView(jn, jl-1) - qcfView(jn, jl-1)) / (1.0 - qView(jn, jl-1)
                        - qclView(jn, jl-1) - qcfView(jn, jl-1));
      vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                           hl_minus_hm1 * vptdrydens;
      vptdrydensHat = - dryrhoView(jn, jl) * dryrhoHatView(jn, jl) * hl_minus_hm1
                      / vptdrydens_intp_times_h_minus_hm1;
      vptdrydensHat_jlm1 = - dryrhoView(jn, jl) * dryrhoHatView(jn, jl) * h_minus_hl /
                           vptdrydens_intp_times_h_minus_hm1;
      ptHatView(jn, jl) += (1.0 + ::mo::constants::c_virtual * qView(jn, jl)
                            - qclView(jn, jl) - qcfView(jn, jl))
                           / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                           * vptdrydensHat;
      qHatView(jn, jl) += (1.0 + ::mo::constants::c_virtual) * ptView(jn, jl)
                           * (1.0 - qclView(jn, jl) - qcfView(jn, jl))
                           / ((1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                               * (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                           * vptdrydensHat;
      qclHatView(jn, jl) += ptView(jn, jl) *  qView(jn, jl) * (1.0 + ::mo::constants::c_virtual)
                            / ((1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                                * (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                            * vptdrydensHat;
      qcfHatView(jn, jl) += ptView(jn, jl) *  qView(jn, jl) * (1.0 + ::mo::constants::c_virtual)
                            / ((1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                                * (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                            * vptdrydensHat;
      ptHatView(jn, jl-1) += (1.0 + ::mo::constants::c_virtual * qView(jn, jl-1)
                              - qclView(jn, jl-1) - qcfView(jn, jl-1))
                             / (1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                             * vptdrydensHat_jlm1;
      qHatView(jn, jl-1) += (1.0 + ::mo::constants::c_virtual) * ptView(jn, jl-1)
                            * (1.0 - qclView(jn, jl-1) - qcfView(jn, jl-1))
                            / ((1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                                * (1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1)))
                            * vptdrydensHat_jlm1;
      qclHatView(jn, jl-1) += ptView(jn, jl-1) *  qView(jn, jl-1)
                              * (1.0 + ::mo::constants::c_virtual)
                              / ((1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                                  * (1.0 - qView(jn, jl-1) - qclView(jn, jl-1)
                                     - qcfView(jn, jl-1)))
                              * vptdrydensHat_jlm1;
      qcfHatView(jn, jl-1) += ptView(jn, jl-1) *  qView(jn, jl-1)
                              * (1.0 + ::mo::constants::c_virtual)
                              / ((1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                                  * (1.0 - qView(jn, jl-1) - qclView(jn, jl-1)
                                     - qcfView(jn, jl-1)))
                              * vptdrydensHat_jlm1;
      dryrhoHatView(jn, jl) = 0.0;
    }

    exnerHatView(jn, 0) += dryrhoView(jn, 0) *
      (1.0 - ::mo::constants::rd_over_cp) * dryrhoHatView(jn, 0) /
      (exnerView(jn, 0) * ::mo::constants::rd_over_cp);

    vptdrydens = ptView(jn, 0) *
      (1.0 + ::mo::constants::c_virtual * qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)) /
      (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0));
    vptdrydensHat = - dryrhoView(jn, 0) * dryrhoHatView(jn, 0) / vptdrydens;
    ptHatView(jn, 0) += (1.0 + ::mo::constants::c_virtual * qView(jn, 0)
                         - qclView(jn, 0) - qcfView(jn, 0))
                        / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)) * vptdrydensHat;
    qHatView(jn, 0) += (1.0 + ::mo::constants::c_virtual) * ptView(jn, 0)
                       * (1.0 - qclView(jn, 0) - qcfView(jn, 0))
                       / ((1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                       * (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))) * vptdrydensHat;
    qclHatView(jn, 0) += ptView(jn, 0) *  qView(jn, 0) * (1.0 + ::mo::constants::c_virtual)
                         / ((1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                         * (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))) * vptdrydensHat;
    qcfHatView(jn, 0) += ptView(jn, 0) *  qView(jn, 0) * (1.0 + ::mo::constants::c_virtual)
                         / ((1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                         * (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))) * vptdrydensHat;
    dryrhoHatView(jn, 0) = 0.0;
  }
  hatFlds["dimensionless_exner_function_levels_minus_one"].set_dirty();
  hatFlds["air_potential_temperature"].set_dirty();
  hatFlds[specific_humidity_mo].set_dirty();
  hatFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"].set_dirty();
  hatFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"].set_dirty();
  hatFlds["dry_air_density_levels_minus_one"].set_dirty();

  oops::Log::trace() << "[eval_dry_air_density_from_exner_levels_minus_one_ad()] ... exit"
                     << std::endl;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryAirDensityFromExnerm1>
  makerDryAirDensityFromExnerm1_("mo_dry_air_density_from_exnerm1");

// -----------------------------------------------------------------------------

DryAirDensityFromExnerm1::DryAirDensityFromExnerm1(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::DryAirDensityFromExnerm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

DryAirDensityFromExnerm1::~DryAirDensityFromExnerm1() {
  oops::Log::trace() << classname() << "::~DryAirDensityFromExnerm1 starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensityFromExnerm1");
  oops::Log::trace() << classname() << "::~DryAirDensityFromExnerm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  eval_dry_air_density_from_exner_levels_minus_one_tl(fset.fieldSet(),
                                                      xb_.fieldSet());
  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  eval_dry_air_density_from_exner_levels_minus_one_ad(fset.fieldSet(),
                                                      xb_.fieldSet());

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  oops::Log::info() << classname() << "::leftInverseMultiply (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  oops::Log::info() << classname() << "::directCalibration (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
