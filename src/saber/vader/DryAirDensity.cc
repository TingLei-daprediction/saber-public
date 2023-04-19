/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryAirDensity.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/eval_dry_air_density.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryAirDensity> makerDryAirDensity_("mo_dry_air_density");

// -----------------------------------------------------------------------------

DryAirDensity::DryAirDensity(const oops::GeometryData & outerGeometryData,
                             const std::vector<size_t> & activeVariableSizes,
                             const oops::Variables & outerVars,
                             const eckit::Configuration & covarConf,
                             const Parameters_ & params,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg)
  : innerGeometryData_(outerGeometryData), innerVars_(outerVars), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::DryAirDensity starting" << std::endl;

  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredStateVariables{ "air_temperature",
                                                   "air_pressure",
                                                   "air_pressure_levels_minus_one",
                                                   "dlsvpdT",
                                                   "dry_air_density_levels_minus_one",
                                                   "exner",
                                                   "exner_levels_minus_one",
                                                   "height",
                                                   "height_levels",
                                                   "m_ci",
                                                   "m_cl",
                                                   "m_r",
                                                   "m_v",
                                                   "m_t",
                                                   "potential_temperature",
                                                   "qsat",
                                                   "specific_humidity",
                                                   "svp",
                                                   "virtual_potential_temperature"};


  // Check that they are allocated (i.e. exist in the state fieldset)
  for (auto & s : requiredStateVariables) {
    if (!xb.has(s)) {
      oops::Log::error() << "::DryAirDensity variable " << s <<
                            "is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  std::vector<std::string> requiredGeometryVariables{"height_levels",
                                                     "height"};
  for (const auto & s : requiredGeometryVariables) {
    if (outerGeometryData.fieldSet().has(s)) {
      augmentedStateFieldSet_.add(outerGeometryData.fieldSet()[s]);
    } else {
      augmentedStateFieldSet_.add(xb[s]);
    }
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::eval_sat_vapour_pressure_nl(params.svp_file, augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalVirtualPotentialTemperature(augmentedStateFieldSet_);
  mo::eval_dry_air_density_nl(augmentedStateFieldSet_);

  augmentedStateFieldSet_.haloExchange();

  oops::Log::trace() << classname() << "::DryAirDensity done" << std::endl;
}

// -----------------------------------------------------------------------------

DryAirDensity::~DryAirDensity() {
  oops::Log::trace() << classname() << "::~DryAirDensity starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensity");
  oops::Log::trace() << classname() << "::~DryAirDensity done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::eval_dry_air_density_tl(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------

void DryAirDensity::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::eval_dry_air_density_ad(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::leftInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
