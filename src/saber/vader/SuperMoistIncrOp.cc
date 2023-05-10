/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/SuperMoistIncrOp.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SuperMoistIncrOp>
  makerSuperMoistIncrOp_("mo_super_mio");

// -----------------------------------------------------------------------------

SuperMoistIncrOp::SuperMoistIncrOp(const oops::GeometryData & outerGeometryData,
                                   const std::vector<size_t> & activeVariableSizes,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConf,
                                   const Parameters_ & params,
                                   const atlas::FieldSet & xb,
                                   const atlas::FieldSet & fg,
                                   const util::DateTime & validTimeOfXbFg)
  : innerGeometryData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars)),
    exnerThetaToTemp_(std::make_unique<AirTemperature>(outerGeometryData,
                                                       activeVariableSizes,
                                                       outerVars,
                                                       covarConf,
                                                       params.airTemperature,
                                                       xb, fg, validTimeOfXbFg)),
    MIO_(std::make_unique<MoistIncrOp>(outerGeometryData,
                                       activeVariableSizes,
                                       outerVars,
                                       covarConf,
                                       params.moistIncrOp, xb, fg, validTimeOfXbFg))
{
  oops::Log::trace() << classname() << "::SuperMoistIncrOp starting" << std::endl;
  oops::Log::trace() << classname() << "::SuperMoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

SuperMoistIncrOp::~SuperMoistIncrOp() {
  oops::Log::trace() << classname() << "::~SuperMoistIncrOp starting" << std::endl;
  util::Timer timer(classname(), "~SuperMoistIncrOp");
  oops::Log::trace() << classname() << "::~SuperMoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::addZeroFieldToFieldSet("air_temperature", "potential_temperature", fset);
  exnerThetaToTemp_->multiply(fset);
  MIO_->multiply(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::addZeroFieldToFieldSet("air_temperature", "potential_temperature", fset);
  MIO_->multiplyAD(fset);
  exnerThetaToTemp_->multiplyAD(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  util::addZeroFieldToFieldSet("air_temperature", "potential_temperature", fset);
  exnerThetaToTemp_->multiply(fset);
  MIO_->leftInverseMultiply(fset);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber