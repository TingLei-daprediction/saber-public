/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/mgbf/covariance/mgbf_Interpolation.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Variables.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/mgbf/covariance/mgbf_Grid.h"

#include "saber/oops/Utilities.h"

namespace saber {
namespace mgbf {

// -------------------------------------------------------------------------------------------------

//clt template <typename T>
static SaberOuterBlockMaker<mgbf_Interpolation> makerInterpolation_("mgbf interpolation to model grid");

// -------------------------------------------------------------------------------------------------

//clt mgbfInterpolation::mgbfInterpolation(const oops::GeometryData & outerGeometryData,
mgbf_Interpolation::mgbf_Interpolation(const oops::GeometryData & outerGeometryData,
                             const oops::Variables & outerVars,
                             const eckit::Configuration & covarConf,
                             const Parameters_ & params,
                             const oops::FieldSet3D & xb,
                             const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()), innerVars_(outerVars)

{
  oops::Log::trace() << classname() << "::mgbf_Interpolation constructor starting" << std::endl;
  util::Timer timer(classname(), "MGBF Interpolation");

  // Grid
  oops::Log::trace()<<"in mgbf interp params "<<params<<std::endl;  
  mgbfGrid grid(outerGeometryData.comm(), params.mgbfgrid.value());

  // Inner geometry and variables
  oops::Log::trace()<<"in mgbf interp after grid "<<std::endl;  
  innerGeometryData_.reset(new oops::GeometryData(grid.functionSpace(),
                                                  outerGeometryData.fieldSet(),
                                                  outerGeometryData.levelsAreTopDown(),
                                                  outerGeometryData.comm()));

  // Active variables
  const oops::Variables activeVars = getActiveVars(params, outerVars);
  std::vector<size_t> activeVariableSizes;
  for (const auto & var : activeVars) {
    activeVariableSizes.push_back(var.getLevels());
  }
  oops::Log::trace()<<"in mgbf interp before interpolator "<<std::endl;  
  interpolator_.reset(new gsi::UnstructuredInterpolation(outerGeometryData.comm(),
                                                    params.toConfiguration(),
                                                    innerGeometryData_->functionSpace(),
                                                    outerGeometryData.functionSpace(),
                                                    activeVariableSizes,
                                                    activeVars));
  oops::Log::trace() << classname() << "mgbf::Interpolator constructor  done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

mgbf_Interpolation::~mgbf_Interpolation() {
  oops::Log::trace() << classname() << "::~mgbfInterpolation starting" << std::endl;
  util::Timer timer(classname(), "~mgbfInterpolation");
  oops::Log::trace() << classname() << "::~mgbfInterpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void mgbf_Interpolation::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting <<fset" << std::endl;
//  std::cout <<"mgbf_Interpoalation::multiply starting <<fset "<<std::endl; 
//  fset.print(std::cout) ;
  util::Timer timer(classname(), "multiply");
  interpolator_->apply(fset.fieldSet());
//  std::cout <<"mgbf_Interpolation::multiply starting 2 fset" << std::endl;
//  fset.print(std::cout) ;
  oops::Log::trace() << "mgbf_Interpolation::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void mgbf_Interpolation::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
//  std::cout<<"thinkdeb mgbf::interpolation::mutliplyAD cout<<fset "<<std::endl;
//  fset.print(std::cout);
  util::Timer timer(classname(), "multiplyAD");
  interpolator_->applyAD(fset.fieldSet());
//  std::cout<<"thinkdeb after mgbf::interpolation::mutliplyAD cout<<fset "<<std::endl;
//  fset.print(std::cout);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
/*
//clt void Interpolation::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "leftInverseMultiply");
  inverseInterpolator_->apply(fset.fieldSet());
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
} */

// -------------------------------------------------------------------------------------------------

void mgbf_Interpolation::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace mgbf
}  // namespace saber
