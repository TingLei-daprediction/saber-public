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

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"


#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/mgbf/covariance/MGBF_Covariance.interface.h"
//clt #include "saber/mgbf/grid/MGBF_Grid.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include <iostream>


using atlas::option::levels;
using atlas::option::name;

namespace oops {
  class Variables;
}

namespace saber {
namespace mgbf {

// -------------------------------------------------------------------------------------------------
class CovarianceParameters: public SaberBlockParametersBase  {
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters,SaberBlockParametersBase)
  public:
    // Mandatory active variables
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -------------------------------------------------------------------------------------------------

//clt template <typename MODEL>
class Covariance : public SaberCentralBlockBase {
//clt  typedef oops::Increment<MODEL>                          Increment_;

 public:
  static const std::string classname() {return "saber::mgbf::Covariance";}
   typedef CovarianceParameters Parameters_;


//cltorg  Covariance(const Geometry_ &, const Parameters_ &, const State_ &, const State_ &);
Covariance(const oops::GeometryData &,
        const std::vector<size_t> &,
        const oops::Variables &,
        const Parameters_ &,
        const atlas::FieldSet &,
        const atlas::FieldSet &,
        const std::vector<atlas::FieldSet> &,
        const size_t &);

  virtual ~Covariance();

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;
 std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override;
  void setReadFields(const std::vector<oops::FieldSet3D> &) override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override;
  void iterativeCalibrationFinal() override;

  void dualResolutionSetup(const oops::GeometryData &) override;

  void write() const override;
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
    override;

//clttodo  size_t ctlVecSize() const override {return bump_->getCvSize();}
  void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const override;


 private:
  void print(std::ostream &) const override;
  // Fortran LinkedList key
  CovarianceKey keySelf_;
  // Variables
  std::vector<std::string> variables_;
  // Function space
  atlas::FunctionSpace mgbfGridFuncSpace_;
  // Grid
//clt  Grid grid_;
};

// -------------------------------------------------------------------------------------------------


Covariance::Covariance(const oops::GeometryData  & geometryData,
        const std::vector<size_t> & activeVariableSizes,
        const oops::Variables & centralVars,
        const Parameters_ & params,
        const atlas::FieldSet & xbg,
        const atlas::FieldSet & xfg,
        const std::vector<atlas::FieldSet> & fsetVec,
        const size_t & timeRank)
  :  variables_() 


//clt Covariance<MODEL>::Covariance(const Geometry_ & geom, const Parameters_ & params,
 //clt                             const State_ & xbg, const State_ & xfg)
//clt  : SaberCentralBlockBase(params), variables_() 
//clt  : SaberCentralBlockBase<MODEL>(params), variables_(), grid_(geom.getComm(), params)
{
  oops::Log::trace() << classname() << "MGBF::Covariance starting" << std::endl;
  util::Timer timer(classname(), "Covariance");
  std::cout<<"thinkdebconfig0 ifhas -1 "<<std::endl;
  eckit::LocalConfiguration mgbf_config = params.toConfiguration();
  std::cout<<"thinkdebconfig0 ifhas "<<mgbf_config<<std::endl;
//  std::cout<<"thinkdebconfig0 ifhas "<<mgbf_config.has("background error")<<std::endl;
//  std::cout<<"thinkdebconfig0 ifhas "<<mgbf_config.has("test")<<std::endl;
//  std::cout<<"thinkdebconfig "<<mgbf_config.getString("test")<<std::endl;
  
  
  
 
  // Assert that there is no variable change in this block
//clt  ASSERT(params.inputVars.value() == params.outputVars.value());
//clt  variables_ = params.inputVars.value().variables();

  // Function space

  // Need to convert background and first guess to Atlas and MGBF grid.

  // Create covariance module
  mgbf_covariance_create_f90(keySelf_, geometryData.comm(), params.toConfiguration(),
                            xbg.get(), xfg.get());

  oops::Log::trace() << classname() << "::Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Covariance::~Covariance() {
  oops::Log::trace() << classname() << "::~Covariance starting" << std::endl;
  util::Timer timer(classname(), "~Covariance");
  mgbf_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // Ignore incoming fields and create new ones based on the block function space
  // ----------------------------------------------------------------------------
  atlas::FieldSet newFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding fields on mgbf grid
  for (auto sabField : fset) {
      // Get the name
      const auto fieldName = name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Create the mgbf grid field and add to Fieldset
//clt      newFields.add(mgbfGridFuncSpace_.createField<double>(fieldName | levels(sabField.levels())));
  }

  // Replace whatever fields are coming in with the mgbf grid fields
  fset = newFields;

  // Call implementation
//clt  MGBF_covariance_randomize_f90(keySelf_, fset.get());
  mgbf_covariance_randomize_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  mgbf_covariance_multiply_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//
//
// -------------------------------------------------------------------------------------------------

void Covariance::iterativeCalibration(const atlas::FieldSet & fset, const size_t & ie) {
  oops::Log::trace() << classname() << "::iterativeCalibration starting" << std::endl;
//clt  bump_->iterativeUpdate(fset, ie);
  oops::Log::trace() << classname() << "::iterativeCalibration done" << std::endl;
}
void Covariance::getOutputFields(const eckit::LocalConfiguration & config , atlas::FieldSet & fset) const  {
  oops::Log::trace() << classname() << "dummy getOutFields" << std::endl;
 };

// -------------------------------------------------------------------------------------------------
void Covariance::finalSetup() {
  oops::Log::trace() << classname() << "::calibration starting" << std::endl;
//clttothink dump
  oops::Log::trace() << classname() << "::calibration done" << std::endl;
}

void Covariance::print(std::ostream & os) const {
  os << classname();
}



}  // namespace mgbf
}  // namespace saber
