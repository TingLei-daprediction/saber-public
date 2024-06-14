/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"
#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.interface.h"
#include <fstream>
#include "mpi.h"
#include <string>

namespace saber {
namespace gsi {

// -----------------------------------------------------------------------------
UnstructuredInterpolation::UnstructuredInterpolation(
  const eckit::mpi::Comm & comm,
  const eckit::Configuration & config,
  const atlas::FunctionSpace & innerFuncSpace,
  const atlas::FunctionSpace & outerFuncSpace,
  const std::vector<size_t> & activeVariableSizes,
  const oops::Variables & activeVars)
  : innerFuncSpace_(innerFuncSpace), outerFuncSpace_(outerFuncSpace),
    activeVariableSizes_(activeVariableSizes), activeVars_(activeVars)
{
  oops::Log::trace()<<"thinkdeb in gsi::UnstracutredInterpolation CTOR begin"<<std::endl;
  oops::Log::trace()<<"thinkdeb in gsi::UnstracutredInterpolation inner*.lonlat "<<innerFuncSpace_.lonlat()<<std::endl;
   int mpirank;
   MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  std::ofstream file("gsi_mgbf_latlon_"+std::to_string(mpirank)+".txt");
  innerFuncSpace_.lonlat().dump(file);
  oops::Log::trace()<<"thinkdeb in gsi::UnstracutredInterpolation outer*.lonlat "<<outerFuncSpace_.lonlat()<<std::endl;
  std::ofstream file2("gsi_model_latlon_"+std::to_string(mpirank)+".txt");
  outerFuncSpace_.lonlat().dump(file2);
  oops::Log::trace()<<"xx         "<<std::endl;
  
  
  saber_unstrc_create_f90(keyUnstructuredInterpolator_, &comm,
                          innerFuncSpace_.lonlat().get(), outerFuncSpace_.lonlat().get(), config);
  oops::Log::trace()<<"thinkdeb in gsi::UnstracutredInterpolation CTOR done"<<std::endl;
}

// -----------------------------------------------------------------------------
int UnstructuredInterpolation::write(const eckit::Configuration & config) {
  saber_unstrc_write_f90(keyUnstructuredInterpolator_, config);
  return 0;
}

// -----------------------------------------------------------------------------
UnstructuredInterpolation::~UnstructuredInterpolation() {
  saber_unstrc_delete_f90(keyUnstructuredInterpolator_);
}
// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply(const atlas::Field & innerField,
                                      atlas::Field & outerField) {
  saber_unstrc_apply_f90(keyUnstructuredInterpolator_, innerField.get(), outerField.get());
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::applyAD(const atlas::Field & outerField,
                                        atlas::Field & innerField) {
  saber_unstrc_apply_ad_f90(keyUnstructuredInterpolator_, outerField.get(), innerField.get());
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply(atlas::FieldSet & fset) {
  // TODO(Someone): check if we can get rid of activeVariableSizes everywhere
  // and use Variables levels instead.
  for (size_t i = 0; i < activeVars_.size(); ++i) {
    atlas::Field outerField = outerFuncSpace_.createField<double>(
      atlas::option::name(activeVars_[i].name()) | atlas::option::levels(activeVariableSizes_[i]));
    this->apply(fset[activeVars_[i].name()], outerField);
    util::removeFieldsFromFieldSet(fset, {activeVars_[i].name()});
    fset.add(outerField);
  }
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::applyAD(atlas::FieldSet & fset) {
  for (size_t i = 0; i < activeVars_.size(); ++i) {
    atlas::Field innerField = innerFuncSpace_.createField<double>(
      atlas::option::name(activeVars_[i].name()) | atlas::option::levels(activeVariableSizes_[i]));
    this->applyAD(fset[activeVars_[i].name()], innerField);
    util::removeFieldsFromFieldSet(fset, {activeVars_[i].name()});
    fset.add(innerField);
  }
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::print(std::ostream & os) const {
  os << " UnstructuredInterpolation: print not implemented yet.";
}
// -----------------------------------------------------------------------------
}  // namespace gsi
}  // namespace saber
