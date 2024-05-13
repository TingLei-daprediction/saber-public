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

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/gsi/grid/Grid.interface.h"
//cltorg modified from gsi/Grid.h

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class GridParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GridParameters, Parameters)
};

// -------------------------------------------------------------------------------------------------

class Grid {
 public:
  static const std::string classname() {return "saber::gsi::Grid";}

  // Constructor & destructor
  Grid(const eckit::mpi::Comm &, const eckit::Configuration &);
  ~Grid();

  // Accessor functions
  int levels() {return gsiLevels_;}
  const atlas::FunctionSpace & functionSpace() const {return mgbfGridFuncSpace_;}
  atlas::FunctionSpace & functionSpace() {return mgbfGridFuncSpace_;}

 private:
  void print(std::ostream &) const;
  // Fortran LinkedList key
  GridKey keySelf_;
  // Function spaces
  atlas::FunctionSpace mgbfGridFuncSpace_;
  // Number of levels
  int mgbfLevels_;
};

Grid::Grid(const eckit::mpi::Comm & comm, const eckit::Configuration & conf)
{
  oops::Log::trace() << classname() << "::Grid starting" << std::endl;
  util::Timer timer(classname(), "Grid");
// set up grid and functionspace based on description in the yaml
   const atlas::StructuredGrid grid(conf);
//clt  how about PointCloud functionspace 
   const atlas::functionspace::StructuredColumns mgbfGridFuncSpace_(grid);


 //clt mgbfGridFuncSpace_ = atlas::functionspace::PointCloud(lonlat);

  oops::Log::trace() << classname() << "::Grid done" << std::endl;
//
}

// -------------------------------------------------------------------------------------------------

Grid::~Grid() {
  oops::Log::trace() << classname() << "::~Grid starting" << std::endl;
  util::Timer timer(classname(), "~Grid");
}

// -------------------------------------------------------------------------------------------------

void Grid::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
//

}  // namespace mgbf 
}  // namespace saber
