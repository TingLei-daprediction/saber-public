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
#include "atlas/grid.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include <fstream> //clt

//cltorg modified from gsi/Grid.h

namespace saber {
namespace mgbf {

// -------------------------------------------------------------------------------------------------

class mgbfGridParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(mgbfGridParameters, Parameters)
};

// -------------------------------------------------------------------------------------------------

class mgbfGrid {
 public:
  static const std::string classname() {return "saber::mgbf::mgbfGrid";}

  // Constructor & destructor
  mgbfGrid(const eckit::mpi::Comm &, const eckit::Configuration &);
  ~mgbfGrid();

  // Accessor functions
//clt  int levels() {return gsiLevels_;}
  const atlas::FunctionSpace & functionSpace() const {return mgbfFuncSpace_;}
  atlas::FunctionSpace & functionSpace() {return mgbfFuncSpace_;}

 private:
  void print(std::ostream &) const;
  // Fortran LinkedList key
//clt  GridKey keySelf_;
  // Function spaces
  atlas::FunctionSpace mgbfFuncSpace_;
  // Number of levels
  int mgbfLevels_;
};

mgbfGrid::mgbfGrid(const eckit::mpi::Comm & comm, const eckit::Configuration & conf)
{
  oops::Log::trace() << classname() << "::Grid starting" << std::endl;
  oops::Log::trace()<<"mgbf config is "<<conf<<std::endl;
  util::Timer timer(classname(), "mgbfGrid");
// set up grid and functionspace based on description in the yaml
   const atlas::StructuredGrid grid(conf);
//clt  how about PointCloud functionspace 
//clt   const atlas::functionspace::StructuredColumns mgbfGridFuncSpace_(grid);
  oops::Log::trace() << "mgbf grid before mgbfFunction space  " << std::endl;
//clt   mgbfFuncSpace_=atlas::functionspace::StructuredColumns(grid); //todo 
   mgbfFuncSpace_=atlas::functionspace::PointCloud(grid); //todo 
  oops::Log::trace() << "mgbfFunction space.lonlat   " <<mgbfFuncSpace_.lonlat()<< std::endl;
  std::ofstream file("latlon.txt");
  mgbfFuncSpace_.lonlat().dump(file) ;
  oops::Log::trace() << "mgbfFunction space  is " <<mgbfFuncSpace_<< std::endl;


 //clt mgbfGridFuncSpace_ = atlas::functionspace::PointCloud(lonlat);

  oops::Log::trace() << classname() << "::Grid done" << std::endl;
//
}

// -------------------------------------------------------------------------------------------------

mgbfGrid::~mgbfGrid() {
  oops::Log::trace() << classname() << "::~Grid starting" << std::endl;
  util::Timer timer(classname(), "~Grid");
}

// -------------------------------------------------------------------------------------------------

void mgbfGrid::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
//

}  // namespace mgbf 
}  // namespace saber
