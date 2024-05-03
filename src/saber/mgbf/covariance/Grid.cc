/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <memory>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"


namespace saber {
namespace mgbf {

// -------------------------------------------------------------------------------------------------

Grid::Grid(const eckit::mpi::Comm & comm, const eckit::Configuration & conf)
{
  oops::Log::trace() << classname() << "::Grid starting" << std::endl;
  util::Timer timer(classname(), "Grid");

  // Create grid
 todo
  // Get number of levels
 todo
  // Create a regional latlon  grid
  //
  //
  //
  //
   // Goal: find number of grid points on a 90deg arc giving the requested resolution or higher.
  //       N will be used to construct an atlas LN latlon grid with (2N+1) lats and 4N lons.
  const double resolDegrees = conf.getDouble("resolution in degrees");
  ASSERT(resolDegrees > 0.0 && resolDegrees < 90.0);
  const double ratio = 90.0 / resolDegrees;
  const auto isNearlyInt = [](const double x) { return fabs(round(x) - x) < 1.e-6; };
//copied from oops:LatLonGridWriter.h
  // Want final results on rank 0 only, so use atlas's serial partitioner to assign all grid points
  // to mpi rank 0.
//clt   const std::string atlasGridName = "L" + std::to_string(gridRes_);
  const atlas::RegularRegionalLonLatGrid grid(atlasGridName); //todo
  eckit::LocalConfiguration serialPartConf{};
  serialPartConf.set("partition", 0);
  const atlas::grid::Partitioner part("serial", serialPartConf);
  targetFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(grid, part));

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

}  // namespace mgbf
}  // namespace saber
