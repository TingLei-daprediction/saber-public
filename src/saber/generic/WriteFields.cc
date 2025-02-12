/*
 * (C) Crown Copyright 2023-2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/WriteFields.h"

#include <sys/stat.h>

#include <sstream>

#include "atlas/functionspace.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator/detail/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/mpi/mpi.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

namespace {
void makeLocalGmshOutput(const std::string& fileName,
                         const atlas::FieldSet& fields) {
  const auto& fs = fields[0].functionspace();
  auto configGmsh =
    atlas::util::Config("coordinates", "xyz") |
    atlas::util::Config("ghost", true) |
    atlas::util::Config("info", true);
  if (fs.type() == "NodeColumns") {
    const auto& mesh = atlas::functionspace::NodeColumns(fs).mesh();
    const auto gmsh = atlas::output::Gmsh(fileName, configGmsh);
    gmsh.write(mesh);
    gmsh.write(fields, fs);
  } else if (fs.type() == "StructuredColumns") {
    if ((fs.distribution() == "ectrans") || ((fs.distribution() == "serial") &&
       (eckit::mpi::comm().size() == 1))) {
      atlas::StructuredMeshGenerator meshgenerator;
      atlas::Mesh mesh = meshgenerator.generate(
        atlas::functionspace::StructuredColumns(fs).grid(),
        atlas::grid::Partitioner(new atlas::grid::detail::partitioner::TransPartitioner()));
      const auto gmsh = atlas::output::Gmsh(fileName, configGmsh);
      gmsh.write(mesh);
      gmsh.write(fields, fs);
    } else {
      oops::Log::info() << "WriteFields::StructuredColumns partitioner name = "
                        << fs.distribution() << " not accounted for in code."
                        << std::endl;
    }
  } else {
    std::cout << fs.type() <<std::endl;
    throw eckit::NotImplemented(
      "functionspace type", Here());
  }
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<WriteFields> makerWriteFields_("write fields");

// -----------------------------------------------------------------------------

void WriteFields::writeToFile(const oops::FieldSet3D & fset,
                              const std::string & description,
                              size_t & count) const {
  // Select the fields to write out.
  // If the `fieldNames` list is empty, write out all fields in `fset`.
  oops::FieldSet3D fsetWrite(fset.validTime(), fset.commGeom());
  oops::Variables variablesToWrite;
  if (params_.fieldNames.value().empty()) {
    variablesToWrite = fset.variables();
  } else {
    variablesToWrite = oops::Variables(params_.fieldNames.value());
  }
  for (const auto & field : fset) {
    if (variablesToWrite.has(field.name())) {
      fsetWrite.add(field);
    }
  }

  // Check there is at least one field to write out.
  if (fsetWrite.size() == 0) {
    oops::Log::warning() << "No field information to output." << std::endl;
    return;
  }

  // Produce a filename and increment the relevant counter.
  std::stringstream filepath;
  filepath << params_.outputPath.value() << "/"
           << description << "_" << count++;
  const std::string filepathnc = filepath.str() + ".nc";

  if (params_.saveNetCDFFile) {
    eckit::LocalConfiguration conf;
    conf.set("filepath", filepath.str());

    // todo: replace the use of `stat` with `std::filesystem::exists`
    // when compilers have been upgraded.
    struct stat buffer;
    if (stat(filepathnc.c_str(), &buffer) == 0) {
      oops::Log::warning() << "File " << filepathnc << " already exists. Overwriting." << std::endl;
    }

    // Write field set to file.
    fsetWrite.write(conf);

    // Output filename to test stream.
    oops::Log::test() << "Wrote file " << filepathnc << std::endl;
  } else {
    // Output filename to test stream.
    oops::Log::test() << "Did not write file " << filepathnc << std::endl;
  }

  if (params_.saveGMSHFile) {
    std::string fileName = filepath.str() + ".gmsh";
    makeLocalGmshOutput(fileName, fsetWrite.fieldSet());
    // Output filename to test stream.
    oops::Log::test() << "Wrote file " << fileName << std::endl;
  }

  // Output basic field information to test output stream.
  oops::Log::test() << "Fields:" << std::endl;
  for (const auto & field : fsetWrite) {
    oops::Log::test() << "  " << field.name() << ": "
                      << util::normField(field, innerGeometryData_.comm()) << std::endl;
  }
  oops::Log::info() << std::endl;
}

// -----------------------------------------------------------------------------

WriteFields::WriteFields(const oops::GeometryData & outerGeometryData,
                         const oops::Variables & outerVars,
                         const eckit::Configuration & covarConfig,
                         const Parameters_ & params,
                         const oops::FieldSet3D & xb,
                         const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    params_(params),
    count_xb_(1),
    count_fg_(1),
    count_multiply_(1),
    count_multiplyad_(1),
    count_leftinversemultiply_(1)
{
  oops::Log::trace() << classname() << "::WriteFields starting" << std::endl;

  if (params_.XbFileName.value() != ::boost::none) {
    oops::FieldSet3D fset(xb.validTime(), outerGeometryData.comm());
    fset.shallowCopy(xb.fieldSet());
    writeToFile(fset, params_.XbFileName.value().value() , count_xb_);
  }

  if (params_.FgFileName.value() != ::boost::none) {
    oops::FieldSet3D fset(fg.validTime(), outerGeometryData.comm());
    fset.shallowCopy(fg.fieldSet());
    writeToFile(fset, params_.FgFileName.value().value(), count_fg_);
  }

  oops::Log::trace() << classname() << "::WriteFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (params_.multiplyFileName.value() != ::boost::none) {
    writeToFile(fset, params_.multiplyFileName.value().value(), count_multiply_);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (params_.multiplyADFileName.value() != ::boost::none) {
    writeToFile(fset, params_.multiplyADFileName.value().value(), count_multiplyad_);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (params_.leftInverseFileName.value() != ::boost::none) {
    writeToFile(fset, params_.leftInverseFileName.value().value(), count_leftinversemultiply_);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
