#include "saber/mgbf/utils/MGBF_GeometryBridge.h"

#include <memory>
#include <stdexcept>

#include "atlas/array.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/log/Log.h"

#include "fckit/config/Configuration.h"
#include "fckit/mpi/Comm.h"

#include "saber/interpolation/Geometry.h"

namespace saber {
namespace mgbf {

namespace {

const char *kInnerGeometryKey = "inner geometry";

const eckit::Configuration &ensureInnerGeometry(const fckit::Configuration &conf,
                                                std::unique_ptr<eckit::Configuration> &holder) {
  if (!conf.has(kInnerGeometryKey)) {
    throw eckit::BadParameter("inner geometry section missing in SABER configuration");
  }
  holder.reset(new eckit::LocalConfiguration(conf.getSubConfiguration(kInnerGeometryKey)));
  return *holder;
}

}  // namespace

extern "C" void saber_mgbf_inner_geom_build(const void *conf_ptr,
                                            const void *comm_ptr,
                                            double **lonlat_out,
                                            int *npts_total_out,
                                            int *npts_owned_out,
                                            int *status_out) {
  if (lonlat_out == nullptr || npts_total_out == nullptr ||
      npts_owned_out == nullptr || status_out == nullptr) {
    if (status_out != nullptr) *status_out = 1;
    return;
  }

  *lonlat_out = nullptr;
  *npts_total_out = 0;
  *npts_owned_out = 0;
  *status_out = 0;

  try {
    const auto *conf_wrapper = reinterpret_cast<const fckit::Configuration *>(conf_ptr);
    const auto *comm_wrapper = reinterpret_cast<const fckit::mpi::Comm *>(comm_ptr);

    if (conf_wrapper == nullptr || comm_wrapper == nullptr) {
      throw eckit::SeriousBug("Null configuration or communicator pointer passed to geometry bridge");
    }

    std::unique_ptr<eckit::Configuration> inner_holder;
    const eckit::Configuration &inner_conf = ensureInnerGeometry(*conf_wrapper, inner_holder);

    saber::interpolation::Geometry geom(inner_conf, comm_wrapper->mpiComm());

    const atlas::FunctionSpace &fs = geom.functionSpace();
    if (fs.type() != "StructuredColumns") {
      throw eckit::BadParameter("Inner geometry must be StructuredColumns for MGBF");
    }

    atlas::functionspace::StructuredColumns structured(fs);
    const atlas::Field lonlatField = structured.lonlat();
    auto lonlatView = atlas::array::make_view<double, 2>(lonlatField);

    const std::size_t npts_total = lonlatView.shape(0);
    const std::size_t ncoords = lonlatView.shape(1);
    if (ncoords != 2) {
      throw eckit::SeriousBug("Unexpected lonlat field rank in geometry bridge");
    }

    std::unique_ptr<double[]> buffer(new double[npts_total * 2]);
    for (std::size_t i = 0; i < npts_total; ++i) {
      buffer[i + 0 * npts_total] = lonlatView(i, 0);
      buffer[i + 1 * npts_total] = lonlatView(i, 1);
    }

    *npts_total_out = static_cast<int>(npts_total);
    *npts_owned_out = static_cast<int>(structured.sizeOwned());
    *lonlat_out = buffer.release();
  } catch (const std::exception &e) {
    *status_out = 1;
    *lonlat_out = nullptr;
    *npts_total_out = 0;
    *npts_owned_out = 0;
    eckit::Log::error() << "saber_mgbf_inner_geom_build: " << e.what() << std::endl;
  }
}

extern "C" void saber_mgbf_inner_geom_free(double *lonlat) {
  delete[] lonlat;
}

}  // namespace mgbf
}  // namespace saber
