#pragma once

#include <cstddef>

namespace saber {
namespace mgbf {

/// Build the inner SABER geometry and return lon/lat coordinates together with
/// the number of owned points.
/// \param[in]  conf_ptr   Pointer to the fckit configuration (C handle)
/// \param[in]  comm_ptr   Pointer to the fckit MPI communicator (C handle)
/// \param[out] lonlat     Newly allocated array of size (npts_total * 2)
/// \param[out] npts_total Total number of grid points (owned + halo)
/// \param[out] npts_owned Number of locally owned grid points
/// \param[out] status     0 on success, non-zero otherwise
extern "C" void saber_mgbf_inner_geom_build(const void *conf_ptr,
                                            const void *comm_ptr,
                                            double **lonlat,
                                            int *npts_total,
                                            int *npts_owned,
                                            int *status);

/// Release the lon/lat array allocated by saber_mgbf_inner_geom_build.
extern "C" void saber_mgbf_inner_geom_free(double *lonlat);

}  // namespace mgbf
}  // namespace saber

