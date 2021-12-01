/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_BUMP_H_
#define SABER_OOPS_BUMP_H_

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/bump/type_bump.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class BUMP_Parameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMP_Parameters, oops::Parameters)

 public:
  // External parameters

  // Universe radius (increment)
  oops::OptionalParameter<eckit::LocalConfiguration> universeRadius{"universe radius", this};
  // Input parameters
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> input{"input", this};
  // Ensemble parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  // Missing value (real)
  oops::OptionalParameter<double> msvalr{"msvalr", this};
  // Grids
  oops::OptionalParameter<eckit::LocalConfiguration> grids{"grids", this};
  // Output parameters
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> output{"output", this};
  // Operators application
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> appConfs{"operators application",
    this};

  // Internal parameters

  // general_param

  // Data directory
  oops::OptionalParameter<std::string> datadir{"datadir", this};
  // Files prefix
  oops::OptionalParameter<std::string> prefix{"prefix", this};
  // Model name ('aro', 'arp', 'fv3', 'gem', 'geos', 'gfs', 'ifs', 'mpas', 'nemo', 'norcpm',
  // 'online', 'qg, 'res' or 'wrf')
  oops::OptionalParameter<std::string> model{"model", this};
  // Verbosity level ('all', 'main' or 'none')
  oops::OptionalParameter<std::string> verbosity{"verbosity", this};
  // Add colors to the log (for display on terminal)
  oops::OptionalParameter<bool> colorlog{"colorlog", this};
  // Default seed for random numbers
  oops::OptionalParameter<bool> default_seed{"default_seed", this};
  // Inter-compilers reproducibility
  oops::OptionalParameter<bool> repro{"repro", this};
  // Reproducibility threshold
  oops::OptionalParameter<double> rth{"rth", this};
  // Parallel NetCDF I/O
  oops::OptionalParameter<bool> parallel_io{"parallel_io", this};
  // Number of I/O processors
  oops::OptionalParameter<int> nprocio{"nprocio", this};
  // Universe radius [in meters]
  oops::OptionalParameter<double> universe_rad{"universe_rad", this};
  // Use CGAL for mesh generation (or STRIPACK instead)
  oops::OptionalParameter<bool> use_cgal{"use_cgal", this};

  // driver_param

  // Localization/hybridization to compute ('cor', 'loc', 'hyb-avg', 'hyb-rnd' or 'dual-ens')
  oops::OptionalParameter<std::string> method{"method", this};
  // Localization strategy ('diag_all', 'common', 'common_weighted', 'specific_univariate' or
  // 'specific_multivariate')
  oops::OptionalParameter<std::string> strategy{"strategy", this};
  // New normality test
  oops::OptionalParameter<bool> new_normality{"new_normality", this};
  // New vertical covariance
  oops::OptionalParameter<bool> new_vbal_cov{"new_vbal_cov", this};
  // Update vertical covariance sequentially
  oops::OptionalParameter<bool> update_vbal_cov{"update_vbal_cov", this};
  // Load local vertical covariance
  oops::OptionalParameter<bool> load_vbal_cov{"load_vbal_cov", this};
  // Write local vertical covariancee
  oops::OptionalParameter<bool> write_vbal_cov{"write_vbal_cov", this};
  // Compute new vertical balance operator
  oops::OptionalParameter<bool> new_vbal{"new_vbal", this};
  // Load local vertical balance operator
  oops::OptionalParameter<bool> load_vbal{"load_vbal", this};
  // Write vertical balance operator
  oops::OptionalParameter<bool> write_vbal{"write_vbal", this};
  // Compute new variance
  oops::OptionalParameter<bool> new_var{"new_var", this};
  // Update variance sequentially
  oops::OptionalParameter<bool> update_var{"update_var", this};
  // Load variance
  oops::OptionalParameter<bool> load_var{"load_var", this};
  // Write variance
  oops::OptionalParameter<bool> write_var{"write_var", this};
  // Compute new sampling moments
  oops::OptionalParameter<bool> new_mom{"new_mom", this};
  // Update sampling moments sequentially
  oops::OptionalParameter<bool> update_mom{"update_mom", this};
  // Load sampling moments
  oops::OptionalParameter<bool> load_mom{"load_mom", this};
  // Write sampling moments
  oops::OptionalParameter<bool> write_mom{"write_mom", this};
  // Compute new HDIAG diagnostics
  oops::OptionalParameter<bool> new_hdiag{"new_hdiag", this};
  // Write HDIAG diagnostics
  oops::OptionalParameter<bool> write_hdiag{"write_hdiag", this};
  // Compute new LCT
  oops::OptionalParameter<bool> new_lct{"new_lct", this};
  // Write LCT
  oops::OptionalParameter<bool> write_lct{"write_lct", this};
  // Load C matrix
  oops::OptionalParameter<bool> load_cmat{"load_cmat", this};
  // Write C matrix
  oops::OptionalParameter<bool> write_cmat{"write_cmat", this};
  // Compute new NICAS parameters
  oops::OptionalParameter<bool> new_nicas{"new_nicas", this};
  // Load local NICAS parameters
  oops::OptionalParameter<bool> load_nicas_local{"load_nicas_local", this};
  // Load global NICAS parameters
  oops::OptionalParameter<bool> load_nicas_global{"load_nicas_global", this};
  // Write local NICAS parameters
  oops::OptionalParameter<bool> write_nicas_local{"write_nicas_local", this};
  // Write global NICAS parameters
  oops::OptionalParameter<bool> write_nicas_global{"write_nicas_global", this};
  // Compute wind transform
  oops::OptionalParameter<bool> new_wind{"new_wind", this};
  // Load local wind transform
  oops::OptionalParameter<bool> load_wind_local{"load_wind_local", this};
  // Write local wind transform
  oops::OptionalParameter<bool> write_wind_local{"write_wind_local", this};
  // Test vertical balance inverse and adjoint
  oops::OptionalParameter<bool> check_vbal{"check_vbal", this};
  // Test NICAS adjoints
  oops::OptionalParameter<bool> check_adjoints{"check_adjoints", this};
  // Test NICAS application on diracs
  oops::OptionalParameter<bool> check_dirac{"check_dirac", this};
  // Test NICAS randomization
  oops::OptionalParameter<bool> check_randomization{"check_randomization", this};
  // Test HDIAG-NICAS consistency
  oops::OptionalParameter<bool> check_consistency{"check_consistency", this};
  // Test HDIAG optimality
  oops::OptionalParameter<bool> check_optimality{"check_optimality", this};
  // Test BUMP with no grid point on the last MPI task
  oops::OptionalParameter<bool> check_no_point_mpi{"check_no_point_mpi", this};
  // Test BUMP with all grid points masked on half of the domain
  oops::OptionalParameter<bool> check_no_point_mask{"check_no_point_mask", this};
  // Test set_parameter interface for correlation
  oops::OptionalParameter<bool> check_set_param_cor{"check_set_param_cor", this};
  // Test set_parameter interface for hybrid case
  oops::OptionalParameter<bool> check_set_param_hyb{"check_set_param_hyb", this};
  // Test set_parameter interface for LCT
  oops::OptionalParameter<bool> check_set_param_lct{"check_set_param_lct", this};
  // Test get_parameter interface for standard-deviation
  oops::OptionalParameter<bool> check_get_param_stddev{"check_get_param_stddev", this};
  // Test get_parameter interface for correlation
  oops::OptionalParameter<bool> check_get_param_cor{"check_get_param_cor", this};
  // Test get_parameter interface for hybrid case
  oops::OptionalParameter<bool> check_get_param_hyb{"check_get_param_hyb", this};
  // Test get_parameter interface for anisotropic localization
  oops::OptionalParameter<bool> check_get_param_Dloc{"check_get_param_Dloc", this};
  // Test get_parameter interface for LCT
  oops::OptionalParameter<bool> check_get_param_lct{"check_get_param_lct", this};
  // Test apply_vbal interfaces
  oops::OptionalParameter<bool> check_apply_vbal{"check_apply_vbal", this};
  // Test apply_stddev interfaces
  oops::OptionalParameter<bool> check_apply_stddev{"check_apply_stddev", this};
  // Test apply_nicas interfaces
  oops::OptionalParameter<bool> check_apply_nicas{"check_apply_nicas", this};

  // files_param

  // Variance files
  oops::OptionalParameter<std::vector<std::string>> fname_var{"fname_var", this};
  // Sampling file
  oops::OptionalParameter<std::string> fname_samp{"fname_samp", this};
  // Vertical covariance files
  oops::OptionalParameter<std::vector<std::string>> fname_vbal_cov{"fname_vbal_cov", this};
  // Vertical balance file
  oops::OptionalParameter<std::string> fname_vbal{"fname_vbal", this};
  // Moments files
  oops::OptionalParameter<std::vector<std::string>> fname_mom{"fname_mom", this};
  // C matrix file
  oops::OptionalParameter<std::string> fname_cmat{"fname_cmat", this};
  // NICAS file
  oops::OptionalParameter<std::string> fname_nicas{"fname_nicas", this};
  // Wind transform file
  oops::OptionalParameter<std::string> fname_wind{"fname_wind", this};

  // model_param

  // Number of levels
  oops::OptionalParameter<int> nl{"nl", this};
  // Levels
  oops::OptionalParameter<std::vector<int>> levs{"levs", this};
  // Level for 2D variables ('first' or 'last')
  oops::OptionalParameter<std::string> lev2d{"lev2d", this};
  // Use pressure logarithm as vertical coordinate (model level if .false.)
  oops::OptionalParameter<bool> logpres{"logpres", this};
  // Number of variables
  oops::OptionalParameter<int> nv{"nv", this};
  // Variables names
  oops::OptionalParameter<std::vector<std::string>> variables{"variables", this};
  // Variable change
  oops::OptionalParameter<std::string> variable_change{"variable_change", this};
  // Do not use geometry mask
  oops::OptionalParameter<bool> nomask{"nomask", this};
  // I/O keys
  oops::OptionalParameter<std::vector<std::string>> io_keys{"io_keys", this};
  // I/O values
  oops::OptionalParameter<std::vector<std::string>> io_values{"io_values", this};
  // Redundant points configuration for the QG model
  oops::OptionalParameter<bool> qg_redundant{"qg_redundant", this};
  // Regional domain configuration for the QG model
  oops::OptionalParameter<bool> qg_regional{"qg_regional", this};
  // Urban domain configuration for the QG model
  oops::OptionalParameter<bool> qg_urban{"qg_urban", this};

  // ens1_param

  // Ensemble 1 size
  oops::OptionalParameter<int> ens1_ne{"ens1_ne", this};
  // Ensemble 1 sub-ensembles number
  oops::OptionalParameter<int> ens1_nsub{"ens1_nsub", this};

  // ens2_param

  // Ensemble 2 size
  oops::OptionalParameter<int> ens2_ne{"ens2_ne", this};
  // Ensemble 2 sub-ensembles number
  oops::OptionalParameter<int> ens2_nsub{"ens2_nsub", this};

  // sampling_param

  // Load local sampling
  oops::OptionalParameter<bool> load_samp_local{"load_samp_local", this};
  // Load global sampling
  oops::OptionalParameter<bool> load_samp_global{"load_samp_global", this};
  // Write local sampling
  oops::OptionalParameter<bool> write_samp_local{"write_samp_local", this};
  // Write global sampling
  oops::OptionalParameter<bool> write_samp_global{"write_samp_global", this};
  // Write sampling grids
  oops::OptionalParameter<bool> write_samp_grids{"write_samp_grids", this};
  // Mask restriction type
  oops::OptionalParameter<std::string> mask_type{"mask_type", this};
  // Mask threshold side ('lower' if mask_th is the lower bound, 'upper' if mask_th is the
  // upper bound)
  oops::OptionalParameter<std::vector<std::string>> mask_lu{"mask_lu", this};
  // Mask threshold
  oops::OptionalParameter<std::vector<double>> mask_th{"mask_th", this};
  // Threshold on vertically contiguous points for sampling mask (0 to skip the test)
  oops::OptionalParameter<int> ncontig_th{"ncontig_th", this};
  // Check that sampling couples and interpolations do not cross mask boundaries
  oops::OptionalParameter<bool> mask_check{"mask_check", this};
  // Sampling draw type ('random_uniform','random_coast' or 'octahedral')
  oops::OptionalParameter<std::string> diag_draw_type{"diag_draw_type", this};
  // Length-scale to increase sampling density along coasts [in meters]
  oops::OptionalParameter<double> Lcoast{"Lcoast", this};
  // Minimum value to increase sampling density along coasts
  oops::OptionalParameter<double> rcoast{"rcoast", this};
  // Number of sampling points
  oops::OptionalParameter<int> nc1{"nc1", this};
  // Number of diagnostic points
  oops::OptionalParameter<int> nc2{"nc2", this};
  // Number of classes
  oops::OptionalParameter<int> nc3{"nc3", this};
  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  oops::OptionalParameter<double> dc{"dc", this};
  // Reduced number of levels for diagnostics
  oops::OptionalParameter<int> nl0r{"nl0r", this};
  // Maximum number of random number draws
  oops::OptionalParameter<int> irmax{"irmax", this};

  // diag_param

  // Ensemble size
  oops::OptionalParameter<int> ne{"ne", this};
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  oops::OptionalParameter<double> gen_kurt_th{"gen_kurt_th", this};
  // Gaussian approximation for asymptotic quantities
  oops::OptionalParameter<bool> gau_approx{"gau_approx", this};
  // Number of bins for averaged statistics histograms
  oops::OptionalParameter<int> avg_nbins{"avg_nbins", this};
  // Activation of vertical balance (ordered line by line in the lower triangular formulation)
  oops::OptionalParameter<std::vector<bool>> vbal_block{"vbal_block", this};
  // Vertical balance diagnostic radius [in meters]
  oops::OptionalParameter<double> vbal_rad{"vbal_rad", this};
  // Vertical balance diagnostic latitude band half-width [in degrees]
  oops::OptionalParameter<double> vbal_dlat{"vbal_dlat", this};
  // Diagonal auto-covariance for the inversion
  oops::OptionalParameter<std::vector<bool>> vbal_diag_auto{"vbal_diag_auto", this};
  // Diagonal regression
  oops::OptionalParameter<std::vector<bool>> vbal_diag_reg{"vbal_diag_reg", this};
  // Pseudo-inverse for auto-covariance
  oops::OptionalParameter<bool> vbal_pseudo_inv{"vbal_pseudo_inv", this};
  // Dominant mode for pseudo-inverse
  oops::OptionalParameter<int> vbal_pseudo_inv_mmax{"vbal_pseudo_inv_mmax", this};
  // Variance threshold to compute the dominant mode for pseudo-inverse
  oops::OptionalParameter<double> vbal_pseudo_inv_var_th{"vbal_pseudo_inv_var_th", this};
  // Force specific variance
  oops::OptionalParameter<bool> forced_var{"forced_var", this};
  // Forced standard-deviation
  oops::OptionalParameter<eckit::LocalConfiguration> stddev{"stddev", this};
  // Filter variance
  oops::OptionalParameter<bool> var_filter{"var_filter", this};
  // Number of iterations for the variance filtering (0 for uniform variance)
  oops::OptionalParameter<int> var_niter{"var_niter", this};
  // Number of passes for the variance filtering (0 for uniform variance)
  oops::OptionalParameter<int> var_npass{"var_npass", this};
  // Variance initial filtering support radius [in meters]
  oops::OptionalParameter<eckit::LocalConfiguration> var_rhflt{"var_rhflt", this};
  // Activate local diagnostics
  oops::OptionalParameter<bool> local_diag{"local_diag", this};
  // Local diagnostics calculation radius [in meters]
  oops::OptionalParameter<double> local_rad{"local_rad", this};
  // Local diagnostics calculation latitude band half-width [in degrees]
  oops::OptionalParameter<double> local_dlat{"local_dlat", this};

  // fit_param

  // Horizontal filtering suport radius [in meters]
  oops::OptionalParameter<double> diag_rhflt{"diag_rhflt", this};
  // Vertical filtering support radius
  oops::OptionalParameter<double> diag_rvflt{"diag_rvflt", this};
  // Number of levels between interpolation levels
  oops::OptionalParameter<int> fit_dl0{"fit_dl0", this};
  // Number of LCT scales
  oops::OptionalParameter<int> lct_nscales{"lct_nscales", this};
  // Factor between diffusion scales
  oops::OptionalParameter<double> lct_scale_ratio{"lct_scale_ratio", this};
  // Minimum relevant correlation for LCT first guess
  oops::OptionalParameter<double> lct_cor_min{"lct_cor_min", this};
  // Diagnostic of diagonal LCT components only
  oops::OptionalParameter<std::vector<bool>> lct_diag{"lct_diag", this};
  // LCT quality control threshold
  oops::OptionalParameter<double> lct_qc_th{"lct_qc_th", this};
  // LCT quality control maximum
  oops::OptionalParameter<double> lct_qc_max{"lct_qc_max", this};
  // Write full correlations
  oops::OptionalParameter<bool> lct_write_cor{"lct_write_cor", this};

  // nicas_param

  // Non-unit diagonal for the NICAS application
  oops::OptionalParameter<bool> nonunit_diag{"nonunit_diag", this};
  // Resolution
  oops::OptionalParameter<double> resol{"resol", this};
  // Maximum size of the Sc1 subset
  oops::OptionalParameter<int> nc1max{"nc1max", this};
  // Subsampling draw type ('random_uniform','random_coast' or 'octahedral')
  oops::OptionalParameter<std::string> nicas_draw_type{"nicas_draw_type", this};
  // Network-base convolution calculation (distance-based if false)
  oops::OptionalParameter<bool> network{"network", this};
  // Force specific support radii
  oops::OptionalParameter<bool> forced_radii{"forced_radii", this};
  // Forced horizontal support radius [in meters]
  oops::OptionalParameter<eckit::LocalConfiguration> rh{"rh", this};
  // Forced vertical support radius
  oops::OptionalParameter<eckit::LocalConfiguration> rv{"rv", this};
  // Minimum level
  oops::OptionalParameter<eckit::LocalConfiguration> min_lev{"min_lev", this};
  // Maximum level
  oops::OptionalParameter<eckit::LocalConfiguration> max_lev{"max_lev", this};
  // Positive-definiteness test
  oops::OptionalParameter<bool> pos_def_test{"pos_def_test", this};
  // Write NICAS fields on model grid (should be written via OOPS if .false.)
  oops::OptionalParameter<bool> write_nicas_c0{"write_nicas_c0", this};
  // Write NICAS grids
  oops::OptionalParameter<bool> write_nicas_grids{"write_nicas_grids", this};

  // dirac_param

  // Number of Diracs
  oops::OptionalParameter<int> ndir{"ndir", this};
  // Diracs longitudes [in degrees]
  oops::OptionalParameter<std::vector<double>> londir{"londir", this};
  // Diracs latitudes [in degrees]
  oops::OptionalParameter<std::vector<double>> latdir{"latdir", this};
  // Diracs level
  oops::OptionalParameter<std::vector<int>> levdir{"levdir", this};
  // Diracs variable indices
  oops::OptionalParameter<std::vector<int>> ivdir{"ivdir", this};

  // output_param

  // Number of neighbors for the full grid smoother
  oops::OptionalParameter<int> full_grid_smoother_nn{"full_grid_smoother_nn", this};
  // Number of local diagnostics profiles to write (for local_diag = .true.)
  oops::OptionalParameter<int> nldwv{"nldwv", this};
  // Index on model grid of the local diagnostics profiles to write
  oops::OptionalParameter<std::vector<int>> img_ldwv{"img_ldwv", this};
  // Longitudes of the local diagnostics profiles to write [in degrees]
  oops::OptionalParameter<std::vector<double>> lon_ldwv{"lon_ldwv", this};
  // Latitudes of the local diagnostics profiles to write [in degrees]
  oops::OptionalParameter<std::vector<double>> lat_ldwv{"lat_ldwv", this};
  // Name of the local diagnostics profiles to write
  oops::OptionalParameter<std::vector<std::string>> name_ldwv{"name_ldwv", this};

  // wind_param

  // Streamfunction variable name
  oops::OptionalParameter<std::string> wind_streamfunction{"wind_streamfunction", this};
  // Velocity potential variable name
  oops::OptionalParameter<std::string> wind_velocity_potential{"wind_velocity_potential", this};
  // Zonal wind variable name
  oops::OptionalParameter<std::string> wind_zonal{"wind_zonal", this};
  // Meridional variable name
  oops::OptionalParameter<std::string> wind_meridional{"wind_meridional", this};
  // Number of longitudes for the regular grid
  oops::OptionalParameter<int> wind_nlon{"wind_nlon", this};
  // Number of latitudes for the regular grid
  oops::OptionalParameter<int> wind_nlat{"wind_nlat", this};
  // Half-width of the Savitzky-Golay to compute derivatives
  oops::OptionalParameter<int> wind_nsg{"wind_nsg", this};
  // Wind inflation to compensate the Savitzky-Golay smoothing
  oops::OptionalParameter<double> wind_inflation{"wind_inflation", this};
};

// -----------------------------------------------------------------------------

template<typename MODEL> class BUMP {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  // Constructors
  BUMP(const Geometry_ &,
       const oops::Variables &,
       const BUMP_Parameters &,
       const EnsemblePtr_ ens1 = NULL,
       const EnsemblePtr_ ens2 = NULL);

  // Copy
  explicit BUMP(BUMP &);

  // Destructor
  ~BUMP();

  // C++ interfaces
  size_t getSize() {return keyBUMP_.size();}
  int getKey(int igrid) const {return keyBUMP_[igrid];}
  void clearKey() {keyBUMP_.clear();}

  // Write / apply operators
  void write() const;
  void apply() const;

  // Fortran interfaces
  void addMember(const atlas::FieldSet &, const int &, const int &) const;
  void updateVbalCov(const Increment_ &, const int &) const;
  void updateVar(const Increment_ &, const int &) const;
  void updateMom(const Increment_ &, const int &) const;
  void runDrivers() const;
  void multiplyVbal(atlas::FieldSet *) const;
  void inverseMultiplyVbal(atlas::FieldSet *) const;
  void multiplyVbalAd(atlas::FieldSet *) const;
  void inverseMultiplyVbalAd(atlas::FieldSet *) const;
  void multiplyStdDev(atlas::FieldSet *) const;
  void inverseMultiplyStdDev(atlas::FieldSet *) const;
  void randomizeNicas(atlas::FieldSet *) const;
  void multiplyNicas(atlas::FieldSet *) const;
  void multiplyPsiChiToUV(atlas::FieldSet *) const;
  void multiplyPsiChiToUVAd(atlas::FieldSet *) const;
  void getParameter(const std::string &, Increment_ &) const;
  void setParameter(const std::string &, const Increment_ &) const;
  void partialDealloc() const;

 private:
  const Geometry_ resol_;
  const oops::Variables activeVars_;
  BUMP_Parameters params_;
  std::vector<int> keyBUMP_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP<MODEL>::BUMP(const Geometry_ & resol,
                  const oops::Variables & activeVars,
                  const BUMP_Parameters & params,
                  const EnsemblePtr_ ens1,
                  const EnsemblePtr_ ens2)
  : resol_(resol), activeVars_(activeVars), params_(params), keyBUMP_() {
  oops::Log::trace() << "BUMP<MODEL>::BUMP construction starting" << std::endl;


  // Get ensemble 1 size if ensemble 1 is available
  int ens1_ne = 0;
  if (ens1) ens1_ne = ens1->size();
  const boost::optional<eckit::LocalConfiguration> &ensembleConfig = params.ensemble.value();
  if (ensembleConfig != boost::none) {
    std::vector<eckit::LocalConfiguration> memberConfig;
    (*ensembleConfig).get("members", memberConfig);
    ens1_ne = memberConfig.size();
  }

  // Get ensemble 2 size if ensemble 2 is available
  int ens2_ne = 0;
  if (ens2) ens2_ne = ens2->size();

  // Read universe size
  oops::Log::info() << "Read universe radius" << std::endl;
  std::unique_ptr<atlas::FieldSet> universe_rad(new atlas::FieldSet());
  const boost::optional<eckit::LocalConfiguration> &universeRadius = params.universeRadius.value();
  if (universeRadius != boost::none) {
    // Setup increment
    util::DateTime time(1977, 5, 25, 0, 0, 0);
    Increment_ dx(resol_, activeVars_, time);
    dx.read(*universeRadius);

    // Get ATLAS fieldset
    dx.toAtlas(universe_rad.get());
  }

  // Add ensemble sizes
  eckit::LocalConfiguration conf(params.toConfiguration());
  conf.set("ens1_ne", ens1_ne);
  conf.set("ens2_ne", ens2_ne);

  // Add missing value
  const double msvalr = util::missingValue(double());
  conf.set("msvalr", msvalr);

  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get global prefix
  std::string prefix;
  if (conf.has("prefix")) {
    conf.get("prefix", prefix);
  } else {
    prefix = "bump";
  }

  // Get the grids configuration from input configuration and complete it
  if (conf.has("grids")) {
    // Get grids from input configuration
    conf.get("grids", grids);
    ASSERT(grids.size() > 0);
  } else {
    // Create one empty configuration
    eckit::LocalConfiguration emptyConf;
    grids.push_back(emptyConf);
  }

  // Loop over grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Add prefix
    if (!grids[jgrid].has("prefix")) {
      std::ostringstream ss;
      ss << std::setw(2) << std::setfill('0') << jgrid;
      grids[jgrid].set("prefix", prefix + "_" + ss.str());
    }

    // Dummy date
    util::DateTime date(1977, 5, 25, 0, 0, 0);

    // Get ATLAS variable names
    Increment_ dx(resol, activeVars_, date);
    std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
    dx.setAtlas(atlasFieldSet.get());
    std::vector<std::string> vars_atlas;
    for (int jvar = 0; jvar < atlasFieldSet->size(); ++jvar) {
      atlas::Field atlasField = atlasFieldSet->field(jvar);
      vars_atlas.push_back(atlasField.name());
    }

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grids[jgrid].has("variables")) {
      grids[jgrid].get("variables", vars_str);
    } else {
      vars_str = vars_atlas;
      grids[jgrid].set("variables", vars_str);
    }
    grids[jgrid].set("nv", vars_str.size());

    // Get the required number of levels add it to the grid configuration
    int nl = 0;
    for (size_t jvar = 0; jvar < vars_str.size(); ++jvar) {
      atlas::Field atlasField = atlasFieldSet->field(vars_str[jvar]);
      nl = std::max(nl, std::max(atlasField.levels(), 1));
    }
    grids[jgrid].set("nl", nl);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }
  }

  // Check grids number
  ASSERT(grids.size() > 0);

  // Print configuration
  oops::Log::info() << "Configuration: " << conf << std::endl;

  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Print configuration for this grid
    oops::Log::info() << "Grid " << jgrid << ": " << grids[jgrid] << std::endl;

    // Create BUMP instance
    int keyBUMP = 0;
    bump_create_f90(keyBUMP, &resol.getComm(),
                    resol.atlasFunctionSpace()->get(),
                    resol.atlasFieldSet()->get(),
                    conf, grids[jgrid], universe_rad->get());
    keyBUMP_.push_back(keyBUMP);
  }

  // Add members of ensemble 1
  if (ens1) {
    oops::Log::info() << "--- Add members of ensemble 1" << std::endl;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      oops::Log::info() << "      Member " << ie+1 << " / " << ens1_ne << std::endl;
      this->addMember((*ens1)[ie].atlas(), ie, 1);
    }
  }

  // Add members of ensemble 2
  if (ens2) {
    oops::Log::info() << "--- Add members of ensemble 2" << std::endl;
    for (int ie = 0; ie < ens2_ne; ++ie) {
      oops::Log::info() << "      Member " << ie+1 << " / " << ens2_ne << std::endl;
      this->addMember((*ens2)[ie].atlas(), ie, 2);
    }
  }

  // Reset parameters
  params_.validateAndDeserialize(conf);

  // Read data from files
  oops::Log::info() << "    Read data from files" << std::endl;
  const boost::optional<std::vector<eckit::LocalConfiguration>> &input = params_.input.value();
  if (input != boost::none) {
    // Set BUMP input parameters
    for (const auto & inputConf : *input) {
      // Get date
      const util::DateTime date(inputConf.getString("date"));

      // Setup increment
      Increment_ dx(resol_, activeVars_, date);
      dx.read(inputConf);

      // Set parameter to BUMP
      std::string param = inputConf.getString("parameter");
      this->setParameter(param, dx);
      oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                        << std::setprecision(3) << dx.norm() << std::endl;
    }
  }

  // Load ensemble members sequentially
  if (ensembleConfig != boost::none) {
    // Get ensemble and members configurations
    std::vector<eckit::LocalConfiguration> memberConfig;
    (*ensembleConfig).get("members", memberConfig);

    // Check what needs to be updated
    const boost::optional<bool> &update_vbal_cov = params_.update_vbal_cov.value();
    const boost::optional<bool> &update_var = params_.update_var.value();
    const boost::optional<bool> &update_mom = params_.update_mom.value();

    // Loop over all ensemble members
    for (int ie = 0; ie < ens1_ne; ++ie) {
      // Get date
      const util::DateTime date(memberConfig[ie].getString("date"));

      // Define increment
      Increment_ incr(resol, activeVars_, date);

      // Read member
      oops::Log::info() <<
      "-------------------------------------------------------------------" << std::endl;
      oops::Log::info() << "--- Load member " << ie+1 << " / " << ens1_ne << std::endl;
      incr.read(memberConfig[ie]);

      if (update_vbal_cov != boost::none) {
        if (*update_vbal_cov) {
          // Update vertical covariance
          this->updateVbalCov(incr, ie);
        }
      }
      if (update_var != boost::none) {
        if (*update_var) {
          // Update variance
          this->updateVar(incr, ie);
        }
      }
      if (update_mom != boost::none) {
        if (*update_mom) {
          // Update moments
          this->updateMom(incr, ie);
        }
      }
    }
  }

  // Run drivers
  this->runDrivers();

  // Partial deallocation
  this->partialDealloc();

  oops::Log::trace() << "BUMP:BUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP<MODEL>::BUMP(BUMP & other) : keyBUMP_() {
  for (unsigned int jgrid = 0; jgrid < other.getSize(); ++jgrid) {
    keyBUMP_.push_back(other.getKey(jgrid));
  }
  other.clearKey();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP<MODEL>::~BUMP() {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    if (keyBUMP_[jgrid] > 0) bump_dealloc_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::write() const {
  oops::Log::trace() << "BUMP::write starting" << std::endl;

  // Write parameters
  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::info() << "--- Write parameters" << std::endl;

  const boost::optional<std::vector<eckit::LocalConfiguration>> &output = params_.output.value();
  if (output != boost::none) {
    for (const auto & outputConf : *output) {
      // Get date
      const util::DateTime date(outputConf.getString("date"));

      // Setup increment
      Increment_ dx(resol_, activeVars_, date);

      // Set increment to zero
      dx.zero();

      // Get parameter from BUMP
      std::string param = outputConf.getString("parameter");
      this->getParameter(param, dx);

      // Write parameter
      dx.write(outputConf);
      oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                        << std::setprecision(3) << dx.norm() << std::endl;
    }
  } else {
    oops::Log::test() << "No output configuration" << std::endl;
  }

  oops::Log::trace() << "BUMP::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::apply() const {
  oops::Log::trace() << "BUMP::apply starting" << std::endl;

  // Aplying operators
  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::info() << "--- Apply operators" << std::endl;

  const boost::optional<std::vector<eckit::LocalConfiguration>>
    &appConfs = params_.appConfs.value();
  if (appConfs != boost::none) {
    if (appConfs->size() > 0) {
      for (const auto & appConf : *appConfs) {
        // Get date
        const util::DateTime date(appConf.getString("date"));

        // Setup increments
        Increment_ dxi(resol_, activeVars_, date);
        Increment_ dxo(resol_, activeVars_, date);

        // Read input file
        eckit::LocalConfiguration inputConf(appConf, "input");
        oops::Log::info() << "       - Input file: " << inputConf << std::endl;
        dxi.read(inputConf);

        // ATLAS transfer
        std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
        dxo.setAtlas(atlasFieldSet.get());
        dxi.toAtlas(atlasFieldSet.get());

        // Apply BUMP operator
        std::vector<std::string> bumpOperators;
        appConf.get("bump operators", bumpOperators);
        for (const auto & bumpOperator : bumpOperators) {
          oops::Log::info() << "         Apply " << bumpOperator << std::endl;
          if (bumpOperator == "multiplyVbal") {
            this->multiplyVbal(atlasFieldSet.get());
          } else if (bumpOperator == "inverseMultiplyVbal") {
            this->inverseMultiplyVbal(atlasFieldSet.get());
          } else if (bumpOperator == "multiplyVbalAd") {
            this->multiplyVbalAd(atlasFieldSet.get());
          } else if (bumpOperator == "inverseMultiplyAd") {
            this->inverseMultiplyVbalAd(atlasFieldSet.get());
          } else if (bumpOperator == "multiplyStdDev") {
            this->multiplyStdDev(atlasFieldSet.get());
          } else if (bumpOperator == "inverseMultiplyStdDev") {
            this->inverseMultiplyStdDev(atlasFieldSet.get());
          } else if (bumpOperator == "multiplyNicas") {
            this->multiplyNicas(atlasFieldSet.get());
          } else {
              ABORT("Wrong bump operator: " + bumpOperator);
          }
        }

        // ATLAS transfer
        dxo.fromAtlas(atlasFieldSet.get());

        // Write file
        eckit::LocalConfiguration outputConf(appConf, "output");
        oops::Log::info() << "         Output file: " << outputConf << std::endl;
        dxo.write(outputConf);
      }
    }
  }

  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::trace() << "BUMP::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::addMember(const atlas::FieldSet & atlasFieldSet, const int & ie,
                              const int & iens) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_add_member_f90(keyBUMP_[jgrid], atlasFieldSet.get(), ie+1, iens);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::updateVbalCov(const Increment_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  dx.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_vbal_cov_f90(keyBUMP_[jgrid], atlasFieldSet->get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::updateVar(const Increment_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  dx.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_var_f90(keyBUMP_[jgrid], atlasFieldSet->get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::updateMom(const Increment_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  dx.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_mom_f90(keyBUMP_[jgrid], atlasFieldSet->get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::runDrivers() const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_run_drivers_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyVbal(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::inverseMultiplyVbal(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyVbalAd(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::inverseMultiplyVbalAd(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_ad_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyStdDev(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::inverseMultiplyStdDev(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::randomizeNicas(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_randomize_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyNicas(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyPsiChiToUV(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyPsiChiToUVAd(atlas::FieldSet * atlasFieldSet) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_ad_f90(keyBUMP_[jgrid], atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::getParameter(const std::string & param, Increment_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_get_parameter_f90(keyBUMP_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
  dx.fromAtlas(atlasFieldSet.get());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::setParameter(const std::string & param, const Increment_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  dx.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_set_parameter_f90(keyBUMP_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::partialDealloc() const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_partial_dealloc_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_BUMP_H_