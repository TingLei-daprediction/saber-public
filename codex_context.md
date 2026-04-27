# codex_context.md

## Scope
- Package: saber (MGBF components)
- Focus: memory usage and filtering paths in mgbf_lib and covariance wrappers.

## Key paths
- MGBF covariance wrapper: src/saber/mgbf/covariance/mgbf_covariance_mod.f90
- MGBF core library: src/saber/mgbf/mgbf_lib/*.f90
- Filtering implementation: src/saber/mgbf/mgbf_lib/mg_filtering.f90
- Parameters: src/saber/mgbf/mgbf_lib/mg_parameter.f90
- Internal state: src/saber/mgbf/mgbf_lib/mg_intstate.f90

## Current configuration assumptions
- l_loc = .false.
- l_vertical_filter = .true.
- l_mgbf_inhomogeneous = .false.
- filtering_fast_bkg(this) is the active filter path
- mgbf_line may be .false. in the fast-bkg flow

## Memory drivers (high-level)
- Large arrays scale with km_all (= km * n_ens), im/jm (plus halo), lm, and count of arrays.
- Additional multiplicative factor when nscale or nvargrp > 1.

## Large allocations (mgbf_lib)
- VALL/HALL: (km_all, im+2*hx, jm+2*hy)
- a_diff_f/a_diff_h/b_diff_f/b_diff_h: (km_all, im+2*hx, jm+2*hy)
- paspx4d/paspy4d and ssx4d/ssy4d: (lm, im+2*hx, jm+2*hy, 2)
- pasp3 (3x3 tensor per grid point): (3,3,im,jm,lm)
- hss3: (im, jm, lm, 6)
- qcols: (0:7, im, jm, lm)
- (localization only) w1_loc..w4_loc
- (setup only; can be scoped) weig_var: (km_all, im+2*hx, jm+2*hy, gm)

## filtering_fast_bkg usage
Used arrays in filtering_fast_bkg:
- paspx4d, paspy4d, ssx4d, ssy4d
- pasp1, ss1 (vertical filter)
- VALL, HALL

Not used by filtering_fast_bkg:
- vpasp2, hss2, dixs, diys
- vpasp3, hss3, qcols, dixs3, diys3, dizs3
- ss3
- pasp2/pasp3 only needed for radial/line paths outside fast-bkg

## Notes on mgbf_covariance_mod multiply
- work_mgbf is a full 3D buffer; work2d_mgbf is a packed 2D buffer.
- work_mgbf2 was removed by user; ensure no leftover deallocate or references.

## Cleanup candidates
- Consider guarding allocation of line-operator arrays when only filtering_fast_bkg is used.
- Consider scoping weig_var allocation to def_mg_weights and deallocating after use.
- Add deallocation of mg_parameter allocatables (ixm/jym/nxy/im0/jm0/Fimax/Fjmax/FimaxL/FjmaxL/zofis/isofz) if reinit happens.

## Jim new beta-function integration
Invocation hint: "continue our work on using Jim's new function".

Current patch status:
- Added MGBF-local copies of Jim's 1D calibrated beta-filter constants and routines in `src/saber/mgbf/mgbf_lib/jp_pbfil.f90`.
- `sres2_jim_new` and `sres3_jim_new` correspond to Jim's `sres2` and `sres3` in `../dr-jim/ybfil.f90`; these tables were corrected to match Jim's source values.
- Added type-bound declarations/interfaces in `src/saber/mgbf/mgbf_lib/mg_parameter.f90` for:
  - `rcalib1_jim_new`
  - `rbeta1_jim_new`, `rbeta3d_1_jim_new`
  - `rbeta1T_jim_new`, `rbeta3d_1T_jim_new`
- Added `paspx4d_jim_new` and `paspy4d_jim_new` in `src/saber/mgbf/mgbf_lib/mg_intstate.f90` with shape `(0:1,lm,x,y,2)`.
- Generation-1 Jim coefficients are computed from `paspx4d(:,:,:,1)` / `paspy4d(:,:,:,1)` by `rcalib1_jim_new`.
- Higher-generation Jim coefficients are transferred from generation 1, following the existing `paspx4d` / `paspy4d` pattern: `boco_2d` then `upsending_normalized`, separately for coefficient index `0` and `1`.
- Added missing type-bound procedure/interface for `filtering_fast_bkg_new_jim` in `mg_intstate.f90`.
- In `src/saber/mgbf/mgbf_lib/mg_filtering.f90`, `case(15)` now calls `filtering_fast_bkg_new_jim`.
- Inside `filtering_fast_bkg_new_jim`, horizontal old `rbeta` / `rbetaT` calls were replaced by Jim versions:
  - adjoint: `rbeta3d_1T_jim_new`
  - direct: `rbeta3d_1_jim_new`
  - using `paspx4d_jim_new` / `paspy4d_jim_new`
- After Jim's boundary-condition guidance, added local `FLIPT`/`FLIP` equivalents:
  - `rflip1T_jim_new`, `rflip3d_1T_jim_new`
  - `rflip1_jim_new`, `rflip3d_1_jim_new`
- `filtering_fast_bkg_new_jim` now calls `rflip3d_1T_jim_new` after the adjoint y/x halo side exchanges and calls `rflip3d_1_jim_new` after the direct y/x halo preparation, immediately before the direct beta filters.
- Boundary control uses the existing `Flwest/Fleast/Flsouth/Flnorth` flags. Boundary aspect inputs are taken from the corresponding edge values of `paspx4d` / `paspy4d`; if a boundary aspect is `0.0`, the flip routines implement Jim's truncate/zero-reflectivity option.
- The original `filtering_fast_bkg` path was intentionally left unchanged.

Verification so far:
- `git diff --check` passes.
- No compile/test run yet; no existing CMake/build tree was found under or above `saber-public`.
