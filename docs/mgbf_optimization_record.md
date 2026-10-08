# MGBF performance optimizations: record for re-application

This file records the MGBF performance work on branch `feature/saber-mgbf-sdl2`
so that it can be carried to a new branch that starts without it. Each entry
says what changed, where, why, and what to check after porting.

Base before the optimizations: `584f5367` (merge of `develop`).

| # | Change | Commit | Files |
|---|---|---|---|
| 1 | Cached filter coefficients | `5aaf0715` | `jp_pbfil.f90`, `mg_filtering.f90`, `mg_parameter.f90` |
| 2 | Precision-independent level check, clearer names | `9b30c728` | same three files |
| 3 | Skip inter-generation transfers when `gm=1` | `73dbcf93` | `mg_filtering.f90` |
| 4 | Namelists read on rank 0 and broadcast | `d3fa1eaa` | `mg_namelist_io.f90` (new), `mg_parameter.f90`, `mgbf_covariance_mod.f90`, `mgbf_lib/CMakeLists.txt` |

To port 1–4 onto another branch:

```bash
git format-patch 584f5367..d3fa1eaa -- src/saber/mgbf
git am 000*.patch          # on the new branch; resolve conflicts if MGBF moved
```

If a change turns out to be wrong, revert its commit (`git revert <hash>`) and
update its entry here.

All paths below are relative to `src/saber/mgbf/`.

## 1. Cached filter coefficients (`5aaf0715`)

**Problem.** The beta-filter kernels recomputed the stencil bounds and weights
(`ceiling`/`floor` of the support, `(1-r^2)^p`) at every grid point on every
call, although they depend only on the filter configuration, which is fixed
for the whole run.

**Change.**
- `mgbf_lib/jp_pbfil.f90`: `vrbeta1`/`vrbeta1T` and `rbeta3d_1`/`rbeta3d_1T`
  take optional cache arguments `rebuild, gx_lo, gx_hi, weights` (and
  `same_levels` for the 3D versions). Pass all of them or none:
  - none: original behaviour;
  - `rebuild=.true.`: compute, store and apply;
  - `rebuild=.false.`: apply the stored stencil only.
- `rbeta3d_1*` also detect when `el` and `ss` are identical on every level and
  then keep a single level of the cache (`reuse_level1_coeffs`).
- `mgbf_lib/mg_filtering.f90`: `sup_vrbeta1_bkg`/`sup_vrbeta1T_bkg` (vertical)
  and `filtering_fast_bkg` (horizontal sweeps) own **thread-private** caches,
  allocated inside each `!$omp parallel` region, one per coefficient
  configuration.
- `mgbf_lib/mg_parameter.f90`: interface blocks updated for the new optional
  arguments.

**Assumptions to keep when porting.**
- A cache is valid for one coefficient configuration
  (`hx, lx, mx, el, ss, p, rmom2_1`). The caller is responsible for
  rebuilding it when any of these changes.
- Surface (`km2`) variables reuse level `lm` of the 3D configuration's cache.
  This is correct only while surface coefficients equal that level's.
- The cached path is meant to reproduce the uncached arithmetic exactly (same
  operations, same order). Re-check bitwise agreement against the uncached
  path after porting.

## 2. Precision-independent level check (`9b30c728`)

- The "same on every level" test compares `el`/`ss` byte-wise with
  `transfer(..., [0_int8])` instead of as `int64` words. This removes the
  64-bit-real requirement of #1 and still tells `+0.0` from `-0.0`.
- Renames: optional argument `same_levels` → `same_over_levels`; kernel-local
  `shared_levels` → `reuse_level1_coeffs`.
- No arithmetic change.

## 3. Skip inter-generation transfers when `gm=1` (`73dbcf93`)

**Problem.** In `filtering_fast_bkg` (`mgbf_proc=5`), `upsending` and
`downsending` perform their g1→g2 and g2→g1 steps outside their generation
loops. With a single generation, the code still ran adjoint interpolation,
boundary exchanges, upsend/downsend MPI messages and direct interpolation for
a generation that does not exist.

**Change.** Both calls are guarded with `gm>1`. For `gm=1`, `HALL` is set to
zero instead of upsending.

**Notes.**
- `gm` is the same on all ranks, so the skipped MPI traffic is skipped
  consistently.
- `gm>1` is unchanged. For `gm=1` the removed path only added an exact zero, so
  results match the baseline except that `-0.0` is no longer turned into `+0.0`.
- Scope: `filtering_fast_bkg` only. The other filtering procedures in
  `mg_filtering.f90` still call `upsending_all`/`downsending_all`
  unconditionally.

## Measured effect of 1–3 (Cactus, 2026-10-06)

| run | `mgbf::Covariance::multiply` before | after | speed-up |
|---|---|---|---|
| Parallel B, 10 groups × 196 ranks, OMP 4 | 396.4 s | 245.4 s | 1.62x |
| Serial ensemble, 1936 ranks | 337.2 s | 235.3 s | 1.43x |

Both over 102 B-multiplies (parallel: cross-rank average; serial: rank 0).
Configuration: `mgbf_line=.true.`, `mgbf_proc=5`, `gm=1`.

## 4. Namelists read on rank 0 and broadcast (`d3fa1eaa`, 2026-10-07)

**Status:** committed before being built or run. Verification below is still
open.

**Problem.** Every MPI rank opened and read the MGBF namelist files from the
file system: the SDL/VDL init namelist in `mgbf_covariance_mod.f90` and each
group namelist in `init_mg_parameter`, once per scale × variable group. With
2940 ranks this sends thousands of simultaneous open requests to the Lustre
metadata servers.

The MGBF constructor time went from 9.3 s (10 groups, 1960 ranks) to 12.2 s
(15 groups, 2940 ranks), cross-rank average. Its minimum stayed at about 4 s
and its imbalance was 80–93%: ranks were waiting, not computing.

**Change.**
- New module `mgbf_lib/mg_namelist_io.f90`, `read_namelist_lines(filename,
  comm, lines)`:
  - rank 0 reads the file into a character array and broadcasts it;
  - every rank then parses the text with `read(lines, nml=...)`, a namelist
    read from an internal file (Fortran 2003);
  - aborts if the file cannot be opened or a line exceeds `nml_line_len`
    (512) characters, so a long line is never truncated silently.
- `mgbf_lib/mg_parameter.f90` (`init_mg_parameter`) and
  `covariance/mgbf_covariance_mod.f90` (`create`) call it in place of
  `open`/`read`/`close`. The communicator is the MGBF communicator
  (`this%mpi_comm_comp`, `self%mp_comm_world`), so in Parallel B each group
  has one reader.
- `mgbf_lib/CMakeLists.txt`: `mg_namelist_io.f90` added to the parent sources,
  before `mg_parameter.f90`.

**Behaviour.** Parsed values are unchanged: same namelist text, same
defaults for variables the file does not set.

**To verify.**
- Build with both GNU and Intel.
- Compare the MGBF constructor timer and the analysis output with a run of the
  previous build.

## Not done yet: other file-system work in the constructor

Found in the same review, not changed:

- `mgbf_lib/mg_intstate.f90`, `def_mg_weights`: every rank opens unit 12 on
  `mgbf_tmpfile_<rank>.txt` and never writes to it. The file is closed only in
  the `l_mgbf_inhomogeneous` branch. In Parallel B, ranks with the same group
  rank in different groups create the same file name at the same moment.
  Candidate for deletion.
- Same routine: rank 0 of each MGBF communicator writes `latlon.txt`, so with
  Parallel B several groups write one file. Candidate for the `debug print`
  flag.
- `init_mg_parameter`: if `dir_coef_normalization` is set, every rank reads
  `profile_subdomain_XXXX.txt`, and in Parallel B each file is read once per
  group. Could be read once and scattered, or merged into one file.

## Related optimizations outside this repository

The Parallel B redistribution was fused in **OOPS**, branch
`feature/ensemble_read2changres_and_cached_interpolator`:
- `26371c13` fuse the per-colour `allToAllv` loop;
- `52c3ca2c` integer-overflow checks when scaling a routing by the level count.

Redistribution went from 348 s to 80 s with a bit-identical answer. A new
branch without recent optimizations needs these too, if it uses Parallel B.
