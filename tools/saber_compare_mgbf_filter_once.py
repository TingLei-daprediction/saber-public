#!/usr/bin/env python3
#
# (C) Copyright 2026 DOC/NOAA
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

"""
Equivalence check for the MGBF localization switch l_loc_filter_once (V2).

Compares the Dirac output files of a switch-off run (reference r) and a
switch-on run (result p) of the same build, value by value over every field:

    |p - r| <= rtol * |r| + atol * S_j

S_j is the absolute scale of output field j:

    S_j = max over points and levels of  sum_c |r_c|

where r_c are the outputs of switch-off runs with a Dirac in one variable only
(--contribution, one per Dirac variable). The output is linear in the Diracs, so
these are the contributions that add up to r, weights included; the bound
therefore stays meaningful where contributions cancel. Without --contribution,
S_j = max |r| over field j.

Values equal to the field's _FillValue must be identical. Any non-finite
reference value, result or bound is a failure.

Only the Python standard library is used: the files are netCDF classic format
(CDF-1 or 64-bit offset CDF-2), as written by oops util::writeFieldSet.

Call as:
  saber_compare_mgbf_filter_once.py --reference DIR --result DIR
                                    [--contribution DIR ...]
                                    [--rtol 1e-12] [--atol 1e-12]

Each DIR holds the run's output files (*.nc); files are paired by name.
Exit status 1 on any failure, 0 on success.
"""

import argparse
import glob
import math
import os
import struct
import sys

# netCDF classic type code -> (struct format, size in bytes)
NC_TYPES = {1: ('b', 1), 2: ('c', 1), 3: ('h', 2), 4: ('i', 4), 5: ('f', 4), 6: ('d', 8)}
NC_DIMENSION, NC_VARIABLE, NC_ATTRIBUTE = 10, 11, 12


class ClassicReader:
    """Minimal reader for netCDF classic files without record variables."""

    def __init__(self, path):
        with open(path, 'rb') as infile:
            self.data = infile.read()
        self.path = path
        self.pos = 0
        magic = self.data[:4]
        if magic[:3] != b'CDF' or magic[3] not in (1, 2):
            raise ValueError('%s: not a netCDF classic (CDF-1/CDF-2) file' % path)
        self.offset_size = 8 if magic[3] == 2 else 4
        self.pos = 4
        self._int()  # numrecs
        self.dims = self._dim_list()
        self._att_list()  # global attributes
        self.vars = self._var_list()

    def _int(self):
        value = struct.unpack_from('>i', self.data, self.pos)[0]
        self.pos += 4
        return value

    def _offset(self):
        fmt = '>q' if self.offset_size == 8 else '>i'
        value = struct.unpack_from(fmt, self.data, self.pos)[0]
        self.pos += self.offset_size
        return value

    def _name(self):
        length = self._int()
        name = self.data[self.pos:self.pos + length].decode('utf-8')
        self.pos += length + (-length % 4)
        return name

    def _list_header(self, tag):
        kind, count = self._int(), self._int()
        if kind == 0 and count == 0:
            return 0
        if kind != tag:
            raise ValueError('%s: unexpected header tag %d' % (self.path, kind))
        return count

    def _dim_list(self):
        return [(self._name(), self._int()) for _ in range(self._list_header(NC_DIMENSION))]

    def _att_list(self):
        atts = {}
        for _ in range(self._list_header(NC_ATTRIBUTE)):
            name = self._name()
            nc_type, count = self._int(), self._int()
            fmt, size = NC_TYPES[nc_type]
            raw = self.data[self.pos:self.pos + count * size]
            self.pos += count * size + (-(count * size) % 4)
            atts[name] = raw if fmt == 'c' else struct.unpack('>%d%s' % (count, fmt), raw)
        return atts

    def _var_list(self):
        variables = {}
        for _ in range(self._list_header(NC_VARIABLE)):
            name = self._name()
            dimids = [self._int() for _ in range(self._int())]
            atts = self._att_list()
            nc_type = self._int()
            self._int()  # vsize
            begin = self._offset()
            if any(self.dims[d][1] == 0 for d in dimids):
                raise ValueError('%s: record variable %s is not supported' % (self.path, name))
            variables[name] = (dimids, atts, nc_type, begin)
        return variables

    def values(self, name):
        dimids, _, nc_type, begin = self.vars[name]
        count = 1
        for d in dimids:
            count *= self.dims[d][1]
        fmt, size = NC_TYPES[nc_type]
        return struct.unpack_from('>%d%s' % (count, fmt), self.data, begin)

    def fill_value(self, name):
        fill = self.vars[name][1].get('_FillValue')
        return fill[0] if fill else None

    def numeric_names(self):
        return [n for n, v in self.vars.items() if NC_TYPES[v[2]][0] != 'c']


def compare_file(name, ref_dir, res_dir, contrib_dirs, rtol, atol):
    """Compare one output file; return the number of failures."""
    ref = ClassicReader(os.path.join(ref_dir, name))
    res = ClassicReader(os.path.join(res_dir, name))
    contribs = [ClassicReader(os.path.join(d, name)) for d in contrib_dirs]
    failures = 0

    if sorted(ref.numeric_names()) != sorted(res.numeric_names()):
        print('%s: different variables %s vs %s'
              % (name, sorted(ref.numeric_names()), sorted(res.numeric_names())))
        return 1

    for var in sorted(ref.numeric_names()):
        r, p = ref.values(var), res.values(var)
        fill = ref.fill_value(var)
        if len(r) != len(p):
            print('%s/%s: different sizes %d vs %d' % (name, var, len(r), len(p)))
            failures += 1
            continue

        # Absolute scale S_j
        if contribs:
            parts = [c.values(var) for c in contribs]
            sums = [sum(abs(part[k]) for part in parts) for k in range(len(r))]
        else:
            sums = [abs(value) for value in r]
        valid = [s for k, s in enumerate(sums) if r[k] != fill]
        # check every term: max() can silently skip a NaN
        if not all(math.isfinite(s) for s in valid):
            print('%s/%s: non-finite value in the absolute scale' % (name, var))
            failures += 1
            continue
        scale = max(valid, default=0.0)

        nbad, worst = 0, 0.0
        for k, (rk, pk) in enumerate(zip(r, p)):
            if fill is not None and (rk == fill or pk == fill):
                if rk != pk:
                    nbad += 1
                continue
            bound = rtol * abs(rk) + atol * scale
            if not (math.isfinite(rk) and math.isfinite(pk) and math.isfinite(bound)):
                nbad += 1
                continue
            diff = abs(pk - rk)
            if diff > bound:
                nbad += 1
            if bound > 0.0:
                worst = max(worst, diff / bound)
        print('%s/%s: %d values, %d failure(s), worst |diff|/bound = %.3e, S = %.6e'
              % (name, var, len(r), nbad, worst, scale))
        failures += nbad
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--reference', required=True, help='switch-off run output directory')
    parser.add_argument('--result', required=True, help='switch-on run output directory')
    parser.add_argument('--contribution', action='append', default=[],
                        help='switch-off single-Dirac run output directory (repeatable)')
    parser.add_argument('--rtol', type=float, default=1.0e-12)
    parser.add_argument('--atol', type=float, default=1.0e-12)
    args = parser.parse_args()

    names = sorted(os.path.basename(f) for f in glob.glob(os.path.join(args.reference, '*.nc')))
    if not names:
        print('No *.nc output files in ' + args.reference)
        return 1

    failures = 0
    for name in names:
        for directory in [args.result] + args.contribution:
            if not os.path.exists(os.path.join(directory, name)):
                print('Missing %s in %s' % (name, directory))
                failures += 1
        if failures:
            continue
        failures += compare_file(name, args.reference, args.result, args.contribution,
                                 args.rtol, args.atol)

    print('Total failures: %d (rtol=%g, atol=%g)' % (failures, args.rtol, args.atol))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
