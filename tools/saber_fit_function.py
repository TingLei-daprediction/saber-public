#!/usr/bin/env python3

import argparse
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("srcdir", help="SABER source directory")
args = parser.parse_args()

# General parameters
nnd = 51
dnd = 1.0/float(nnd-1)
nd = np.linspace(0, (nnd-1)*dnd, nnd)
epsabs_hor = 1.0e-2
epsabs_ver = 1.0e-4

# Parameters
axis_invmin = 0.2
axis_invmax = 0.9
daxis_inv = 0.1
naxis_inv = int((axis_invmax-axis_invmin)/daxis_inv+1.0e-6)+1
axis_inv = np.linspace(axis_invmin, axis_invmax, naxis_inv)
run_horizontal = True
run_vertical = True

# Distance
axis = np.zeros(nnd)
for ind in range(0,nnd):
   axis[ind] = float(ind)/float(nnd-1)

# Functions
def S(r):
   if np.abs(r) <= 0.5:
      return 1.0-(2.0*np.abs(r))
   else:
      return 0.0

def S_hor(x,y):
   r = np.sqrt(x**2+y**2)
   return S(r)

def S_ver(z):
   return S(z)

# Initialize arrays
func_sqrt_hor = np.zeros((nnd))
func_hor = np.zeros((nnd))
func_inv_hor = np.zeros((naxis_inv))
func_sqrt_ver = np.zeros((nnd))
func_ver = np.zeros((nnd))
func_inv_ver = np.zeros((naxis_inv))
scaled_axis = np.zeros((nnd))

if run_horizontal:
   for ind in range(0, nnd):
      # Square-root function
      func_sqrt_hor[ind] = S(axis[ind])

      # Horizontal integration (2D)
      f = lambda  y, x: S_hor(x,y)*S_hor(axis[ind]-x,y)
      fint = integrate.dblquad(f, -0.5, 0.5, lambda x: -0.5, lambda x: 0.5, epsabs = epsabs_hor)
      func_hor[ind] = fint[0]
      if ind == 0:
         norm = func_hor[ind]
      func_hor[ind] = func_hor[ind]/norm

   # Inverse function
   for iaxis_inv in range(0, naxis_inv):
      func_inv_hor[iaxis_inv] = 1.0
      for ind in range(0, nnd-1):
         if func_hor[ind]>axis_inv[iaxis_inv] and func_hor[ind+1]<axis_inv[iaxis_inv]:
            A = (func_hor[ind]-func_hor[ind+1])/(axis[ind]-axis[ind+1])
            B = func_hor[ind]-A*axis[ind]
            func_inv_hor[iaxis_inv] = (axis_inv[iaxis_inv]-B)/A
            break

   if True:
      # Plot curves
      fig, ax = plt.subplots(ncols=2, figsize=(14,7))
      ax[0].set_xlim([0,1.0])
      ax[0].set_ylim([0,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(axis, func_sqrt_hor, 'k')
      ax[1].set_xlim([0,1.0])
      ax[1].set_ylim([0,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(axis, func_hor)
      plt.savefig("fit_hor.jpg", format="jpg", dpi=300)
      plt.close()

if run_vertical:
   for ind in range(0, nnd):
      # Square-root function
      func_sqrt_ver[ind] = S_ver(axis[ind])/S_ver(0)

      # Vertical integration (1D)
      f = lambda  z: S_ver(z)*S_ver(axis[ind]-z)
      fint = integrate.quad(f, -0.5, 0.5, epsabs = epsabs_ver)
      func_ver[ind] = fint[0]
      if ind == 0:
         norm = func_ver[ind]
      func_ver[ind] = func_ver[ind]/norm

   # Inverse function
   for iaxis_inv in range(0, naxis_inv):
      func_inv_ver[iaxis_inv] = 1.0
      for ind in range(0, nnd-1):
         if func_ver[ind]>axis_inv[iaxis_inv] and func_ver[ind+1]<axis_inv[iaxis_inv]:
            A = (func_ver[ind]-func_ver[ind+1])/(axis[ind]-axis[ind+1])
            B = func_ver[ind]-A*axis[ind]
            func_inv_ver[iaxis_inv] = (axis_inv[iaxis_inv]-B)/A
            break

   if True:
      # Plot curves
      fig, ax = plt.subplots(ncols=2, figsize=(14,7))
      ax[0].set_xlim([0,1.0])
      ax[0].set_ylim([0,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(axis, func_sqrt_ver, 'k')
      ax[1].set_xlim([0,1.0])
      ax[1].set_ylim([0,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(axis, func_ver)
      plt.savefig("fit_ver.jpg", format="jpg", dpi=300)
      plt.close()

if run_horizontal and run_vertical:
   # Open file
   file = open(args.srcdir + "/src/saber/bump/tools_gc99.fypp", "w")

   # Write file
   file.write("#:include 'instrumentation.fypp'\n")
   file.write("#:include 'generics.fypp'\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Module: tools_gc99\n")
   file.write("!> Gaspari and Cohn (1999)-inspired functions and their square-roots\n")
   file.write("! Author: Benjamin Menetrier\n")
   file.write("! Licensing: this code is distributed under the CeCILL-C license\n")
   file.write("! Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT\n")
   file.write("! WARNING: this module is generated by the python script\n")
   file.write("!            tools/saber_fit_function.py\n")
   file.write("!          to modify this module, update and rerun the python script\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("module tools_gc99\n")
   file.write("\n")
   file.write("use tools_const, only: zero,half,one,two\n")
   file.write("use tools_kinds, only: kind_real\n")
   file.write("use tools_repro, only: eq,inf,infeq\n")
   file.write("use type_mpl, only: mpl_type\n")
   file.write("@:use_probe()\n")
   file.write("\n")
   file.write("implicit none\n")
   file.write("\n")
   file.write("! Public parameters\n")
   file.write("integer,parameter :: nnd = " + str(nnd) + "\n")
   file.write("integer,parameter :: naxis_inv = " + str(naxis_inv) + "\n")
   file.write("real(kind_real),parameter :: ndmin = %.8f_kind_real\n" % (min(nd)))
   file.write("real(kind_real),parameter :: ndmax = %.8f_kind_real\n" % (max(nd)))
   file.write("real(kind_real),parameter :: dnd = %.8f_kind_real\n" % (dnd))
   file.write("real(kind_real),parameter :: axis_invmin = %.8f_kind_real\n" % (axis_invmin))
   file.write("real(kind_real),parameter :: axis_invmax = %.8f_kind_real\n" % (axis_invmax))
   file.write("real(kind_real),parameter :: axis_inv(naxis_inv) = (/ &\n")
   for iaxis_inv in range(0, naxis_inv):
      if iaxis_inv != naxis_inv-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (axis_inv[iaxis_inv]) + suffix + "\n")
   file.write("real(kind_real),parameter :: func_inv_hor(naxis_inv) = (/ &\n")
   for iaxis_inv in range(0, naxis_inv):
      if iaxis_inv != naxis_inv-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (func_inv_hor[iaxis_inv]) + suffix + "\n")
   file.write("real(kind_real),parameter :: func_hor(nnd) = (/ &\n")
   for ind in range(0, nnd):
      if ind != nnd-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (func_hor[ind]) + suffix + "\n")
   file.write("real(kind_real),parameter :: func_inv_ver(naxis_inv) = (/ &\n")
   for iaxis_inv in range(0, naxis_inv):
      if iaxis_inv != naxis_inv-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (func_inv_ver[iaxis_inv]) + suffix + "\n")
   file.write("real(kind_real),parameter :: func_ver(nnd) = (/ &\n")
   for ind in range(0, nnd):
      if ind != nnd-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (func_ver[ind]) + suffix + "\n")
   file.write("\n")
   file.write("interface fit_func\n")
   file.write("   module procedure gc99_fit_func\n")
   file.write("end interface\n")
   file.write("interface fit_func_sqrt\n")
   file.write("   module procedure gc99_fit_func_sqrt\n")
   file.write("end interface\n")
   file.write("\n")
   file.write("private\n")
   file.write("public :: naxis_inv,axis_inv,axis_invmin,axis_invmax\n")
   file.write("public :: func_inv_hor,func_inv_ver\n")
   file.write("public :: fit_func,fit_func_sqrt\n")
   file.write("\n")
   file.write("contains\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func\n")
   file.write("!> Fit function\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func(mpl,dir,nd) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("character(len=*),intent(in) :: dir  !< Direction\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: indm,indp\n")
   file.write("real(kind_real) :: rndm,rndp\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (inf(nd,zero)) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("\n")
   file.write("! Initialization\n")
   file.write("value = zero\n")
   file.write("\n")
   file.write("if (eq(nd,zero)) then\n")
   file.write("   ! Origin\n")
   file.write("   value = one\n")
   file.write("elseif (infeq(nd,one)) then\n")
   file.write("   ! Inside support\n")
   file.write("\n")
   file.write("   ! Indices\n")
   file.write("   indm = floor(nd/dnd)+1\n")
   file.write("   if (indm==nnd) then\n")
   file.write("      indp = indm\n")
   file.write("   else\n")
   file.write("      indp = indm+1\n")
   file.write("   end if\n")
   file.write("\n")
   file.write("   ! Coefficients\n")
   file.write("   if (indm==nnd) then\n")
   file.write("      rndm = one\n")
   file.write("   else\n")
   file.write("      rndm = real(indp-1,kind_real)-nd/dnd\n")
   file.write("   end if\n")
   file.write("   rndp = (one-rndm)\n")
   file.write("\n")
   file.write("   ! Interpolated value\n")
   file.write("   if (dir=='hor') then\n")
   file.write("      ! Horizontal fit function\n")
   file.write("      value = rndm*func_hor(indm)+rndp*func_hor(indp)\n")
   file.write("   elseif (dir=='ver') then\n")
   file.write("      ! Vertical fit function\n")
   file.write("      value = rndm*func_ver(indm)+rndp*func_ver(indp)\n")
   file.write("   else\n")
   file.write("      call mpl%abort('${subr}$','wrong direction: '//dir)\n")
   file.write("   end if\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func_sqrt\n")
   file.write("!> Fit function function square-root\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func_sqrt(mpl,nd) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func_sqrt)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (inf(nd,zero)) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("\n")
   file.write("! Initialization\n")
   file.write("value = zero\n")
   file.write("\n")
   file.write("if (eq(nd,zero)) then\n")
   file.write("   ! Origin\n")
   file.write("   value = one\n")
   file.write("elseif (infeq(nd,half)) then\n")
   file.write("   ! Inside support\n")
   file.write("   value = one-(two*nd)\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_sqrt\n")
   file.write("\n")
   file.write("end module tools_gc99\n")

   # Close file
   file.close()
