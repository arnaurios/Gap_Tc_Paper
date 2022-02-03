##############################################################################
# THIS CODE FITS THE DATA FOR BCS, BCS+HF OR SRC GAPS AND TCs FROM
# THE APPROPIATE ../method/potential/neumat/gaps/gap_tc_data.dat FILE
#
# IT OUTPUTS THE FIT PARAMETERS ON SCREEN AS WELL AS fit_gap_* and fit_tc_*
# FILES THAT HAVE 2 COLUMNS: kf AND fit data.
##############################################################################
# coding: utf-8
# manage data and fit
#import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from fitfunction import *

# first part with least squares
from scipy.optimize import curve_fit

# second part about ODR
from scipy.odr import ODR, Model, Data, RealData

from scipy.stats import t

method="SRC";
#method="BCS_HF"
#method="BCS"
#potential="N3LO";
potential="N3LO500_TBF_regfull";
#potential="CDBONN"
file="../FiniteT_" + method + "/" + potential + "/neumat/gaps/gap_tc_data.dat"

# STARTING PARAMETERS FOR FIT VALUES
beta00=[0.05,1,1.5,1.2,20]
print(method)
print(potential)
print(file)

# READING DATA
kf,gap = np.loadtxt(file,usecols=(0,1),unpack=True);

# ERRORS FOR GAPS
if(method=="BCS") : ergap=0.001
if(method=="BCS_HF") : ergap=0.001
if(method=="SRC") : ergap=0.1
egap = ergap + np.zeros_like(gap)

# MODEL PROPERTIES
num_pars = 5
num_data = len( kf )

# MESH OF FERMI MOMENTA WHERE FIT IS EVALUATED FOR FIGURES
xkf=np.linspace(0,2,200)

# CURVE FIT
popt, pcov = curve_fit(
    f=fitfunction,       # model function
    xdata=kf,   # DATA X
    ydata=gap,   # DATA Y
    sigma=egap,   # ERROR IN Y
    p0=beta00,    # INITIAL VALUES OF FITS
    maxfev=10000,
    bounds=( 0.,+np.inf ), # WE BOUND THE PARAMETERS TO BE POSITIVE
    absolute_sigma=True
)

# FIT GAP USING FIT CURVE
k0,k1,k2,k3,d0=popt
perr = np.sqrt(np.diag(pcov))
ek0,ek1,ek2,ek3,ed0= perr
print("Curve FIT: GAP")
print("d0 = %6.2f +/- %4.3f" % (d0, ed0))
print("k0 = %6.2f +/- %4.3f" % (k0, ek0))
print("k1 = %6.2f +/- %4.3f" % (k1, ek1))
print("k2 = %6.2f +/- %4.3f" % (k2, ek2))
print("k3 = %6.2f +/- %4.3f" % (k3, ek3))

print("k1/k3 = %6.2f" % (k1/k3))

gap_model=fitfunction(kf,k0,k1,k2,k3,d0)
rmax=np.abs( max(gap_model-gap))
print("r_max = %6.2f" % (rmax))

gap_model_extended=fitfunction(xkf,k0,k1,k2,k3,d0)
plt.plot(xkf,gap_model_extended)
plt.plot(kf,gap,'o')
plt.show()

################################################################################################################
# READING DATA FOR TC
kf,tc,etc = np.loadtxt(file,usecols=(0,2,3),unpack=True);

# FIT TC (WITH Y ERRORS)
popt, pcov = curve_fit(
    f=fitfunction,       # model function
    xdata=kf,   # DATA X
    ydata=tc,   # DATA Y
    sigma=etc,   # ERROR IN Y
    p0=beta00,    # INITIAL VALUES OF FITS
    maxfev=100000,
    bounds=( [0.,0.0,0.0,0.0,0.],+np.inf ),
    absolute_sigma=True
)

k0,k1,k2,k3,d0=popt
perr = np.sqrt(np.diag(pcov))
ek0, ek1, ek2, ek3, ed0 = perr
print("Curve FIT: TC")
print("d0 = %6.2f +/- %4.3f" % (d0, ed0))
print("k0 = %6.2f +/- %4.3f" % (k0, ek0))
print("k1 = %6.2f +/- %4.3f" % (k1, ek1))
print("k2 = %6.2f +/- %4.3f" % (k2, ek2))
print("k3 = %6.2f +/- %4.3f" % (k3, ek3))

print("k1/k3 = %6.2f" % (k1/k3))

tc_model=fitfunction(kf,k0,k1,k2,k3,d0)
rmax=np.abs( max(tc_model-tc))
print("r_max = %6.2f" % (rmax))


tc_model_extended=fitfunction(xkf,k0,k1,k2,k3,d0)

plt.plot(xkf,tc_model_extended)
plt.errorbar(kf,tc,yerr=etc)
plt.show()

# GAP FIT
outputfile="fit_gap_" + potential + "_FiniteT_" + method + ".dat"
header="#".ljust(6) + "kf [fm-1]".ljust(12) + "Gap Fit [MeV]".ljust(14)
with open(outputfile,"w+") as file_id :
    np.savetxt(file_id,[header],fmt="%s")
    data_to_write =np.column_stack( [ xkf,gap_model_extended ] )
    np.savetxt(file_id,data_to_write,fmt="%14.6E")

# TC FIT
outputfile="fit_tc_" + potential + "_FiniteT_" + method + ".dat"
header="#".ljust(6) + "kf [fm-1]".ljust(12) + "Tc Fit [MeV]".ljust(14)
with open(outputfile,"w+") as file_id :
    np.savetxt(file_id,[header],fmt="%s")
    data_to_write =np.column_stack( [ xkf,tc_model_extended ] )
    np.savetxt(file_id,data_to_write,fmt="%14.6E")
