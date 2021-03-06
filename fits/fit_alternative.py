# coding: utf-8

##############################################################################
# THIS CODE FITS THE DATA FOR BCS, BCS+HF OR SRC GAPS AND TCs FROM
# THE APPROPIATE ../method/potential/neumat/gaps/gap_tc_data.dat FILE
#
# THIS FITS THE FUNCTION IN THE PAPER
#
# IT OUTPUTS THE FIT PARAMETERS ON SCREEN AS WELL AS fit_gap_* and fit_tc_*
# FILES THAT HAVE 2 COLUMNS: kf AND fit data.
##############################################################################
import numpy as np

import matplotlib.pyplot as plt
from fithaensel import *

# first part with least squares
from scipy.optimize import curve_fit

method="SRC";
#method="BCS_HF"
#method="BCS"
potential="N3LO";
#potential="N3LO500_TBF_regfull";
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
    f=fithaensel,       # model function
    xdata=kf,   # DATA X
    ydata=gap,   # DATA Y
    sigma=egap,   # ERROR IN Y
    p0=beta00,    # INITIAL VALUES OF FITS
    maxfev=10000,
    bounds=( 0.,+np.inf ), # WE BOUND THE PARAMETERS TO BE POSITIVE
    absolute_sigma=True
)

# FIT GAP USING FIT CURVE
k0,k1,aa,bb,d0=popt
perr = np.sqrt(np.diag(pcov))
ek0,ek1,eaa,ebb,ed0= perr
print("Curve FIT: GAP")
print("d0 = %6.2f +/- %4.3f" % (d0, ed0))
print("k0 = %6.2f +/- %4.3f" % (k0, ek0))
print("k1 = %6.2f +/- %4.3f" % (k1, ek1))
print("aa = %6.2f +/- %4.3f" % (aa, eaa))
print("bb = %6.2f +/- %4.3f" % (bb, ebb))

gap_model=fithaensel(kf,k0,k1,aa,bb,d0)
rmax=np.abs( max(gap_model-gap))
print("r_max = %6.2f" % (rmax))

gap_model_extended=fithaensel(xkf,k0,k1,aa,bb,d0)
plt.plot(xkf,gap_model_extended)
plt.plot(kf,gap,'o')
plt.show()

################################################################################################################
# READING DATA FOR TC
kf,tc,etc = np.loadtxt(file,usecols=(0,2,3),unpack=True);

# FIT TC (WITH Y ERRORS)
popt, pcov = curve_fit(
    f=fithaensel,       # model function
    xdata=kf,   # DATA X
    ydata=tc,   # DATA Y
    sigma=etc,   # ERROR IN Y
    p0=beta00,    # INITIAL VALUES OF FITS
    maxfev=100000,
    bounds=( [0.,0.0,0.0,0.0,0.],+np.inf ),
    absolute_sigma=True
)

k0,k1,aa,bb,d0=popt
perr = np.sqrt(np.diag(pcov))
ek0,ek1,eaa,ebb,ed0= perr
print("Curve FIT: GAP")
print("d0 = %6.2f +/- %4.3f" % (d0, ed0))
print("k0 = %6.2f +/- %4.3f" % (k0, ek0))
print("k1 = %6.2f +/- %4.3f" % (k1, ek1))
print("aa = %6.2f +/- %4.3f" % (aa, eaa))
print("bb = %6.2f +/- %4.3f" % (bb, ebb))

tc_model=fithaensel(kf,k0,k1,aa,bb,d0)
rmax=np.abs( max(tc_model-tc))
print("r_max = %6.2f" % (rmax))

tc_model_extended=fithaensel(xkf,k0,k1,aa,bb,d0)

plt.plot(xkf,tc_model_extended)
plt.errorbar(kf,tc,yerr=etc)
plt.show()

# GAP FIT
outputfile="fit2_gap_" + potential + "_FiniteT_" + method + ".dat"
header="#".ljust(6) + "kf [fm-1]".ljust(12) + "gap".ljust(14)
with open(outputfile,"w+") as file_id :
    np.savetxt(file_id,[header],fmt="%s")
    data_to_write =np.column_stack( [ xkf,gap_model_extended ] )
    np.savetxt(file_id,data_to_write,fmt="%14.6E")

# TC FIT
outputfile="fit2_tc_" + potential + "_FiniteT_" + method + ".dat"
header="#".ljust(6) + "kf [fm-1]".ljust(12) + "gap".ljust(14)
with open(outputfile,"w+") as file_id :
    np.savetxt(file_id,[header],fmt="%s")
    data_to_write =np.column_stack( [ xkf,tc_model_extended ] )
    np.savetxt(file_id,data_to_write,fmt="%14.6E")
