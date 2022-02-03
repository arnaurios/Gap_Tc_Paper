# coding: utf-8
# manage data and fit
#import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from fithaensel import *

# first part with least squares
from scipy.optimize import curve_fit

# second part about ODR
from scipy.odr import ODR, Model, Data, RealData

from scipy.stats import t

method="SRC";
#method="BCS_HF"
method="BCS"
potential="N3LO";
#potential="N3LO500_TBF_regfull";
#potential="CDBONN"
file="../FiniteT_" + method + "/" + potential + "/neumat/gaps/gap_tc_data.dat"

beta00=[0.05,1.5,3,2,2]
print(method)
print(potential)
print(file)

# READING DATA
kf,gap = np.loadtxt(file,usecols=(0,1),unpack=True);
if(method=="BCS") : ergap=0.001
if(method=="BCS_HF") : ergap=0.001
if(method=="SRC") : ergap=0.2
egap = ergap + np.zeros_like(gap)
print(kf,gap,egap)

# MODEL PROPERTIES
num_pars = 5
num_data = len( kf )
print(num_pars)
print(num_data)
# convert provided value to a significance level (e.g. 95. -> 0.05), then calculate alpha
siglevel=68.
al = 1 - (1 - siglevel / 100) / 2
# student’s t test
tval = t.ppf(al, num_data - num_pars)

xkf=np.linspace(0,2,200)

#print(fithaensel(0.5,0.05,1.5,3,2,2))
# CURVE FIT
popt, pcov = curve_fit(
    f=fithaensel,       # model function
    xdata=kf,   # x data
     ydata=gap,   # y data
     sigma=egap,   # y data
     p0=beta00,    # initial value of the parameters
     maxfev=10000,
     bounds=( 0.,+np.inf ),
     absolute_sigma=True
)

# FIT GAP (WITH NO Y ERRORS) USING FIT CURVE
k1,k2,a,b,d0=popt
perr = np.sqrt(np.diag(pcov))
ek1, ek2, ea, eb, ed0 = perr
print("Curve FIT")
print("k1 = %6.3f +/- %4.3f" % (k1, ek1))
print("k2 = %6.3f +/- %4.3f" % (k2, ek2))
print("a = %6.3f +/- %4.3f" % (a, ea))
print("b = %6.3f +/- %4.3f" % (b, eb))
print("d0 = %6.3f +/- %4.3f" % (d0, ed0))

dfdp = derfithaensel(kf,k1,k2,a,b,d0)
#dfdp[:,0:4] =np.zeros_like(kf) # derfithaensel(kf,k0,k1,k2,k3,d0)
# process covarianvce matrix and derivatives
d=np.zeros_like(kf)
for ij in range(num_pars) :
    for ik in range(num_pars) :
        d += dfdp[:,ij] * dfdp[:,ik] * pcov[ij,ik]
tval=0
f_68int_model=tval*np.sqrt(d+egap)
gap_model=fithaensel(kf,k1,k2,a,b,d0)

gap_model_extended=fithaensel(xkf,k1,k2,a,b,d0)

#plt.fill_between(xkf,gap_model_extended-f_68int_model_extended,gap_model_extended+f_68int_model_extended)
plt.fill_between(kf,gap_model-f_68int_model,gap_model+f_68int_model)
plt.plot(xkf,gap_model_extended)
plt.plot(kf,gap,'o')
plt.show()

################################################################################################################
# READING DATA
kf,tc,etc = np.loadtxt(file,usecols=(0,2,3),unpack=True);
print(kf,tc,etc)

# FIT TC (WITH Y ERRORS)
popt, pcov = curve_fit(
    f=fithaensel,       # model function
    xdata=kf,   # x data
     ydata=tc,   # y data
     sigma=etc,   # y data
     p0=beta00,    # initial value of the parameters
     maxfev=100000,
     bounds=( [0.,0.1,0.5,0.1,0.],+np.inf ),
     absolute_sigma=True
)

k1,k2,a,b,d0=popt
perr = np.sqrt(np.diag(pcov))
ek1, ek2, ea, eb, ed0 = perr
print("Curve FIT")
print("k1 = %6.3f +/- %4.3f" % (k1, ek1))
print("k2 = %6.3f +/- %4.3f" % (k2, ek2))
print("a = %6.3f +/- %4.3f" % (a, ea))
print("b = %6.3f +/- %4.3f" % (b, eb))
print("d0 = %6.3f +/- %4.3f" % (d0, ed0))

dfdp = derfithaensel(kf,k1,k2,a,b,d0)
# process covarianvce matrix and derivatives
d=np.zeros_like(kf)
for ij in range(num_pars) :
    for ik in range(num_pars) :
        d += dfdp[:,ij] * dfdp[:,ik] * pcov[ij,ik]
tval=0
f_68int=tval*np.sqrt(etc+d)
tc_model=fithaensel(kf,k1,k2,a,b,d0)
tc_model_extended=fithaensel(xkf,k1,k2,a,b,d0)

plt.fill_between(kf,tc_model-f_68int,tc_model+f_68int)
plt.plot(xkf,tc_model_extended)
plt.errorbar(kf,tc,yerr=etc)
plt.show()
#print(f_68int)
#print(pcov)

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
