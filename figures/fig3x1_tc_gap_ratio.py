# coding: utf-8
# manage data and fit
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, FixedLocator)
from matplotlib.font_manager import FontProperties
import math

pi=math.pi
print(pi)

params = {'axes.linewidth': 1.4,
         'axes.labelsize': 16,
         'axes.titlesize': 14,
         'lines.markeredgecolor': "black",
         'lines.linewidth': 1.5,
         'xtick.labelsize': 14,
         'ytick.labelsize': 14,
         "text.usetex": True,
         "font.family": "serif",
         "font.serif": ["Palatino"]
         }

plt.rcParams.update(params)


numcols=range(3)
print(len(numcols))
#fig, axs = plt.subplots(figsize=(3,3*(len(numcols))), nrows=len(numcols), ncols=1, sharex=True)
fig, axs = plt.subplots(figsize=(4,10), nrows=len(numcols), ncols=1, sharex=True)


markers=['o','^','s']
linedash=["--", "-.", ":"]

# CHOOSE BETWEEN TWO METHODS TO GENERATE TWO DIFFERENT PLOTS
method="SRC"
#method="BCSHF"
if(method == "BCSHF") :
     meth="FiniteT_BCS_HF"
     colors=['#92c5de','#0571b0']
     mcolors=colors
     pots=["N3LO","N3LO500_TBF_regfull"]
     pot_label=["EM","EM+3NF"]
     #ymaxi=[2.8,1.6,3];
     ymaxi=[3.0,1.7,3];
     egap=0
axis_spacing=[0.5,0.2,0.5];

if(method == "SRC") :
    meth="FiniteT_SRC"
    colors=['#92c5de','#0571b0','#1a9850']
    mcolors=["none","none","none"]
    pots=["N3LO","N3LO500_TBF_regfull","CDBONN"]
    pots=["N3LO","N3LO500_TBF_regfull"]
    pot_label=["EM","EM+3NF","CD-Bonn"]
    #ymaxi=[2.8,1.2,5];
    ymaxi=[3.0,1.7,5];
    egap=0.1

print(ymaxi)

# READ DATA
filename="gap_tc_data.dat"
filename_BCS_reference="../../BCS_free/N3LO/neumat/gaps/" + filename
for ix,potential in enumerate(pots) :
    print(potential)
    fname='../' + meth + '/' + potential + '/neumat/gaps/' + filename
    print(fname)
    kf, gap, tc, error_tc = np.loadtxt(fname,usecols=(0,1,2,3),unpack=True)
    print(kf,gap)
    error_gap=egap+np.zeros_like(gap)


    axs[0].plot(kf,gap,linewidth=0,marker=markers[ix],color=colors[ix],markeredgecolor="black",markerfacecolor=colors[ix],label=pot_label[ix],zorder=4)
    axs[0].errorbar(kf,gap,error_gap,linewidth=0,elinewidth=1.5,color=colors[ix],zorder=3)

    fname_fit_gap="../fits/fit_gap_" + potential + "_" + meth + ".dat"
    f_kf, f_gap = np.loadtxt(fname_fit_gap,usecols=(0,1),unpack=True)
    axs[0].plot(f_kf,f_gap,"-",color=colors[ix],zorder=2)

    axs[1].plot(kf,tc,linewidth=0,marker=markers[ix],color=colors[ix],markeredgecolor="black",markerfacecolor=colors[ix],label=pot_label[ix],zorder=4)
    axs[1].errorbar(kf,tc,error_tc,linewidth=0,elinewidth=1.5,color=colors[ix],zorder=3)
    fname_fit_tc="../fits/fit_tc_" + potential + "_" + meth + ".dat"
    if(os.path.exists(fname_fit_tc)) :
        f_kf, f_tc = np.loadtxt(fname_fit_tc,usecols=(0,1),unpack=True)
        axs[1].plot(f_kf,f_tc,"-",color=colors[ix],zorder=2)

    ratio=gap/tc
    error_ratio=ratio*np.sqrt( np.power(error_gap/gap,2) + np.power(error_tc/tc,2))
    print("ratio:",ratio)
    print("av_ratio",np.average(ratio[kf<1]))
    print("e_ratio:",error_ratio)
    print("max_e_ratio:",np.max(error_ratio[kf<1]))
    axs[2].plot(kf,ratio,linewidth=0,marker=markers[ix],color=mcolors[ix],label=pot_label[ix],markerfacecolor=colors[ix],zorder=4)
    axs[2].fill_between(kf,ratio - error_ratio,ratio + error_ratio,linewidth=0,color=colors[ix],alpha=0.6,zorder=3)

fname_fit_gap_bcs="../fits/fit_gap_N3LO_FiniteT_BCS.dat"
f_kf_bcs, f_gap_bcs = np.loadtxt(fname_fit_gap_bcs,usecols=(0,1),unpack=True)
axs[0].plot(f_kf_bcs,f_gap_bcs,"--",color="#984ea3",alpha=0.6,zorder=1,label="EM (BCS)")

fname_fit_tc_bcs="../fits/fit_tc_N3LO_FiniteT_BCS.dat"
f_kf_bcs, f_tc_bcs = np.loadtxt(fname_fit_tc_bcs,usecols=(0,1),unpack=True)
axs[1].plot(f_kf_bcs,f_tc_bcs,"--",color="#984ea3",alpha=0.6,zorder=1,label="EM (BCS)")


xxx=[0,0.5,1,10]
bcs_ratio=1.764*np.ones_like(xxx)
axs[2].plot(xxx,bcs_ratio,"--",color="#984ea3",zorder=1)
axs[2].text(0.08,1.82,"BCS",color="#984ea3",fontsize=14)

if(method == "SRC") :
    # THIS ADDS GAP VALUES OF REF 31 IN TOP PANEL
    fname="../Gandolfi2022/afdmc_swave.dat"
    kfmc, gapmc, egapmc = np.loadtxt(fname,usecols=(0,1,2),unpack=True)
    axs[0].errorbar(kfmc,gapmc,egapmc,linewidth=0,elinewidth=1.5,color="#fdbb84",zorder=-10)
    axs[0].plot(kfmc,gapmc,linewidth=0,marker="s",color="#fdbb84",markeredgecolor="black",markerfacecolor="#fdbb84",label="Ref. [30]",zorder=-10,markersize=4)

    # THIS ADDS VALUES OF RATIOS AS BANDS IN LOWEST PLOT
    w_ratio=2.7*np.ones_like(xxx)
    w_ratio_error=0.3*np.ones_like(w_ratio)
    axs[2].plot(xxx,w_ratio,linestyle="--",linewidth=1,color="#868686",zorder=1)
    axs[2].fill_between(xxx,w_ratio - w_ratio_error,w_ratio + w_ratio_error,linewidth=0,color="#868686",alpha=0.3,zorder=0)

    axs[2].text(0.08,3.1,"Ref. [41]",color="#868686",fontsize=14)

axs[0].set_ylabel(r'Gap, $\Delta_0$ [MeV]')
axs[1].set_ylabel(r'Critical temperature, $T_c$ [MeV]')
axs[2].set_ylabel(r'Ratio, $\Delta_0/T_c$')

axs[2].set_xlabel(r'Fermi momentum, $k_F$ [fm$^{-1}$]')

def kf2den(x) :
    return np.power(x,3)/3/np.power(pi,2)

def den2kf(x) :
    return np.power(3*np.power(pi,2)*x,1./3.)

for ix in numcols :
    axs[ix].set_ylim([0,ymaxi[ix]])
    axs[ix].set_xlim([0,1.5])
    ttt=np.arange(0.,1.5,0.2)
    axs[ix].set_xticks(ttt)
    #axs[ix].legend(prop={'size': 14},frameon=False,loc=2,handletextpad=0,borderaxespad=0.2,labelspacing=0.1,borderpad=0.0)
    #axs[ix].legend(prop={'size': 14},frameon=False,bbox_to_anchor=(0, 0.8),handletextpad=0,borderaxespad=0.2,labelspacing=0.1,borderpad=0.0)
    axs[ix].legend(prop={'size': 12},frameon=False,loc='lower right', bbox_to_anchor=(0.65, -0.02),handletextpad=0.1)

    #axs[ix].xaxis.set_tick_params(width=1.5)
    #axs[ix].yaxis.set_tick_params(width=1.5)

    axs[ix].tick_params(which='both',direction='in',top=False,right=True)
    axs[ix].xaxis.set_minor_locator(MultipleLocator(0.1))
    axs[ix].yaxis.set_major_locator(MultipleLocator(axis_spacing[ix]))
    axs[ix].yaxis.set_minor_locator(MultipleLocator(0.1))

maj_tic=[]
min_tic=[]
for decade in range(-4,1) :
    dec=np.power(10.,decade)
    maj_tic.append( dec )
    for unit in range(1,10) :
        min_tic.append( unit*dec )

print(maj_tic)
#print(min_tic)

for ix in numcols :
    axx = axs[ix]
    secax = axx.secondary_xaxis("top", functions=(kf2den,den2kf))
    secax.set_xticks(maj_tic,labels=["","0.001","0.01","0.1","1"])

    secax.tick_params(which='both',direction='in',top=True,right=True)
    secax.xaxis.set_minor_locator(FixedLocator(min_tic))
    #secax.ticklabel_format(axis='x',style='scientific',useMathText=True)

    if( ix== 0 ) :
        secax.set_xlabel(r'Density, $\rho$ [fm$^{-3}$]',labelpad=10)
    else :
        secax.xaxis.set_ticklabels([])



plt.tight_layout(pad=0.5)
fileout="gap_tc_panel_" + meth + ".pdf"
print(fileout)
fig.savefig(fileout)
