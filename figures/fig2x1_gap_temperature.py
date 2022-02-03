 # coding: utf-8
# manage data and fit
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.font_manager import FontProperties

params = {'axes.linewidth': 1.4,
          'axes.labelsize': 16,
          'axes.titlesize': 14,
          'lines.markeredgecolor': "black",
          'lines.linewidth': 1.5,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          "text.usetex": True,
          "font.family": "serif",
          "font.serif": ["Palatino"],
          "axes.axisbelow": True
         }

plt.rcParams.update(params)




numcols=range(2)
#fig, axs = plt.subplots(figsize=(3,3*(len(numcols))), nrows=len(numcols), ncols=1, sharex=True)
fig, axs = plt.subplots(figsize=(7,4), nrows=1, ncols=2, sharey=True)

# CHOOSE BETWEEN TWO METHODS TO GENERATE TWO DIFFERENT PLOTS
colors=['#92c5de','#0571b0']
mcolors=colors
pots=["N3LO","N3LO500_TBF_regfull"]
pot_label=["EM","EM+3NF"]

markers=['o','^','s']
linedash=["-", "--"]

# READ DATA
rho="0.04"
kf="1.058019"

file="gap_kf_1S0.dat"
dir1_BCS="../FiniteT_BCS_HF/"
dir2_BCS="/neumat/gaps/kf_" + kf + "/"

dir1_SRC="../FiniteT_SRC/"
dir2_SRC="/neumat/gaps/" + rho + "/"



for ix,potential in enumerate(pots) :
    filename=dir1_BCS + potential + dir2_BCS + file
    print(filename)
    rrr, temp, kkf, gap = np.loadtxt(filename,usecols=(0,1,2,3),unpack=True)
    ttt=temp[gap>2e-6]
    ggg=gap[gap>2e-6]
    axs[0].plot(ttt,ggg,linewidth=0.,marker=markers[ix],color=colors[ix],markeredgecolor="black",markerfacecolor=colors[ix],label=pot_label[ix],zorder=10,clip_on=False)
    axs[0].plot(ttt,ggg,linedash[ix],color=colors[ix],zorder=9,clip_on=False)

    axs[0].set_title(r"BCS+HF, $\rho=$" + rho + " fm$^{-3}$")

    #axs[0].fill_between(kf,gap-error_gap,gap+error_gap,linewidth=0,color=colors[ix],alpha=0.3,zorder=3)
    #axs[0].errorbar(kf,gap,error_gap,linewidth=0,elinewidth=1.5,color=colors[ix],zorder=3)

    #fname_fit_gap="../fits/fit_gap_" + potential + "_" + meth + ".dat"
    #f_kf, f_gap = np.loadtxt(fname_fit_gap,usecols=(0,1),unpack=True)
    #axs[0].plot(f_kf,f_gap,"-",color=colors[ix],zorder=2)

    #fname_fit_gap2="../fits/fit2_gap_" + potential + "_" + meth + ".dat"
    #if(os.path.exists(fname_fit_gap2)) :
    #    f2_kf, f2_gap = np.loadtxt(fname_fit_gap2,usecols=(0,1),unpack=True)
    #    axs[0].plot(f2_kf,f2_gap,"--",color=colors[ix],zorder=2)

    filename=dir1_SRC + potential + dir2_SRC + file
    print(filename)
    kkf, temp, gap = np.loadtxt(filename,usecols=(0,1,2),unpack=True)
    glim=[2e-6,1e-6]
    ttt=temp[gap>glim[ix]]
    ggg=gap[gap>glim[ix]]

    axs[1].plot(ttt,ggg,linewidth=0.,marker=markers[ix],color=colors[ix],markeredgecolor="black",markerfacecolor=colors[ix],label=pot_label[ix],zorder=10,clip_on=False)
    axs[1].plot(ttt,ggg,linedash[ix],color=colors[ix],zorder=9,clip_on=False)
    axs[1].set_axisbelow(True)
    axs[1].set_title(r"SRC, $\rho=$" + rho + " fm$^{-3}$")

axs[0].set_ylabel(r'Gap, $\Delta(T)$ [MeV]')
for ix in numcols :
    axs[ix].set_xlabel(r'Temperature, $T$ [MeV]')
    axs[ix].set_xlim([0,1.3])
    axs[ix].set_ylim([0,2.1])
    axs[ix].legend(prop={'size': 12},frameon=False,loc=1,handletextpad=0.)

    axs[ix].tick_params(which='both',direction='in',top=True,right=True)
    axs[ix].xaxis.set_minor_locator(AutoMinorLocator(2))
    axs[ix].yaxis.set_minor_locator(AutoMinorLocator(5))
    ttt=np.arange(0.,1.3,0.2)
    axs[ix].set_xticks(ttt)

#plt.tight_layout(pad=0.1)
plt.tight_layout()
fileout='gap_temperature.pdf'
fig.savefig(fileout)
