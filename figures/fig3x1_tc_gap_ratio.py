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
         "font.serif": ["Palatino"]
         }

plt.rcParams.update(params)


numcols=range(3)
print(len(numcols))
#fig, axs = plt.subplots(figsize=(3,3*(len(numcols))), nrows=len(numcols), ncols=1, sharex=True)
fig, axs = plt.subplots(figsize=(4,10), nrows=len(numcols), ncols=1, sharex=True)

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


markers=['o','^','s']
linedash=["--", "-.", ":"]
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
    #axs[0].fill_between(kf,gap-error_gap,gap+error_gap,linewidth=0,color=colors[ix],alpha=0.3,zorder=3)
    axs[0].errorbar(kf,gap,error_gap,linewidth=0,elinewidth=1.5,color=colors[ix],zorder=3)

    fname_fit_gap="../fits/fit_gap_" + potential + "_" + meth + ".dat"
    f_kf, f_gap = np.loadtxt(fname_fit_gap,usecols=(0,1),unpack=True)
    axs[0].plot(f_kf,f_gap,"-",color=colors[ix],zorder=2)

    #fname_fit_gap2="../fits/fit2_gap_" + potential + "_" + meth + ".dat"
    #if(os.path.exists(fname_fit_gap2)) :
    #    f2_kf, f2_gap = np.loadtxt(fname_fit_gap2,usecols=(0,1),unpack=True)
    #    axs[0].plot(f2_kf,f2_gap,"--",color=colors[ix],zorder=2)


    axs[1].plot(kf,tc,linewidth=0,marker=markers[ix],color=colors[ix],markeredgecolor="black",markerfacecolor=colors[ix],label=pot_label[ix],zorder=4)
    #axs[1].fill_between(kf,tc-error_tc,tc+error_tc,linewidth=0,color=colors[ix],alpha=0.3,zorder=3)
    axs[1].errorbar(kf,tc,error_tc,linewidth=0,elinewidth=1.5,color=colors[ix],zorder=3)
    fname_fit_tc="../fits/fit_tc_" + potential + "_" + meth + ".dat"
    if(os.path.exists(fname_fit_tc)) :
        f_kf, f_tc = np.loadtxt(fname_fit_tc,usecols=(0,1),unpack=True)
        axs[1].plot(f_kf,f_tc,"-",color=colors[ix],zorder=2)

    #fname_fit_tc2="../fits/fit2_tc_" + potential + "_" + meth + ".dat"
    #if(os.path.exists(fname_fit_tc2)) :
    #    f2_kf, f2_tc = np.loadtxt(fname_fit_tc2,usecols=(0,1),unpack=True)
    #    axs[1].plot(f2_kf,f2_tc,"--",color=colors[ix],zorder=2)

    ratio=gap/tc
    error_ratio=ratio*np.sqrt( np.power(error_gap/gap,2) + np.power(error_tc/tc,2))
    print("ratio:",ratio)
    print("av_ratio",np.average(ratio[kf<1]))
    print("e_ratio:",error_ratio)
    print("max_e_ratio:",np.max(error_ratio[kf<1]))
    axs[2].plot(kf,ratio,linewidth=0,marker=markers[ix],color=mcolors[ix],label=pot_label[ix],markerfacecolor=colors[ix],zorder=4)
    #axs[2].errorbar(kf,ratio,error_ratio,linestyle=linedash[ix],color=colors[ix],zorder=3)
    axs[2].fill_between(kf,ratio - error_ratio,ratio + error_ratio,linewidth=0,color=colors[ix],alpha=0.6,zorder=3)
    #axs[2].plot(f_kf,f_gap/f_tc,linewidth=1,color=colors[ix])

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
    w_ratio=2.7*np.ones_like(xxx)
    w_ratio_error=0.3*np.ones_like(w_ratio)
    axs[2].plot(xxx,w_ratio,linestyle="--",linewidth=1,color="#868686",zorder=1)
    axs[2].fill_between(xxx,w_ratio - w_ratio_error,w_ratio + w_ratio_error,linewidth=0,color="#868686",alpha=0.3,zorder=0)

    axs[2].text(0.08,3.1,"Ref. [35]",color="#868686",fontsize=14)



axs[0].set_ylabel(r'Gap, $\Delta_0$ [MeV]')
axs[1].set_ylabel(r'Critical temperature, $T_c$ [MeV]')
axs[2].set_ylabel(r'Ratio, $\Delta_0/T_c$')

axs[2].set_xlabel(r'Fermi momentum, $k_F$ [fm$^{-1}$]')

for ix in numcols :
    axs[ix].set_ylim([0,ymaxi[ix]])
    axs[ix].set_xlim([0,1.5])
    ttt=np.arange(0.,1.5,0.2)
    axs[ix].set_xticks(ttt)
    #axs[ix].legend(prop={'size': 14},frameon=False,loc=2,handletextpad=0,borderaxespad=0.2,labelspacing=0.1,borderpad=0.0)
    #axs[ix].legend(prop={'size': 14},frameon=False,bbox_to_anchor=(0, 0.8),handletextpad=0,borderaxespad=0.2,labelspacing=0.1,borderpad=0.0)
    axs[ix].legend(prop={'size': 12},frameon=False,loc=8,handletextpad=0.1)

    #axs[ix].xaxis.set_tick_params(width=1.5)
    #axs[ix].yaxis.set_tick_params(width=1.5)

    axs[ix].tick_params(which='both',direction='in',top=True,right=True)
    axs[ix].xaxis.set_minor_locator(MultipleLocator(0.1))
    axs[ix].yaxis.set_major_locator(MultipleLocator(axis_spacing[ix]))
    axs[ix].yaxis.set_minor_locator(MultipleLocator(0.1))

#axs[1].yaxis.set_minor_locator(AutoMinorLocator(5))

plt.tight_layout()
fileout="gap_tc_panel_" + meth + ".pdf"
print(fileout)
fig.savefig(fileout)
