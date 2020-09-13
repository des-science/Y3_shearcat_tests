import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import numpy as np
import itertools
import matplotlib.patches as mpatches
import matplotlib.colors as colors


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    parser.add_argument('--taus1',
                        default='/data/catalogs/rhos-taus/marco/Y3_03_31_20/tau__JK_Y3_1.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--taus2',
                        #default='/data/catalogs/rhos-taus/marco/Y3_03_31_20/tau__sw__JK_Y3_1.fits',
                        default='/data/catalogs/rhos-taus/TAUS_Y3_03-31-20-mod_SPONLY_non-tomographic.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--plotspath',
                        default='/data/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args

def pretty_tau(ax,  meanr, rho, sig,  legend=None, legsz=24, color='black', marker='o', marsz=12,  xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho(\theta)$', labsz=24 , title=None,  xlim=None,  ylim=None):    
    if sig is None: sig =  np.zeros(len(rho))
    ax.plot(meanr, rho, color=color, label=legend)
    ax.plot(meanr, -rho, color=color, ls=':')
    #rigth quadrants
    ax.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker,  capsize=2,  markersize=marsz)
    ax.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker,  capsize=2,  markersize=marsz)
    #leftquadrants
    ax.errorbar( -meanr, rho, yerr=sig, color=color,  marker='^',  capsize=2,  markersize=marsz)
    ax.errorbar( -meanr,-rho, yerr=sig, color=color,  marker='^', ls=':', capsize=2,  markersize=marsz)
    if legend is not None: ax.legend(loc='best', fontsize=legsz)
    if ylim is not None: ax.set_ylim( ylim )
    if xlim is not None: ax.set_xlim(xlim)
    ax.tick_params(axis='both', which='major', labelsize=10)
    if xlabel is not None: ax.set_xlabel(r'$\theta$ (arcmin)', fontsize=labsz)
    if ylabel is not None: ax.set_ylabel(ylabel, fontsize=labsz)
    ax.set_xscale('log')
    ax.set_yscale('log', nonposy='clip')
    if title is not None: ax.set_title(title)


def main():
    import sys; sys.path.append(".")
    from src.readfits import  read_taus_plots

 
    #legends=['DES Y1', 'DES Y3'];  colors=['gray', 'red'];  markers=['o', '^']
    #legends=[r'e$_{\textrm{model}}$', r'e$_{\textrm{obs}}$'];  colors=['blue', 'red'];  markers=['o', '^']
    #legends=['Mean-subtracted', r'Non-mean-subtracted'];  colors=['blue', 'red'];  markers=['o', '^']
    #legends=['DES Y3', r'DES Y3 with weights'];  colors=['blue', 'red'];  markers=['o', '^']
    #legends=['DES Y3 (w=epiff*dT, e=eobs)', r'DES Y3 (w=eobs*dT, e=epiff)'];  colors=['blue', 'red'];  markers=['o', '^']
    legends=['DES Y3 (e= eobs, w=epiff*dT, e=eobs)', r'DES Y3 (e=epiff,w=eobs*dT, e=epiff)'];  colors=['blue', 'red'];  markers=['o', '^']


    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    #One fig for each correlation
    xlim = [0.1, 250.]
    figs = [];  axs = []; filenames = []
    ylabels  = np.concatenate([ np.array([r'$\tau_{%i+}(\theta)$'%(i), r'$\tau_{%i-}(\theta)$'%(i)]) for i in [0, 2, 5] ])
    names =  np.concatenate([ np.array(['Y1-Y3_tau%ip.png'%(i),'Y1-Y3_tau%im.png'%(i) ]) for i in [0, 2, 5] ] ) 

    for i in range(len(ylabels)):
        figaux, axaux = plt.subplots()
        figs.append(figaux); axs.append(axaux)
        filenames.append(os.path.join(plotspath,names[i]))

    
    for j, fil in enumerate([args.taus1, args.taus2]):
        meanr, taus,  cov_taus =  read_taus_plots(fil)
        for i in range(len(taus)):
            sig_tau = np.sqrt(np.diag(cov_taus[i]))
            pretty_tau(axs[i], meanr, taus[i], sig_tau, legend=legends[j], color=colors[j], marker=markers[j], ylabel=ylabels[i], xlim=xlim,  ylim=None)
    
    for i, fig in enumerate(figs):
        fig.tight_layout()
        fig.savefig(filenames[i], dpi=200)
        plt.close(fig)
        print(filenames[i], 'Printed!')

    ##SINGLE FIG
    
    ylims = [ [2.e-7, 5.e-5] , [1.e-10, 5.e-5] ,[1.e-10, 5.e-6] ,[1.e-11, 5.e-6] ,[1.e-12, 5.e-7] ,[1.e-12, 2.e-6]]
    #ylims =  [None]*6
    nbinsx = 1; nbinsy = 3
    fig_tp, ax_tp = plt.subplots(nbinsx, nbinsy, figsize=(10*nbinsx, 1.0*nbinsy), sharey=False, sharex=False)
    fig_tm, ax_tm = plt.subplots(nbinsx, nbinsy, figsize=(10*nbinsx, 1.0*nbinsy), sharey=False, sharex=False)
    a=[i for i in range(nbinsx)]
    b=[j for j in range(nbinsy)]
    bin_pairs=[]
    for p in itertools.product(a, b): bin_pairs.append(p)

    
    for q, fil in enumerate([args.taus1, args.taus2]):
        meanr, taus,  cov_taus =  read_taus_plots(fil)
        for k, [i,j] in enumerate(bin_pairs):
            sig_tau = np.sqrt(np.diag(cov_taus[2*k]))
            pretty_tau(ax_tp[j], meanr, taus[2*k], sig_tau, legend=None, color=colors[q], marker=markers[q], ylabel=ylabels[2*k], xlim=xlim,  ylim=ylims[2*k],  legsz=10,  labsz=10,  marsz=2)
            sig_taum = np.sqrt(np.diag(cov_taus[2*k + 1]))
            pretty_tau(ax_tm[j], meanr, taus[2*k + 1], sig_taum, legend=None, color=colors[q], marker=markers[q], ylabel=ylabels[2*k + 1], xlim=xlim,  ylim=ylims[2*k + 1],  legsz=10,  labsz=10,  marsz=2)

    line0, = plt.plot([1, 2, 3], label=legends[0], color= colors[0])
    line1, = plt.plot([1, 2, 3], label=legends[1], color= colors[1])
    fig_tp.legend(handles=[line0, line1], bbox_to_anchor=(0.5,-0.02), loc='lower center',  ncol=2,  fontsize= 10)
    fig_tm.legend(handles=[line0, line1], bbox_to_anchor=(0.5,-0.02), loc='lower center',  ncol=2,  fontsize= 10)
    
    
    filenames = [os.path.join(plotspath, s ) for s in ['Y3-Y3w_alltausp.png', 'Y3-Y3w_alltausm.png']]
    
    for i, fig in enumerate([fig_tp, fig_tm]):
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        fig.savefig(filenames[i], dpi=200)
        plt.close(fig)
        print(filenames[i], 'Printed!')
    

 

if __name__ == "__main__":
    main()
