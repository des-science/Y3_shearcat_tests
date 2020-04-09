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
    parser.add_argument('--rhos1',
                        default='/data/catalogs/rhos-taus/32bins/Y3binning/nomasks/RHOS_Y1-obs.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--rhos2',
                        default='/data/catalogs/rhos-taus/32bins/Y3binning/RHOS_Y3_7-24-19-obs.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args

def pretty_rho(ax,  meanr, rho, sig,  legend=None, legsz=24, color='black', marker='o', marsz=12,  xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho(\theta)$', labsz=24 , title=None,  xlim=None,  ylim=None):    
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
    from src.readfits import  read_rhos_plots

 
    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    #One fig for each correlation
    xlim = [0.5, 250.]
    figs = [];  axs = []; filenames = []
    ylabels  = np.concatenate([ np.array([r'$\rho_{%i+}(\theta)$'%(i), r'$\rho_{%i-}(\theta)$'%(i)]) for i in range(0, 6) ])
    names =  np.concatenate([ np.array(['Y1-Y3_rho%ip.png'%(i),'Y1-Y3_rho%im.png'%(i) ]) for i in range(0, 6) ] ) 

    for i in range(len(ylabels)):
        figaux, axaux = plt.subplots()
        figs.append(figaux); axs.append(axaux)
        filenames.append(os.path.join(plotspath,names[i]))

    legends=['DES Y1', 'DES Y3'];  colors=['gray', 'red'];  markers=['o', '^']
    for j, fil in enumerate([args.rhos1, args.rhos2]):
        meanr, rhos,  cov_rhos =  read_rhos_plots(fil)
        for i in range(len(rhos)):
            sig_rho = np.sqrt(np.diag(cov_rhos[i]))
            pretty_rho(axs[i], meanr, rhos[i], sig_rho, legend=legends[j], color=colors[j], marker=markers[j], ylabel=ylabels[i], xlim=xlim,  ylim=None)
    
    for i, fig in enumerate(figs):
        fig.tight_layout()
        fig.savefig(filenames[i], dpi=200)
        plt.close(fig)
        print(filenames[i], 'Printed!')

    ##SINGLE FIG
    ylims = [[8.e-6, 5.e-4] , None ,[1.e-9, 5.e-6] ,None ,[1.e-7, 2.e-6] ,None , [1.e-11, 1.e-8] ,None  ,[1.e-12, 5.e-8] ,None ,[1.e-8, 8.e-7] ,[1.e-11 ,4.e-8 ] ] 
    nbinsx = 2; nbinsy = 3
    fig_tp, ax_tp = plt.subplots(nbinsx, nbinsy, figsize=(4.8*nbinsx, 1.6*nbinsy), sharey=False, sharex=False)
    fig_tm, ax_tm = plt.subplots(nbinsx, nbinsy, figsize=(4.8*nbinsx, 1.6*nbinsy), sharey=False, sharex=False)
    a=[i for i in range(nbinsx)]
    b=[j for j in range(nbinsy)]
    bin_pairs=[]
    for p in itertools.product(a, b): bin_pairs.append(p)
    for q, fil in enumerate([args.rhos1, args.rhos2]):
        meanr, rhos,  cov_rhos =  read_rhos_plots(fil)
        for k, [i,j] in enumerate(bin_pairs):
            sig_rho = np.sqrt(np.diag(cov_rhos[2*k]))
            pretty_rho(ax_tp[i][j], meanr, rhos[2*k], sig_rho, legend=None, color=colors[q], marker=markers[q], ylabel=ylabels[2*k], xlim=xlim,  ylim=ylims[2*k],  legsz=10,  labsz=10,  marsz=2)
            sig_rhom = np.sqrt(np.diag(cov_rhos[2*k + 1]))
            pretty_rho(ax_tm[i][j], meanr, rhos[2*k + 1], sig_rhom, legend=None, color=colors[q], marker=markers[q], ylabel=ylabels[2*k + 1], xlim=xlim,  ylim=ylims[2*k + 1],  legsz=10,  labsz=10,  marsz=2)

    line0, = plt.plot([1, 2, 3], label='DES Y1', color= colors[0])
    line1, = plt.plot([1, 2, 3], label='DES Y3', color= colors[1])
    fig_tp.legend(handles=[line0, line1], bbox_to_anchor=(0.4,-0.03), loc='lower left',  ncol=2,  fontsize= 10)
    fig_tm.legend(handles=[line0, line1], bbox_to_anchor=(0.4,-0.03), loc='lower left',  ncol=2,  fontsize= 10)
    
    
    filenames = [os.path.join(plotspath, s ) for s in ['Y1-Y3_allrhosp.png', 'Y1-Y3_allrhosm.png']]
    
    for i, fig in enumerate([fig_tp, fig_tm]):
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        fig.savefig(filenames[i], dpi=200)
        plt.close(fig)
        print(filenames[i], 'Printed!')

 

if __name__ == "__main__":
    main()
