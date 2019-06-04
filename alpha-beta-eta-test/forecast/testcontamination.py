import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('/home/dfa/sobreira/alsina/alpha-beta-gamma/code/SVA1StyleSheet.mplstyle')
import matplotlib.patches as mpatches
import matplotlib.colors as colors


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--fiducial',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/pipeline/2pt_sim_1110_baseline_Y3cov.fits',
                        help='fit file with fiducial data vectors, covariance matrix and so on.')
    parser.add_argument('--contaminant_marg',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/marg_ab_dxi.fits',
                        help='fit file with contamination data vector, covariance matrix, marginalized best fit')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/ab_dxi.fits',
                        help='fit file with contamination data vector, overall best fit')
    parser.add_argument('--contaminated',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/pipeline/2pt_sim_1110_baseline_Y3cov_contaminated.fits', 
                        help='fit file with contamination data vector, covariance matrix')
    parser.add_argument('--plotpath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='path where output will be send') 
    
    args = parser.parse_args()
    return args

def pretty_rho(meanr, rho, sig,  legend=None, lfontsize=24, color='black', marker='o', ylabel=r'$\rho(\theta)$',title=None,  xlim=None,  ylim=None):
    plt.plot(meanr, rho, color=color, label=legend)
    plt.plot(meanr, -rho, color=color, ls=':')
    #rigth quadrants
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker,  capsize=2)
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker,  capsize=2)
    #leftquadrants
    plt.errorbar( -meanr, rho, yerr=sig, color=color,  marker='^',  capsize=2)
    plt.errorbar( -meanr,-rho, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
    plt.legend(loc='best', fontsize=lfontsize)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    plt.tight_layout()
def plotcorrmat(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    plt.colorbar()
def plotfiducial(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    #PLOTING FIDUCIAL INFO
    fiducialfit = fitfile
    covmatrixfit=fitsio.read(fiducialfit,ext=1)
    xipfit=fitsio.read(fiducialfit,ext=2)
    ximfit=fitsio.read(fiducialfit,ext=3)
    gammatfit=fitsio.read(fiducialfit,ext=4)
    wthetafit=fitsio.read(fiducialfit,ext=5)
    nz_sourcefit=fitsio.read(fiducialfit,ext=6)
    nz_lensfit=fitsio.read(fiducialfit,ext=7)
    lengths = [len(xipfit),len(ximfit),len(gammatfit),len(wthetafit)]

    ##Covariance Matrix.
    plotcorrmat(covmatrixfit)
    plt.title(r'$\xi_{+}(\theta) \mid \xi_{-}(\theta) \mid \gamma_{t}(\theta) \mid \omega(\theta)$')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = out + 'CovariancematrixFiducial.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    ylim = [5.e-9, 5.e-5]
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin]
        xip=xipfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
    filename = out + 'xip_fiducial.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')
    
    ##xim
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{-}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (ximfit['BIN1']==i)&(ximfit['BIN2']==j)
        theta=ximfit['ANG'][bin]
        xim=ximfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr, xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
    filename = out + 'xim_fiducial.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')

def plotcontaminant_marg(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    contaminantfit = fitfile
    covmatrixfit=fitsio.read(contaminantfit,ext=1)
    xipfit=fitsio.read(contaminantfit,ext=2)
    
    plt.clf()
    plotcorrmat(covmatrixfit)
    plt.title(r'$\xi_{+}(\theta)$')
    filename = out + 'Covariancematrix_contaminant.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\delta \xi_{+}(\theta)$'
    nbins=4
    xlim = [2., 300.]
    plt.clf()
    meanr = xipfit['ANG']
    dxi = xipfit['VALUE']
    pretty_rho(meanr, dxi, np.sqrt(np.diag(covmatrixfit)) , legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
    filename = out + 'xi_contaminant.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')

def plotcontaminant(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    contaminantfit = fitfile
    covmatrixfit=fitsio.read(contaminantfit,ext=1)
    xipfit=fitsio.read(contaminantfit,ext=2)
    
    plt.clf()
    plotcorrmat(covmatrixfit)
    plt.title(r'$\xi_{+}(\theta)$')
    filename = out + 'Covariancematrix_contaminant.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\delta \xi_{+}(\theta)$'
    nbins=4
    xlim = [2., 300.]
    plt.clf()
    meanr = xipfit['ANG']
    dxi = xipfit['VALUE']
    pretty_rho(meanr, dxi, np.sqrt(np.diag(covmatrixfit)) , legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
    filename = out + 'xi_contaminant.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')
    
def plotcontaminantandfiducial_old(contaminant, fiducial, out):
    import fitsio
    import itertools
    import numpy as np
    contaminantfit = contaminant
    covmatrixfit_cont=fitsio.read(contaminantfit,ext=1)
    xifit_cont=fitsio.read(contaminantfit,ext=2)
    
    fiducialfit = fiducial
    covmatrixfit_fid=fitsio.read(fiducialfit,ext=1)
    xipfit_fid=fitsio.read(fiducialfit,ext=2)
    ximfit_fid=fitsio.read(fiducialfit,ext=3)
    lengths_fid = [len(xipfit_fid), len(ximfit_fid)]

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    ylim = [5.e-9, 5.e-5]
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit_fid['BIN1']==i)&(xipfit_fid['BIN2']==j)
        theta=xipfit_fid['ANG'][bin]
        xip=xipfit_fid['VALUE'][bin]
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
        plot_tomograpically_bin(ax, i, j, xifit_cont['ANG'],
                                xifit_cont['VALUE'],
                                yerr=np.sqrt(np.diag(covmatrixfit_cont)),
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
        
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{+}^{sim}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    fig.tight_layout()
    filename = out + 'xicont_and_xipfid.png'
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')

    
    ##xim
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{-}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (ximfit_fid['BIN1']==i)&(ximfit_fid['BIN2']==j)
        theta=ximfit_fid['ANG'][bin]
        xim=ximfit_fid['VALUE'][bin]
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
        plot_tomograpically_bin(ax, i, j, xifit_cont['ANG'],
                                xifit_cont['VALUE'],
                                yerr=np.sqrt(np.diag(covmatrixfit_cont)),
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
    
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{-}^{sim}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    fig.tight_layout()
    filename = out + 'xicont_and_ximfid.png'
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')

def plotcontaminantandfiducial(contaminant, fiducial, out):
    import fitsio
    import itertools
    import numpy as np
    contaminantfit = contaminant
    xipfit_cont=fitsio.read(contaminantfit,ext=1)
    ximfit_cont=fitsio.read(contaminantfit,ext=2)
    
    fiducialfit = fiducial
    covmatrixfit_fid=fitsio.read(fiducialfit,ext=1)
    xipfit_fid=fitsio.read(fiducialfit,ext=2)
    ximfit_fid=fitsio.read(fiducialfit,ext=3)
    lengths_fid = [len(xipfit_fid), len(ximfit_fid)]

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    ylim = [5.e-9, 5.e-5]
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit_fid['BIN1']==i)&(xipfit_fid['BIN2']==j)
        theta=xipfit_fid['ANG'][bin]
        xip=xipfit_fid['VALUE'][bin]
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
        bin2 = (xipfit_cont['BIN1']==i)&(xipfit_cont['BIN2']==j)
        theta_cont=xipfit_cont['ANG'][bin2]
        xip_cont=xipfit_cont['VALUE'][bin2]
        yerr_cont=None
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xip_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
        
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{+}^{sim}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    fig.tight_layout()
    filename = out + 'xicont_and_xipfid.png'
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')

    
    ##xim
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{-}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (ximfit_fid['BIN1']==i)&(ximfit_fid['BIN2']==j)
        theta=ximfit_fid['ANG'][bin]
        xim=ximfit_fid['VALUE'][bin]
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
        bin2 = (ximfit_cont['BIN1']==i)&(ximfit_cont['BIN2']==j)
        theta_cont=ximfit_cont['ANG'][bin2]
        xim_cont=ximfit_cont['VALUE'][bin2]
        yerr_cont=None
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xim_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
    
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{-}^{sim}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    fig.tight_layout()
    filename = out + 'xicont_and_ximfid.png'
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')
 
def plotcontaminantandfiducial_marg(contaminant, fiducial, out):
    import fitsio
    import itertools
    import numpy as np
    contaminantfit = contaminant
    covmatrixfit_cont=fitsio.read(contaminantfit,ext=1)
    xipfit_cont=fitsio.read(contaminantfit,ext=2)
    ximfit_cont=fitsio.read(contaminantfit,ext=3)
    lengths_cont = [len(xipfit_cont), len(ximfit_cont)]
    
    fiducialfit = fiducial
    covmatrixfit_fid=fitsio.read(fiducialfit,ext=1)
    xipfit_fid=fitsio.read(fiducialfit,ext=2)
    ximfit_fid=fitsio.read(fiducialfit,ext=3)
    lengths_fid = [len(xipfit_fid), len(ximfit_fid)]

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    ylim = [5.e-9, 5.e-5]
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit_fid['BIN1']==i)&(xipfit_fid['BIN2']==j)
        theta=xipfit_fid['ANG'][bin]
        xip=xipfit_fid['VALUE'][bin]
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
        bin = (xipfit_cont['BIN1']==i)&(xipfit_cont['BIN2']==j)
        theta_cont=xipfit_cont['ANG'][bin]
        xip_cont=xipfit_cont['VALUE'][bin]
        yerr_cont=get_error_cont(covmatrixfit_cont, lengths_cont, 'delta_xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xip_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
        
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{+}^{sim}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    fig.tight_layout()
    filename = out + 'xicont_and_xipfid_marg.png'
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')

    
    ##xim
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{-}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (ximfit_fid['BIN1']==i)&(ximfit_fid['BIN2']==j)
        theta=ximfit_fid['ANG'][bin]
        xim=ximfit_fid['VALUE'][bin]
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
        bin = (ximfit_cont['BIN1']==i)&(ximfit_cont['BIN2']==j)
        theta_cont=ximfit_cont['ANG'][bin]
        xim_cont=ximfit_cont['VALUE'][bin]
        yerr_cont=get_error_cont(covmatrixfit_cont, lengths_cont, 'delta_xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xim_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
    
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{-}^{sim}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    fig.tight_layout()
    filename = out + 'xicont_and_ximfid_marg.png'
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')
    
def plotcontaminated(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    #PLOTING FIDUCIAL INFO
    fiducialfit = fitfile
    covmatrixfit=fitsio.read(fiducialfit,ext=1)
    xipfit=fitsio.read(fiducialfit,ext=2)
    ximfit=fitsio.read(fiducialfit,ext=3)
    gammatfit=fitsio.read(fiducialfit,ext=4)
    wthetafit=fitsio.read(fiducialfit,ext=5)
    nz_sourcefit=fitsio.read(fiducialfit,ext=6)
    nz_lensfit=fitsio.read(fiducialfit,ext=7)
    lengths = [len(xipfit),len(ximfit),len(gammatfit),len(wthetafit)]

    
    ##Covariance Matrix.
    plt.clf()
    plotcorrmat(covmatrixfit)
    plt.title(r'$\xi_{+}(\theta) \mid \xi_{-}(\theta) \mid \gamma_{t}(\theta) \mid \omega(\theta)$')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = out + 'Covariancematrixcontaminated.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')
    
    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    ylim = [5.e-9, 5.e-5]
    plt.clf()
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin]
        xi=xipfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xi, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
    filename = out + 'xi_contaminated.png'
    fig.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')
def checkcontamination(contaminantedfile, fiducial,  out):
    import fitsio
    import itertools
    import numpy as np

    fiducialfit = fiducial
    covmatrixfit=fitsio.read(fiducialfit,ext=1)
    xipfit=fitsio.read(fiducialfit,ext=2)
    ximfit=fitsio.read(fiducialfit,ext=3)
    gammatfit=fitsio.read(fiducialfit,ext=4)
    wthetafit=fitsio.read(fiducialfit,ext=5)
    nz_sourcefit=fitsio.read(fiducialfit,ext=6)
    nz_lensfit=fitsio.read(fiducialfit,ext=7)
    lengths = [len(xipfit),len(ximfit),len(gammatfit),len(wthetafit)]

    contaminatedfit = contaminantedfile
    covmatrixfit_cont=fitsio.read(contaminatedfit,ext=1)
    xifit_cont=fitsio.read(contaminatedfit,ext=2)

    #RESIDUAL COVARIANCE MATRIX
    diff_covmatrix = covmatrixfit_cont - covmatrixfit
    if (np.linalg.det(diff_covmatrix) != 0): 
        plotcorrmat(diff_covmatrix)
        plt.title(r'$\xi_{+}(\theta) \mid \xi_{-}(\theta) \mid \gamma_{t}(\theta) \mid \omega(\theta)$')
        pos_lines = [0]
        for i in range(len(lengths)):
            pos_lines.append(pos_lines[i] + lengths[i])
        pos_lines = pos_lines[1:-1]
        for line in pos_lines:
            plt.axvline(x=line, c='k', lw=1, ls='-')
            plt.axhline(y=line, c='k', lw=1, ls='-')
        plt.tight_layout()
        filename = out + 'contaminated-fiducial_covmat.png'
        plt.savefig(filename, dpi=500)
        print(filename, 'Printed!')


    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    plt.clf()
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    ylim = [5.e-9, 5.e-5]
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin1 = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin1]
        xi=xipfit['VALUE'][bin1]
        bin2 = (xifit_cont['BIN1']==i)&(xifit_cont['BIN2']==j)
        theta_cont=xifit_cont['ANG'][bin2]
        xi_cont=xifit_cont['VALUE'][bin2]
        yerr=get_error(diff_covmatrix, lengths, 'xip')[bin2]
        plot_tomograpically_bin(ax, i, j, theta, xi_cont - xi, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue',  ylim=ylim)
    filename = out + 'xicont-xipfid.png'
    fig.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')
    
def plot_tomograpically_bin(ax, i, j, x, y, xerr=None,
                            yerr=None,xlabel='', ylabel='',
                            nbins=None, nbins1=None, nbins2=None,
                            color='blue',label=None, ylim=None):
    import numpy as np
    if(nbins):
        nbins1=nbins;nbins2=nbins
    elif(nbins1 and nbins2):
        print("Entering in not symmetric mode")
    else:
        print('Error, defining symultaneosly nbins and nbins1 or nbins2')
        raise
        
    if xerr is None: xerr =  np.zeros(len(x))
    if yerr is None: yerr =  np.zeros(len(y))

    if y.size !=0 :
        ax[j-1][i-1].plot(x, y, color=color, label=label, linewidth=0.4)
        ax[j-1][i-1].plot(x,-y, color=color,  ls=':', linewidth=0.4)
        ax[j-1][i-1].errorbar(x[y>0], y[y>0],xerr=xerr[y>0],yerr=yerr[y>0], fmt='o' ,capsize=1,markersize=2, color=color, mec=color, elinewidth=0.3, ls='')
        ax[j-1][i-1].errorbar(x[y<0], -y[y<0],xerr=xerr[y<0],yerr=yerr[y<0], fmt='o' ,capsize=1,markersize=2, color=color, mec=color, elinewidth=0.3, ls='')
        ax[j-1][i-1].errorbar( -x, y,xerr=xerr,yerr=yerr, fmt='^' ,capsize=1,markersize=2, color=color, mec=color, elinewidth=0.3, ls='')
        ax[j-1][i-1].errorbar( -x,-y,xerr=xerr,yerr=yerr, fmt='^' ,capsize=1,markersize=2, color=color, mec=color, elinewidth=0.3, ls=':')
        ax[j-1][i-1].text(0.85, 0.85, "{},{}".format(i, j), horizontalalignment='center', verticalalignment='center', transform=ax[j-1][i-1].transAxes, fontsize=12)
        ax[j-1][i-1].set_xscale('log', nonposx='clip')
        ax[j-1][i-1].set_yscale('log', nonposy='clip')
        if(label):
            ax[j-1][i-1].legend(loc='best')
        if (j == nbins2):
            ax[j-1][i-1].set_xlabel(xlabel)
        if (i == 1):
            ax[j-1][i-1].set_ylabel(ylabel)
        if ylim is not None:
            ax[j-1][i-1].set_ylim(ylim)
    else:
        ax[j-1][i-1].set_visible(False)
def get_error(covmatrix, lengths, name):
    import numpy as np
    if name is not None:
        if (name=='xip'):
            start = 0
            end =start + lengths[0]
        elif (name=='xim'):
            start = lengths[0]
            end =start + lengths[1]
        elif (name=='gammat'):
            start = lengths[0] + lengths[1]
            end =start + lengths[2]
        elif (name=='wtheta'):
            start = lengths[0] + lengths[1]+ lengths[2]
            end =start + lengths[3]
        return np.diagonal(covmatrix)[start:end]**0.5
    else:
        print("Correlation function no defined")
        return None

def get_error_cont(covmatrix, lengths, name):
    import numpy as np
    if name is not None:
        if (name=='delta_xip'):
            start = 0
            end =start + lengths[0]
        elif (name=='delta_xim'):
            start = lengths[0]
            end =start + lengths[1]
        elif (name=='delta_gammat'):
            start = lengths[0] + lengths[1]
            end =start + lengths[2]
        elif (name=='delta_wtheta'):
            start = lengths[0] + lengths[1]+ lengths[2]
            end =start + lengths[3]
        return np.diagonal(covmatrix)[start:end]**0.5
    else:
        print("Correlation function no defined")
        return None

##Ploting covariance matrix might be not than ilustrative as
##correlations matrix


def main():
    
    import numpy as np
    args = parse_args()
    out = os.path.expanduser(args.plotpath)
    try:
        if not os.path.isdir(out):
            os.makedirs(out)
    except OSError as e:
        print "Ignore OSError from makedirs(work):"
        print e
        pass
    
    #plotfiducial(args.fiducial,  out)
    #plotcontaminant(args.contaminant, out)
    plotcontaminantandfiducial(args.contaminant, args.fiducial, out)
    plotcontaminantandfiducial_marg(args.contaminant_marg, args.fiducial, out)
    #plotcontaminated(args.contaminated, out)
    #checkcontamination(args.contaminated,args.fiducial,  out)
    


if __name__ == "__main__":
    main()
