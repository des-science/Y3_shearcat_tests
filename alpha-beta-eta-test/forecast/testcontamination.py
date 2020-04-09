import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import matplotlib.patches as mpatches
import matplotlib.colors as colors


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--fiducial',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/pipeline/2pt_sim_1110_baseline_Y3cov.fits',
                        help='fit file with fiducial data vectors, covariance matrix and so on.')
    parser.add_argument('--contaminant_marg',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/topropagatecosmology/marg_abe_dxi_eq_4_1sigma.fits',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/marg_abe_eq4_contaminant.fits', 
                        help='fit file with contamination data vector, covariance matrix, marginalized best fit')
    parser.add_argument('--contaminant_over',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/overall_ab_dxi_eq_4_.fits',
                        help='fit file with contamination data vector, overall best fit')
    parser.add_argument('--contaminated',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/2pt_sim_1110_baseline_Y3cov_contaminated_sup_2sig.fits', 
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

def plot_covmat(fitfile, out,n2pts=4, filename='CovarianceMatrix'):
    import fitsio
    covmatrixfit=fitsio.read(fitfile,ext=1)
    lengths = [len(fitsio.read(fitfile,ext=i)) for i in range(2, 2 + n2pts)]
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
    filename = os.path.join(out,filename)
    plt.savefig(filename, dpi=150)
    print(filename, 'Printed!')

def plot_tomotwopoint(fitfile, out, n2pts=4, overall=False,  xlabels=[r'$\theta$ [arcmin]',r'$\theta$ [arcmin]'],  ylabels=[r'$\chi_{+}$ [arcmin]',r'$\theta$ [arcmin]'],  filenames=[r'$\theta$ [arcmin]',r'$\theta$ [arcmin]']):
    import fitsio
    import itertools
    import numpy as np
    #PLOTING FIDUCIAL INFO
    if overall:
        xipfit=fitsio.read(fitfile,ext=1)
        ximfit=fitsio.read(fitfile,ext=2)
        lengths = [len(fitsio.read(fitfile,ext=i)) for i in range(1, 1+n2pts)]
    else:
        covmatrixfit=fitsio.read(fitfile,ext=1)
        xipfit=fitsio.read(fitfile,ext=2)
        ximfit=fitsio.read(fitfile,ext=3)
        lengths = [len(fitsio.read(fitfile,ext=i)) for i in range(2, 2+n2pts)]
        
    #xip
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
        if overall: yerr = None
        else: yerr=get_error(covmatrixfit, lengths, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabels[0], ylabel=ylabels[0], nbins=4,
                                color='blue', ylim=ylim)
    filename = os.path.join(out,filenames[0])
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')
    
    ##xim
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
        if overall: yerr = None
        else: yerr=get_error(covmatrixfit, lengths, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr, xlabel=xlabels[1], ylabel=ylabels[1], nbins=4,
                                color='blue')
    filename = os.path.join(out,filenames[1])
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')


def plotcontaminantandfiducial(contaminant, fiducial, out,overall=False, title=None, filenames=None,  nsig=1):
    import fitsio
    import itertools
    import numpy as np

    if overall:
        xipfit_cont=fitsio.read(contaminant,ext=1)
        ximfit_cont=fitsio.read(contaminant,ext=2)
    else:
        covmatrixfit_cont=fitsio.read(contaminant,ext=1)
        covmatrixfit_cont *= (nsig**2)
        xipfit_cont=fitsio.read(contaminant,ext=2)
        ximfit_cont=fitsio.read(contaminant,ext=3)
        
    lengths_cont = [len(xipfit_cont), len(ximfit_cont)]
    
    covmatrixfit_fid=fitsio.read(fiducial,ext=1)
    xipfit_fid=fitsio.read(fiducial,ext=2)
    ximfit_fid=fitsio.read(fiducial,ext=3)
    lengths_fid = [len(xipfit_fid), len(ximfit_fid)]

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    ylim = [1.e-10, 1.e-4]
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
        if overall: yerr = None
        else: yerr=get_error(covmatrixfit_fid, lengths_fid, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
        bin = (xipfit_cont['BIN1']==i)&(xipfit_cont['BIN2']==j)
        theta_cont=xipfit_cont['ANG'][bin]
        xip_cont=xipfit_cont['VALUE'][bin]
        if overall: yerr_cont = None
        else: yerr_cont=get_error(covmatrixfit_cont, lengths_cont, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xip_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
        
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{+}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    if title is not None: fig.suptitle(title)
    fig.tight_layout()
    filename = os.path.join(out,filenames[0])
    plt.savefig(filename, dpi=200)
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
        if overall: yerr = None
        else: yerr=get_error(covmatrixfit_fid, lengths_fid, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
        bin = (ximfit_cont['BIN1']==i)&(ximfit_cont['BIN2']==j)
        theta_cont=ximfit_cont['ANG'][bin]
        xim_cont=ximfit_cont['VALUE'][bin]
        if overall: yerr_cont = None
        else: yerr_cont=get_error(covmatrixfit_cont, lengths_cont, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xim_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
    
    red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{-}$')
    fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    if title is not None: fig.suptitle(title)
    fig.tight_layout()
    filename = os.path.join(out,filenames[1])
    plt.savefig(filename, dpi=200)
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
    xipfit_cont=fitsio.read(contaminatedfit,ext=2)
    ximfit_cont=fitsio.read(contaminatedfit,ext=3)

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
        filename = os.path.join(out,'contaminated-fiducial_covmat.png')
        plt.savefig(filename, dpi=500)
        print(filename, 'Printed!')
    else:
        print("COVMAT of contaminated and fiducial is the same")


    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    plt.clf()
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    ylim = [1.e-10, 1.e-4]
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin1 = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin1]
        xi=xipfit['VALUE'][bin1]
        bin2 = (xipfit_cont['BIN1']==i)&(xipfit_cont['BIN2']==j)
        theta_cont=xipfit_cont['ANG'][bin2]
        xi_cont=xipfit_cont['VALUE'][bin2]
        yerr=get_error(diff_covmatrix, lengths, 'xip')[bin2]
        plot_tomograpically_bin(ax, i, j, theta, xi_cont - xi, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue',  ylim=ylim)
    filename = os.path.join(out, 'xipcont-xipfid.png')
    fig.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(filename, 'Printed!')

    #xim
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{-}(\theta)$'
    nbins=4
    plt.clf()
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    ylim = [1.e-10, 1.e-4]
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin1 = (ximfit['BIN1']==i)&(ximfit['BIN2']==j)
        theta=ximfit['ANG'][bin1]
        xi=ximfit['VALUE'][bin1]
        bin2 = (ximfit_cont['BIN1']==i)&(ximfit_cont['BIN2']==j)
        theta_cont=ximfit_cont['ANG'][bin2]
        xi_cont=ximfit_cont['VALUE'][bin2]
        yerr=get_error(diff_covmatrix, lengths, 'xim')[bin2]
        plot_tomograpically_bin(ax, i, j, theta, xi_cont - xi, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue',  ylim=ylim)
    filename = os.path.join(out, 'ximcont-ximfid.png')
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
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass

    #COVARIANCES MATRIX
    #plot_covmat(args.fiducial, out, filename='CovariancematrixFiducial.png')
    #plot_covmat(args.contaminated, out, filename='CovariancematrixContaminated.png')

    ##TWO POINT STATS
    
    #plot_tomotwopoint(args.fiducial,out,ylabels=[r'$\xi_{+}(\theta)$',r'$\xi_{-}(\theta)$'],
    #                  filenames=['xip_fiducial.png','xim_fiducial.png'])
    #plot_tomotwopoint(args.contaminated,out,ylabels=[r'$\xi_{+}(\theta)$',r'$\xi_{-}(\theta)$'],
    #                  filenames=['xip_contaminated.png','xim_contaminated.png'])
    #plot_tomotwopoint(args.contaminant_marg,out,n2pts=2,
    #                  ylabels=[r'$\delta \xi_{+}(\theta)$',r'$\delta \xi_{-}(\theta)$'],
    #                  filenames=['xip_contaminant.png','xim_contaminant.png'])
    #plot_tomotwopoint(args.contaminant_over,out,n2pts=2,overall=True,
    #                  ylabels=[r'$\delta \xi_{+}(\theta)$',r'$\delta\xi_{-}(\theta)$'],
    #                  filenames=['xip_contaminant_over.png','xim_contaminant_over.png'])
    
    
    
    

    plotcontaminantandfiducial(args.contaminant_marg, args.fiducial, out, title='Alpha-Beta-eta', filenames=['xipcont_xipfid_abe_1sig.png','ximcont_ximfid_abe_1sig.png'], nsig=2 )
    #plotcontaminantandfiducial(args.contaminant_marg, args.fiducial, out, title='Alpha-Beta', filenames=['xipcont_xipfid_ab_1sig.png','ximcont_ximfid_ab_1sig.png'] )
    #plotcontaminantandfiducial(args.contaminant_over, args.fiducial, out, overall=True, title='Alpha-Beta', filenames=['xipcontover_xipfid_abe.png','ximcontover_ximfid_abe.png'] )

    

    #Checking contamination
    #plotting contaminated minus fiducial
    #checkcontamination(args.contaminated,args.fiducial,  out)
    


if __name__ == "__main__":
    main()
