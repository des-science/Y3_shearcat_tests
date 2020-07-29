import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import numpy as np

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='samples total PSF bias using rhos cosmology and previous samples of abe')
    parser.add_argument('--samplesabe', default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_ab_margin_Y3_03-31-20_JKmarco/sample_ab_Y3_03-31-20_JKmarco.fits')
    parser.add_argument('--rhoscosmo', default='/data/catalogs/rhos-taus/RHOS_Y3_7-24-19-mod_cosmo.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--contaminant', default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_ab_margin_Y3_03-31-20_JKmarco/2pcfpsfbias_ab_Y3_03-31-20_JKmarco.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--fiducial',
                        default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/forecast/pipeline/2pt_sim_1110_baseline_Y3cov.fits',
                        help='fit file with fiducial data vectors, covariance matrix and so on.')
    parser.add_argument('--burn_frac', default=0.0, type=float, 
                        help='Burn frac of samples')
    parser.add_argument('--nsig', default=1,  type=int, 
                        help='Burn frac of samples')
    parser.add_argument('--corner', default=False,
                        action='store_const', const=True, help='plot samples abe contours and posterios using corner instead getdist')
    parser.add_argument('--plotspath',
                        default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_ab_margin_Y3_03-31-20_JKmarco/plots', 
                        help='location of the plots.')


    args = parser.parse_args()

    return args
                        
def corner_plot(samples, labels, filename, title=None):
    import corner
    import matplotlib.ticker as ticker
    #burn = 5000
    plt.clf()
    fig = corner.corner(np.c_[samples].T, labels=labels, 
                        quantiles=[0.16, 0.5, 0.84],  #-1sigma,0sigma,1sigma
                        levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-9./2)), #1sigma, 2sigma and 3sigma contours
                        show_titles=True, title_kwargs={"fontsize": 16}, title_fmt= '.4f', 
                        smooth1d=None, plot_contours=True, 
                        no_fill_contours=False, plot_density=True, use_math_text=True)
    for i in range(len(fig.axes)):
        fig.axes[i].locator_params(axis='x', nbins=2)
        fig.axes[i].locator_params(axis='y', nbins=2)
        fig.axes[i].tick_params(axis='x', rotation =0, labelsize=16)
        fig.axes[i].tick_params(axis='y', rotation =90, labelsize=16)
        fig.axes[i].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
        fig.axes[i].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
  
    if title is not None:
        plt.suptitle(title,  fontsize=24,  color='blue', x=0.8  )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename,  dpi=150)
    plt.close(fig)
    print(filename, "Printed")
def plotcorrmat(cov):
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    clb = plt.colorbar()
    
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

def plot_tomograpically_bin(ax, i, j, x, y, xerr=None,
                            yerr=None,xlabel='', ylabel='',
                            nbins=None, nbins1=None, nbins2=None,
                            color='blue',label=None, xlim=None,  ylim=None):
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
        if xlim is not None:
            ax[j-1][i-1].set_xlim(xlim)
    else:
        ax[j-1][i-1].set_visible(False)
def plotcontaminantandfiducial(contaminant, fiducial, out, title=None, filenames=None,  nsig=1):
    import fitsio
    import itertools
    import numpy as np

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
        
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue', ylim=ylim)
        bin = (xipfit_cont['BIN1']==i)&(xipfit_cont['BIN2']==j)
        theta_cont=xipfit_cont['ANG'][bin]
        xip_cont=xipfit_cont['VALUE'][bin]
        yerr_cont=get_error(covmatrixfit_cont, lengths_cont, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xip_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)
        
    #red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    #blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{+}$')
    #fig.legend(handles=[red_patch, blue_patch], fontsize=20)
    red_patch  = plt.errorbar(0, 0,  yerr=0.5, fmt='o' ,capsize=1,markersize=7,  color='red', label=r'$\delta \xi_{+}^{PSF} $')
    blue_patch = plt.errorbar(0, 0, yerr=0.5, fmt='o' ,capsize=1,markersize=7,  color='blue', label=r'$\xi_{+}$')
    fig.legend(handles=[ red_patch, blue_patch], fontsize=20,  loc='upper right') #bbox_to_anchor=(0.6,0.0),
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
        yerr=get_error(covmatrixfit_fid, lengths_fid, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
        bin = (ximfit_cont['BIN1']==i)&(ximfit_cont['BIN2']==j)
        theta_cont=ximfit_cont['ANG'][bin]
        xim_cont=ximfit_cont['VALUE'][bin]
        yerr_cont=get_error(covmatrixfit_cont, lengths_cont, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta_cont,
                                xim_cont,
                                yerr=yerr_cont,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red', ylim=ylim)

    #red_patch = mpatches.Patch(color='red', label=r'$\delta \xi $')
    #blue_patch = mpatches.Patch(color='blue', label=r'$\xi_{-}$')
    #fig.legend(handles=[red_patch, blue_patch], fontsize=20)

    red_patch  = plt.errorbar(0, 0,  yerr=0.5, fmt='o' ,capsize=1,markersize=7,  color='red', label=r'$\delta \xi_{-}^{PSF} $')
    blue_patch = plt.errorbar(0, 0, yerr=0.5, fmt='o' ,capsize=1,markersize=7,  color='blue', label=r'$\xi_{-}$')
    fig.legend(handles=[ red_patch, blue_patch], fontsize=20,  loc='upper right') #bbox_to_anchor=(0.6,0.0),

    if title is not None: fig.suptitle(title)
    fig.tight_layout()
    filename = os.path.join(out,filenames[1])
    plt.savefig(filename, dpi=200)
    plt.close(fig)
    print(filename, 'Printed!')
    
def plotrhoscontributions(abesamps, rhos_cosmo_file, plotspath, mflags):
    from src.readfits import  read_rhos
    import itertools
    from src.maxlikelihood import percentiles, bestparameters
    
    meanr, rhos,  covmat = read_rhos(rhos_cosmo_file)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m = rhos

    if ('a1' in abesamps.dtype.names): alist = [abesamps['a1'], abesamps['a2'], abesamps['a3'], abesamps['a4']]
    else: alist = [[0]*len(abesamps) ]*4; print('ALPHA is not in samples')
    if ('b1' in abesamps.dtype.names): blist = [abesamps['b1'], abesamps['b2'], abesamps['b3'], abesamps['b4']]
    else: blist = [[0]*len(abesamps) ]*4; print('BETA is not in samples')
    if ('e1' in abesamps.dtype.names): elist = [abesamps['e1'], abesamps['e2'], abesamps['e3'], abesamps['e4']]
    else: elist = [[0]*len(abesamps) ]*4; print('ETA is not in samples')

    nbins=4
    a=[i for i in range(nbins)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2): bin_pairs.append(p)

    veclistrho0 = [];    veclistrho1 = []
    veclistrho2 = [];    veclistrho3 = []
    veclistrho4 = [];    veclistrho5 = []
    for z in range(len(abesamps['a1'])):       
        drho0p = [alist[i][z]*alist[j][z]*rho0p for i,j in bin_pairs]
        drho0m = [alist[i][z]*alist[j][z]*rho0m for i,j in bin_pairs]
        drho0 = drho0p + drho0m
        veclistrho0.append(np.concatenate(np.c_[drho0]))

        drho1p = [blist[i][z]*blist[j][z]*rho1p for i,j in bin_pairs]
        drho1m = [blist[i][z]*blist[j][z]*rho1m for i,j in bin_pairs]
        drho1 = drho1p + drho1m
        veclistrho1.append(np.concatenate(np.c_[drho1]))

        drho2p = [(blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2p for i,j in bin_pairs]
        drho2m = [(blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2m for i,j in bin_pairs]
        drho2 = drho2p + drho2m
        veclistrho2.append(np.concatenate(np.c_[drho2]))

        drho3p = [elist[i][z]*elist[j][z]*rho3p for i,j in bin_pairs]
        drho3m = [elist[i][z]*elist[j][z]*rho3m for i,j in bin_pairs]
        drho3 = drho3p + drho3m
        veclistrho3.append(np.concatenate(np.c_[drho3]))

        drho4p = [(blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4p for i,j in bin_pairs]
        drho4m = [(blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4m for i,j in bin_pairs]
        drho4 = drho4p + drho4m
        veclistrho4.append(np.concatenate(np.c_[drho4]))

        drho5p = [(elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5p for i,j in bin_pairs]
        drho5m = [(elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5m for i,j in bin_pairs]
        drho5 = drho5p + drho5m
        veclistrho5.append(np.concatenate(np.c_[drho5]))
   
    covmatrho0 = np.cov(np.c_[veclistrho0].T)
    covmatrho1 = np.cov(np.c_[veclistrho1].T)
    covmatrho2 = np.cov(np.c_[veclistrho2].T)
    covmatrho3 = np.cov(np.c_[veclistrho3].T)
    covmatrho4 = np.cov(np.c_[veclistrho4].T)
    covmatrho5 = np.cov(np.c_[veclistrho5].T)
 


    al = bestparameters(alist)
    bl = bestparameters(blist)
    el = bestparameters(elist)

    xlabel=r'$\theta$ [arcmin]'
    ylim1 = [1.e-12, 1.e-6]
    ylim2 = [1.e-14, 1.e-6]
    ylim3 = [1.e-12, 1.e-6]
    ylim4 = [1.e-14, 1.e-6]
    fig1, ax1 = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    fig2, ax2 = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    fig3, ax3 = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    fig4, ax4 = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)

    for k,[i,j] in enumerate(bin_pairs):
         drho0p = al[i]*al[j]*rho0p; drho0m = al[i]*al[j]*rho0m
         drho1p = bl[i]*bl[j]*rho1p; drho1m = bl[i]*bl[j]*rho1m
         drho2p = (bl[i]*al[j] + bl[j]*al[i])*rho2p; drho2m = (bl[i]*al[j] + bl[j]*al[i])*rho2m
         drho3p = el[i]*el[j]*rho3p; drho3m = el[i]*el[j]*rho3m
         drho4p = (bl[i]*el[j] + bl[j]*el[i])*rho4p; drho4m = (bl[i]*el[j] + bl[j]*el[i])*rho4m
         drho5p = (el[i]*al[j] +el[j]*al[i])*rho5p; drho5m = (el[i]*al[j] +el[j]*al[i])*rho5m

         infp, supp = [2*k*len(meanr),(2*k+1)*len(meanr)]
         infm, supm = [(2*k + 1)*len(meanr),2*(k+1)*len(meanr)]
         plot_tomograpically_bin(ax1, i + 1, j + 1, meanr,
                                drho0p,
                                yerr=(np.diagonal(covmatrho0)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                color='red', ylim=ylim1)
         plot_tomograpically_bin(ax1, i + 1, j + 1, meanr,
                                drho1p,
                                yerr=(np.diagonal(covmatrho1)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                color='blue', ylim=ylim3)
         plot_tomograpically_bin(ax1, i + 1, j + 1, meanr,
                                drho2p,
                                yerr=(np.diagonal(covmatrho2)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                color='green', ylim=ylim1)
         plot_tomograpically_bin(ax2, i + 1, j + 1, meanr,
                                drho3p,
                                yerr=(np.diagonal(covmatrho3)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                color='red', ylim=ylim3)
         plot_tomograpically_bin(ax2, i+1 , j + 1, meanr,
                                drho4p,
                                yerr=(np.diagonal(covmatrho4)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                color='blue', ylim=ylim3)
         plot_tomograpically_bin(ax2, i + 1, j + 1, meanr,
                                drho5p,
                                yerr=(np.diagonal(covmatrho5)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                color='green', ylim=ylim1)
         
         plot_tomograpically_bin(ax3, i + 1, j + 1, meanr,
                                drho0m,
                                yerr=(np.diagonal(covmatrho0)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color='red', ylim=ylim2)
         plot_tomograpically_bin(ax3, i + 1, j + 1, meanr,
                                drho1m,
                                yerr=(np.diagonal(covmatrho1)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color='blue', ylim=ylim3)
         plot_tomograpically_bin(ax3, i + 1, j + 1, meanr,
                                drho2m,
                                yerr=(np.diagonal(covmatrho2)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color='green', ylim=ylim2)      
         plot_tomograpically_bin(ax4, i + 1, j + 1, meanr,
                                drho3m,
                                yerr=(np.diagonal(covmatrho3)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color='red', ylim=ylim4)
         plot_tomograpically_bin(ax4, i + 1, j + 1, meanr,
                                drho4m,
                                yerr=(np.diagonal(covmatrho4)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color='blue', ylim=ylim4)
         plot_tomograpically_bin(ax4, i + 1, j + 1, meanr,
                                drho5m,
                                yerr=(np.diagonal(covmatrho5)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color='green', ylim=ylim2)


    for ax in [ax1, ax2, ax3, ax4]:
        ax[0][1].set_visible(False); ax[0][2].set_visible(False)
        ax[1][2].set_visible(False); ax[0][3].set_visible(False)
        ax[2][3].set_visible(False); ax[1][3].set_visible(False)
        
    labels1 = [r'$(\alpha^{i}\eta^{j}+\eta^{i}\alpha^{j})\rho_{0+} $', r'$\alpha^{i}\alpha^{j}\rho_{1+}$', r'$(\beta^{i}\alpha^{j}+\alpha^{i}\beta^{j})\rho_{2+}$']
    labels2 = [r'$(\beta^{i}\eta^{j}+\eta^{i}\beta^{j})\rho_{3+} $', r'$\beta^{i}\beta^{j}\rho_{4+}$', r'$\eta^{i}\eta^{j}\rho_{5+}$']
    labels3 = [r'$(\alpha^{i}\eta^{j}+\eta^{i}\alpha^{j})\rho_{0-} $', r'$\alpha^{i}\alpha^{j}\rho_{1-}$', r'$(\beta^{i}\alpha^{j}+\alpha^{i}\beta^{j})\rho_{2-}$']
    labels4 = [r'$(\beta^{i}\eta^{j}+\eta^{i}\beta^{j})\rho_{3-} $', r'$\beta^{i}\beta^{j}\rho_{4-}$', r'$\eta^{i}\eta^{j}\rho_{5-}$']
    labels = [labels1, labels2, labels3, labels4]
    filenames = ['rhosp_a.png', 'rhosm_a.png', 'rhosp_b.png', 'rhosm_b.png']

    for i, fig in enumerate([fig1, fig2, fig3, fig4]):
        red_patch = mpatches.Patch(color='red', label=labels[i][0])
        blue_patch = mpatches.Patch(color='blue', label=labels[i][1])
        green_patch = mpatches.Patch(color='green', label=labels[i][2])
        fig.legend(handles=[ blue_patch, green_patch, red_patch], fontsize=20)
        fig.tight_layout()
        filename = os.path.join(plotspath, filenames[i] )
        print('Printing', filename)
        fig.savefig(filename, dpi=200)
        plt.close(fig)

def plotrhoscontributions2(contaminant, abesamps, rhos_cosmo_file, plotspath, mflags,  nsig=1):
    from src.readfits import  read_rhos
    import itertools
    import fitsio
    from src.maxlikelihood import percentiles, bestparameters

    covmatrixfit_cont=fitsio.read(contaminant,ext=1)
    covmatrixfit_cont *= (nsig**2)
    xipfit_cont=fitsio.read(contaminant,ext=2)
    ximfit_cont=fitsio.read(contaminant,ext=3)
    

    meanr, rhos,  covmat = read_rhos(rhos_cosmo_file)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m = rhos

    if ('a1' in abesamps.dtype.names): alist = [abesamps['a1'], abesamps['a2'], abesamps['a3'], abesamps['a4']]
    else: alist = [[0]*len(abesamps) ]*4; print('ALPHA is not in samples')
    if ('b1' in abesamps.dtype.names): blist = [abesamps['b1'], abesamps['b2'], abesamps['b3'], abesamps['b4']]
    else: blist = [[0]*len(abesamps) ]*4; print('BETA is not in samples')
    if ('e1' in abesamps.dtype.names): elist = [abesamps['e1'], abesamps['e2'], abesamps['e3'], abesamps['e4']]
    else: elist = [[0]*len(abesamps) ]*4; print('ETA is not in samples')

    nbins=4
    a=[i for i in range(nbins)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2): bin_pairs.append(p)

    veclistrho0 = [];    veclistrho1 = []
    veclistrho2 = [];    veclistrho3 = []
    veclistrho4 = [];    veclistrho5 = []
    for z in range(len(abesamps['a1'])):       
        drho0p = [alist[i][z]*alist[j][z]*rho0p for i,j in bin_pairs]
        drho0m = [alist[i][z]*alist[j][z]*rho0m for i,j in bin_pairs]
        drho0 = drho0p + drho0m
        veclistrho0.append(np.concatenate(np.c_[drho0]))

        drho1p = [blist[i][z]*blist[j][z]*rho1p for i,j in bin_pairs]
        drho1m = [blist[i][z]*blist[j][z]*rho1m for i,j in bin_pairs]
        drho1 = drho1p + drho1m
        veclistrho1.append(np.concatenate(np.c_[drho1]))

        drho2p = [(blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2p for i,j in bin_pairs]
        drho2m = [(blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2m for i,j in bin_pairs]
        drho2 = drho2p + drho2m
        veclistrho2.append(np.concatenate(np.c_[drho2]))

        drho3p = [elist[i][z]*elist[j][z]*rho3p for i,j in bin_pairs]
        drho3m = [elist[i][z]*elist[j][z]*rho3m for i,j in bin_pairs]
        drho3 = drho3p + drho3m
        veclistrho3.append(np.concatenate(np.c_[drho3]))

        drho4p = [(blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4p for i,j in bin_pairs]
        drho4m = [(blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4m for i,j in bin_pairs]
        drho4 = drho4p + drho4m
        veclistrho4.append(np.concatenate(np.c_[drho4]))

        drho5p = [(elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5p for i,j in bin_pairs]
        drho5m = [(elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5m for i,j in bin_pairs]
        drho5 = drho5p + drho5m
        veclistrho5.append(np.concatenate(np.c_[drho5]))
   
    covmatrho0 = np.cov(np.c_[veclistrho0].T)
    covmatrho1 = np.cov(np.c_[veclistrho1].T)
    covmatrho2 = np.cov(np.c_[veclistrho2].T)
    covmatrho3 = np.cov(np.c_[veclistrho3].T)
    covmatrho4 = np.cov(np.c_[veclistrho4].T)
    covmatrho5 = np.cov(np.c_[veclistrho5].T)

    covmatrhos = [covmatrho0, covmatrho1, covmatrho2, covmatrho3, covmatrho4, covmatrho5]
 

    
    al = bestparameters(alist)
    bl = bestparameters(blist)
    el = bestparameters(elist)

    af, bf, ef = mflags
        
    labels1 = [ r'$\alpha^{i}\alpha^{j}\rho_{0+}$'*(af), r'$\beta^{i}\beta^{j}\rho_{1+}$'*(bf) ]
    labels2 = [ r'$(\beta^{i}\alpha^{j}+\alpha^{i}\beta^{j})\rho_{2+}$'*(af&bf),r'$\eta^{i}\eta^{j}\rho_{3+}$'*(ef) ]
    labels3 = [ r'$(\beta^{i}\eta^{j}+\eta^{i}\beta^{j})\rho_{4+} $'*(bf&ef),  r'$(\alpha^{i}\eta^{j}+\eta^{i}\alpha^{j})\rho_{5+} $'*(af&ef)]
    labels4 = [ r'$\alpha^{i}\alpha^{j}\rho_{0-}$'*(af), r'$\beta^{i}\beta^{j}\rho_{1-}$'*(bf) ]
    labels5 = [ r'$(\beta^{i}\alpha^{j}+\alpha^{i}\beta^{j})\rho_{2-}$'*(af&bf),r'$\eta^{i}\eta^{j}\rho_{3-}$'*(ef) ]
    labels6 = [ r'$(\beta^{i}\eta^{j}+\eta^{i}\beta^{j})\rho_{4-} $'*(bf&ef),  r'$(\alpha^{i}\eta^{j}+\eta^{i}\alpha^{j})\rho_{5-} $'*(af&ef)]
    labels = [labels1, labels2, labels3, labels4, labels5, labels6]
    filenames = ['contribution_rhosp_a.png', 'contribution_rhosp_b.png','contribution_rhosp_c.png', 'contribution_rhosm_a.png', 'contribution_rhosm_b.png',   'contribution_rhosm_c.png']

    colors = ['red', 'blue', 'green', 'gray', 'fuchsia', 'darkslategray' ]

    xlabel=r'$\theta$ [arcmin]'
    figs = [];  axs = [];
    ylimsp = [[1.e-10, 5.e-7]] * 6
    ylimsm = [[1.e-13, 5.e-7]] * 6
    xlim = [2.5, 250.]
    #ylims = [[1.e-12, 1.e-6],  [1.e-14, 1.e-6],  [1.e-12, 1.e-6],  [1.e-14, 1.e-6], [1.e-12, 1.e-6],  [1.e-14, 1.e-6]]
    for i in range(len(labels)):
        figaux, axaux = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
        figs.append(figaux); axs.append(axaux)
    
 


    lengths_cont = [len(fitsio.read(contaminant,ext=2)), len(fitsio.read(contaminant,ext=3))]
   
    for ax in axs[:3]: 
        for i,j in bin_pairs:
            bin = (xipfit_cont['BIN1']==i+1)&(xipfit_cont['BIN2']==j+1)
            theta_cont=xipfit_cont['ANG'][bin]
            xip_cont=xipfit_cont['VALUE'][bin]
            yerr_cont=get_error(covmatrixfit_cont, lengths_cont, 'xip')[bin]
            plot_tomograpically_bin(ax, i + 1, j + 1, theta_cont, xip_cont,
                                    yerr=yerr_cont, xlabel=xlabel,
                                    ylabel=r'$\xi_{+}(\theta)$', nbins=4, color='black',
                                    xlim=xlim, ylim=None)
    for ax in axs[3: ]: 
        for i,j in bin_pairs:
            bin = (ximfit_cont['BIN1']==i+1)&(ximfit_cont['BIN2']==j+1)
            theta_cont=ximfit_cont['ANG'][bin]
            xim_cont=ximfit_cont['VALUE'][bin]
            yerr_cont=get_error(covmatrixfit_cont, lengths_cont, 'xim')[bin]
            plot_tomograpically_bin(ax, i + 1, j + 1, theta_cont, xim_cont,
                                    yerr=yerr_cont, xlabel=xlabel,
                                    ylabel=r'$\xi_{-}(\theta)$', nbins=4, color='black',
                                    xlim=xlim, ylim=None)
    

    
    for k,[i,j] in enumerate(bin_pairs):
         drho0p = al[i]*al[j]*rho0p; drho0m = al[i]*al[j]*rho0m
         drho1p = bl[i]*bl[j]*rho1p; drho1m = bl[i]*bl[j]*rho1m
         drho2p = (bl[i]*al[j] + bl[j]*al[i])*rho2p; drho2m = (bl[i]*al[j] + bl[j]*al[i])*rho2m
         drho3p = el[i]*el[j]*rho3p; drho3m = el[i]*el[j]*rho3m
         drho4p = (bl[i]*el[j] + bl[j]*el[i])*rho4p; drho4m = (bl[i]*el[j] + bl[j]*el[i])*rho4m
         drho5p = (el[i]*al[j] +el[j]*al[i])*rho5p; drho5m = (el[i]*al[j] +el[j]*al[i])*rho5m

         drhosp = [drho0p, drho1p, drho2p, drho3p, drho4p, drho5p]
         drhosm = [drho0m, drho1m, drho2m, drho3m, drho4m, drho5m]

         infp, supp = [2*k*len(meanr),(2*k+1)*len(meanr)]
         infm, supm = [(2*k + 1)*len(meanr),2*(k+1)*len(meanr)]

         idx = [0, 0, 1, 1, 2, 2]
         for z, covmat in enumerate(covmatrhos):
             plot_tomograpically_bin(axs[idx[z]], i + 1, j + 1, meanr,
                                drhosp[z],
                                yerr=(np.diagonal(covmat)**0.5)[infp:supp],
                                xlabel=xlabel, ylabel=r'$\xi_{+}(\theta)$', nbins=4,
                                     color=colors[z], xlim=xlim, ylim=ylimsp[z])
             plot_tomograpically_bin(axs[idx[z]+3], i + 1, j + 1, meanr,
                                drhosm[z],
                                yerr=(np.diagonal(covmat)**0.5)[infm:supm],
                                xlabel=xlabel, ylabel=r'$\xi_{-}(\theta)$', nbins=4,
                                color=colors[z], xlim=xlim, ylim=ylimsm[z])

    for ax in axs:
        ax[0][1].set_visible(False); ax[0][2].set_visible(False)
        ax[1][2].set_visible(False); ax[0][3].set_visible(False)
        ax[2][3].set_visible(False); ax[1][3].set_visible(False)

   
    for i, fig in enumerate(figs):
        black_patch  = plt.errorbar(0, 0,   yerr=0.5, fmt='o' ,capsize=1,markersize=7, color='black', label=r'$\delta \xi $')
        red_patch  = plt.errorbar(0, 0,  yerr=0.5, fmt='o' ,capsize=1,markersize=7,  color=(colors + colors)[2*i], label=labels[i][0])
        blue_patch = plt.errorbar(0, 0, yerr=0.5, fmt='o' ,capsize=1,markersize=7,  color=(colors + colors)[2*i + 1] , label=labels[i][1])
        fig.legend(handles=[ red_patch, blue_patch, black_patch], fontsize=20,  loc='upper right') #bbox_to_anchor=(0.6,0.0),
        fig.tight_layout()
        #fig.tight_layout(rect=[0, 0, 1, 0.95])
        filename = os.path.join(plotspath, filenames[i] )
        print('Printing', filename)
        fig.savefig(filename, dpi=200)
        plt.close(fig)
        
def main():
    import fitsio
    import getdist
    from getdist import plots, MCSamples
    from astropy.io import fits
    import itertools
    from src.maxlikelihood import percentiles, bestparameters
    


    args = parse_args()

    #Make directory where the ouput data will be
    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise

    abesamps = fitsio.read(args.samplesabe)
    #Burning just to be sure
    frac= args.burn_frac
    abesamps =abesamps[int(frac*len(abesamps) ): ]

    af, bf, ef = [name in abesamps.dtype.names for name in ['a1', 'b1', 'e1'] ]
    ext = '%s%s%s'%('a'*af,'b'*bf, 'e'*ef)

  
    title_covmatsample = r'$ %s %s %s %s %s $'%(r'\alpha'*af,r'\mid'*(af&bf),r'\beta'*bf, r'\mid'*(ef&bf), r'\eta'*ef)
    filename_covmatsample = 'Covariancematrix_samples_%s.png'%(ext)
    filename_contoursabe =  'contours_samples_%s.png'%(ext)
    labels_contourabe = np.concatenate([ np.array([r'\alpha^{(%d)}'%(i),r'\beta^{(%d)}'%(i),r'\eta^{(%d)}'%(i)])[[af, bf, ef]]  for i in range(1, 5)])

    filename_covmat2pcfpsfbias = 'Covariancematrix_totalbiases_%s.png'%(ext)
    filename_covmat2pcfpsfxipbias = 'Covariancematrix_totalbiasxip_%s.png'%(ext)
    filename_covmat2pscfpsfbiasxip11 = 'Covariancematrix_totalbiasxip11_%s.png'%(ext)

    title_contfid = None
    #title_contfid = '%s %s %s'%('Alpha'*af , 'Beta'*bf,'Eta'*ef)
    filename_xip_contfid = 'xipcont_xipfid_%s_%dsig.png'%(ext, args.nsig)
    filename_xim_contfid = 'ximcont_ximfid_%s_%dsig.png'%(ext, args.nsig)


    veclist = [abesamps[name] for name in abesamps.dtype.names]
    covmat = np.cov(veclist)
    
    #PLOT COVMAT
    

    plt.clf()
    lengths = [4, 4, 4]
    plotcorrmat(covmat)
    plt.title(title_covmatsample)
    pos_lines = [ -0.5]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = os.path.join(plotspath, filename_covmatsample)
    print('Printing', filename)
    plt.savefig(filename, dpi=200)
    print('Printed', filename)
        
    #PLOT marginalized and contours
    if(args.corner):
        filename = os.path.join(plotspath,'corner_%s'%(filename_contoursabe))
        #corner_plot(veclist, ['a1', 'b1', 'e1'], filename, title=None)
        corner_plot(veclist, abesamps.dtype.names, filename, title=None)
    else:
        
        #mcmcpars = percentiles(veclist, nsig=1) 
        #print( ' mcmc parameters xi+',  'nsig=', 1, ' percentiles: ',  mcmcpars)
        samples = MCSamples(samples=veclist, names=abesamps.dtype.names, labels=labels_contourabe)
      
        g = plots.getSubplotPlotter()
        g.settings.plot_meanlikes = False
        g.settings.alpha_factor_contour_lines = True
        #g.settings.axis_marker_lw = 5
        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add=0.35
        g.settings.title_limit_fontsize = 16
        g.settings.figure_legend_loc = 'best'
        g.settings.rcSizes(axes_fontsize = 12, lab_fontsize=20, legend_fontsize =40)
        g.triangle_plot([samples], filled_compare=[True], 
                    line_args=[{'ls':'solid', 'lw':2, 'color':'green'}],
                    contour_colors=['green'], title_limit=1)
        #g.add_legend(legend_labels=[legend_name], fontsize=36, legend_loc=(-3.5,7))
        filename = os.path.join(plotspath,filename_contoursabe)
        print('Printing', filename)
        plt.tight_layout()
        plt.savefig(filename, dpi=200)
        print('Printed', filename)

    ##PSF BIAS COV MAT
    covmatrixfit_cont=fitsio.read(args.contaminant,ext=1)
    covmatrixfit_cont *= (args.nsig**2)        
    lengths = [len(fitsio.read(args.contaminant,ext=2)), len(fitsio.read(args.contaminant,ext=3))]

    covmat = covmatrixfit_cont 
    plt.clf()
    plt.cla()
    plt.close()

    plotcorrmat(covmat)
    plt.title(r'$\delta \xi_{+}^{PSF} \mid \delta \xi_{-}^{PSF} $')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    
    filename = os.path.join(plotspath,filename_covmat2pcfpsfbias)
    print('Printing', filename)
    plt.savefig(filename, dpi=200)

    plt.clf()
    plt.cla()
    plt.close()
    plotcorrmat(covmat[:200, :200])
    plt.title(r'$\delta \xi_{+}^{PSF}$')
    plt.tight_layout()
    filename = os.path.join(plotspath, filename_covmat2pcfpsfxipbias)
    print('Printing', filename)
    plt.savefig(filename, dpi=200)

    plt.clf()
    plt.cla()
    plt.close()
    plotcorrmat(covmat[:20, :20])
    plt.title(r'$\delta \xi_{+}^{PSF} (1,1) $')
    plt.tight_layout()
    
    filename = os.path.join(plotspath,filename_covmat2pscfpsfbiasxip11)
    print('Printing', filename)
    plt.savefig(filename, dpi=200)
    
    ##CONTAMINATION VC FIDUCIAL
    plotcontaminantandfiducial(args.contaminant, args.fiducial, plotspath, title=title_contfid, filenames=[filename_xip_contfid,filename_xim_contfid], nsig=args.nsig )
    

    ##RHOS CONTRIBUTION IN THE TOTAL PSF BIAS
    plotrhoscontributions2(args.contaminant, abesamps, args.rhoscosmo, plotspath, [af, bf, ef], nsig=args.nsig)

   
  
    
if __name__ == "__main__":
    main()



        
