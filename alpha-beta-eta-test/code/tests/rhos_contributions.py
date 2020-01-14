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
    parser.add_argument('--samplesabe', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/samples_abe.fits')
    parser.add_argument('--rhoscosmo', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS_Y3.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--plots', default=False,
                        action='store_const', const=True, help='Plot correlations functions')
    
    parser.add_argument('--outpath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the final contaminant')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    parser.add_argument('--filename',
                        default='marg_abe_eq4_contaminant.fits',
                        help='Filename of the final contamination')

    args = parser.parse_args()

    return args
                        
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
                        
def main():
    import sys; sys.path.append(".")
    import fitsio
    import getdist
    from getdist import plots, MCSamples
    from astropy.io import fits
    import itertools
    from src.maxlikelihood import percentiles, bestparameters
    from src.readfits import  read_rhos
    from src.chi2 import chi2nu
    import numpy as np

    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise

    abesamps = fitsio.read(args.samplesabe)
    #Burning just to be sure
    frac= 0.0
    abesamps =abesamps[int(frac*len(abesamps) ): ]
  
        
        
    ##SAMPLING BIAS
    meanr, rhos,  covmat = read_rhos(args.rhoscosmo)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m = rhos

    alist = [abesamps['a1'], abesamps['a2'], abesamps['a3'], abesamps['a4']]
    blist = [abesamps['b1'], abesamps['b2'], abesamps['b3'], abesamps['b4']]
    elist = [abesamps['e1'], abesamps['e2'], abesamps['e3'], abesamps['e4']]

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



         

   
  
    
if __name__ == "__main__":
    main()



        
