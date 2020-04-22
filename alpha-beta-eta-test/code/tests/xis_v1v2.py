import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Comparing xis for different flips')
    parser.add_argument('--xisflask1',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/xis_noflip/',
                        help='Full Path to the taus measurement with flask version 1')
    parser.add_argument('--xisflask2',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/xis_g1flip/',
                        help='Full Path to the taus measurement with flask version 2')
    parser.add_argument('--zbin', default=4 , type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args

def plotflask(axs, zbin, xisflask, plotspath, color, label):
    from src.readfits import read_xis
    import numpy as np

    ax1, ax2 = axs
    
    veclist = []
    count = 0
   
    for seed in range(1, 401 ):
        for ck in range(1, 2):
            name = os.path.join(xisflask, 'xis_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            exist =  os.path.isfile(name)
            if exist:
                meanr, xis, covxis = read_xis(name)
                if (np.count_nonzero(xis) == 0):
                    print("Warning, weird measurement, skipping", name)
                else:
                    veclist.append(np.concatenate(np.c_[xis]))
                    count +=1

    print(count, "FLASK catalogs were read")
    meanvec = np.mean(veclist, axis=0)
    
    nrows = len(meanr)
    xipmean =  meanvec[0:nrows]
    ximmean =  meanvec[nrows:2*nrows]

    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
    print('matrix covariance shape', covmat.shape)
    sig_xip = np.sqrt(np.diag(covmat[0:nrows, 0:nrows]))
    sig_xim = np.sqrt(np.diag(covmat[nrows:2*nrows, nrows:2*nrows]))

    ximeans = [xipmean,ximmean ]
    sig_xis = [sig_xip,sig_xim ]

    ylabels = [r'$\xi_{+}$', r'$\xi_{-}$']; xlabels = [r'$\theta$', r'$\theta$']
    for i, ax in enumerate(axs):
        ax.errorbar(meanr,ximeans[i],yerr=sig_xis[i],color=color, ls='', marker='.', capsize=2, label=label)
        ax.legend(loc='best', fontsize=18, frameon=True)
        ax.set_ylabel(ylabels[i]); ax.set_xlabel(xlabels[i])
        ax.set_xscale('log')
        #ax.set_yscale('log')
        #ax.set_ylim([ -2.e-6,2.e-6 ])

def main():
    import sys; sys.path.append(".")
    from src.readfits import read_xis
    import numpy as np
    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise
        
    if False:
        ck = 1; seed = 1; zbin = args.zbin
        name1 = os.path.join(args.xisflask1, 'xis_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
        name2 = os.path.join(args.xisflask2, 'xis_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
        meanr, xis1, covxis1 = read_xis(name1)
        meanr, xis2, covxis2 = read_xis(name2)
        nrows = len(meanr)
        sig_xip1 = np.sqrt(np.diag(covxis1[0:nrows, 0:nrows]))
        sig_xim1 = np.sqrt(np.diag(covxis1[nrows:2*nrows, nrows:2*nrows]))
        sig_xip2 = np.sqrt(np.diag(covxis2[0:nrows, 0:nrows]))
        sig_xim2 = np.sqrt(np.diag(covxis2[nrows:2*nrows, nrows:2*nrows]))
        
        plt.errorbar(meanr,xis1[1],yerr=sig_xim1,color='red', ls='', marker='.', capsize=2, label='No flip')
        plt.errorbar(meanr,xis2[1],yerr=sig_xim2,color='blue', ls='', marker='.', capsize=2, label='One component flip')
        plt.legend(loc='best', fontsize=25)
        plt.ylabel(r'$\xi_{-}$'); plt.xlabel(r'$\theta$')
        plt.xscale('log')
        filename = os.path.join(args.plotspath, 'xism_flask_zbin%d_onefile.png'%(args.zbin))
        plt.savefig(filename, dpi=200)
   
    
    

    # Mean of samples 
    figs = [];  axs = []; filenames = []
    names = ['xisp', 'xism']
    for i in range(len(names)):
        figaux, axaux = plt.subplots()
        figs.append(figaux); axs.append(axaux)
        filenames.append(os.path.join(plotspath,'%s_flask_zbin%d%s'%(names[i], args.zbin, '.png') ))

        
    plotflask(axs, args.zbin, args.xisflask1, args.plotspath, 'red', 'No flip')
    plotflask(axs, args.zbin, args.xisflask2, args.plotspath, 'blue', 'One component flip')

    for i, fig in enumerate(figs):
        fig.tight_layout()
        fig.savefig(filenames[i], dpi=200)
        plt.close(fig)
        print(filenames[i], 'Printed!')
    
if __name__ == "__main__":
    main()
