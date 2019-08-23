import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    parser.add_argument('--tausflask1',
                        default='/home2/dfa/sobreira/alsina/catalogs/flask/taus/',
                        help='Full Path to the taus measurement with flask version 1')
    parser.add_argument('--tausflask2',
                        default='/home2/dfa/sobreira/alsina/catalogs/flask/taus_v2/',
                        help='Full Path to the taus measurement with flask version 2')
    parser.add_argument('--zbin', default=4 , type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args

def plotflask(axs, zbin, tausflask, plotspath, color, label):
    from src.readfits import read_taus
    import numpy as np

    ax1, ax2, ax3, ax4, ax5, ax6 = axs
    
    veclist = []
    count = 0
   
    for seed in range(1, 401 ):
        for ck in range(1, 2):
            name = os.path.join(tausflask, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            exist =  os.path.isfile(name)
            if exist:
                meanr, taus, covtaus = read_taus(name)
                if (np.count_nonzero(taus) == 0):
                    print("Warning, weird measurement, skipping", name)
                else:
                    veclist.append(np.concatenate(np.c_[taus]))
                    count +=1           
    print(count, "FLASK catalogs were read")
    meanvec = np.mean(veclist, axis=0)
    nrows = len(meanr)
    tau0pmean =  meanvec[0:nrows]
    tau0mmean =  meanvec[nrows:2*nrows]
    tau2pmean =  meanvec[2*nrows:3*nrows]
    tau2mmean =  meanvec[3*nrows:4*nrows]
    tau5pmean =  meanvec[4*nrows:5*nrows]
    tau5mmean =  meanvec[5*nrows:6*nrows]

    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
    print('matrix covariance shape', covmat.shape)
    sig_tau0p = np.sqrt(np.diag(covmat[0:nrows, 0:nrows]))
    sig_tau0m = np.sqrt(np.diag(covmat[nrows:2*nrows, nrows:2*nrows]))
    sig_tau2p = np.sqrt(np.diag(covmat[2*nrows:3*nrows, 2*nrows:3*nrows]))
    sig_tau2m = np.sqrt(np.diag(covmat[3*nrows:4*nrows, 3*nrows:4*nrows]))
    sig_tau5p = np.sqrt(np.diag(covmat[4*nrows:5*nrows, 4*nrows:5*nrows]))
    sig_tau5m = np.sqrt(np.diag(covmat[5*nrows:6*nrows, 5*nrows:6*nrows]))

    taumeans = [tau0pmean,tau0mmean,tau2pmean,tau2mmean,tau5pmean,tau5mmean ]
    sig_taus = [sig_tau0p,sig_tau0m,sig_tau2p,sig_tau2m,sig_tau5p,sig_tau5m  ]

    ylabels = [r'$\tau_{0+}$', r'$\tau_{0-}$', r'$\tau_{2+}$', r'$\tau_{2-}$', r'$\tau_{5+}$', r'$\tau_{5-}$']
    for i, ax in enumerate(axs):
        ax.errorbar(meanr,taumeans[i],yerr=sig_taus[i],color=color, ls='', marker='.', capsize=2, label=label)
        ax.legend(loc='best', fontsize=10)
        ax.set_ylabel(ylabels[i]); ax1.set_xlabel(r'$\theta$')
        ax.set_xscale('log')

def main():
    from src.readfits import read_taus
    import numpy as np
    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise
        

    figs = [];  axs = []; filenames = []
    names = ['taus0p', 'taus0m', 'taus2p', 'taus2m' , 'taus5p' , 'taus5m']
    for i in range(6):
        figaux, axaux = plt.subplots()
        figs.append(figaux); axs.append(axaux)
        filenames.append(os.path.join(plotspath,'%s_flask_zbin%d%s'%(names[i], args.zbin, '.png') ))

        
    plotflask(axs, args.zbin, args.tausflask1, args.plotspath, 'red', 'Taus flask v1')
    plotflask(axs, args.zbin, args.tausflask2, args.plotspath, 'blue', 'Taus flask v2')

    for i, fig in enumerate(figs):
        fig.tight_layout()
        fig.savefig(filenames[i], dpi=500)
        plt.close(fig)
        print(filenames[i], 'Printed!')
    
if __name__ == "__main__":
    main()
