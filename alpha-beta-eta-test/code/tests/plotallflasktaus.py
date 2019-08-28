import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    parser.add_argument('--tausmetacal',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_4.fits',
                        help='Full Path to the taus measurement of metcal')
    parser.add_argument('--tausflask',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/taus_g1g2flip/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--zbin', default=4 , type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--plots', default=True,
                        action='store_const', const=True, help='Plot correlations functions')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args
def plotcorrmat(cov, zbin):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    clb = plt.colorbar()
    clb.ax.set_title("zbin" + str(zbin))
    plt.tight_layout()
 
def plotflask(tausmetacal, zbin, tausflask, plotspath):
    from src.readfits import read_taus
    import numpy as np
    meanr, taus,  covtau =  read_taus(tausmetacal)
    nrows = len(meanr)
    tau0p, tau0m, tau2p, tau2m, tau5p, tau5m = taus
    tau0pvar_meta =  np.diag(covtau[0:nrows, 0:nrows])
    tau0mvar_meta =  np.diag(covtau[nrows:2*nrows , nrows:2*nrows])
    tau2pvar_meta =  np.diag(covtau[2*nrows:3*nrows, 2*nrows:3*nrows])
    tau2mvar_meta =  np.diag(covtau[3*nrows:4*nrows, 3*nrows:4*nrows])
    tau5pvar_meta =  np.diag(covtau[4*nrows:5*nrows, 4*nrows:5*nrows])
    tau5mvar_meta =  np.diag(covtau[5*nrows:6*nrows, 5*nrows:6*nrows])

    veclist = []
    count = 0;
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()
    
    for seed in range(1, 500 ):
        for ck in range(1, 2):
            name = os.path.join(tausflask, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            exist =  os.path.isfile(name)
            if exist:
                meanr, taus, covtaus = read_taus(name)
                if (np.count_nonzero(taus) == 0):
                    print("Warning, weird measurement, skipping", name)
                else:
                    ax1.plot(meanr, taus[0], 'k',  color='blue', alpha=0.15)
                    ax1.set_ylabel(r'$\tau_{0+}$')
                    ax1.set_xscale('log')
                    ax2.plot(meanr, taus[1], 'k',  color='blue', alpha=0.15)
                    ax2.set_ylabel(r'$\tau_{0-}$')
                    ax2.set_xscale('log')
                    ax3.plot(meanr, taus[2], 'k',  color='blue', alpha=0.15)
                    ax3.set_ylabel(r'$\tau_{2+}$')
                    ax3.set_xscale('log')
                    ax4.plot(meanr, taus[3], 'k',  color='blue', alpha=0.15)
                    ax4.set_ylabel(r'$\tau_{2-}$')
                    ax4.set_xscale('log')    
                    ax5.plot(meanr, taus[4], 'k',  color='blue', alpha=0.15)
                    ax5.set_ylabel(r'$\tau_{5+}$')
                    ax5.set_xscale('log')
                    ax6.plot(meanr, taus[5], 'k',  color='blue', alpha=0.15)
                    ax6.set_ylabel(r'$\tau_{5-}$')
                    ax6.set_xscale('log')

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
    sig_tau0p = np.sqrt(np.diag(covmat[0:nrows, 0:nrows]))
    sig_tau0m = np.sqrt(np.diag(covmat[nrows:2*nrows, nrows:2*nrows]))
    sig_tau2p = np.sqrt(np.diag(covmat[2*nrows:3*nrows, 2*nrows:3*nrows]))
    sig_tau2m = np.sqrt(np.diag(covmat[3*nrows:4*nrows, 3*nrows:4*nrows]))
    sig_tau5p = np.sqrt(np.diag(covmat[4*nrows:5*nrows, 4*nrows:5*nrows]))
    sig_tau5m = np.sqrt(np.diag(covmat[5*nrows:6*nrows, 5*nrows:6*nrows]))

    ax1.errorbar(meanr, tau0pmean, yerr=sig_tau0p, color='red', ls='', marker='.',  capsize=2, label='Mean tau flask')
    ax1.errorbar(meanr, tau0p, yerr=np.sqrt(tau0pvar_meta), color='green', ls='', marker='.',  capsize=2, label='tau metacal')
    ax1.legend(loc='best', fontsize=10) 
    fig1.tight_layout()
    filename = plotspath + 'taus0pflask_zbin%d%s'%(zbin, '.png')
    fig1.savefig(filename, dpi=500)
    plt.close(fig1)
    print(filename, 'Printed!')
    
    ax2.errorbar(meanr, tau0mmean, yerr=sig_tau0m, color='red', ls='', marker='.',  capsize=2, label='Mean tau flask')
    ax2.errorbar(meanr, tau0m, yerr=np.sqrt(tau0mvar_meta), color='green', ls='', marker='.',  capsize=2, label='tau metacal')
    ax2.legend(loc='best', fontsize=10) 
    fig2.tight_layout()
    filename = plotspath + 'taus0mflask_zbin%d%s'%(zbin, '.png')
    fig2.savefig(filename, dpi=500)
    plt.close(fig2)
    print(filename, 'Printed!')

    ax3.errorbar(meanr, tau2pmean, yerr=sig_tau2p, color='red', ls='', marker='.',  capsize=2, label='Mean tau flask')
    ax3.errorbar(meanr, tau2p, yerr=np.sqrt(tau2pvar_meta), color='green', ls='', marker='.',  capsize=2, label='tau metacal')
    ax3.legend(loc='best', fontsize=10) 
    fig3.tight_layout()
    filename = plotspath + 'taus2pflask_zbin%d%s'%(zbin, '.png')
    fig3.savefig(filename, dpi=500)
    plt.close(fig3)
    print(filename, 'Printed!')
    
    ax4.errorbar(meanr, tau2mmean, yerr=sig_tau2m, color='red', ls='', marker='.',  capsize=2, label='Mean tau flask')
    ax4.errorbar(meanr, tau2m, yerr=np.sqrt(tau2mvar_meta), color='green', ls='', marker='.',  capsize=2, label='tau metacal')
    ax4.legend(loc='best', fontsize=10) 
    fig4.tight_layout()
    filename = plotspath + 'taus2mflask_zbin%d%s'%(zbin, '.png')
    fig4.savefig(filename, dpi=500)
    plt.close(fig4)
    print(filename, 'Printed!')

    ax5.errorbar(meanr, tau5pmean, yerr=sig_tau5p, color='red', ls='', marker='.',  capsize=2, label='Mean tau flask')
    ax5.errorbar(meanr, tau5p, yerr=np.sqrt(tau5pvar_meta), color='green', ls='', marker='.',  capsize=2, label='tau metacal')
    ax5.legend(loc='best', fontsize=10) 
    fig5.tight_layout()
    filename = plotspath + 'taus5pflask_zbin%d%s'%(zbin, '.png')
    fig5.savefig(filename, dpi=500)
    plt.close(fig5)
    print(filename, 'Printed!')

    ax6.errorbar(meanr, tau5mmean, yerr=sig_tau5m, color='red', ls='', marker='.',  capsize=2, label='Mean tau flask')
    ax6.errorbar(meanr, tau5m, yerr=np.sqrt(tau5mvar_meta), color='green', ls='', marker='.',  capsize=2, label='tau metacal')
    ax6.legend(loc='best', fontsize=10) 
    fig6.tight_layout()
    filename = plotspath + 'taus5mflask_zbin%d%s'%(zbin, '.png')
    fig6.savefig(filename, dpi=500)
    plt.close(fig6)
    print(filename, 'Printed!')
    
def main():
    import sys; sys.path.append(".")
    from src.readfits import read_taus
    from src.plot_stats import pretty_rho
    import fitsio
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    import itertools
    import numpy as np
        
    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    plotflask(args.tausmetacal, args.zbin, args.tausflask, plotspath )
    
    veclist = [];     count = 0; zbin = args.zbin;
    for seed in range(1, 500 ):
        for ck in range(1, 2):
            name = os.path.join(args.tausflask, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            exist =  os.path.isfile(name)
            if exist:
                meanr, taus, covtaus = read_taus(name)
                if (np.count_nonzero(taus) == 0):
                    print("Warning, weird measurement, skipping", name)
                else:
                    veclist.append(np.concatenate(np.c_[taus]))
                    count +=1
    print(count, "FLASK catalogs were read")
    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
    

    lengths = [len(taus[0]), len(taus[1]), len(taus[2]), len(taus[3]),
               len(taus[4]), len(taus[5])]

    plt.clf()
    plt.cla()
    plotcorrmat(covmat, zbin)
    plt.title(r'$\tau_{0+}(\theta) \mid \tau_{0-}(\theta) \mid \tau_{2+}(\theta) \mid \tau_{2-}(\theta) \mid \tau_{5+}(\theta) \mid \tau_{5-}(\theta) $')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = plotspath + 'CovariancematrixTausFlask_%d%s'%(args.zbin, '.png')
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')



    
if __name__ == "__main__":
    main()

