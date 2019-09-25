import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    parser.add_argument('--taus',
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/taus_zbin1_large_scales.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/taus_zbin2_large_scales.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/taus_zbin3_large_scales.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/taus_zbin4_large_scales.fits'], 
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args

def pretty_tau_ax(ax,  meanr, rho, sig,  legend=None, lfontsize=24, color='black', marker='o', ylabel=r'$\rho(\theta)$',title=None,  xlim=None,  ylim=None):
    import numpy as np
     
    if sig is None: sig =  np.zeros(len(rho))
    ax.plot(meanr, rho, color=color, label=legend)
    ax.plot(meanr, -rho, color=color, ls=':')
    #rigth quadrants
    ax.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker,  capsize=2)
    ax.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker,  capsize=2)
    #leftquadrants
    ax.errorbar( -meanr, rho, yerr=sig, color=color,  marker='^',  capsize=2)
    ax.errorbar( -meanr,-rho, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
    ax.legend(loc='best', fontsize=lfontsize)
    if ylim is not None: ax.ylim( ylim )
    if xlim is not None: ax.xlim(xlim)
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.set_xlabel(r'$\theta$ (arcmin)', fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_xscale('log')
    ax.set_yscale('log', nonposy='clip')
    if title is not None: ax.set_title(title)
  

def main():
    import sys; sys.path.append(".")
    from src.readfits import read_taus
    from src.plot_stats import pretty_rho

    import numpy as np
        
    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()

    colors=['red','green','blue','black']
    for i, taubin in enumerate(args.taus):
       meanr, taus,  covtau =  read_taus(taubin)
       nrows = len(meanr)
       tau0p, tau0m, tau2p, tau2m, tau5p, tau5m = taus
       tau0psig =  np.sqrt(np.diag(covtau[0:nrows, 0:nrows]))
       tau0msig =  np.sqrt(np.diag(covtau[nrows:2*nrows , nrows:2*nrows]))
       tau2psig =  np.sqrt(np.diag(covtau[2*nrows:3*nrows, 2*nrows:3*nrows]))
       tau2msig =  np.sqrt(np.diag(covtau[3*nrows:4*nrows, 3*nrows:4*nrows]))
       tau5psig =  np.sqrt(np.diag(covtau[4*nrows:5*nrows, 4*nrows:5*nrows]))
       tau5msig =  np.sqrt(np.diag(covtau[5*nrows:6*nrows, 5*nrows:6*nrows]))
       pretty_tau_ax(ax1,  meanr, tau0p, tau0psig, legend='Bin %d'%(i + 1), lfontsize=14, color=colors[i], marker='o', ylabel=r'$\tau_{0+}$')
       pretty_tau_ax(ax2,  meanr, tau0m, tau0msig, legend='Bin %d'%(i + 1), lfontsize=14, color=colors[i], marker='o', ylabel=r'$\tau_{0-}$')
       pretty_tau_ax(ax3,  meanr, tau2p, tau2psig, legend='Bin %d'%(i + 1), lfontsize=14, color=colors[i], marker='o', ylabel=r'$\tau_{2+}$')
       pretty_tau_ax(ax4,  meanr, tau2m, tau2msig, legend='Bin %d'%(i + 1), lfontsize=14, color=colors[i], marker='o', ylabel=r'$\tau_{2-}$')
       pretty_tau_ax(ax5,  meanr, tau5p, tau5psig, legend='Bin %d'%(i + 1), lfontsize=14, color=colors[i], marker='o', ylabel=r'$\tau_{5+}$')
       pretty_tau_ax(ax6,  meanr, tau5m, tau5msig, legend='Bin %d'%(i + 1), lfontsize=14, color=colors[i], marker='o', ylabel=r'$\tau_{5-}$')
    exts = ['tau0p', 'tau0m', 'tau2p', 'tau2m', 'tau5p', 'tau5m']
    for i, fig  in enumerate([fig1, fig2, fig3, fig4, fig5, fig6]):
        fig.tight_layout()
        filename = os.path.join(plotspath,  '%s_metacal.png'%(exts[i]) )
        fig.savefig(filename, dpi=500)
        plt.close(fig)
        print(filename, 'Printed!')
    


if __name__ == "__main__":
    main()
