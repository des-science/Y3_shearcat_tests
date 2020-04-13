import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Comparing errors of data with error of model i.é rho errors versus tau errors')
    parser.add_argument('--rhos',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/rhos_large_scales.fits',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    parser.add_argument('--taus',
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_1.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_2.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_3.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_4.fits'],
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    args = parser.parse_args()

    return args

def main():
    import numpy as np
    import sys; sys.path.append(".")
    from src.readfits import read_taus, read_rhos
    
    args = parse_args()

    #Make directory where the ouput data will be
    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    
    zbin = 1
    meanr, taus,  covtau =  read_taus(args.taus[zbin-1])
    nrows = len(meanr)
    tau0p, tau0m, tau2p, tau2m, tau5p, tau5m = taus
    tau0psig =  np.sqrt(np.diag(covtau[0:nrows, 0:nrows]))
    tau0msig =  np.sqrt(np.diag(covtau[nrows:2*nrows , nrows:2*nrows]))
    tau2psig =  np.sqrt(np.diag(covtau[2*nrows:3*nrows, 2*nrows:3*nrows]))
    tau2msig =  np.sqrt(np.diag(covtau[3*nrows:4*nrows, 3*nrows:4*nrows]))
    tau5psig =  np.sqrt(np.diag(covtau[4*nrows:5*nrows, 4*nrows:5*nrows]))
    tau5msig =  np.sqrt(np.diag(covtau[5*nrows:6*nrows, 5*nrows:6*nrows]))

    meanr, rhos,  covrho =  read_rhos(args.rhos)
    nrows = len(meanr)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m =  rhos
    sig_rhos = [np.sqrt(np.diag(covrho[i*nrows:(i + 1)*nrows, i*nrows:(i + 1)*nrows])) for i in range(len(rhos)) ]
    sig_rho0p, sig_rho0m, sig_rho1p, sig_rho1m, sig_rho2p, sig_rho2m, sig_rho3p, sig_rho3m, sig_rho4p, sig_rho4m, sig_rho5p, sig_rho5m =  sig_rhos 


    #First equation contrast
    plt.clf()
    plt.plot(meanr, tau0psig, color='blue', marker='o', label=r'VAR($\tau_{0+}$)')
    plt.plot(meanr, sig_rho0p, color='red', marker='o', label=r'VAR($\rho_{0+}$)')
    plt.plot(meanr, sig_rho2p, color='green', marker='o', label=r'VAR($\rho_{2+}$)')
    plt.plot(meanr, sig_rho5p, color='gray', marker='o', label=r'VAR($\rho_{5+}$)')
    plt.ylim([1.e-11, 2.e-5])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)  
    plt.ylabel('Error', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.legend(loc='best', fontsize=18)
    plt.tight_layout()
    filename = os.path.join(plotspath, 'sigmas_eq1p.png')
    plt.savefig(filename)
    print(filename, 'printed')

    plt.clf()
    plt.plot(meanr, tau0msig, color='blue', marker='o', label=r'VAR($\tau_{0-}$)')
    plt.plot(meanr, sig_rho0m, color='red', marker='o', label=r'VAR($\rho_{0-}$)')
    plt.plot(meanr, sig_rho2m, color='green', marker='o', label=r'VAR($\rho_{2-}$)')
    plt.plot(meanr, sig_rho5m, color='gray', marker='o', label=r'VAR($\rho_{5-}$)')
    plt.ylim([1.e-11, 2.e-5])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)  
    plt.ylabel('Error', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.legend(loc='best', fontsize=18)
    plt.tight_layout()
    filename = os.path.join(plotspath, 'sigmas_eq1m.png')
    plt.savefig(filename)
    print(filename, 'printed')

    #Second equation contrast
    plt.clf()
    plt.plot(meanr, tau2psig, color='blue', marker='o', label=r'VAR($\tau_{2+}$)')
    plt.plot(meanr, sig_rho2p, color='red', marker='o', label=r'VAR($\rho_{2+}$)')
    plt.plot(meanr, sig_rho1p, color='green', marker='o', label=r'VAR($\rho_{1+}$)')
    plt.plot(meanr, sig_rho4p, color='gray', marker='o', label=r'VAR($\rho_{4+}$)')
    plt.ylim([1.e-12, 3.e-5])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)  
    plt.ylabel('Error', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.legend(loc='best', fontsize=18)
    plt.tight_layout()
    filename = os.path.join(plotspath, 'sigmas_eq2p.png')
    plt.savefig(filename)
    print(filename, 'printed')

    plt.clf()
    plt.plot(meanr, tau2msig, color='blue', marker='o', label=r'VAR($\tau_{2-}$)')
    plt.plot(meanr, sig_rho2m, color='red', marker='o', label=r'VAR($\rho_{2-}$)')
    plt.plot(meanr, sig_rho1m, color='green', marker='o', label=r'VAR($\rho_{1-}$)')
    plt.plot(meanr, sig_rho4m, color='gray', marker='o', label=r'VAR($\rho_{4-}$)')
    plt.ylim([1.e-12, 3.e-5])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)  
    plt.ylabel('Error', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.legend(loc='best', fontsize=18)
    plt.tight_layout()
    filename = os.path.join(plotspath, 'sigmas_eq2m.png')
    plt.savefig(filename)
    print(filename, 'printed')

    #Third equation contrast
    plt.clf()
    plt.plot(meanr, tau5psig, color='blue', marker='o', label=r'VAR($\tau_{5+}$)')
    plt.plot(meanr, sig_rho5p, color='red', marker='o', label=r'VAR($\rho_{5+}$)')
    plt.plot(meanr, sig_rho4p, color='green', marker='o', label=r'VAR($\rho_{4+}$)')
    plt.plot(meanr, sig_rho3p, color='gray', marker='o', label=r'VAR($\rho_{3+}$)')
    plt.ylim([1.e-13, 2.e-6])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)  
    plt.ylabel('Error', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.legend(loc='best', fontsize=18)
    plt.tight_layout()
    filename = os.path.join(plotspath, 'sigmas_eq3p.png')
    plt.savefig(filename)
    print(filename, 'printed')

    plt.clf()
    plt.plot(meanr, tau5msig, color='blue', marker='o', label=r'VAR($\tau_{5-}$)')
    plt.plot(meanr, sig_rho5m, color='red', marker='o', label=r'VAR($\rho_{5-}$)')
    plt.plot(meanr, sig_rho4m, color='green', marker='o', label=r'VAR($\rho_{4-}$)')
    plt.plot(meanr, sig_rho3m, color='gray', marker='o', label=r'VAR($\rho_{3-}$)')
    plt.ylim([1.e-13, 2.e-6])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)  
    plt.ylabel('Error', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.legend(loc='best', fontsize=18)
    plt.tight_layout()
    filename = os.path.join(plotspath, 'sigmas_eq3m.png')
    plt.savefig(filename)
    print(filename, 'printed')

    

if __name__ == "__main__":
    main()
