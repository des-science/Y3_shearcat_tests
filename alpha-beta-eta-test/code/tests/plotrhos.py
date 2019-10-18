import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    parser.add_argument('--rhos',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/rhos_large_scales.fits',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args


def main():
    import sys; sys.path.append(".")
    from src.plot_stats import plotallrhosfits
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

    xlim = [0.1, 250.]
    #xlim = [None, None]
    ylims = [[1.e-12,1.e-5 ],[1.e-13,1.e-6 ],[5.e-9 ,3.e-4 ],[1.e-11 ,4.e-3 ]]
    #ylims = [[None,None ],[None,None ],[None , None ],[None ,None]]
    rhostitle = ''
    plotallrhosfits(args.rhos, outpath=plotspath, title=rhostitle, xlim=xlim, ylims=ylims)
    

if __name__ == "__main__":
    main()
