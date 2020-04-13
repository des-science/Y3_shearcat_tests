import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import numpy as np

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='samples total PSF bias using rhos cosmology and previous samples of abe')
    parser.add_argument('--samplesabe', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/samples_abe.fits')
    parser.add_argument('--taus',
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    parser.add_argument('--rhos',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')

    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--minscale', default=None, type=float, 
                        help='Limit the analysis to certain minimum scale')
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
                        
   
                        
def main():
    import fitsio
    import getdist
    from getdist import plots, MCSamples
    from astropy.io import fits
    import itertools
    from src.maxlikelihood import percentiles, bestparameters
    from src.readfits import  read_rhos,  read_taus
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
    samplesbin1 = [abesamps['a1'], abesamps['b1'], abesamps['e1']]
    samplesbin2 = [abesamps['a2'], abesamps['b2'], abesamps['e2']]
    samplesbin3 = [abesamps['a3'], abesamps['b3'], abesamps['e3']]
    samplesbin4 = [abesamps['a4'], abesamps['b4'], abesamps['e4']]
    samples = [samplesbin1, samplesbin2, samplesbin3, samplesbin4]

    nbins = 100
    data = {}
    data['rhos'] = read_rhos(args.rhos, minscale=args.minscale, maxscale=args.maxscale)[1]
    data['cov_rhos'] = read_rhos(args.rhos, minscale=args.minscale, maxscale=args.maxscale)[2]
    for i,  taufile in enumerate(args.taus):
        data['taus'] = read_taus(taufile, minscale=args.minscale, maxscale=args.maxscale)[1]
        data['cov_taus'] = read_taus(taufile, minscale=args.minscale, maxscale=args.maxscale)[2]
        chis = [ chi2nu( np.c_[samples[i]].T[j] ,data, eq=4, mflags=[True, True, True], xip=True, xim=True) for j in range(len(abesamps['a1'])) ]
        plt.close()
        plt.clf()
       
        plt.hist(chis, bins=np.linspace(0.7, 1.5, nbins), label='zbin%d'%(i + 1) )
        plt.ylabel('Counts'); plt.xlabel(r'$\chi^{2}_{\nu}$')
        #plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        filename =  os.path.join(plotspath,'chis_sample_zbin%d.png'%(i + 1))
        print("Printing File", filename)
        plt.savefig(filename)
  
    
if __name__ == "__main__":
    main()



        
