import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='shape noise flask vs metacal')
    
    parser.add_argument('--npairstaus_metacal',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/npairstaus_metacal.fits',
                        help='Full Path to the Metacalibration npairs file')
    parser.add_argument('--npairstaus_flask',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/npairstausflask_seed1.fits',
                        help='Full Path to the flask npairs file')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()
    return args


def main():
    import sys; sys.path.append(".")
    import numpy as np
    import fitsio

    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise


    for zbin in [1, 2, 3, 4]:
        metacal_pairs_data =  fitsio.read(args.npairstaus_metacal, ext=zbin)
        meanr_metacal =  metacal_pairs_data['ANG']
        npairs_metacal = metacal_pairs_data['NPAIRS']
        flask_pairs_data =  fitsio.read(args.npairstaus_flask, ext=zbin)
        meanr_flask =  flask_pairs_data['ANG']
        npairs_flask = flask_pairs_data['NPAIRS']
        plt.clf()
        plt.plot( meanr_metacal, npairs_metacal, color='red', label='Metacal zbin%d'%(zbin), marker='o')
        plt.plot( meanr_flask, npairs_flask, color='blue', label='Flask zbin%d'%(zbin), marker='o')
        plt.legend(loc='best', fontsize=24)
        plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
        plt.ylabel('N pairs', fontsize=24)
        plt.xscale('log')
        plt.yscale('log', nonposy='clip')
        plt.tight_layout()
        filename = os.path.join(plotspath, 'npairs_zbin%d.png'%(zbin))
        print('Printing :',  filename)
        plt.savefig(filename)
        plt.clf()
        plt.plot( meanr_metacal, npairs_flask/npairs_metacal, color='red', label='zbin%d'%(zbin), marker='o')
        plt.ylim([0.5, 2])
        plt.legend(loc='best', fontsize=24)
        plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
        plt.ylabel('Npairs ratio', fontsize=24)
        plt.xscale('log')
        #plt.yscale('log', nonposy='clip')
        plt.tight_layout()
        filename = os.path.join(plotspath, 'ratio_npairs_zbin%d.png'%(zbin))
        print('Printing :',  filename)
        plt.savefig(filename)
      
        

        
        
if __name__ == "__main__":
    main()
