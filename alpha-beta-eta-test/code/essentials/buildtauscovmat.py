import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--tausflask',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/taus_flask_01-250/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--tausjk',
                        default='/home/dfa/sobreira/alsina/catalogs/JK/taus_jk_01-250_V2/',
                        help='Full Path to the taus measurement of JK catalogs')
    parser.add_argument('--input_tau',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS.fits',
                        help='Fit file with the taus correlations, and which covariance matrix will be replaced and writen in filename')
    parser.add_argument('--filename',
                        default='TAUS_FLASK.fits',
                        help='Fit file based on inputfile but now with Flask Covariance matrix')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the files')
    parser.add_argument('--zbin', default=None,
                        help='int indicating the number of the redshift bin')
    parser.add_argument('--plots', default=True,
                        action='store_const', const=True, help='Plot correlations functions')
    parser.add_argument('--jk', default=False,
                        action='store_const', const=True, help='flag to switch betwen flask and JK')
    parser.add_argument('--njk', default=1000 , type=int,
                         help='number of jk patches, only used if args.jk is True')
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
 
def main():
    import fitsio
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    import itertools
    import numpy as np
    from src.readfits import read_taus

 
    
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

    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f4',  'f4']
    dtype = dict(names = names, formats=forms)
    #nrows = len(tau0marr)
    #outdata = np.recarray((nrows, ), dtype=dtype)
    veclist = []
    count = 0; 
    for seed in range(1, 1000 ): #version2
        for ck in range(1,2):
            if args.zbin is not None: 
                zmin = args.zbin; zmax=args.zbin+1
            else:
                zmin=1; zmax=5
            for zbin in range(zmin,zmax):
                if args.jk: name = os.path.join(args.tausjk, 'taus_src-cat_jk%d_z%d.fits'%(seed,zbin ))
                else: name = os.path.join(args.tausflask, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
                exist =  os.path.isfile(name)
                if exist:
                    try: meanr, taus, covtaus = read_taus(name)
                    except: print(name, 'could not be open')
                    if (np.count_nonzero(taus) == 0):
                        print("Warning, weird measurement, skipping", name)
                    else:
                        veclist.append(np.concatenate(np.c_[taus]))
                        count +=1
                else:
                    print(name, 'Does not exist')
    print(count, "catalogs were read")
                
    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
    if args.jk: covmat *= (args.njk - 1)   
    
    if(args.plots):
        lengths = [len(taus[0]), len(taus[1]), len(taus[2]), len(taus[3]),
               len(taus[4]), len(taus[5])]
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
        if args.zbin is not None:
            if args.jk: filename = os.path.join(plotspath,'CovariancematrixTausJK_zbin%d.png'%(args.zbin))
            else: filename = os.path.join(plotspath,'CovariancematrixTausFlask_zbin%d.png'%(args.zbin))
        else:
            if args.jk: filename = os.path.join(plotspath,'CovariancematrixTausJK.png')
            else: filename = os.path.join(plotspath,'CovariancematrixTausFlask.png')
        plt.savefig(filename, dpi=500)
        print(filename, 'Printed!')

    hdulist = fits.open(args.input_tau)
    oldheaders =  [hdulist[1].header]
    hdulist.pop(index=1);
    
    
    covmathdu= fits.ImageHDU(covmat)
    hdulist.insert(1, covmathdu)
    hdulist[1].header = oldheaders[0]
    #print(hdulist)
    print("writting, ", outpath + args.filename)
    hdulist.writeto(outpath + args.filename, overwrite=True)
    
if __name__ == "__main__":
    main()

