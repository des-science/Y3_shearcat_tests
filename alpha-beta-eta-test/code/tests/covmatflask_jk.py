import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--tausflask',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/taus_flask_01-1/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--tausjk',
                        default='/home/dfa/sobreira/alsina/catalogs/JK/taus_jk1000/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--zbin', default=1 , type=int,
                        help='seed used, useful to run parallel')
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

def buildcovarianceflask(tauspath, zbin,  ndraws):
    import numpy as np
    from src.readfits import read_taus
    veclist = []
    count = 0 ;seed = 1;ck = 1
    while(count <ndraws):
        name = os.path.join(tauspath, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
        exist =  os.path.isfile(name)
        if exist:
            meanr, taus, covtaus = read_taus(name)
            if (np.count_nonzero(taus) == 0):
                print("Warning, weird measurement, skipping", name)
            else:
                veclist.append(np.concatenate(np.c_[taus]))
                count +=1
        else:
            print(name, 'Does not exist')          
        seed += 1
    print(count, "FLASK catalogs were read")
                
    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
  
    lengths = [len(taus[0]), len(taus[1]), len(taus[2]), len(taus[3]),
               len(taus[4]), len(taus[5])]
    return lengths, covmat

def buildcovarianceJK(jkpath, zbin,  ndraws):
    import numpy as np
    from src.readfits import read_taus
    veclist = []
    count = 0 ; jkid = 1;ck = 1
    while(count <ndraws):
        name = os.path.join(jkpath, 'taus_src-cat_jk%d_z%d.fits'%(jkid,zbin ))
        exist =  os.path.isfile(name)
        if exist:
            meanr, taus, covtaus = read_taus(name)
            if (np.count_nonzero(taus) == 0):
                print("Warning, weird measurement, skipping", name)
            else:
                veclist.append(np.concatenate(np.c_[taus]))
                count +=1
        else:
            print(name, 'Does not exist')          
        jkid += 1
    print(count, "JK catalogs were read")
                
    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
  
    lengths = [len(taus[0]), len(taus[1]), len(taus[2]), len(taus[3]),
               len(taus[4]), len(taus[5])]
    return lengths, covmat
 
def plotcov(covmat, lengths,zbin,  filename):
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
    
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

def main():
    import sys; sys.path.append(".")
    import fitsio
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    import itertools
    import numpy as np
    
    args = parse_args()

    #Make directory where the ouput data will be
    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    zbin = args.zbin
    lengths, covmat =  buildcovarianceflask(args.tausflask,zbin,  45)
    filename = os.path.join(plotspath, 'CovariancematrixTausFlask_zbin%d.png'%(args.zbin))
    plotcov(covmat, lengths, zbin,  filename)

    lengths, covmatjk=  buildcovarianceJK(args.tausjk,zbin,  45)
    filename = os.path.join(plotspath, 'CovariancematrixTausJK_zbin%d.png'%(args.zbin))
    plt.clf()
    plotcov(covmatjk, lengths, zbin,  filename)

    filename = os.path.join(plotspath, 'CovariancematrixTausJK-FLASK_zbin%d.png'%(args.zbin))
    plt.clf()
    plotcov(covmatjk - covmat, lengths, zbin,  filename)

    print("COVMAT FLASK")
    print(np.diag(covmat))
    print("COVMAT JK")
    print(np.diag(covmatjk))
    print("COVMAT DIFF")
    print(np.diag(covmat - covmatjk))

    
if __name__ == "__main__":
    main()

