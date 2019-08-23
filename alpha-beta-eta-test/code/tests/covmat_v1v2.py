def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Comparing covariance matrix of two flask catalogs')
    parser.add_argument('--file1',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                        help='Full Path to the taus measurement with flask version 1')
    parser.add_argument('--file2',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1_v2.fits',
                        help='Full Path to the taus measurement with flask version 2')
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
def main():
    from src.readfits import read_taus
    import numpy as np
    args = parse_args()
    meanr1, taus1,  covtau1 =  read_taus(args.file1)
    meanr2, taus2,  covtau2 =  read_taus(args.file2)
    covtau1_inv = np.linalg.pinv(covtau1,  rcond=1e-10)
    covtau2_inv = np.linalg.pinv(covtau2,  rcond=1e-10)
    
    print(covtau1)
    print("SPACE SPACE SPACE")
    print(covtau2)
    print("SPACE SPACE SPACE")
    
    print(covtau1_inv)
    print("SPACE SPACE SPACE")
    print(covtau2_inv)
    print("SPACE SPACE SPACE")

    print(np.linalg.det(covtau1_inv))
    print("SPACE SPACE SPACE")
    print(np.linalg.det(covtau2_inv))
    
if __name__ == "__main__":
    main()
