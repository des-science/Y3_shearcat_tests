import os
import numpy as np

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='samples total PSF bias using rhos cosmology and previous samples of abe')
    parser.add_argument('--samplesabe', default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_ab_margin_Y3_03-31-20_JKmarco/sample_ab_Y3_03-31-20_JKmarco.fits')
    parser.add_argument('--rhoscosmo', default='/data/catalogs/rhos-taus/RHOS_Y3_7-24-19-mod_cosmo.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--burn_frac', default=0.0, type=float, 
                        help='Burn frac of samples')
    parser.add_argument('--outpath',
                        default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_abe_margin_Y3_03-31-20_JKmarco',
                        help='location of the output of the final contaminant')
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
    from src.readfits import  read_rhos
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

    abesamps = fitsio.read(args.samplesabe)
    #Burning just to be sure
    frac= args.burn_frac
    abesamps =abesamps[int(frac*len(abesamps) ): ]

    ##SAMPLING BIAS
    meanr, rhos,  covmat = read_rhos(args.rhoscosmo)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m = rhos

    
    if ('a1' in abesamps.dtype.names): alist = [abesamps['a1'], abesamps['a2'], abesamps['a3'], abesamps['a4']]
    else: alist = [[0]*len(abesamps) ]*4; print('ALPHA is not in samples')
    if ('b1' in abesamps.dtype.names): blist = [abesamps['b1'], abesamps['b2'], abesamps['b3'], abesamps['b4']]
    else: blist = [[0]*len(abesamps) ]*4; print('BETA is not in samples')
    if ('e1' in abesamps.dtype.names): elist = [abesamps['e1'], abesamps['e2'], abesamps['e3'], abesamps['e4']]
    else: elist = [[0]*len(abesamps) ]*4; print('ETA is not in samples')

    nbins=4
    a=[i for i in range(nbins)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2): bin_pairs.append(p)
    veclist = []
    
    for z in range(len(abesamps)):
        dxip = [alist[i][z]*alist[j][z]*rho0p + blist[i][z]*blist[j][z]*rho1p + elist[i][z]*elist[j][z]*rho3p + (blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2p + (blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4p + (elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5p for i,j in bin_pairs]
        dxim = [alist[i][z]*alist[j][z]*rho0m + blist[i][z]*blist[j][z]*rho1m + elist[i][z]*elist[j][z]*rho3m + (blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2m + (blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4m + (elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5m for i,j in bin_pairs]
        dxi = dxip + dxim
        veclist.append(np.concatenate(np.c_[dxi]))
    covmat = np.cov(np.c_[veclist].T)





    
    #Writting final contamination
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])

    dxip_list = []; dxim_list = []; bin1_list = []; bin2_list = []; angbin_list = []; ang_list = []
    covxip_list = []; covxim_list = []
    
    nbins=4
    a=[i for i in range(nbins)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2): bin_pairs.append(p)
    al = bestparameters(alist)
    bl = bestparameters(blist)
    el = bestparameters(elist)
    for i,j in bin_pairs:
         dxip = al[i]*al[j]*rho0p + bl[i]*bl[j]*rho1p + el[i]*el[j]*rho3p + (bl[i]*al[j] + bl[j]*al[i])*rho2p + (bl[i]*el[j] + bl[j]*el[i])*rho4p + (el[i]*al[j] +el[j]*al[i])*rho5p
         dxim = al[i]*al[j]*rho0m + bl[i]*bl[j]*rho1m + el[i]*el[j]*rho3m + (bl[i]*al[j] + bl[j]*al[i])*rho2m + (bl[i]*el[j] + bl[j]*el[i])*rho4m + (el[i]*al[j] +el[j]*al[i])*rho5m

        
         ang_list.append(meanr)
         bin1_list.append(np.array( [i + 1]*len(meanr)))
         bin2_list.append(np.array( [j + 1]*len(meanr)))
         angbin_list.append(np.arange(len(meanr)))
         dxip_list.append(dxip)
         dxim_list.append(dxim)

    hdul.insert(1, fits.ImageHDU(covmat, name='COVMAT'))
    bin1array = np.concatenate(bin1_list)
    bin2array = np.concatenate(bin2_list)
    angbinarray = np.concatenate(angbin_list)
    valuearray = np.concatenate(dxip_list)
    angarray = np.concatenate(ang_list)
        
    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f4',  'f4']
    dtype = dict(names = names, formats=forms)
    nrows = len(angarray)
    outdata = np.recarray((nrows, ), dtype=dtype)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xip')
    hdul.insert(2, corrhdu)
    valuearray = np.concatenate(dxim_list)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xim')
    hdul.insert(3, corrhdu)
    outname = os.path.join(outpath, args.filename)
    hdul.writeto(outname, overwrite=True)
    print(outname,'Written!')

   
  
    
if __name__ == "__main__":
    main()



        
