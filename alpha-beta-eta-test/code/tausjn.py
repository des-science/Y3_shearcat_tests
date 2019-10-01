import numpy as np
import os
import kmeans_radec
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--metacal_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/Y3_mastercat_7_24_19.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--piff_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/code/ally3.grizY',
                        #default='/home/dfa/sobreira/alsina/DESWL/psf/testexp',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--use_reserved', default=True,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')
    parser.add_argument('--mod', default=True,
                        action='store_const', const=True,
                        help='If true it substracts the mean to each field before calculate correlations')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/catalogs/JK/tausjn_njk4',
                        help='location of the output of the files')
    parser.add_argument('--bin_config', default=None,
                        help='bin_config file for running taus')
    parser.add_argument('--njks', default=4,type=int, 
                        help='Number of patches to divide footprint')
    parser.add_argument('--zbin', default=None,type=int, 
                        help='Run particular tomobin')
    parser.add_argument('--nz_source',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')    
    
    args = parser.parse_args()

    return args

def jk_kmeans(ra_sam, dec_sam, ra,dec,njk,plot=False):
    '''
    Function that takes RA and Dec from a given catalog and computes JK regions using kmeans_radec module by Erin Sheldon.
    
    Parameters
    ----------
    ra, dec : numpy arrays of RA and Dec. len(ra) = len(dec) = number of galaxies.
    njk : number of JK regions.
    
    Returns
    -------
    jk = JK region for each galaxy: integer ranging from 0 to njk-1. It is numpy array with the same length as ra and dec. 

    '''
    print("Running kmeans")
    from astropy.coordinates import SkyCoord, Angle
    from astropy import units
    radec = np.zeros((len(ra),2)); radec_sam = np.zeros((len(ra_sam),2))
    radec[:,0] = ra; radec_sam[:,0] = ra_sam
    radec[:,1] = dec; radec_sam[:,1] = dec_sam
    km = kmeans_radec.kmeans_sample(radec_sam,njk,maxiter=500,tol=1e-05)
    jk = km.find_nearest(radec)
    if not km.converged:
        print('k means did not converge')
    if plot:
        coords = SkyCoord(ra=ra, dec=dec, unit='degree')
        ra = coords.ra.wrap_at(180 * units.deg)
        dec = coords.dec
        plt.figure()
        plt.scatter(ra,dec,c=jk,lw=0,cmap='Paired',rasterized=True)
        plt.xlabel(r'RA',fontsize=12)
        plt.ylabel(r'Dec',fontsize=12)
        plt.tight_layout()
        plt.savefig('jk_kmeans.png')
    print("kmeans finished!")
    return jk


def main():
    import treecorr
    from src.read_cats import read_data_stars, toList, read_metacal
    from src.runcorr import measure_tau
    from astropy.io import fits
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        
    names=['JKR','ANGBIN','THETA','TAU0P','TAU0M','VAR_TAU0','TAU2P','TAU2M','VAR_TAU2','TAU5P','TAU5M', 'VAR_TAU5']
    forms = ['i4', 'i4', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)
    
  
    #Reading Mike stars catalog
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T', 'mag']
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
 
    exps = toList(args.exps_file)
    data_sam = read_data_stars(exps, args.piff_cat , keys,
                                          limit_bands=args.bands,
                                          use_reserved=args.use_reserved, frac=0.01)
    data_stars = read_data_stars(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    if args.bin_config is not None:
        print("Using external bin config")
        bin_config = treecorr.read_config(args.bin_config)
        print(bin_config)
    else:
        bin_config = dict( sep_units = 'arcmin', min_sep = 1.0, max_sep = 250, nbins = 20,)
    
    if args.zbin is not None:
        print('Starting measurement for zbin', args.zbin)
        
        data_gal = read_metacal(args.metacal_cat,  galkeys,  zbin=args.zbin,  nz_source_file=args.nz_source)
        njk = args.njks
        ##TODO generate km first an later finnearest,
        jkindexes_gals = jk_kmeans(data_sam['ra'], data_sam['dec'], data_gal['ra'], data_gal['dec'],njk)
        
        for jkidx in range(njk):
            print("running jackkniffe region",  jkidx)
            booljk = [jkindexes_gals!=jkidx ] 
            tau0, tau2, tau5= measure_tau(data_stars, data_gal[booljk], bin_config, mod=args.mod)

            tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
            tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
            taus = [tau0parr, tau0marr, tau2parr, tau2marr, tau5parr, tau5marr]
            taus_names = ['TAU0P', 'TAU0M','TAU2P','TAU2M', 'TAU5P', 'TAU5M']

            ##Format of the fit file output
            names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
            forms = ['i4', 'i4', 'i4',  'f4',  'f4']
            dtype = dict(names = names, formats=forms)
            nrows = len(tau0marr)
            outdata = np.recarray((nrows, ), dtype=dtype)

     
            covmat = np.diag(np.concatenate( (tau0.varxip, tau0.varxim, tau2.varxip, tau2.varxim, tau5.varxip, tau5.varxim ) ))
            hdu = fits.PrimaryHDU()
            hdul = fits.HDUList([hdu])
            covmathdu = fits.ImageHDU(covmat, name='COVMAT')
            hdul.insert(1, covmathdu)
        
            bin1array = np.array([ -999]*nrows)
            bin2array = np.array([ -999]*nrows)
            angbinarray = np.arange(nrows)
            angarray = np.exp(tau0.meanlogr)
        
            for j, nam in enumerate(taus_names):
                array_list = [bin1array, bin2array, angbinarray,np.array(taus[j]),  angarray ]
                for array, name in zip(array_list, names): outdata[name] = array 
                corrhdu = fits.BinTableHDU(outdata, name=nam)
                hdul.insert(j+2, corrhdu)

            hdul[1].header['COVDATA'] = True
            hdul[1].header['EXTNAME'] =  'COVMAT'
            hdul[1].header['NAME_0'] =  'TAU0P'
            hdul[1].header['STRT_0'] =  0
            hdul[1].header['LEN_0'] = nrows
            hdul[1].header['NAME_1'] =  'TAU0M'
            hdul[1].header['STRT_1'] =  nrows
            hdul[1].header['LEN_1'] = nrows
            hdul[1].header['NAME_2'] =  'TAU2P'
            hdul[1].header['STRT_2'] =  2*nrows
            hdul[1].header['LEN_2'] = nrows
            hdul[1].header['NAME_3'] =  'TAU2M'
            hdul[1].header['STRT_3'] =  3*nrows
            hdul[1].header['LEN_3'] = nrows
            hdul[1].header['NAME_4'] =  'TAU5P'
            hdul[1].header['STRT_4'] =  4*nrows
            hdul[1].header['LEN_4'] = nrows
            hdul[1].header['NAME_5'] =  'TAU5M'
            hdul[1].header['STRT_5'] =  5*nrows
            hdul[1].header['LEN_5'] = nrows
            
            hdul[2].header['QUANT1'] = 'GeR'; hdul[3].header['QUANT1'] = 'GeR'
            hdul[2].header['QUANT2'] = 'PeR'; hdul[3].header['QUANT2'] = 'PeR'
            hdul[4].header['QUANT1'] = 'GeR'; hdul[5].header['QUANT1'] = 'GeR'
            hdul[4].header['QUANT2'] = 'PqR'; hdul[5].header['QUANT2'] = 'PqR'
            hdul[6].header['QUANT1'] = 'GeR'; hdul[7].header['QUANT1'] = 'GeR'
            hdul[6].header['QUANT2'] = 'PwR'; hdul[7].header['QUANT2'] = 'PwR'

            filename = os.path.join(outpath, 'taus_jk%d_z%d.fits'%(jkidx, args.zbin))
            
            print("Printing file:", filename)
            hdul.writeto(filename, overwrite=True)


            
 
    
if __name__ == "__main__":
    main()
