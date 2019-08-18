import os
today = '02-06-19_'
import numpy as np
import treecorr

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce Tau correlations, i.e correlation among galaxies and reserved stars')
    
    parser.add_argument('--metacal_cat',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3fullmaster/Y3_mastercat_v2_6_20_18.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.grizY',
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
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='TAUZ_zbin_n.fits', type=str,
                        help='filename of the tau output file')
    parser.add_argument('--zbin', default=None,type=int, 
                        help='Run particular tomobin')
    parser.add_argument('--nz_source',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        help='Indexes catalog to select galaxies in a particular redshift bin in Metacal')
    args = parser.parse_args()
    return args

def main():
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
  
    #Reading Mike stars catalog
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T', 'mag']
    galkeys = ['ra','dec','e_1','e_2','R11','R22']

    data_stars = read_data_stars(toList(args.exps_file),args.piff_cat, keys,limit_bands=args.bands,use_reserved=args.use_reserved)
    
   

    bin_config = dict( sep_units = 'arcmin', min_sep = 1.0, max_sep = 250, nbins = 20,)
    #bin_config = dict(sep_units = 'arcmin' , bin_slop = 0.1, min_sep = 0.1, max_sep = 300, bin_size = 0.2)
    
    

    if args.zbin is not None:
        print('Starting Tomography!, measuring tau for zbin=', args.zbin)
        data_galaxies = read_metacal(args.metacal_cat, galkeys, zbin=args.zbin,nz_source_file=args.nz_source)
        tau0, tau2, tau5= measure_tau( data_stars , data_galaxies,
                                       bin_config, mod=args.mod)
        tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
        tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
        taus = [tau0parr, tau0marr, tau2parr, tau2marr, tau5parr, tau5marr]
        taus_names = ['TAU0P', 'TAU0M','TAU2P','TAU2M', 'TAU5P', 'TAU5M']

        ##Format of the fit file output
        names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
        forms = ['i4', 'i4', 'i4',  'f8',  'f8']
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

        print("Printin file:", outpath + args.filename)
        hdul.writeto(outpath + args.filename, overwrite=True)
                
    else:
        data_galaxies = read_metacal(args.metacal_cat,  galkeys )
        tau0, tau2, tau5= measure_tau(data_stars, data_galaxies,
                                      bin_config, mod=args.mod)
        tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
        tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
        taus = [tau0parr, tau0marr, tau2parr, tau2marr, tau5parr, tau5marr]
        taus_names = ['TAU0P', 'TAU0M','TAU2P','TAU2M', 'TAU5P', 'TAU5M']

        ##Format of the fit file output
        names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
        forms = ['i4', 'i4', 'i4',  'f8',  'f8']
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

        print("Printin file:", outpath + 'TAUS.fits' )
        hdul.writeto(outpath + 'TAUS.fits', overwrite=True)
    
if __name__ == "__main__":
    main()
