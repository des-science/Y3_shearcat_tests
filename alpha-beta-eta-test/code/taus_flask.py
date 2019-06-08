today = '06-06-19_'
import numpy as np
import os
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--flask_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/flask/desy3_6x2pt_lognormal-maskedcats/', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.grizY',
                        #default='/home/dfa/sobreira/alsina/DESWL/psf/testexp',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--seed', type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--zbin', default=1 , type=int,
                        help='seed used, useful to run parallel')
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
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/flask/taus/',
                        help='location of the output of the files')
    parser.add_argument('--tomo', default=False,
                        action='store_const', const=True,
                        help='Run all tomographic correlations')
    
    args = parser.parse_args()

    return args

def main():
    from astropy.io import fits
    import numpy as np
    from src.runcorr import measure_tau
    from src.read_cats import read_data, toList, read_flask
    import h5py as h
    
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
 
    exps = toList(args.exps_file)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data_stars))
    data_stars = data_stars[data_stars['mag']<20]
    print("Objects with magnitude <20",  len(data_stars))
 
    bin_config = dict( sep_units = 'arcmin', min_sep = 1.0, max_sep = 250, nbins = 20,)
    #bin_config = dict(sep_units = 'arcmin' , bin_slop = 0.1, min_sep = 0.1, max_sep = 300, bin_size = 0.2)
    for seed in range (args.seed, 401):
        for zbin in range(args.zbin, 5):
            for ck in range(1, 3):
                data_galaxies =  read_flask(args.flask_cat, seed, zbin, ck)
                print("Total objects in catalog:", len(data_galaxies))
                tau0, tau2, tau5= measure_tau( data_stars , data_galaxies,
                                       bin_config, mod=args.mod)
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
                
                bin1array = np.array([ zbin]*nrows)
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
                hdul[1].header['NAME_0'] =  'TAU0'
                hdul[1].header['STRT_0'] =  0
                hdul[1].header['LEN_0'] = nrows
                hdul[1].header['NAME_1'] =  'TAU2'
                hdul[1].header['STRT_1'] =  nrows
                hdul[1].header['LEN_1'] = nrows
                hdul[1].header['NAME_2'] =  'TAU5'
                hdul[1].header['STRT_2'] =  2*nrows
                hdul[1].header['LEN_2'] = nrows
            
                hdul[2].header['QUANT1'] = 'GeR'; hdul[3].header['QUANT1'] = 'GeR'
                hdul[2].header['QUANT2'] = 'PeR'; hdul[3].header['QUANT2'] = 'PeR'
                hdul[4].header['QUANT1'] = 'GeR'; hdul[5].header['QUANT1'] = 'GeR'
                hdul[4].header['QUANT2'] = 'PqR'; hdul[5].header['QUANT2'] = 'PqR'
                hdul[6].header['QUANT1'] = 'GeR'; hdul[7].header['QUANT1'] = 'GeR'
                hdul[6].header['QUANT2'] = 'PwR'; hdul[7].header['QUANT2'] = 'PwR'
            
                outname = os.path.join(outpath, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
                hdul.writeto(outname, overwrite=True)

            args.zbin = 1
if __name__ == "__main__":
    main()
