today = '08-19-19_'
import numpy as np
import os
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce xip and xim correlations, i.e correlation among galaxies for all Flask realisations')
    
    parser.add_argument('--flask_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/desy3_6x2pt_lognormal-maskedcats/', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--mod', default=False,
                        action='store_const', const=True,
                        help='If true it substracts the mean to each field before calculate correlations')
    parser.add_argument('--seed', default=1,  type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--zbin', default=1 , type=int,
                        help='zbin used, useful to run parallel')
    parser.add_argument('--cookie', default=1 , type=int,
                        help='cookie used, useful to run parallel')
    parser.add_argument('--g1flip', default=False,
                        action='store_const', const=True,
                        help='If true invert the sig of g1')
    parser.add_argument('--g2flip', default=False,
                        action='store_const', const=True,
                        help='If true invert the sig of g2')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/catalogs/FLASK/xis_noflip/',
                        help='location of the output of the files')    
    args = parser.parse_args()

    return args

def main():
    import sys; sys.path.append(".")
    from astropy.io import fits
    import numpy as np
    from src.runcorr import measure_xi
    from src.read_cats import read_data, toList, read_flask
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
    
 
    bin_config = dict( sep_units = 'arcmin', min_sep = 1.0, max_sep = 250, nbins = 20,)
    #bin_config = dict(sep_units = 'arcmin' , bin_slop = 0.1, min_sep = 0.1, max_sep = 300, bin_size = 0.2)
    flip = [args.g1flip, args.g2flip]
    ck=args.cookie
    for seed in range (args.seed, 401):
        for zbin in range(args.zbin, 5):
            #SKIP ALREADY PROCESS OR missing flask cats
            outname = os.path.join(outpath, 'xi_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            inname = os.path.join(args.flask_cat, 'src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            if os.path.isfile(outname):
                print(outname, "Already exist. Skipping")
                continue
            if not os.path.isfile(inname):
                print(inname, "does not exist. Skipping")
                continue
            data_galaxies =  read_flask(args.flask_cat, seed, zbin, ck, flip=flip)
            print("Total objects in catalog:", len(data_galaxies))
            xi = measure_xi( data_galaxies, bin_config, mod=args.mod)[0]
            ximarr = xi.xim
            xiparr = xi.xip
            xis = [xiparr, ximarr]
            xis_names = ['XIP', 'XIM']

            
            ##Format of the fit file output
            names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
            forms = ['i4', 'i4', 'i4',  'f4',  'f4']
            dtype = dict(names = names, formats=forms)
            nrows = len(ximarr)
            outdata = np.recarray((nrows, ), dtype=dtype)

            covmat = np.diag(np.concatenate( (xi.varxip, xi.varxim ) ))
            hdu = fits.PrimaryHDU()
            hdul = fits.HDUList([hdu])
            covmathdu = fits.ImageHDU(covmat, name='COVMAT')
            hdul.insert(1, covmathdu)
                
            bin1array = np.array([ zbin]*nrows)
            bin2array = np.array([ -999]*nrows)
            angbinarray = np.arange(nrows)
            angarray = np.exp(xi.meanlogr)
            for j, nam in enumerate(xis_names):
                array_list = [bin1array, bin2array, angbinarray,np.array(xis[j]),  angarray ]
                for array, name in zip(array_list, names): outdata[name] = array 
                corrhdu = fits.BinTableHDU(outdata, name=nam)
                hdul.insert(j+2, corrhdu)

            hdul[1].header['COVDATA'] = True
            hdul[1].header['EXTNAME'] =  'COVMAT'
            hdul[1].header['NAME_0'] =  'XIP'
            hdul[1].header['STRT_0'] =  0
            hdul[1].header['LEN_0'] = nrows
            hdul[1].header['NAME_1'] =  'XIM'
            hdul[1].header['STRT_1'] =  nrows
            hdul[1].header['LEN_1'] = nrows
          
            
            hdul[2].header['QUANT1'] = 'GeR'; hdul[3].header['QUANT1'] = 'GeR'
            hdul[2].header['QUANT2'] = 'GeR'; hdul[3].header['QUANT2'] = 'GeR'

            
            outname = os.path.join(outpath, 'xis_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            hdul.writeto(outname, overwrite=True)
            
        args.zbin = 1
if __name__ == "__main__":
    main()
