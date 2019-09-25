today = '08-19-19_'
import numpy as np
import os
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce Tau correlations, i.e correlation among galaxies and reserved stars for all Flask realisations')
    
    parser.add_argument('--flask_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/flask/desy3_6x2pt_lognormal-maskedcats_v3/', 
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
    parser.add_argument('--seed', default=1,  type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--cookie', default=1 , type=int,
                        help='cookie used, useful to run parallel')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/catalogs/flask/taus_v3/',
                        help='location of the output of the files')    
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
    ck=args.cookie
    for seed in range (args.seed, 701):
        for zbin in range(1, 5):
            #SKIP ALREADY PROCESS OR missing flask cats
            outname = os.path.join(outpath, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            inname = os.path.join(args.flask_cat, 'src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            if os.path.isfile(outname):
                print(outname, "Already exist. Skipping")
                continue
            if not os.path.isfile(inname):
                print(inname, "does not exist. Skipping")
                continue
            data_galaxies =  read_flask(args.flask_cat, seed, zbin, ck)
            results = measure_npairs( data_stars , data_galaxies, bin_config, mod=args.mod)[0]

            ##Format of the fit file output
            names=['ANG', 'NPAIRS']
            forms = ['f8',  'i4']
            dtype = dict(names = names, formats=forms)
            nrows = len(results.npairs)
            outdata = np.recarray((nrows, ), dtype=dtype)
      
            angarray = np.exp(tau0.meanlogr)
            array_list = [angarray, np.array(results.npairs)  ]
            for array, name in zip(array_list, names): outdata[name] = array
            corrhdu = fits.BinTableHDU(outdata, name='NPAIRS%d'%(zbin))
            hdul.insert(zbin, corrhdu)
        
        filename = os.path.join(outpath, 'npairstausflask_seed%d.fits'%(seed) )
        print("Printing file:", filename)
        hdul.writeto(filename, overwrite=True)
            

if __name__ == "__main__":
    main()
