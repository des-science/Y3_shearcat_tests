import os
today = '02-06-19_'
import numpy as np
import treecorr

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce Tau correlations, i.e correlation among galaxies and reserved stars')
    
    parser.add_argument('--metacal_cat',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3fullmaster/Y3_mastercat_v2_6_20_18.h5',
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
    parser.add_argument('--frac', default=1.,
                        type=float,
                        help='Choose a random fraction of the input stars')
    parser.add_argument('--mod', default=True,
                        action='store_const', const=True,
                        help='If true it substracts the mean to each field before calculate correlations')
    parser.add_argument('--bin_config', default=None,
                        help='bin_config file for running taus')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='npairstaus_metacal.fits', type=str,
                        help='filename of the tau output file')
    parser.add_argument('--nz_source',
                        #default='/home/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Indexes catalog to select galaxies in a particular redshift bin in Metacal')
    args = parser.parse_args()
    return args

def main():
    from src.read_cats import read_data_stars, toList, read_metacal
    from src.runcorr import measure_npairs
    from astropy.io import fits
    import treecorr
 
    
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
    
   
    if args.bin_config is not None:
        print("Using external bin config")
        bin_config = treecorr.read_config(args.bin_config)
        print(bin_config)
    else:
        bin_config = dict( sep_units = 'arcmin', min_sep = 1.0, max_sep = 250, nbins = 20,)
        #bin_config = dict(sep_units = 'arcmin' , bin_slop = 0.1, min_sep = 0.1, max_sep = 300, bin_size = 0.2)
    
    
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
     
    for zbin in range(1, 5):
        print('Starting Tomography!, measuring tau for zbin=', zbin)
        data_galaxies = read_metacal(args.metacal_cat, galkeys, zbin=zbin,nz_source_file=args.nz_source)
        results = measure_npairs( data_stars , data_galaxies, bin_config, mod=args.mod)[0]

        ##Format of the fit file output
        names=['ANG', 'NPAIRS']
        forms = ['f8',  'i4']
        dtype = dict(names = names, formats=forms)
        nrows = len(results.npairs)
        outdata = np.recarray((nrows, ), dtype=dtype)

        print(results.npairs)
      
        angarray = np.exp(results.meanlogr)
        array_list = [angarray, np.array(results.npairs)  ]
        for array, name in zip(array_list, names): outdata[name] = array
        corrhdu = fits.BinTableHDU(outdata, name='NPAIRS%d'%(zbin))
        hdul.insert(zbin, corrhdu)
        
    filename = os.path.join(outpath, args.filename)
    print("Printin file:", filename)
    hdul.writeto(filename, overwrite=True)
                
   
       
    
if __name__ == "__main__":
    main()
