import os
today = '27-04-19_'

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
    parser.add_argument('--tomo', default=False,
                        action='store_const', const=True,
                        help='Run all tomographic correlations')
    parser.add_argument('--nz_source',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')
    args = parser.parse_args()
    return args
    
def measure_tau(data_stars, data_galaxies, min_sep = 0.1,  max_sep=300, bin_size=0.2, sep_units='arcmin', prefix='piff', mod=True):
    """Compute the tau statistics
    """
    import treecorr
    import numpy as np
    import gc
    
    e1 = data_stars['obs_e1']
    e2 = data_stars['obs_e2']
    p_e1 = data_stars[prefix+'_e1']
    p_e2 = data_stars[prefix+'_e2']
    T = data_stars['obs_T']
    p_T = data_stars[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T

    #w1 = p_e1*dt
    #w2 = p_e2*dt
    w1 = e1*dt
    w2 = e2*dt

    e1gal = data_galaxies['e_1']
    e2gal = data_galaxies['e_2']
    
    #Modified ellipticities reserved stars and galaxies
    if(mod):
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
        e1gal = e1gal - np.array(np.mean(e1gal))
        e2gal = e2gal - np.array(np.mean(e2gal))

        
    ra = data_stars['ra']
    dec = data_stars['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    ragal = data_galaxies['ra']
    decgal = data_galaxies['dec']
    print('ragal = ',ragal)
    print('decgal = ',decgal)
    
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
    wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    egal_cat.name = 'egal_cat'


    #del data_stars, data_galaxies,  ra, dec, ragal, decgal, p_e1, p_e2, de1, de2, w1, w2, e1gal, e2gal, e1, e2, T, p_T, dt
    del data_stars, data_galaxies, e1, e2, T, p_T, dt
    gc.collect()
    
    #bin_config = dict( sep_units = sep_units, min_sep = 2.5, max_sep = 250, nbins = 20,)
    #bin_config = dict(sep_units = sep_units, bin_slop = 0.1, min_sep = 0.5,  max_sep=300, bin_size=0.2)

    bin_config = dict(sep_units = sep_units , bin_slop = 0.1, min_sep = min_sep, max_sep = max_sep, bin_size = bin_size)
    
    results = []

    for (cat1, cat2) in [(egal_cat, ecat), 
                         (egal_cat, decat),
                          (egal_cat, wcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=3)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)
    print('All correlations done sucessfully')
    return results

def main():
    import numpy as np
    from src.read_cats import read_data, toList, read_metacal
    from astropy.io import fits
    import gc

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

    del bands, tilings, exps, keys
    gc.collect()

    min_sep = 0.1;  max_sep=300; bin_size=0.2
    galkeys = ['ra','dec','e_1','e_2','R11','R22']

    if(args.tomo):
        print('Starting Tomography!')
        nbins = 4
    
        for bin_c in range(nbins):
            print('Starting bin!',  bin_c)
            
            data_gal = read_metacal(args.metacal_cat,  galkeys,  zbin=bin_c,  nz_source_file=args.nz_source)
            tau0, tau2, tau5= measure_tau(data_stars, data_gal, min_sep = min_sep, max_sep = max_sep, bin_size = bin_size,  mod=args.mod)
            tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
            tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
            vartau0arr = tau0.varxi; vartau2arr = tau2.varxi; vartau5arr = tau5.varxi;
            taus = [tau0parr, tau0marr, tau2parr, tau2marr, tau5parr, tau5marr]
            taus_names = ['TAU0P', 'TAU0M','TAU2P','TAU2M', 'TAU5P', 'TAU5M']
            vares = [vartau0arr, vartau2arr, vartau5arr]

            ##Format of the fit file output
            names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
            forms = ['i4', 'i4', 'i4',  'f8',  'f8']
            dtype = dict(names = names, formats=forms)
            nrows = len(tau0marr)
            outdata = np.recarray((nrows, ), dtype=dtype)
        
            covmat = np.diag(np.concatenate((vares[0], vares[0], vares[1], vares[1],  vares[2], vares[2])))
            hdu = fits.PrimaryHDU()
            hdul = fits.HDUList([hdu])
            covmathdu = fits.ImageHDU(covmat, name='COVMAT')
            hdul.insert(1, covmathdu)
   
            bin1array = np.array([bin_c]*nrows)
            bin2array = np.array([-999]*nrows)
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

            
            print("Printin file:", outpath + 'TAUS' +'_bin_' + str(bin_c) +  '.fits' )
            hdul.writeto(outpath + 'TAUS' +'_bin_' + str(bin_c) +  '.fits', overwrite=True)
                
    else:
        data_galaxies =  read_metacal(args.metacal_cat,  galkeys )
        print("Total objects in catalog:", len(data_galaxies))
        tau0, tau2, tau5= measure_tau(data_stars, data_galaxies, min_sep = min_sep, max_sep = max_sep, bin_size = bin_size, mod=args.mod)
        tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
        tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
        vartau0arr = tau0.varxi; vartau2arr= tau2.varxi; vartau5arr = tau5.varxi;

        taus = [tau0parr, tau0marr, tau2parr, tau2marr, tau5parr, tau5marr]
        taus_names = ['TAU0P', 'TAU0M','TAU2P','TAU2M', 'TAU5P', 'TAU5M']
        vares = [vartau0arr, vartau2arr, vartau5arr]

        ##Format of the fit file output
        names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
        forms = ['i4', 'i4', 'i4',  'f8',  'f8']
        dtype = dict(names = names, formats=forms)
        nrows = len(tau0marr)
        outdata = np.recarray((nrows, ), dtype=dtype)

     
        covmat = np.diag(np.concatenate((vares[0], vares[0], vares[1], vares[1],  vares[2], vares[2])))
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
