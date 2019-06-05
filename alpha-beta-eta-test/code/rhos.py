import os
today = '22-04-19_'

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
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
    parser.add_argument('--obs', default=False,
                        action='store_const', const=True,
                        help='Use e_obs instead of e_piff to calculate modified rho stats')
    parser.add_argument('--cosmobin', default=False,
                        action='store_const', const=True,
                        help='Use Y3 cosmology binning . Useful to calculate the bias of xip.')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args

def measure_rho(data, min_sep = 0.1,  max_sep=300, bin_size=0.2,  sep_units='arcmin',  tag=None, prefix='piff', mod=True,  obs=False,  cosmobin=False ):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    T = data['obs_T']
    p_T = data[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    w1 = p_e1*dt
    w2 = p_e2*dt
    w1obs = e1*dt 
    w2obs = e2*dt 
        
    
    #Modified ellipticities
    if(mod):
        e1 = e1 - np.array(np.mean(e1))
        e2 = e2 - np.array(np.mean(e2))       
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
        w1obs = w1obs - np.array(np.mean(w1obs))
        w2obs = w2obs - np.array(np.mean(w2obs))
        
    ra = data['ra']
    dec = data['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    if(obs):
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    else:
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        #wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
            cat.name = tag + ":"  + cat.name

    
    #bin_config = dict(sep_units = sep_units, bin_slop = 0.1, min_sep = 0.5,  max_sep=300, bin_size=0.2)
    if cosmobin:
        #BINNING USED TO PROPAGATE IN COSMOLOGY
        bin_config = dict(sep_units = sep_units, nbins = 20, min_sep = 2.5, max_sep = 250,)
    else:
        bin_config = dict(sep_units = sep_units, nbins = 20, min_sep = 1.0, max_sep = 250,)
        #bin_config = dict(sep_units = sep_units, bin_slop = 0.1, min_sep = min_sep, max_sep = max_sep, bin_size = bin_size)
 
    

    results = []
    for (cat1, cat2) in [(ecat, ecat), 
                         (decat, decat),
                          (decat, ecat),
                          (wcat, wcat),
                          (decat, wcat),
                          (ecat, wcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    return results

def main():
    import numpy as np
    from src.read_cats import read_data,  toList
    from astropy.io import fits
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

    #STATISTIC USING ONLY RESERVED STARS
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T',  'mag']
 
    exps = toList(args.exps_file)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data))
    data = data[data['mag']<20]
    print("Objects with magnitude <20",  len(data))

   
    
    #BINING FOR ESTIMATING ABE
    min_sep = 0.1;  max_sep=300; bin_size=0.2

    rho0, rho1, rho2, rho3, rho4, rho5 = measure_rho(data,min_sep = min_sep, max_sep = max_sep,
                                                     bin_size = bin_size, mod=args.mod,
                                                     obs=args.obs, cosmobin=args.cosmobin)
    rho0parr = rho0.xip; rho1parr = rho1.xip; rho2parr = rho2.xip
    rho3parr = rho3.xip; rho4parr = rho4.xip; rho5parr = rho5.xip
    rho0marr = rho0.xim; rho1marr = rho1.xim; rho2marr = rho2.xim
    rho3marr = rho3.xim; rho4marr = rho4.xim; rho5marr = rho5.xim

    rhos = [rho0parr, rho0marr, rho1parr, rho1marr, rho2parr,
            rho2marr, rho3parr, rho3marr, rho4parr, rho4marr,
            rho5parr, rho5marr]
    rhos_names = ['RHO0P', 'RHO0M', 'RHO1P', 'RHO1M', 'RHO2P','RHO2M', 'RHO3P', 'RHO3M', 'RHO4P', 'RHO4M', 'RHO5P', 'RHO5M']

    
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = len(rho0parr)
    outdata = np.recarray((nrows, ), dtype=dtype)

  
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    
    covmat = np.diag(np.concatenate((rho0.varxip, rho0.varxim, rho1.varxip,
                                     rho1.varxim, rho2.varxip, rho2.varxim,
                                     rho3.varxip, rho3.varxim, rho4.varxip,
                                     rho4.varxim, rho5.varxip, rho5.varxim)))
    covmathdu = fits.ImageHDU(covmat, name='COVMAT')
    hdul.insert(1, covmathdu)
    
    bin1array = np.array([ -999]*nrows)
    bin2array = np.array([ -999]*nrows)
    angbinarray = np.arange(nrows)
    angarray = np.exp(rho0.meanlogr)

    for j, nam in enumerate(rhos_names):
        array_list = [bin1array, bin2array, angbinarray, np.array(rhos[j]),  angarray ]
        for array, name in zip(array_list, names): outdata[name] = array 
        corrhdu = fits.BinTableHDU(outdata, name=nam)
        hdul.insert(j+2, corrhdu)
        hdul[j+2].header['2PTDATA'] = 'T'
        hdul[j+2].header['NANGLE'] = nrows
        hdul[j+2].header['SIMULATD'] = 'T'
        hdul[j+2].header['BLINDED'] = 'F'


    hdul[1].header['COVDATA'] = True
    hdul[1].header['EXTNAME'] =  'COVMAT'
    hdul[1].header['NAME_0'] =  'RHO0P'
    hdul[1].header['STRT_0'] =  0
    hdul[1].header['LEN_0'] = nrows
    hdul[1].header['NAME_1'] =  'RHO0M'
    hdul[1].header['STRT_1'] =  nrows
    hdul[1].header['LEN_1'] = nrows
    hdul[1].header['NAME_2'] =  'RHO1P'
    hdul[1].header['STRT_2'] =  2*nrows
    hdul[1].header['LEN_2'] = nrows
    hdul[1].header['NAME_3'] =  'RHO1M'
    hdul[1].header['STRT_3'] =  3*nrows
    hdul[1].header['LEN_3'] = nrows
    hdul[1].header['NAME_4'] =  'RHO2P'
    hdul[1].header['STRT_4'] =  4*nrows
    hdul[1].header['LEN_4'] = nrows
    hdul[1].header['NAME_5'] =  'RHO2M'
    hdul[1].header['STRT_5'] =  5*nrows
    hdul[1].header['LEN_5'] = nrows
    hdul[1].header['NAME_6'] =  'RHO3P'
    hdul[1].header['STRT_6'] =  6*nrows
    hdul[1].header['LEN_6'] = nrows
    hdul[1].header['NAME_7'] =  'RHO3M'
    hdul[1].header['STRT_7'] =  7*nrows
    hdul[1].header['LEN_7'] = nrows
    hdul[1].header['NAME_8'] =  'RHO4P'
    hdul[1].header['STRT_8'] =  8*nrows
    hdul[1].header['LEN_8'] = nrows
    hdul[1].header['NAME_9'] =  'RHO4M'
    hdul[1].header['STRT_9'] =  9*nrows
    hdul[1].header['LEN_9'] = nrows
    hdul[1].header['NAME_10'] =  'RHO5P'
    hdul[1].header['STRT_10'] =  10*nrows
    hdul[1].header['LEN_10'] = nrows
    hdul[1].header['NAME_11'] =  'RHO5M'
    hdul[1].header['STRT_11'] =  11*nrows
    hdul[1].header['LEN_11'] = nrows
  

    hdul[2].header['QUANT1'] = 'PeR'; hdul[3].header['QUANT1'] = 'PeR'
    hdul[2].header['QUANT2'] = 'PeR'; hdul[3].header['QUANT2'] = 'PeR'
    hdul[4].header['QUANT1'] = 'PqR'; hdul[5].header['QUANT1'] = 'PqR'
    hdul[4].header['QUANT2'] = 'PqR'; hdul[5].header['QUANT2'] = 'PqR'
    hdul[6].header['QUANT1'] = 'PqR'; hdul[7].header['QUANT1'] = 'PqR'
    hdul[6].header['QUANT2'] = 'PeR'; hdul[7].header['QUANT2'] = 'PeR'
    hdul[8].header['QUANT1'] = 'PwR'; hdul[9].header['QUANT1'] = 'PwR'
    hdul[8].header['QUANT2'] = 'PwR'; hdul[9].header['QUANT2'] = 'PwR'
    hdul[10].header['QUANT1'] = 'PqR'; hdul[10].header['QUANT1'] = 'PqR'
    hdul[10].header['QUANT2'] = 'PwR'; hdul[10].header['QUANT2'] = 'PwR'
    hdul[11].header['QUANT1'] = 'PeR'; hdul[11].header['QUANT1'] = 'PeR'
    hdul[11].header['QUANT2'] = 'PwR'; hdul[11].header['QUANT2'] = 'PwR'

    
    if args.cosmobin:
        hdul.writeto(outpath + 'RHOS_Y3.fits', overwrite=True)
    else:
        hdul.writeto(outpath + 'RHOS.fits', overwrite=True)

if __name__ == "__main__":
    main()
