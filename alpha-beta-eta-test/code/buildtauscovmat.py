import os 
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--tausjnfile',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/alltausp_jk.fits',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--xim', default=False,
                        action='store_const', const=True,
                        help='trecorr return xim instead of xip')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args
def main():
    import sys
    import fitsio
    from astropy.io import fits
    import itertools
    import numpy as np
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')

 
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    namesout=[ 'TAU0','TAU2','TAU5']
    if args.xim:
        args.tausjnfile = '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/alltausm_jk.fits'
    alldata =  fitsio.read(args.tausjnfile, ext=1)

    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    a=[i for i in range(0,nrows)]
    b=[j for j in range(0,nrows)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)

    for nam in namesout:
        corr = []
        covmat = np.zeros(shape=(nrows, nrows))
        for i,j in bin_pairs:
            bin1 =  (alldata['ANGBIN']==i)
            bin2 =  (alldata['ANGBIN']==j)
            a = alldata[nam][bin1]; b = alldata[nam][bin2];
            covmat[i,j] = np.cov(a,b)[0][1]
            if (i==j):
                corr.append(np.mean(alldata[nam][bin1]))

        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        covmathdu = fits.ImageHDU(covmat, name='COVMAT')
        hdul.insert(1, covmathdu)
            

        angarray = alldata['THETA'][alldata['JKR']==0]
        valuearray =  np.array(corr)
        bin1array = np.array([ -999]*nrows)
        bin2array = np.array([ -999]*nrows)
        angbinarray = np.arange(nrows)
        array_list = [bin1array, bin2array, angbinarray, valuearray,  angarray ]
        for array, name in zip(array_list, names): outdata[name] = array 

        corrhdu = fits.BinTableHDU(outdata, name=nam)
        hdul.insert(2, corrhdu)
        if args.xim:
            hdul.writeto(outpath + nam + 'm.fits', clobber=True)
        else:
            hdul.writeto(outpath + nam + 'p.fits', clobber=True)
        
 
        

    
if __name__ == "__main__":
    main()

