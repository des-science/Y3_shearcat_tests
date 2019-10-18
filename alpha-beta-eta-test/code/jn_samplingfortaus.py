import numpy as np
import os
import kmeans_radec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='When the catalog is too big we need to assing objects splitting the catalogs. This code does the splitting for the assignation of JK regions to galaxies ')
    
    parser.add_argument('--metacal_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/Y3_mastercat_7_24_19.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--nz_source',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--piff_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/code/ally3.grizY',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--use_reserved', default=True,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/catalogs/JK/jn_njk1000_last',
                        help='location of the output of the files')
    parser.add_argument('--njks', default=1000,type=int, 
                        help='Number of patches to divide footprint')
    parser.add_argument('--plot', default=False,
                        action='store_const', const=True, help='Plot jkmeans')
       
    
    args = parser.parse_args()

    return args

def find_regions(ra_sam, dec_sam, njk,plot=False):
    print("Finding Region")
    from astropy.coordinates import SkyCoord, Angle
    from astropy import units
    radec_sam = np.zeros((len(ra_sam),2))
    radec_sam[:,0] = ra_sam
    radec_sam[:,1] = dec_sam
    km = kmeans_radec.kmeans_sample(radec_sam,njk,maxiter=100,tol=1e-05)
    if plot:
        jk_sam = km.find_nearest(radec_sam)
        coords = SkyCoord(ra=ra_sam, dec=dec_sam, unit='degree')
        ra_sam = coords.ra.wrap_at(180 * units.deg)
        dec_sam = coords.dec
        plt.figure()
        plt.scatter(ra_sam,dec_sam,c=jk_sam,lw=0,cmap='Paired',rasterized=True)
        plt.xlabel(r'RA',fontsize=12)
        plt.ylabel(r'Dec',fontsize=12)
        plt.tight_layout()
        plt.savefig('jk_kmeans.png', dpi=200)
    print(njk, 'regions were produced')
    return km
def find_index(km, ra, dec, plot=False):
    print("Finding Region")
    from astropy.coordinates import SkyCoord, Angle
    from astropy import units
    radec = np.zeros((len(ra),2))
    radec[:,0] = ra
    radec[:,1] = dec
    jk = km.find_nearest(radec)
    print(len(ra), 'indexes were assigned')
    if plot:
        jk = km.find_nearest(radec)
        coords = SkyCoord(ra=ra, dec=dec, unit='degree')
        name = ra[0]
        ra = coords.ra.wrap_at(180 * units.deg)
        dec = coords.dec
        plt.clf()
        plt.scatter(ra,dec,c=jk,lw=0,cmap='Paired',rasterized=True)
        plt.xlabel(r'RA',fontsize=12)
        plt.ylabel(r'Dec',fontsize=12)
        plt.ylim([ - 70, 10])
        plt.xlim([ - 70, 110])
        plt.tight_layout()
        plt.savefig('jk_kmeans%.3f.png'%(name), dpi=200)
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
        
    

    ##Format of the fit file output
    names=['RA', 'DEC', 'GAMMA1', 'GAMMA2', 'JKID']
    forms = ['f4', 'f4', 'f4',  'f4', 'i4']
    dtype = dict(names = names, formats=forms)
 

    ##CREATING JK REGIONS USING full metacal
    nsamples = 10000
    nguests = 1000000
    nobjs = len(read_metacal(args.metacal_cat,  ['ra']))
    mask = np.array( [True]*nobjs)
    mask = np.where(mask)[0]
    mask = np.random.choice(mask, nsamples , replace=False)
    data_sam=read_metacal(args.metacal_cat,  ['ra', 'dec'])[mask]
   
    print('sampling jk regions number of objs', len(data_sam))
    njk = args.njks
    km = find_regions(data_sam['ra'], data_sam['dec'], njk,plot=args.plot)


    ##Assign glaxies to the regions for each redshift bin
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    for zbin in range(1, 5):
        print('Starting measurement for zbin', zbin)
        
        data_gal = read_metacal(args.metacal_cat,  galkeys,  zbin=zbin,  nz_source_file=args.nz_source)
        ngals = len(data_gal['ra'])
        
        jkindexes_gals = []
        for i in range(ngals//nguests + 1):
            inf = i*nguests
            sup = (i + 1)*nguests
            if sup == ngals//nguests: sup = None
            jkindexes_gals += find_index( km, data_gal['ra'][inf:sup], data_gal['dec'][inf:sup]).tolist()


        print(ngals, len(jkindexes_gals))
        nrows = len(data_gal['ra'])
        outdata = np.recarray((nrows, ), dtype=dtype)
        array_list = [data_gal['ra'], data_gal['dec'], data_gal['e_1'], data_gal['e_2'], jkindexes_gals]
        for array, name in zip(array_list, names): outdata[name] = array 

        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        corrhdu = fits.BinTableHDU(outdata, name='METACAL_JK_z%d'%(zbin) )
        hdul.insert(1, corrhdu)

           
        filename = os.path.join(outpath, 'metacal_JK_z%d.fits'%(zbin))
            
        print("Printing file:", filename)
        hdul.writeto(filename, overwrite=True)

        
    ##Assign stars to the regions
    exps = toList(args.exps_file)
    starskeys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T', 'piff_e1', 'piff_e2', 'piff_T', 'mag']
 
    data_stars = read_data_stars(exps, args.piff_cat , starskeys, limit_bands=args.bands, use_reserved=args.use_reserved)
    nstars = len(data_stars['ra'])

    jkindexes_stars = []
    for i in range(nstars//nguests + 1):
        inf = i*nguests
        sup = (i + 1)*nguests
        if sup == nstars//nguests: sup = None
        jkindexes_stars += find_index( km, data_stars['ra'][inf:sup], data_stars['dec'][inf:sup]).tolist()

    print(nstars, len(jkindexes_stars))
    nrows = len(data_stars['ra'])
    forms = ['f4']*len(starskeys) + ['i4']
    names = starskeys + ['JKID']
    dtype = dict(names = names, formats=forms)
    outdata = np.recarray((nrows, ), dtype=dtype)
    array_list = [data_stars[key] for key in starskeys ]+ [jkindexes_stars]
    for array, name in zip(array_list, names): outdata[name] = array 

    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    corrhdu = fits.BinTableHDU(outdata, name='Y3STAR_JK')
    hdul.insert(1, corrhdu)
    
    filename = os.path.join(outpath, 'Y3_reserved_STAR_JK.fits')
            
    print("Printing file:", filename)
    hdul.writeto(filename, overwrite=True)


            
 
    
if __name__ == "__main__":
    main()
