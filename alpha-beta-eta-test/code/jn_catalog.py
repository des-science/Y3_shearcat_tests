import numpy as np
import os
import kmeans_radec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--metacal_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/Y3_mastercat_7_24_19.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--nz_source',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog') 
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/catalogs/JK/jn_njk1000',
                        help='location of the output of the files')
    parser.add_argument('--njks', default=1000,type=int, 
                        help='Number of patches to divide footprint')
    parser.add_argument('--plot', default=False,
                        action='store_const', const=True, help='Plot jkmeans')
       
    
    args = parser.parse_args()

    return args

def jk_kmeans(ra_sam, dec_sam, ra,dec,njk,plot=False):
 
    print("Running kmeans")
    from astropy.coordinates import SkyCoord, Angle
    from astropy import units
    radec = np.zeros((len(ra),2)); radec_sam = np.zeros((len(ra_sam),2))
    radec[:,0] = ra; radec_sam[:,0] = ra_sam
    radec[:,1] = dec; radec_sam[:,1] = dec_sam
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
        '''
        coords = SkyCoord(ra=ra, dec=dec, unit='degree')
        ra = coords.ra.wrap_at(180 * units.deg)
        dec = coords.dec
        plt.figure()
        plt.scatter(ra,dec,c=jk,lw=0,cmap='Paired',rasterized=True)
        plt.xlabel(r'RA',fontsize=12)
        plt.ylabel(r'Dec',fontsize=12)
        plt.tight_layout()
        plt.savefig('jk_kmeans.png')
        '''
    
    jk = km.find_nearest(radec)
    if not km.converged:
        print('k means did not converge')

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
        
    

    ##Format of the fit file output
    names=['RA', 'DEC', 'GAMMA1', 'GAMMA2', 'JKID']
    forms = ['f4', 'f4', 'f4',  'f4', 'i4']
    dtype = dict(names = names, formats=forms)
 
    
    
    nobjs = len(read_metacal(args.metacal_cat,  ['ra']))
    mask = np.array( [True]*nobjs)
    mask = np.where(mask)[0]
    mask = np.random.choice(mask, 100000 , replace=False)
    data_sam=read_metacal(args.metacal_cat,  ['ra', 'dec'])[mask]
   
    print('sampling jk regions number of objs', len(data_sam))

    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    for zbin in range(1, 5):
        print('Starting measurement for zbin', zbin)
        
        data_gal = read_metacal(args.metacal_cat,  galkeys,  zbin=zbin,  nz_source_file=args.nz_source)
        njk = args.njks
        ##TODO generate km first an later finnearest,
        jkindexes_gals = jk_kmeans(data_sam['ra'], data_sam['dec'], data_gal['ra'][:1000000], data_gal['dec'][: 1000000],njk,  plot=args.plot)
        
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


            
 
    
if __name__ == "__main__":
    main()
