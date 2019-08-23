#run abn test and save files arrays of parameters going out to a maximum bin. 
#plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model.
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/code/ally3.grizY',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/code/testexp.txt', 
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--bandcombo', default=False,
                        action='store_const', const=True,
                        help='run rho2 for all combination of bands, if false run particular combination defined in band')
    parser.add_argument('--use_reserved', default=True,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--frac', default=0.01, type=float,
                        help='Choose a random fraction of the input stars')
    parser.add_argument('--flask',
                        default='/home2/dfa/sobreira/alsina/catalogs/flask/desy3_6x2pt_lognormal-maskedcats/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--plotpath', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/code/tests/overlapcats/',
                        help='location of the output of the files')
    args = parser.parse_args()
    return args        
        
def main(): 
    from src.read_cats import read_data, toList, read_flask
    from astropy import units
    from astropy.coordinates import SkyCoord
    args = parse_args()
    #Make directory where the ouput data will be
    plotpath = os.path.expanduser(args.plotpath)
    try:
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
    except OSError:
        if not os.path.exists(plotpath): raise

 
    #Reading Mike stars catalog
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T', 'mag']
 
    exps = toList(args.exps_file)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved, frac=args.frac)
    print("Objects",  len(data_stars))
    data_stars = data_stars[data_stars['mag']<20]
    print("Objects with magnitude <20",  len(data_stars))
    coord_stars= SkyCoord(ra=data_stars['ra'], dec=data_stars['dec'], unit='degree')
    ra_stars=coord_stars.ra.wrap_at(180*units.deg)
    for seed in range (201, 401):
        for zbin in range(1, 5):
            for ck in range(1, 3):
                data_galaxies =  read_flask(args.flask, seed, zbin, ck)
                coord_galaxies= SkyCoord(ra=data_galaxies['ra'], dec=data_galaxies['dec'], unit='degree')
                ra_gal=coord_galaxies.ra.wrap_at(180*units.deg)
                plt.clf()
                plt.scatter(ra_gal, data_galaxies['dec'], c='blue', label='flask gals')
                plt.scatter(ra_stars, data_stars['dec'], c='red', label='stars')
                plt.legend(loc='best', fontsize=14)
                plt.tight_layout() 
                print("Printing File")
                #plt.savefig('aja.png')
                plt.savefig(os.path.join(plotpath,'overlap_s%d_z%d_ck%d'%(seed, zbin, ck)))
    
    


    
if __name__ == "__main__":
    main()
