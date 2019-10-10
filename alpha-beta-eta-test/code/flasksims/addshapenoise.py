today = '08-19-19_'
import numpy as np
import os
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Check which files are missing in flask catalog')
    
    parser.add_argument('--flask_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/desy3_6x2pt_lognormal-maskedcats_v3/', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--metacal_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/Y3_mastercat_7_24_19.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/catalogs/FLASK/desy3_mysp2',
                        help='location of the output of the files')
    parser.add_argument('--sigma_e1', default=0, type=float,
                        help='shapenoise first component')
    parser.add_argument('--sigma_e2', default=0, type=float,
                        help='shapenoise second component')
    parser.add_argument('--seed', default=1,  type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--zbin', default=1 , type=int,
                        help='zbin used, useful to run parallel')
    parser.add_argument('--cookie', default=1 , type=int,
                        help='cookie used, useful to run parallel')
    parser.add_argument('--nz_source',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')  
    args = parser.parse_args()

    return args

def GetSigmasDensities(metacal, nz_source):
    from src.read_cats import read_data_stars, toList, read_metacal
    galkeys = ['e_1','e_2','R11', 'R22']
    sigmas1 = []; sigmas2 = []; densities = []
    for zbin in range(1, 5):
        data_gal = read_metacal(metacal,  galkeys,  zbin=zbin,  nz_source_file=nz_source)
        sigmas1.append(np.std(data_gal['e_1']))
        sigmas2.append(np.std(data_gal['e_2']))
        densities.append(len(data_gal['e_1']))
    return sigmas1, sigmas2, densities

def main():
    import sys; sys.path.append(".")
    from astropy.io import fits
    import numpy as np
    from src.read_cats import read_flask
    
    args = parse_args()

    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    #CALCULATE SIGMA (shape noise) FROM METACAL
    sigma1, sigma2, densities = GetSigmasDensities(args.metacal_cat, args.nz_source)        
    
    names=['RA', 'DEC','GAMMA1', 'GAMMA2']
    forms = ['f4', 'f4', 'f4',  'f4']
    dtype = dict(names = names, formats=forms)
        
    flaskeys = ['RA','DEC','GAMMA1_TRUE', 'GAMMA2_TRUE']
    
    ck=args.cookie
    for seed in range (575, 1000):
        for zbin in range(args.zbin, 5):
            #SKIP ALREADY PROCESS OR missing flask cats
            inname = os.path.join(args.flask_cat, 'src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            if not os.path.isfile(inname):
                #print(inname, "does not exist. Skipping")
                print('src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ), "does not exist. Skipping")
                continue
            else:
                data_flask = read_flask(args.flask_cat,  args.seed, zbin=zbin, ck=1,  keys=flaskeys)
                gamma1sp =  np.random.normal(data_flask['GAMMA1_TRUE'], sigma1[zbin-1])
                gamma2sp =  np.random.normal(data_flask['GAMMA2_TRUE'], sigma2[zbin-1])
                nobjs = len(data_flask['RA'])
                
                mask = np.array( [True]*nobjs)
                if(nobjs > densities[zbin-1]):
                    mask = np.where(mask)[0]
                    mask = np.random.choice(mask, densities[zbin-1], replace=False)
                    print('selecting', len(mask), 'galaxies randomly')
                else:
                    print("Number of objects in sim is smaller than data. Percentual difference: %.3f percent"%(1-(nobjs/densities[zbin-1])) )
                outdata = np.recarray((len(data_flask['RA'][mask]), ), dtype=dtype)
                array_list = [data_flask['RA'][mask], data_flask['DEC'][mask], gamma1sp[mask], gamma2sp[mask]]
                for array, name in zip(array_list, names): outdata[name] = array
                    

                hdu = fits.PrimaryHDU()
                hdul = fits.HDUList([hdu])
                hdul.insert(1, fits.BinTableHDU(outdata, name='FlaskShapeNoise') )
                filename = os.path.join(outpath, 'src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck))
                print("Printing file:", filename)
                hdul.writeto(filename, overwrite=True)
                
        args.zbin = 1
if __name__ == "__main__":
    main()
