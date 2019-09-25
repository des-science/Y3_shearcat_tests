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
    parser.add_argument('--seed', default=1,  type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--zbin', default=1 , type=int,
                        help='zbin used, useful to run parallel')
    parser.add_argument('--cookie', default=1 , type=int,
                        help='cookie used, useful to run parallel')
    args = parser.parse_args()

    return args

def main():
    import sys; sys.path.append(".")
    from astropy.io import fits
    import numpy as np
    from src.runcorr import measure_tau
    from src.read_cats import read_data, toList, read_flask
    import h5py as h
    
    args = parse_args()

    ck=args.cookie
    for seed in range (args.seed, 301):
        for zbin in range(args.zbin, 5):
            #SKIP ALREADY PROCESS OR missing flask cats
            inname = os.path.join(args.flask_cat, 'src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
            if not os.path.isfile(inname):
                #print(inname, "does not exist. Skipping")
                print('src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
                continue            
        args.zbin = 1
if __name__ == "__main__":
    main()
