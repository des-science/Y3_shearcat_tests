import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce Tau correlations, i.e correlation among galaxies and reserved stars')
    
    parser.add_argument('--metacal_cat',
                        default='/home/dfa/sobreira/alsina/catalogs/y3_master/Y3fullmaster/Y3_mastercat_v2_6_20_18.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--nz_source',
                        default='/home/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        help='Indexes catalog to select galaxies in a particular redshift bin in Metacal')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()
    return args

def read_metacal(filename,  keys,  zbin=None,  nz_source_file=None):
    import h5py as h
    import numpy as np

    dgamma = 2*0.01
    
    f = h.File(filename, 'r')
    cat =  f['catalog/metacal/unsheared']
    print('catalog/metacal/unsheared' + ' keys:',cat.keys())
    nrows = len(np.array( cat['ra'] ))
    formats = []
    for key in keys:
        if key == 'ccd' or key == 'tiling':
            formats.append('i2')
        elif key == 'exp' or 'flag' in key:
            formats.append('i4')
        elif key == 'band':
            formats.append('a1')
        else:
            formats.append('f4')
    data = np.recarray(shape=(nrows,), formats=formats, names=keys)
    for key in keys:  data[key] = np.array(cat[key])    
    print('made recarray')

    select = np.array(f['index/select'])
    select_1p = np.array(f['index/select_1p'])
    select_1m = np.array(f['index/select_1m'])
    select_2p = np.array(f['index/select_2p'])
    select_2m = np.array(f['index/select_2m']) 
    if zbin is None:
        R11s = (data['e_1'][select_1p].mean() - data['e_1'][select_1m].mean() )/dgamma
        R22s = (data['e_2'][select_2p].mean() - data['e_2'][select_2m].mean() )/dgamma
        data = data[select]
        
    else:
        print("Reading only data from bin",  zbin)
        n = h.File(nz_source_file, 'r')
        zbin_array = np.array(n['nofz/zbin'])
        ind = np.where( zbin_array==zbin )[0]
        ind_1p = np.where(np.array(n['nofz/zbin_1p'])==zbin)
        ind_1m = np.where(np.array(n['nofz/zbin_1m'])==zbin)
        ind_2p = np.where(np.array(n['nofz/zbin_2p'])==zbin)
        ind_2m = np.where(np.array(n['nofz/zbin_2m'])==zbin)
        R11s = (data['e_1'][select_1p][ind_1p].mean() - data['e_1'][select_1m][ind_1m].mean() )/dgamma
        R22s = (data['e_2'][select_2p][ind_2p].mean() - data['e_2'][select_2m][ind_2m].mean() )/dgamma
        data = data[select][ind]
        
    data['e_1'] = data['e_1']/(R11s + np.mean(data['R11']))
    data['e_2'] = data['e_2']/(R22s + np.mean(data['R22']))
    print('Metal read sucesfully',  len(data),  'objects')

    return data

def main():
    import sys; sys.path.append(".")
    import numpy as np
    #from src.read_cats import read_metacal

    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise

    galkeys = ['ra','dec','e_1','e_2','R11','R22','T', 'psf_T']
    for zbin in range(1, 5):
        print('Plotting histogram for zbin=', zbin)
        data_galaxies = read_metacal(args.metacal_cat, galkeys, zbin=zbin,nz_source_file=args.nz_source)
        Tgal = data_galaxies['T']
        print('Min Tgal:',  min(Tgal), 'Max Tgal', max(Tgal) )
        plt.clf()
        nbins = 20
        plt.hist(Tgal, bins=np.linspace(min(Tgal), max(Tgal), nbins) )
        plt.tight_layout()
        filename = os.path.join(plotspath,'Tgal_zbin%d.png'%(zbin))
        print("Printing File", filename)
        plt.savefig(filename)

        Tpsf = data_galaxies['psf_T']
        print('Min Tpsf:',  min(Tpsf), 'Max Tpsf', max(Tpsf) )
        
        plt.hist(Tpsf, bins=np.linspace(min(Tpsf), max(Tpsf), nbins) )
        plt.tight_layout()
        filename = os.path.join(plotspath,'Tpsf_zbin%d.png'%(zbin))
        print("Printing File", filename)
        plt.savefig(filename)

        ratio = Tgal/Tpsf
        print('Min Tpsf:',  min(ratio), 'Max Tpsf', max(ratio) )
        plt.clf()
        filename = os.path.join(plotspath,'ratio_zbin%d.png'%(zbin))
        print("Printing File", filename)
        plt.savefig(filename)
        
if __name__ == "__main__":
    main()
