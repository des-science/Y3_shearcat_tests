import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce Tau correlations, i.e correlation among galaxies and reserved stars')
    
    parser.add_argument('--metacal_cat',
                        #default='/home/dfa/sobreira/alsina/catalogs/y3_master/Y3fullmaster/Y3_mastercat_v2_6_20_18.h5',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/Y3_mastercat_7_24_19.h5',
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--nz_source',
                        #default='/home/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Indexes catalog to select galaxies in a particular redshift bin in Metacal')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()
    return args

def read_metacal(filename,  keys,  zbin=None,  nz_source_file=None,  size_cut=False):
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
    print('Metacal read sucesfully',  len(data),  'objects. zbin:', zbin)

    if size_cut:
        flag = (data['T']/data['psf_T'] >  0.5)
        data = data[flag]
        print('Metacal have',  len(data),  'objects after Tgal/Tpsf>.5 cut. zbin:', zbin)
    
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
    size_cut = False
    #for zbin in [1, 2, 3, 4, None]:
    for zbin in [None]:
        print('Plotting histogram for zbin=', zbin)
        data_galaxies = read_metacal(args.metacal_cat, galkeys, zbin=zbin,nz_source_file=args.nz_source, size_cut=size_cut)

        Tgal = data_galaxies['T']
        Tpsf = data_galaxies['psf_T']
        ratio = Tgal/Tpsf
        
        print('Len Tgal', len(Tgal), 'Min Tgal:',  min(Tgal), 'Max Tgal', max(Tgal) )
        #Tgal = Tgal[Tgal<100]
        print('Len Tgal', len(Tgal), 'after selecting only Tgal<4')
        plt.clf()
        nbins = 10000000
        if zbin is not None: label = 'zbin:%d'%(zbin)
        else: label = None
        plt.hist(Tgal, bins=np.linspace(min(Tgal), max(Tgal), nbins), label=label )
        plt.ylabel('Counts'); plt.xlabel('Tgal')
        plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'Tgal_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'Tgal.png')
        print("Printing File", filename)
        plt.savefig(filename)

        plt.clf()
        print('Min Tpsf:',  min(Tpsf), 'Max Tpsf', max(Tpsf) )
        plt.hist(Tpsf, bins=np.linspace(min(Tpsf), max(Tpsf), nbins), label=label)
        plt.ylabel('Counts'); plt.xlabel('Tpsf')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'Tpsf_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'Tpsf.png')
        print("Printing File", filename)
        plt.savefig(filename)

        plt.clf()
        
        print('Len ratio', len(ratio), 'Min ratio:',  min(ratio), 'Max ratio', max(ratio) )
        #ratio = ratio[ratio<200]
        print('Len ratio', len(ratio), 'after cutting ratio>10')
        plt.clf()
        #plt.title('ratio zbin:%d'%(zbin))
        plt.hist(ratio, bins=np.linspace(min(ratio), max(ratio), nbins), label=label)
        plt.ylabel('Counts'); plt.xlabel(r'$\frac{Tgal}{Tpsf}$')
        plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'ratio_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'ratio.png')
        print("Printing File", filename)
        plt.savefig(filename)
        
if __name__ == "__main__":
    main()
