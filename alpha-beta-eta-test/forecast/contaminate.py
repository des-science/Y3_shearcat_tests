import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to contaminate a fiducial cosmology using , a dxip contaminant coming from PSF modelling')
    
    parser.add_argument('--original',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/pipeline/2pt_sim_1110_baseline_Y3cov.fits',
                        help='File containing xip to be modified')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/marg__ab_dxi_eq_4_.fits',
                        help='Path for the outputs of this code')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/pipeline/', help='location of the output of the files')
    parser.add_argument('--sup', default=False, 
                        action='store_const', const=True, help='Use the superior limit of the error bar of deltaxip to contaminate')
    parser.add_argument('--inf', default=True, 
                        action='store_const', const=True, help='Use the inferior limit of the error bar of deltaxip to contaminate')
    parser.add_argument('--filename',
                        default='2pt_sim_1110_baseline_Y3cov_contaminated_inf.fits',
                        help='Path for the outputs of this code')  
    
    args = parser.parse_args()
    return args


def get_error(covmatrix, lengths, name):
    import numpy as np
    if name is not None:
        if (name=='xip'):
            start = 0
            end =start + lengths[0]
        elif (name=='xim'):
            start = lengths[0]
            end =start + lengths[1]
        elif (name=='gammat'):
            start = lengths[0] + lengths[1]
            end =start + lengths[2]
        elif (name=='wtheta'):
            start = lengths[0] + lengths[1]+ lengths[2]
            end =start + lengths[3]
        return np.diagonal(covmatrix)[start:end]**0.5
    else:
        print("Correlation function no defined")
        return None


def main():
    import fitsio
    import itertools
    import numpy as np
    import itertools
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        
    
    covmatrixfit_ori=fitsio.read(args.original,ext=1)
    xipfit_ori=fitsio.read(args.original,ext=2)
    ximfit_ori=fitsio.read(args.original,ext=3)

    covmatrixfit_cont=fitsio.read(args.contaminant,ext=1)
    xipfit_cont=fitsio.read(args.contaminant,ext=2)
    ximfit_cont=fitsio.read(args.contaminant,ext=3)
    dxipbin = xipfit_cont['VALUE']
    dximbin = ximfit_cont['VALUE']

    lengths = [len(xipfit_cont), len(ximfit_cont)]
    if args.sup:
        dxipbin += get_error(covmatrixfit_cont, lengths, 'xip')
        dximbin += get_error(covmatrixfit_cont, lengths, 'xim')
    if args.inf:
        dxipbin -= get_error(covmatrixfit_cont, lengths, 'xip')
        dximbin -= get_error(covmatrixfit_cont, lengths, 'xim')

    nrows = 20
    nbins=4
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        binp = (xipfit_ori['BIN1']==i)&(xipfit_ori['BIN2']==j)
        binm = (ximfit_ori['BIN1']==i)&(ximfit_ori['BIN2']==j)
        idxbinsp =  list(itertools.compress(range(len(binp)),  binp))
        idxbinsm =  list(itertools.compress(range(len(binm)),  binm))
        if (len(idxbinsp)!=0): xipfit_ori['VALUE'][binp] -=dxipbin[binp]
        if (len(idxbinsm)!=0): ximfit_ori['VALUE'][binm] -=dxipbin[binm]

    hdulist = fits.open(args.original)
    #delete all xip but saving header
    oldheaders =  [hdulist[2].header, hdulist[3].header]
    hdulist.pop(index=2);
    hdulist.pop(index=2);
 
    xiphdu = fits.BinTableHDU(xipfit_ori)
    ximhdu = fits.BinTableHDU(ximfit_ori)
    hdulist.insert(2, xiphdu)
    hdulist.insert(3, ximhdu)
    hdulist[2].header = oldheaders[0]
    hdulist[3].header = oldheaders[1]
    filename = os.path.join(outpath, args.filename) 
    hdulist.writeto(filename, overwrite=True)
    print(filename, 'written')
    
if __name__ == "__main__":
    main()
