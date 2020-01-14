import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import numpy as np

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='samples total PSF bias using rhos cosmology and previous samples of abe')
    parser.add_argument('--samplesabe', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/samples_abe.fits')
    parser.add_argument('--rhoscosmo', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS_Y3.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--plots', default=False,
                        action='store_const', const=True, help='Plot correlations functions')
    
    parser.add_argument('--outpath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the final contaminant')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    parser.add_argument('--filename',
                        default='marg_abe_eq4_contaminant.fits',
                        help='Filename of the final contamination')

    args = parser.parse_args()

    return args
                        
def corner_plot(samples, labels, filename, title=None):
    import corner
    import matplotlib.ticker as ticker
    #burn = 5000
    plt.clf()
    fig = corner.corner(np.c_[samples].T, labels=labels, 
                        quantiles=[0.16, 0.5, 0.84],  #-1sigma,0sigma,1sigma
                        levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-9./2)), #1sigma, 2sigma and 3sigma contours
                        show_titles=True, title_kwargs={"fontsize": 16}, title_fmt= '.4f', 
                        smooth1d=None, plot_contours=True, 
                        no_fill_contours=False, plot_density=True, use_math_text=True)
    for i in range(len(fig.axes)):
        fig.axes[i].locator_params(axis='x', nbins=2)
        fig.axes[i].locator_params(axis='y', nbins=2)
        fig.axes[i].tick_params(axis='x', rotation =0, labelsize=16)
        fig.axes[i].tick_params(axis='y', rotation =90, labelsize=16)
        fig.axes[i].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
        fig.axes[i].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
  
    if title is not None:
        plt.suptitle(title,  fontsize=24,  color='blue', x=0.8  )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename,  dpi=150)
    plt.close(fig)
    print(filename, "Printed")
def plotcorrmat(cov):
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    clb = plt.colorbar()
    plt.tight_layout()
    
                        
def main():
    import fitsio
    import getdist
    from getdist import plots, MCSamples
    from astropy.io import fits
    import itertools
    from src.maxlikelihood import percentiles, bestparameters
    from src.readfits import  read_rhos
    from src.chi2 import chi2nu
    import numpy as np

    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise

    abesamps = fitsio.read(args.samplesabe)
    #Burning just to be sure
    frac= 0.0
    abesamps =abesamps[int(frac*len(abesamps) ): ]

    veclist = [abesamps['a1'], abesamps['a2'], abesamps['a3'], abesamps['a4'],  abesamps['b1'], abesamps['b2'], abesamps['b3'], abesamps['b4'], abesamps['e1'], abesamps['e2'], abesamps['e3'], abesamps['e4']]
    covmat = np.cov(veclist)

    #PLOT COVMAT
    if(args.plots):
        plt.clf()
        lengths = [4, 4, 4]
        plotcorrmat(covmat)
        plt.title(r'$\alpha \mid \beta \mid \eta $')
        pos_lines = [ -0.5]
        for i in range(len(lengths)):
            pos_lines.append(pos_lines[i] + lengths[i])
        pos_lines = pos_lines[1:-1]
        for line in pos_lines:
            plt.axvline(x=line, c='k', lw=1, ls='-')
            plt.axhline(y=line, c='k', lw=1, ls='-')
        plt.tight_layout()
        filename = os.path.join(plotspath,'Covariancematrix_samplesabe.png')
        print('Printing', filename)
        plt.savefig(filename, dpi=200)
        
    #PLOT marginalized
    if(args.plots): 
        #veclist = [abesamps['a1'], abesamps['b1'], abesamps['e1'], abesamps['a2'], abesamps['b2'], abesamps['e2'], abesamps['a3'], abesamps['b3'], abesamps['e3'], abesamps['a4'], abesamps['b4'], abesamps['e4']]
        veclist = [abesamps['a2'], abesamps['b2'], abesamps['e2']]
        #mcmcpars = percentiles(veclist, nsig=1) 
        #print( ' mcmc parameters xi+',  'nsig=', 1, ' percentiles: ',  mcmcpars)
        filename = os.path.join(plotspath,'cornerabe.png')
        #corner_plot(veclist, ['a1', 'b1', 'e1'], filename, title=None)

        #samples = MCSamples(samples=veclist, names=['a1', 'b1', 'e1', 'a2', 'b2', 'e2', 'a3', 'b3', 'e3', 'a4', 'b4', 'e4'],
        #                    labels=[r'\alpha^{1}',r'\beta^{1}',r'\eta^{1}',r'\alpha^{2}',r'\beta^{2}',r'\eta^{2}',r'\alpha^{3}',r'\beta^{3}',r'\eta^{3}',r'\alpha^{4}',r'\beta^{4}',r'\eta^{4}'])
        samples = MCSamples(samples=veclist, names=['a2', 'b2', 'e2'], labels=[r'\alpha^{2}',r'\beta^{2}',r'\eta^{2}'])
        g = plots.getSubplotPlotter()
        g.settings.plot_meanlikes = False
        g.settings.alpha_factor_contour_lines = True
        #g.settings.axis_marker_lw = 5
        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add=0.35
        g.settings.title_limit_fontsize = 16
        g.settings.figure_legend_loc = 'best'
        g.settings.rcSizes(axes_fontsize = 12, lab_fontsize=20, legend_fontsize =40)
        g.triangle_plot([samples], filled_compare=[True], 
                    line_args=[{'ls':'solid', 'lw':2, 'color':'green'}],
                    contour_colors=['green'], title_limit=1)
        #g.add_legend(legend_labels=[legend_name], fontsize=36, legend_loc=(-3.5,7))
        filename = os.path.join(plotspath,'getdistabe.png')
        print('Printing', filename)
        plt.tight_layout()
        plt.savefig(filename, dpi=200)


        
    ##SAMPLING BIAS
    meanr, rhos,  covmat = read_rhos(args.rhoscosmo)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m = rhos

    alist = [abesamps['a1'], abesamps['a2'], abesamps['a3'], abesamps['a4']]
    blist = [abesamps['b1'], abesamps['b2'], abesamps['b3'], abesamps['b4']]
    elist = [abesamps['e1'], abesamps['e2'], abesamps['e3'], abesamps['e4']]

    nbins=4
    a=[i for i in range(nbins)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2): bin_pairs.append(p)
    veclist = []
    
    for z in range(len(abesamps['a1'])):
        dxip = [alist[i][z]*alist[j][z]*rho0p + blist[i][z]*blist[j][z]*rho1p + elist[i][z]*elist[j][z]*rho3p + (blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2p + (blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4p + (elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5p for i,j in bin_pairs]
        dxim = [alist[i][z]*alist[j][z]*rho0m + blist[i][z]*blist[j][z]*rho1m + elist[i][z]*elist[j][z]*rho3m + (blist[i][z]*alist[j][z] + blist[j][z]*alist[i][z])*rho2m + (blist[i][z]*elist[j][z] + blist[j][z]*elist[i][z])*rho4m + (elist[i][z]*alist[j][z] +elist[j][z]*alist[i][z])*rho5m for i,j in bin_pairs]
        dxi = dxip + dxim
        veclist.append(np.concatenate(np.c_[dxi]))
    covmat = np.cov(np.c_[veclist].T)

    #PLOT COVMAT
    if(args.plots):
        plt.clf()
        plt.cla()
        plt.close()
        lengths = [200, 200]
        plotcorrmat(covmat)
        plt.title(r'$\delta \xi_{+}^{PSF} \mid \delta \xi_{-}^{PSF} $')
        pos_lines = [0]
        for i in range(len(lengths)):
            pos_lines.append(pos_lines[i] + lengths[i])
        pos_lines = pos_lines[1:-1]
        for line in pos_lines:
            plt.axvline(x=line, c='k', lw=1, ls='-')
            plt.axhline(y=line, c='k', lw=1, ls='-')
        plt.tight_layout()
        filename = os.path.join(plotspath,'Covariancematrix_totalbiases.png')
        print('Printing', filename)
        plt.savefig(filename, dpi=200)

        plt.clf()
        plt.cla()
        plt.close()
        plotcorrmat(covmat[:200, :200])
        plt.title(r'$\delta \xi_{+}^{PSF}$')
        plt.tight_layout()
        filename = os.path.join(plotspath,'Covariancematrix_totalbiasxip.png')
        print('Printing', filename)
        plt.savefig(filename, dpi=200)

        plt.clf()
        plt.cla()
        plt.close()
        plotcorrmat(covmat[:20, :20])
        plt.title(r'$\delta \xi_{+}^{PSF} (1,1) $')
        plt.tight_layout()
        filename = os.path.join(plotspath,'Covariancematrix_totalbiasxip11.png')
        print('Printing', filename)
        plt.savefig(filename, dpi=200)
    

    #Writting final contamination
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])

    dxip_list = []; dxim_list = []; bin1_list = []; bin2_list = []; angbin_list = []; ang_list = []
    covxip_list = []; covxim_list = []
    
    nbins=4
    a=[i for i in range(nbins)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2): bin_pairs.append(p)
    al = bestparameters(alist)
    bl = bestparameters(blist)
    el = bestparameters(elist)
    for i,j in bin_pairs:
         dxip = al[i]*al[j]*rho0p + bl[i]*bl[j]*rho1p + el[i]*el[j]*rho3p + (bl[i]*al[j] + bl[j]*al[i])*rho2p + (bl[i]*el[j] + bl[j]*el[i])*rho4p + (el[i]*al[j] +el[j]*al[i])*rho5p
         dxim = al[i]*al[j]*rho0m + bl[i]*bl[j]*rho1m + el[i]*el[j]*rho3m + (bl[i]*al[j] + bl[j]*al[i])*rho2m + (bl[i]*el[j] + bl[j]*el[i])*rho4m + (el[i]*al[j] +el[j]*al[i])*rho5m

        
         ang_list.append(meanr)
         bin1_list.append(np.array( [i + 1]*len(meanr)))
         bin2_list.append(np.array( [j + 1]*len(meanr)))
         angbin_list.append(np.arange(len(meanr)))
         dxip_list.append(dxip)
         dxim_list.append(dxim)

    hdul.insert(1, fits.ImageHDU(covmat, name='COVMAT'))
    bin1array = np.concatenate(bin1_list)
    bin2array = np.concatenate(bin2_list)
    angbinarray = np.concatenate(angbin_list)
    valuearray = np.concatenate(dxip_list)
    angarray = np.concatenate(ang_list)
        
    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f4',  'f4']
    dtype = dict(names = names, formats=forms)
    nrows = len(angarray)
    outdata = np.recarray((nrows, ), dtype=dtype)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xip')
    hdul.insert(2, corrhdu)
    valuearray = np.concatenate(dxim_list)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xim')
    hdul.insert(3, corrhdu)
    outname = os.path.join(outpath, args.filename)
    hdul.writeto(outname, overwrite=True)
    print(outname,'Written!')

   
  
    
if __name__ == "__main__":
    main()



        
