import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the fitting problem of system ofequatiosn, plotting correlations and final correlation function withbias')
    
    parser.add_argument('--taus',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS.fits',
                         #default='/home2/dfa/sobreira/alsina/catalogs/flask/taus/taus_src-cat_s280_z1_ck1.fits',
                        help='Fits file containing all tau stats used to estimate abe')
    parser.add_argument('--rhos', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--rhoscosmo', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS_Y3.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--samecont', default=False,
                        action='store_const', const=True,
                        help='Include xim correlation functions in the vectors. Then use only one contamination for both xip and xim ')
    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--minscale', default=None, type=float, 
                        help='Limit the analysis to certain minimum scale')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Run all tomographic correlations')
    parser.add_argument('--plots', default=False,
                        action='store_const', const=True, help='Plot correlations functions')
    parser.add_argument('--outpath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the final contaminant')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    parser.add_argument('--eq', default=4, type=int, 
                        help='Select equations to be used for istance --eq=0, 4 represent the whole system of equations')
    parser.add_argument('--abe', default=False,
                        action='store_const', const=True, help='Run alpha, beta and eta test.')
    parser.add_argument('--ab', default=False,
                        action='store_const', const=True, help='Run alpha and beta test.')
    parser.add_argument('--ae', default=False,
                        action='store_const', const=True, help='Run alpha and eta test.')
    parser.add_argument('--be', default=False,
                        action='store_const', const=True, help='Run beta and eta test.')
    parser.add_argument('--a', default=False,
                        action='store_const', const=True, help='Run only alpha test.')
    parser.add_argument('--b', default=False,
                        action='store_const', const=True, help='Run only beta test.')
    parser.add_argument('--e', default=False,
                        action='store_const', const=True, help='Run only eta test.')
    
    
    
    args = parser.parse_args()

    return args
def corrmatrix(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    return corr
def writexibias(samples, datarhos,  nsig=1, plots=False, nameterms='terms_dxi.png',dxiname='dxi.png',namecovmat='covm_pars.png',filename='dxi.fits'):
    from src.readfits import read_rhos
    from src.maxlikelihood import bestparameters, percentiles
    from src.plot_stats import pretty_rho
    from astropy.io import fits
    import numpy as np

    mcmcpars = percentiles(samples, nsig=nsig) 
    print( ' mcmc parameters xip',  'nsig=', nsig, ' percentiles: ',  mcmcpars)

    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    #plot covariance matrix of parameters alpha, beta and eta.
    if plots:
        par_matcov = np.cov(samples)
        corr=corrmatrix(par_matcov)
        cov_vmin=np.min(corr)
        plt.clf()
        plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
                   aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
        plt.colorbar()
        plt.title(r'$\alpha \mid \beta \mid \eta $')
        plt.tight_layout()
        print('Printing', namecovmat)
        plt.savefig(namecovmat, dpi=500)
        

    a = b = n = 0; vara =  varb =  varn = 0; covab = covan = covbn = 0
    bestpar = bestparameters(samples)
    print("Best pars", bestpar)
    par_matcov = np.cov(samples) 
    if (par_matcov.size==1 ): variances = par_matcov
    else: variances = np.diagonal(par_matcov)
    covariances = sum( (par_matcov[i,i+1: ].tolist() for i in range(len(samples) - 1)) , [] )
    if(len(samples)==3):
        a, b, n = bestpar
        vara, varb, varn =  variances
        covab, covan, covbn =  covariances
    elif(len(samples)==2):
        a, b = bestpar
        vara, varb =  variances
        covab =  covariances[0]
    elif(len(samples)==1):
        a =  bestpar[0]
        vara =  variances
    else:
        print("Warning, test type not defined")

    meanr, rhos, covrhos =  datarhos
    rho0, rho1, rho2, rho3, rho4, rho5 = rhos
    cov_rho0, cov_rho1, cov_rho2, cov_rho3, cov_rho4, cov_rho5 =  covrhos
    sig_rho0 =  np.sqrt(np.diag(cov_rho0))
    sig_rho1 =  np.sqrt(np.diag(cov_rho1))
    sig_rho2 =  np.sqrt(np.diag(cov_rho2))
    sig_rho3 =  np.sqrt(np.diag(cov_rho3))
    sig_rho4 =  np.sqrt(np.diag(cov_rho4))
    sig_rho5 =  np.sqrt(np.diag(cov_rho5))

    #Ploting each term of the bias
    if(plots):
        xlim = [2., 300.]
        #supposing that a,b and n are idependent of rhos(scale independent)
        var0 = ((2*a*rho0)**2)*vara +  (a**2)*(sig_rho0**2)
        var1 = ((2*b*rho1)**2)*varb +  (b**2)*(sig_rho1**2)
        var2 = ((2*n*rho3)**2)*varn +  (n**2)*(sig_rho3**2)
        varab =  vara*(b**2) + varb*(a**2) + 2*covab*(a*b)
        #varab = ((a*b)**2)*( (vara/((a)**2)) + (varb/((b)**2)) + 2*covab/(a*b) )
        var3 = 4*( (rho2**2)*varab + (sig_rho2**2)*((a*b)**2)  )
        #var3 = 4*((a*b*rho2p)**2)*( varab/((a*b)**2) + (sig_rho2/rho2p)**2 )
        varbn =  varn*(b**2) + varb*(n**2) + 2*covbn*(b*n)
        #varbn = ((n*b)**2)*( (varn/((n)**2)) + (varb/((b)**2)) + 2*covbn/(b*n) ) 
        var4 = 4*( (rho4**2)*varbn + (sig_rho4**2)*((n*b)**2)  )
        #var4 = 4*((n*b*rho4p)**2)*(varbn/((b*n)**2) + (sig_rho4/rho4p)**2)
        varan = varn*(a**2) + vara*(n**2) + 2*covan*(a*n)
        #varan = ((n*a)**2)*( (varn/((n)**2)) + (vara/((a)**2)) + 2*covan/(a*n) ) 
        var5 = 4*( (rho5**2)*varan + (sig_rho5**2)*((n*a)**2)  )
        #var5 = 4*((n*a*rho5p)**2)*(varan/((a*n)**2) + (sig_rho5/rho5p)**2) 
        plt.clf()
        lfontsize = 7
        if (len(samples)==3):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(np.diag(cov_rho0)), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (n**2)*rho3, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*n)*rho4, np.sqrt(var4), legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*n*a)*rho5, np.sqrt(var5), legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (len(samples)==2):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (len(samples)==1):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
    
    #supposing that a,b and n are idependent of rhos(scale independent)
    dxi = (a**2)*rho0 + (b**2)*rho1 + (n**2)*rho3 + (2*a*b)*rho2 + (2*b*n)*rho4 + (2*n*a)*rho5
    f1 = 2*(a*rho0 + b*rho2 + n*rho5)     
    f2 = 2*(b*rho1 + a*rho2 + n*rho4)
    f3 = 2*(n*rho3 + b*rho4 + a*rho5)
    f4 = a**2 ; f5 = b**2; f6 = 2*a*b
    f7 = n**2 ; f8 = 2*b*n; f9 = 2*n*a 
    covmat_dxi = np.diag( (f1**2)*vara + (f2**2)*varb + (f3**2)*varn + + 2*(f1*f2*covab + f1*f3*covan + f2*f3*covbn) ) \
    + (f4**2)*(cov_rho0) + (f5**2)*(cov_rho1) + (f6**2)*(cov_rho2) + (f7**2)*(cov_rho3) +(f8**2)*(cov_rho4) + (f9**2)*(cov_rho5) 

    if(plots):
        plt.clf()
        pretty_rho(meanr, dxi, np.sqrt(np.diag(covmat_dxi)) , legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
        print('Printing',  dxiname)
        plt.savefig(dxiname, dpi=150)

    nrows = len(dxi)
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    covmathdu = fits.ImageHDU(covmat_dxi, name='COVMAT')
    hdul.insert(1, covmathdu)
    angarray = meanr
    valuearray =  np.array(dxi)
    bin1array = np.array([ -999]*nrows)
    bin2array = np.array([ -999]*nrows)
    angbinarray = np.arange(nrows)
    array_list = [bin1array, bin2array, angbinarray, valuearray,  angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='xi')
    hdul.insert(2, corrhdu)
    hdul.writeto(filename, overwrite=True)
    print(filename,'Written!')

def RUNtest(i_guess, data, nwalkers, nsteps, eq='All', mflags=[True, True, True], xip=True, xim=False, moderr=False, uwmprior=False, minimize=True ):
    from src.chi2 import minimizeCHI2
    from src.maxlikelihood import MCMC
    import numpy as np
    if (uwmprior or (not minimize)):
        iguess =  np.array(i_guess)
        chisq = np.inf
    if(minimize):
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            mflags=mflags, xip=xip,
                                            xim=xim, moderr=moderr)
        dof = len(data['rhosp'][0])
        print("reduced Chi2:", chisq/dof )
        print("Found parameters" , fitted_params)
        iguess = fitted_params
    
    samples, chains = MCMC(i_guess,data, nwalkers, nsteps, eq=eq,
                           mflags=mflags, moderr=moderr, uwmprior=uwmprior)
   
    return samples,  chains
    
def main():
    from src.readfits import read_rhos, read_taus
    from src.plot_stats import plotallrhosfits, plotalltausfits, plot_samplesdist
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

 
    if (args.plots):
        xlim = [0.1, 300.]
        ylims = [[1.e-11,5.e-6 ],[1.e-12,1.e-6 ],[3.e-9 ,3.e-4 ],[3.e-9 ,3.e-6 ]]
        rhostitle = ''
        plotallrhosfits(args.rhos, outpath=plotspath, title=rhostitle, xlim=xlim, ylims=ylims)
        plotalltausfits(args.taus, outpath=plotspath, xlim=xlim)
 
    meanr, rhos,  covrhos =  read_rhos(args.rhos, minscale=args.minscale, maxscale=args.maxscale)
    meanr, taus,  covtaus =  read_taus(args.taus, minscale=args.minscale, maxscale=args.maxscale)
    data = {}
    data['rhosp'] = [rhos[2*i] for i in range(6)]; data['rhosm'] = [rhos[2*i+1] for i in range(6)]
    data['covrhosp'] = [ covrhos[2*i] for i in range(6)]; data['covrhosm'] = [ covrhos[2*i + 1] for i in range(6)]
    data['tausp'] = [taus[2*i] for i in range(3)]; data['tausm'] = [taus[2*i + 1] for i in range(3)];
    data['covtausp'] = [covtaus[2*i] for i in range(3)]; data['covtausm'] = [covtaus[2*i + 1] for i in range(3)]; 
    

    #Finding best alpha beta gamma
    nwalkers,  nsteps = 100,  1000
    moderr = False
    minimize = True
    nsig = 1
    eq = args.eq
    print("Using equations: ", eq)
    i_guess0 = [ -0.01, 1, -1 ] #fiducial values

    if not (args.abe or args.ab or args.a): args.abe = True
    
    ## ALPHA-BETA-ETA
    if(args.abe):
        print("### Runing alpha, beta and eta test ### ")
        mflags = [True, True, True] ##alpha,beta,eta
        namemc = 'mcmc_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha-beta-eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_abe_' + str(eq) + '_.png'
        filename =  'abe_dxi.fits'            
    ## ALPHA-BETA
    if(args.ab):
        print("### Runing alpha and beta test ### ")
        mflags = [True, True, False] ##alpha,beta,eta
        namemc = 'mcmc_alpha-beta_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha-beta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha-beta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha-beta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_ab_' + str(eq) + '_.png'
        filename =  'ab_dxi.fits'
    ## ALPHA-ETA
    if(args.ae):
        print("### Runing alpha and eta test ### ")
        mflags = [True, False, True] ##alpha,eta,eta
        namemc = 'mcmc_alpha-eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha-eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha-eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha-eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_ae_' + str(eq) + '_.png'
        filename =  'ae_dxi.fits'
    ## BETA-ETA
    if(args.be):
        print("### Runing beta and eta test ### ")
        mflags = [True, False, True] ##beta,eta,eta
        namemc = 'mcmc_beta-eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_beta-eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_beta-eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_beta-eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_be_' + str(eq) + '_.png'
        filename =  'be_dxi.fits' 
    ## ALPHA
    if(args.a):
        print("### Runing alpha test ### ")
        mflags = [True, False, False] ##alpha,beta,eta
        namemc = 'mcmc_alpha_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_a_' + str(eq) + '_.png'
        filename =  'a_dxi.fits'
    ## Beta
    if(args.b):
        print("### Runing beta test ### ")
        mflags = [False, True, False] ##alpha,beta,eta
        namemc = 'mcmc_beta_eq_' + str(eq) + '_.png'
        namecont = 'contours_beta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_beta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_beta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_b_' + str(eq) + '_.png'
        filename =  'b_dxi.fits'
    ## Eta
    if(args.e):
        print("### Runing eta test ### ")
        mflags = [False, False, True] ##alpha,eta,eta
        namemc = 'mcmc_eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_e_' + str(eq) + '_.png'
        filename =  'e_dxi.fits'

    meanr, rhos,  covrhos =  read_rhos(args.rhoscosmo)
    datarhosp =  [meanr, [rhos[2*i] for i in range(6)],  [ covrhos[2*i] for i in range(6)]  ]
    datarhosm =  [meanr, [rhos[2*i+1] for i in range(6)],  [ covrhos[2*i + 1] for i in range(6)] ]

    i_guess = np.array(i_guess0)[np.array(mflags)].tolist()
    if args.samecont:
        samples, chains = RUNtest(i_guess, data, nwalkers, nsteps,
                                  eq=eq, mflags=mflags, xip=True ,
                                  xim=True , moderr=moderr,
                                  uwmprior=args.uwmprior, minimize=
                                  minimize)
        if(args.plots):
            plot_samplesdist(samples, chains, mflags, nwalkers, nsteps,  namemc, namecont )

        writexibias(samples, datarhosp, plots=args.plots,
                    nameterms=plotspath + 'p_' + nameterms,
                    dxiname=plotspath +'p_' + namedxip,
                    namecovmat=plotspath +'p_' + namecovmat,
                    filename=outpath + 'p_' + filename )
        writexibias(samples, datarhosm, plots=args.plots,
                    nameterms=plotspath +'m_' + nameterms,
                    dxiname=plotspath +'m_' + namedxip,
                    namecovmat=plotspath +'m_' + namecovmat,
                    filename=outpath + 'm_' + filename )
        

    else:
        samplesp, chainsp = RUNtest(i_guess, data, nwalkers, nsteps,
                                    eq=eq, mflags=mflags, xip=True ,
                                    xim=False , moderr=moderr,
                                    uwmprior=args.uwmprior, minimize=
                                    minimize)

        samplesm, chainsm = RUNtest(i_guess, data, nwalkers, nsteps,
                                    eq=eq, mflags=mflags, xip=False ,
                                    xim=True , moderr=moderr,
                                    uwmprior=args.uwmprior, minimize=
                                    minimize)

        if(args.plots):
            plot_samplesdist(samplesp, chainsp, mflags, nwalkers, nsteps, plotspath +'p_' + namemc,plotspath + 'p_' +namecont )
            plot_samplesdist(samplesm, chainsm, mflags, nwalkers, nsteps, plotspath +'m_' + namemc,plotspath + 'm_' +namecont )

        writexibias(samplesp, datarhosp, plots=args.plots,
                    nameterms=plotspath + 'p_' + nameterms,
                    dxiname=plotspath +'p_' + namedxip,
                    namecovmat=plotspath +'p_' + namecovmat,
                    filename=outpath + 'p_' + filename )
        writexibias(samplesm, datarhosm, plots=args.plots,
                    nameterms=plotspath +'m_' + nameterms,
                    dxiname=plotspath +'m_' + namedxip,
                    namecovmat=plotspath +'m_' + namecovmat,
                    filename=outpath + 'm_' + filename )
    
  
    
if __name__ == "__main__":
    main()

