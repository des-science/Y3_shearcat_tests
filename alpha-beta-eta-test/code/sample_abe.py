import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='save samples of abe during MCMC to get the full marginalized posteriors')
    parser.add_argument('--taus',
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tau_newrun_1.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tau_newrun_2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tau_newrun_3.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tau_newrun_4.fits'],
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_JK_zbin1.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_JK_zbin2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_JK_zbin3.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_JK_zbin4.fits'],
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK.fits'],
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_1.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_3.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_zbin_4.fits'],
                        nargs='+',
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    parser.add_argument('--rhos',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS_1-250.fits',
                        help='Fits file containing all rho stats used to estimate abe. This is rhos cosmo')
    parser.add_argument('--splitxipxim', default=False,
                        action='store_const', const=True,
                        help='Instead of use only one set of contamination parameter for both xip and xim, treat xip and xim independently')
    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--minscale', default=None, type=float, 
                        help='Limit the analysis to certain minimum scale')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Use Unweighted moments prior')
    parser.add_argument('--nsig', default=1, type=int, 
                        help='How many sigman for the marginalized confidence interval')
    parser.add_argument('--nsteps', default=1000, type=int, 
                        help='nsteps of MCMC')
    parser.add_argument('--nwalkers', default=100, type=int, 
                        help='nwalkers of MCMC')
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
    parser.add_argument('--plots', default=False,
                        action='store_const', const=True, help='Plot correlations functions')
    parser.add_argument('--outpath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/',
                        help='location of the output of the final contaminant')
    parser.add_argument('--filename',
                        default='samples_abe.fits',
                        help='all abe final samples, from the we get posteriors ')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')

    args = parser.parse_args()

    return args
                           
def getflagsnames(models_combo):
    eq,  abe, ab,  ae, be, a, b, e =  models_combo    
    ## ALPHA-BETA-ETA
    if(abe):
        print("### Runing alpha, beta and eta test ### ")
        mflags = [True, True, True]
        namemc = 'mcmc_alpha-beta-eta_eq_%i'%(eq)
        namecont = 'contours_alpha-beta-eta_eq_%i'%(eq)
        nameterms = 'termsdxi_alpha-beta-eta_eq_%i'%(eq)
        namecovmat = 'covmatrix_alpha-beta-eta_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_alpha-beta-eta_eq_%i'%(eq)
        namedxip = 'xibias_abe_eq_%i'%(eq)
        filename =  'abe_dxi_eq_%i'%(eq)          
    ## ALPHA-BETA
    if(ab):
        print("### Runing alpha and beta test ### ")
        mflags = [True, True, False] ##alpha,beta,eta
        namemc = 'mcmc_alpha-beta_eq_%i'%(eq)
        namecont = 'contours_alpha-beta_eq_%i'%(eq)
        nameterms = 'termsdxi_alpha-beta_eq_%i'%(eq)
        namecovmat = 'covmatrix_alpha-beta_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_alpha-beta_eq_%i'%(eq)
        namedxip = 'xibias_ab_eq%i'%(eq)
        filename =  'ab_dxi_eq_%i'%(eq)
    ## ALPHA-ETA
    if(ae):
        print("### Runing alpha and eta test ### ")
        mflags = [True, False, True]
        namemc = 'mcmc_alpha-eta_eq_%i'%(eq)
        namecont = 'contours_alpha-eta_eq_%i'%(eq)
        nameterms = 'termsdxi_alpha-eta_eq_%i'%(eq)
        namecovmat = 'covmatrix_alpha-eta_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_alpha-eta_eq_%i'%(eq)
        namedxip = 'xibias_ae_eq%i'%(eq)
        filename =  'ae_dxi_eq_%i'%(eq)
    ## BETA-ETA
    if(be):
        print("### Runing beta and eta test ### ")
        mflags = [False, True, True]
        namemc = 'mcmc_beta-eta_eq_%i'%(eq)
        namecont = 'contours_beta-eta_eq_%i'%(eq)
        nameterms = 'termsdxi_beta-eta_eq_%i'%(eq)
        namecovmat = 'covmatrix_beta-eta_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_beta-eta_eq_%i'%(eq)
        namedxip = 'xibias_be_eq%i'%(eq)
        filename =  'be_dxi_eq_%i'%(eq)
    ## ALPHA
    if(a):
        print("### Runing alpha test ### ")
        mflags = [True, False, False]
        namemc = 'mcmc_alpha_eq_%i'%(eq)
        namecont = 'contours_alpha_eq_%i'%(eq)
        nameterms = 'termsdxi_alpha_eq_%i'%(eq)
        namecovmat = 'covmatrix_alpha_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_alpha_eq_%i'%(eq)
        namedxip = 'xibias_a_eq%i'%(eq)
        filename =  'a_dxi_eq_%i'%(eq)
    ## Beta
    if(b):
        print("### Runing beta test ### ")
        mflags = [False, True, False] 
        namemc = 'mcmc_beta_eq_%i'%(eq)
        namecont = 'contours_beta_eq_%i'%(eq)
        nameterms = 'termsdxi_beta_eq_%i'%(eq)
        namecovmat = 'covmatrix_beta_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_beta_eq_%i'%(eq)
        namedxip = 'xibias_b_eq%i'%(eq)
        filename =  'b_dxi_eq_%i'%(eq)
    ## Eta
    if(e):
        print("### Runing eta test ### ")
        mflags = [False, False, True]
        namemc = 'mcmc_eta_eq_%i'%(eq)
        namecont = 'contours_eta_eq_%i'%(eq)
        nameterms = 'termsdxi_eta_eq_%i'%(eq)
        namecovmat = 'covmatrix_eta_eq_%i'%(eq)
        namebfres = 'Bestfitresiduals_eta_eq_%i'%(eq)
        namedxip = 'xibias_e_eq%i'%(eq)
        filename =  'e_dxi_eq_%i'%(eq)
    return [ mflags, namemc, namecont, namecovmat, namebfres, nameterms, namedxip, filename]


def RUNTEST(i_guess, data, nwalkers, nsteps, eq='All', mflags=[True, True, True], xip=True, xim=False, moderr=False, uwmprior=False, minimize=True, margin=True, overall=False ):
    from src.chi2 import minimizeCHI2
    from src.maxlikelihood import MCMC
    import numpy as np
    print("Getting best fit of marginalized likelihood")
    if (uwmprior or (not minimize)):
        iguess =  np.array(i_guess)
        chisq = np.inf
    if(minimize):
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            mflags=mflags, xip=xip,
                                            xim=xim, moderr=moderr)
        i_guess = fitted_params

    print("Used initial_guess", i_guess)
    samples, chains = MCMC(i_guess,data, nwalkers, nsteps, eq=eq,
                           mflags=mflags, xip=xip, xim=xim,
                           moderr=moderr, uwmprior=uwmprior)

    
    return samples, chains


def RUNTEST_PERTAU(rhofile, taufile, minscale, maxscale, models_combo, nwalkers, nsteps,  uwmprior, splitxipxim, margin, overall,  plots,  plotspath, zbin=0, axs=None):
    import numpy as np
    from src.readfits import  read_rhos, read_taus
    from src.plot_stats import plotallrhosfits, plotalltausfits, plot_samplesdist, plotbestfit, plotbestfitresiduals,  plotcovmat

    
    meanr, rhos,  covrho =  read_rhos(rhofile, minscale=minscale, maxscale=maxscale)
    meanr, taus,  covtau =  read_taus(taufile, minscale=minscale, maxscale=maxscale)
    data = {}
    data['rhos'] = rhos
    data['cov_rhos'] = covrho
    data['taus'] = taus
    data['cov_taus'] = covtau
    
    mflags, namemc, namecont, namecovmat, namebfres = getflagsnames(models_combo)[:5]

    #Finding best alpha beta gamma 
    moderr = False
    minimize = True
    eq = models_combo[0]
    print("Using equations: ", eq)
    print('Tomobin', zbin)
    i_guess0 = [ 0, 1, -1 ] #fiducial values
    i_guess = np.array(i_guess0)[np.array(mflags)].tolist()
    
    if splitxipxim:
        print('Using independet contamination parameter for xip and xim')
        print('RUNNING XIP')
        auxp1, auxp2 = RUNTEST(i_guess, data, nwalkers, nsteps, eq=eq,
                               mflags=mflags, xip=True ,
                               xim=False , moderr=moderr,
                               uwmprior=uwmprior, minimize=
                               minimize, margin=margin,
                               overall=overall)
        print('RUNNING XIM')
        auxm1, auxm2 = RUNTEST(i_guess, data, nwalkers, nsteps, eq=eq,
                               mflags=mflags, xip=False ,
                               xim=True , moderr=moderr,
                               uwmprior=uwmprior, minimize=
                               minimize, margin=margin,
                               overall=overall)
    else:
        print('Using same contamination parameter for xip and xim')
        aux1, aux2 = RUNTEST(i_guess, data, nwalkers, nsteps, eq=eq,
                             mflags=mflags, xip=True, xim=True ,
                             moderr=moderr, uwmprior=uwmprior,
                             minimize= minimize, margin=margin,
                             overall=overall)
        
        auxp1 = auxm1 =  aux1; auxp2 = auxm2 = aux2                
        
    #samplesp, samplesm
    return auxp1, auxm1 




    
                        
def main():
    from astropy.io import fits
    from src.maxlikelihood import percentiles, bestparameters
    from src.readfits import  read_rhos, read_taus
    from src.chi2 import chi2nu
    import numpy as np
    import itertools

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
        
    if not (args.abe or args.ab or args.ae or args.be or args.a or args.b or args.e): args.abe = True
    models_combo = [args.eq, args.abe, args.ab, args.ae, args.be, args.a, args.b, args.e]
    eq = args.eq
    nsig = args.nsig
    nwalkers,  nsteps = args.nwalkers, args.nsteps

    mflags = getflagsnames(models_combo)[0]
    ndim =  len(list(itertools.compress(range(len(mflags)),  mflags)))

    ##Format of the fit file output
    nbins = len(args.taus)

    del_flag = [ not a for a in mflags] + [False] #all tables print chi2
    names = np.concatenate(np.array([ np.delete(['a%i'%(i),'b%i'%(i),'e%i'%(i),'chi2nu%i'%i], del_flag) for i in range(1,nbins + 1)]).T)

    forms = ['f4']*len(names) 
    dtype = dict(names = names, formats=forms)
    print("LOOK HERR \n",  names)

    print("STARTING TOMOGRAPHIC ANALYSIS")
    data = {}
    data['rhos'] = read_rhos(args.rhos, minscale=args.minscale, maxscale=args.maxscale)[1]
    data['cov_rhos'] = read_rhos(args.rhos, minscale=args.minscale, maxscale=args.maxscale)[2]

    aplist = []; bplist = []; eplist = []; chi2pnulist = [];
    amlist = []; bmlist = []; emlist = []; chi2mnulist = [];

    
    for i,  taufile in enumerate(args.taus):
        samplesp, samplesm=RUNTEST_PERTAU(args.rhos,taufile,args.minscale, args.maxscale,
                                                  models_combo ,nwalkers,nsteps, args.uwmprior, args.splitxipxim,
                                                  True, False, args.plots, plotspath,  zbin=(i + 1))
        
        data['taus'] =  read_taus(taufile, minscale=args.minscale, maxscale=args.maxscale)[1]
        data['cov_taus'] =  read_taus(taufile, minscale=args.minscale, maxscale=args.maxscale)[2]
        if args.splitxipxim:
            chi2pnulist.append(np.array([ chi2nu(pars, data, eq=eq , mflags=mflags, xip=True, xim=False,  moderr=False) for pars in samplesp.T ]))
            chi2mnulist.append(np.array([ chi2nu(pars, data, eq=eq , mflags=mflags, xip=False, xim=True,  moderr=False) for pars in samplesm.T ]))
        else:
            chi2pnulist.append(np.array( [ chi2nu(pars, data, eq=eq , mflags=mflags, xip=True, xim=True,  moderr=False) for pars in samplesp.T ]))
            chi2mnulist.append(np.array( [ chi2nu(pars, data, eq=eq , mflags=mflags, xip=True, xim=True,  moderr=False) for pars in samplesm.T ]))
        
        if (ndim == 3):
            apl, bpl, epl = samplesp
            aml, bml, eml = samplesm
            aplist.append(apl); bplist.append(bpl); eplist.append(epl)
            amlist.append(aml); bmlist.append(bml); emlist.append(eml)
            print("All this numbers should be the same", len(apl), len(bpl), len(epl) )
            print("All this numbers should be the same", len(aml), len(bml), len(eml) )

        if (ndim == 2):
            apl, bpl = samplesp
            aml, bml = samplesm
            aplist.append(apl); bplist.append(bpl)
            amlist.append(aml); bmlist.append(bml)

            print("All this numbers should be the same", len(apl), len(bpl) )
            print("All this numbers should be the same", len(aml), len(bml) )

        if (ndim == 1):
            apl = samplesp
            aml = samplesm
            aplist.append(apl)
            amlist.append(aml)

        
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
        
    nrows = len(aplist[0])
    
    outdata1 = np.recarray((nrows, ), dtype=dtype)
    if (ndim == 3): array_list = aplist + bplist + eplist + chi2pnulist
    if (ndim == 2): array_list = aplist + bplist + chi2pnulist
    if (ndim == 1): array_list = aplist + chi2pnulist
    for array, name in zip(array_list, names): outdata1[name] = array
    table = fits.BinTableHDU(outdata1, name='abe xip' )
    hdul.insert(1, table)
    
    outdata2 = np.recarray((nrows, ), dtype=dtype)
    if (ndim == 3): array_list = amlist + bmlist + emlist + chi2mnulist
    if (ndim == 2): array_list = amlist + bmlist + chi2mnulist
    if (ndim == 1): array_list = amlist + chi2mnulist
    for array, name in zip(array_list, names): outdata2[name] = array 
    table = fits.BinTableHDU(outdata2, name='abe xim' )
    hdul.insert(2, table)

           
    filename = os.path.join(outpath, args.filename)
    
    print("Printing file:", filename)
    hdul.writeto(filename, overwrite=True)

    

  
    
if __name__ == "__main__":
    main()



        
