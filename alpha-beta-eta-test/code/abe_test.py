import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the fitting problem of system ofequatiosn, plotting correlations and final correlation function withbias')
    parser.add_argument('--taus',
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK.fits'],
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    parser.add_argument('--singletau',
                        default=None, 
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS.fits',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK.fits',
                        #default='/home2/dfa/sobreira/alsina/catalogs/flask/taus/taus_src-cat_s280_z1_ck1.fits',
                        help='Fits file containing all tau stats used to estimate abe')
    parser.add_argument('--rhos', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--rhoscosmo', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS_Y3.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--samecontaminant', default=False,
                        action='store_const', const=True,
                        help='Include xim correlation functions in the vectors. Then use only one set of contamination parameter for both xip and xim ')
    parser.add_argument('--maxscale', default=250, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--minscale', default=2.5, type=float, 
                        help='Limit the analysis to certain minimum scale')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Use Unweighted moments prior')
    parser.add_argument('--nsig', default=1, type=int, 
                        help='How many sigmans for the marginalized confidence interval')
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
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    
    args = parser.parse_args()

    return args
                        
def corrmatrix(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    return corr
def plotcovpars(samples,namecovmat='covmat_pars.png'):
    import numpy as np
    #plot covariance matrix of parameters alpha, beta and eta.
    par_matcov = np.cov(samples)
    corr=corrmatrix(par_matcov)
    cov_vmin=np.min(corr)
    plt.clf()
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin,
               vmax=1.)
    plt.colorbar()
    plt.title(r'$\alpha \mid \beta \mid \eta $')
    plt.tight_layout()
    print('Printing', namecovmat)
    plt.savefig(namecovmat, dpi=500)
def getxibias_margin(samples, datarhos,  nsig=1, plots=False, nameterms='terms_dxi.png',dxiname='dxi.png'):
    from src.readfits import read_rhos
    from src.maxlikelihood import bestparameters
    from src.plot_stats import pretty_rho
    import numpy as np

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

    return meanr, dxi, covmat_dxi

def getxibiastomo_margin(samples_i, samples_j, datarhos, models_combo, nsig=1,  plots=False, nameterms='terms_dxi.png',dxiname='dxi.png'):
    from src.readfits import read_rhos
    from src.maxlikelihood import bestparameters
    from src.plot_stats import pretty_rho
    import numpy as np
    a_i = b_i = e_i = a_j = b_j = e_j =  0; 
    vara_i =  varb_i =  vare_i = 0; 
    vara_j =  varb_j =  vare_j = 0;
    cova_ib_i = cova_ie_i = covb_ie_i = 0
    cova_jb_j = cova_je_j = covb_je_j = 0
    cova_ia_j = cova_ib_j = cova_ie_j = 0
    cova_jb_i = cova_je_i = 0
    covb_ib_j = covb_ie_j = 0
    covb_je_i = cove_ie_j =0

    pars_i = bestparameters(samples_i)
    pars_j = bestparameters(samples_j)
    samples =  np.concatenate([samples_i, samples_j])
    
    par_matcov = np.cov(samples) 
    if (par_matcov.size==1 ): variances = par_matcov
    else: variances = np.diagonal(par_matcov)
    covariances = sum( (par_matcov[i,i+1: ].tolist() for i in range(len(samples) - 1)) , [] )

    abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
    if not (abe_bool or ab_bool or a_bool): abe = True
    if(abe_bool):
        a_i, b_i, e_i = pars_i
        a_j, b_j, e_j = pars_j
        vara_i ,  varb_i ,   vare_i ,  vara_j ,   varb_j ,   vare_j = variances
        cova_ib_i , cova_ie_i , cova_ia_j , cova_ib_j , cova_ie_j , covb_ie_i , cova_jb_i , covb_ib_j , covb_ie_j , cova_je_i , covb_je_i , cove_ie_j , cova_jb_j , cova_je_j , covb_je_j =  covariances     
    if(ab_bool):
        a_i, b_i = pars_i
        a_j, b_j = pars_j
        vara_i ,  varb_i , vara_j ,  varb_j =  variances
        cova_ib_i , cova_ia_j , cova_ib_j , cova_jb_i , covb_ib_j ,cova_jb_j =   covariances
    if(ae_bool):
        a_i, e_i = pars_i
        a_j, e_j = pars_j
        vara_i ,  vare_i , vara_j ,  vare_j = variances
        cova_ie_i , cova_ia_j , cova_ie_j , cova_je_i , cove_ie_j , cova_je_j =  covariances
    if(be_bool):
        b_i, e_i = pars_i
        b_j, e_j = pars_j
        varb_i ,  vare_i ,  varb_j ,  vare_j = variances
        covb_ie_i , covb_ib_j , covb_ie_j , covb_je_i , cove_ie_j , covb_je_j =  covariances
    if(a_bool):
        a_i= pars_i
        a_j= pars_j
        vara_i ,  vara_j  = variances
        cova_ia_j = covariances[0]
    if(b_bool):
        b_i = pars_i
        b_j = pars_j
        varb_i ,  varb_j  = variances
        covb_ib_j = covariances[0]
    if(e_bool):
        e_i = pars_i
        e_j = pars_j
        vare_i ,  vare_j  , variances
        cove_ie_j = covariances[0]

    meanr, rhos, covrhos =  datarhos
    rho0, rho1, rho2, rho3, rho4, rho5 = rhos
    cov_rho0, cov_rho1, cov_rho2, cov_rho3, cov_rho4, cov_rho5 =  covrhos

    #supposing that a,b and n are idependent of rhos(scale independent)
    dxi = a_i*a_j*rho0 + b_i*b_j*rho1 + e_i*e_j*rho3 + (b_i*a_j + b_j*a_i)*rho2 + (b_i*e_j + b_j*e_i)*rho4 + (e_i*a_j +e_j*a_i)*rho5
    f1_j = a_i*rho0 + a_i*rho2 +a_i*rho5; f1_i = a_j*rho0 + a_j*rho2 +a_j*rho5 #partiala
    f2_j = b_i*rho1 + a_i*rho2 +e_i*rho4; f2_i = b_j*rho1 + a_j*rho2 +e_j*rho4 #partialb
    f3_j = e_i*rho3 + b_i*rho4 + a_i*rho5; f3_i = e_j*rho3 + b_j*rho4 + a_j*rho5 #partiale
    f4 = a_i*a_j ; f5 = b_i*b_j; f6 =b_i*a_j + b_j*a_i
    f7 = e_i*e_j ; f8 = (b_i*e_j + b_j*e_i); f9 = (e_i*a_j +e_j*a_i)
    covmat_dxi =np.diag( (f1_i**2)*vara_i + (f1_j**2)*vara_j +(f2_i**2)*varb_i + (f2_j**2)*varb_j + (f3_j**2)*vare_j+(f3_i**2)*vare_i +
                         2*(f1_i*f1_j*cova_ia_j + f1_i*f2_i*cova_ib_i + f1_i*f2_j*cova_ib_j + f1_i*f3_i*cova_ie_i + f1_i*f3_j*cova_ie_j +
                                                  f1_j*f2_i*cova_jb_i + f1_j*f2_j*cova_jb_j + f1_j*f3_i*cova_je_i + f1_j*f3_j*cova_je_j +
                                                                        f2_i*f2_j*covb_ib_j + f2_i*f3_i*covb_ie_i + f2_i*f3_j*covb_ie_j +
                                                                                              f2_j*f3_i*covb_je_i + f2_j*f3_j*covb_je_j +
                                                                                                                    f3_i*f3_j*cove_ie_j) )+ (f4**2)*(cov_rho0) + (f5**2)*(cov_rho1) + (f6**2)*(cov_rho2) +(f7**2)*(cov_rho3) +(f8**2)*(cov_rho4) + (f9**2)*(cov_rho5) 


    if(plots):
        plt.clf()
        pretty_rho(meanr, dxi, np.sqrt(np.diag(covmat_dxi)) , legend=r"$\delta \xi$",  ylabel=r"$\delta \xi$",  xlim=xlim)
        print('Printing',  dxiname)
        plt.savefig(dxiname, dpi=150)
    return meanr, dxi, covmat_dxi
 
def getxibiastomo_overall(pars_i, pars_j , datarhos, models_combo, plots=False, nameterms='terms_dxi.png',dxiname='dxi.png'):
    from src.readfits import read_rhos
    from src.maxlikelihood import bestparameters
    from src.plot_stats import pretty_rho
    import numpy as np
    a_i = b_i = e_i = a_j = b_j = e_j =  0; 
    abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    meanr, rhos, covrhos =  datarhos
    rho0, rho1, rho2, rho3, rho4, rho5 = rhos
    
    if not (abe_bool or ab_bool or a_bool): abe = True
    if(abe_bool):
        a_i, b_i, e_i = pars_i
        a_j, b_j, e_j = pars_j
    if(ab_bool):
        a_i, b_i = pars_i
        a_j, b_j = pars_j    
    if(ae_bool):
        a_i, e_i = pars_i
        a_j, e_j = pars_j
    if(be_bool):
        b_i, e_i = pars_i
        b_j, e_j = pars_j
    if(a_bool):
        a_i= pars_i
        a_j= pars_j
    if(b_bool):
        b_i = pars_i
        b_j = pars_j
    if(e_bool):
        e_i = pars_i
        e_j = pars_j
        
    dxi = a_i*a_j*rho0 + b_i*b_j*rho1 + e_i*e_j*rho3 + (b_i*a_j + b_j*a_i)*rho2 + (b_i*e_j + b_j*e_i)*rho4 + (e_i*a_j +e_j*a_i)*rho5
    
    return meanr, dxi
 
def getmflags(models_combo):
    abe, ab,  ae, be, a, b, e =  models_combo    
    if not (abe or ab or a): abe = True
    if(abe): mflags = [True, True, True] ##alpha,beta,eta
    if(ab): mflags = [True, True, False] ##alpha,beta,eta
    if(ae): mflags = [True, False, True] ##alpha,eta,eta
    if(be): mflags = [True, False, True] ##beta,eta,eta
    if(a): mflags = [True, False, False] ##alpha,beta,eta
    if(b): mflags = [False, True, False] ##alpha,beta,eta
    if(e): mflags = [False, False, True] ##alpha,eta,eta
    return mflags
def getnames(models_combo, eq):
    abe, ab,  ae, be, a, b, e =  models_combo    
    if not (abe or ab or a): abe = True
    ## ALPHA-BETA-ETA
    if(abe):
        print("### Runing alpha, beta and eta test ### ")
        namemc = 'mcmc_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha-beta-eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_abe_' + str(eq) + '_.png'
        filename =  'abe_dxi.fits'            
    ## ALPHA-BETA
    if(ab):
        print("### Runing alpha and beta test ### ")
        mflags = [True, True, False] ##alpha,beta,eta
        namemc = 'mcmc_alpha-beta_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha-beta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha-beta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha-beta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_ab_' + str(eq) + '_.png'
        filename =  'ab_dxi.fits'
    ## ALPHA-ETA
    if(ae):
        print("### Runing alpha and eta test ### ")
        namemc = 'mcmc_alpha-eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha-eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha-eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha-eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_ae_' + str(eq) + '_.png'
        filename =  'ae_dxi.fits'
    ## BETA-ETA
    if(be):
        print("### Runing beta and eta test ### ")
        namemc = 'mcmc_beta-eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_beta-eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_beta-eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_beta-eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_be_' + str(eq) + '_.png'
        filename =  'be_dxi.fits' 
    ## ALPHA
    if(a):
        print("### Runing alpha test ### ")
        namemc = 'mcmc_alpha_eq_' + str(eq) + '_.png'
        namecont = 'contours_alpha_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_alpha_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_alpha_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_a_' + str(eq) + '_.png'
        filename =  'a_dxi.fits'
    ## Beta
    if(b):
        print("### Runing beta test ### ")
        namemc = 'mcmc_beta_eq_' + str(eq) + '_.png'
        namecont = 'contours_beta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_beta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_beta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_b_' + str(eq) + '_.png'
        filename =  'b_dxi.fits'
    ## Eta
    if(e):
        print("### Runing eta test ### ")
        namemc = 'mcmc_eta_eq_' + str(eq) + '_.png'
        namecont = 'contours_eta_eq_' + str(eq) + '_.png'
        nameterms = 'termsdxi_eta_eq_' + str(eq) + '_.png'
        namecovmat = 'covmatrix_eta_eq_' + str(eq) + '_.png'
        namedxip = 'xibias_e_' + str(eq) + '_.png'
        filename =  'e_dxi.fits'
    return [ namemc, namecont, nameterms, namecovmat, namedxip, filename]

def RUNtest(i_guess, data, nwalkers, nsteps, eq='All', mflags=[True, True, True], xip=True, xim=False, moderr=False, uwmprior=False, minimize=True ):
    from src.chi2 import minimizeCHI2
    from src.maxlikelihood import MCMC,  percentiles
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
                           mflags=mflags, xip=xip, xim=xim,
                           moderr=moderr, uwmprior=uwmprior)

    mcmcpars = percentiles(samples, nsig=1) 
    print( ' mcmc parameters xi',  'nsig=', 1, ' percentiles: ',  mcmcpars)
   
    return fitted_params, samples,  chains
    
def write_singlexip_margin( samplesp, samplesm, rhoscosmo,  models_combo,  eq, nsig, plots,  outpath,  plotspath):
    import numpy as np
    from src.readfits import read_rhos
    from astropy.io import fits
    mflags =  getmflags(models_combo)
    namemc, namecont, nameterms, namecovmat, namedxip, filename = getnames(models_combo, eq)
    meanr, rhos,  covrhos =  read_rhos(rhoscosmo)
    datarhosp =  [meanr, [rhos[2*i] for i in range(6)],  [ covrhos[2*i] for i in range(6)]  ]
    datarhosm =  [meanr, [rhos[2*i+1] for i in range(6)],  [ covrhos[2*i + 1] for i in range(6)] ]

    meanr, xip, covxip = getxibias_margin(samplesp, datarhosp, nsig=
                                          nsig, plots=plots,
                                          nameterms=plotspath + 'p_' +
                                          nameterms, dxiname=plotspath +'p_'
                                          + namedxip)
    meanr, xim, covxim = getxibias_margin(samplesm, datarhosm,
                                          nsig=nsig, plots=plots,
                                          nameterms=plotspath +'m_' +
                                          nameterms, dxiname=plotspath +'m_'
                                          + namedxip)
    
    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = len(meanr)
    outdata = np.recarray((nrows, ), dtype=dtype)
    
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])

    hdul.insert(1, fits.ImageHDU(covxip, name='COVMAT_XIP'))
    hdul.insert(2, fits.ImageHDU(covxip, name='COVMAT_XIM'))
    
    bin1array = np.array([ -999]*nrows)
    bin2array = np.array([ -999]*nrows)
    angbinarray = np.arange(nrows)
    array_list = [bin1array, bin2array, angbinarray, xip,  meanr ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='xip')
    hdul.insert(3, corrhdu)
    array_list = [bin1array, bin2array, angbinarray, xim,  meanr ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='xim')
    hdul.insert(4, corrhdu)
    hdul.writeto(filename, overwrite=True)
    print(filename,'Written!')
def write_tomoxip_margin(samplesp_list, samplesm_list, rhoscosmo, models_combo,  eq, nsig,  plots,  outpath , plotspath):
    import itertools
    from src.readfits import read_rhos
    from astropy.io import fits
    import numpy as np
    mflags = getmflags(models_combo)
    namemc, namecont, nameterms, namecovmat, namedxip, filename = getnames(models_combo, eq)
    meanr, rhos,  covrhos =  read_rhos(rhoscosmo)
    datarhosp =  [meanr, [rhos[2*i] for i in range(6)],  [ covrhos[2*i] for i in range(6)]  ]
    datarhosm =  [meanr, [rhos[2*i+1] for i in range(6)],  [ covrhos[2*i + 1] for i in range(6)] ]

    
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])

    xip_list = []; xim_list = []; bin1_list = []; bin2_list = []; angbin_list = []; ang_list = []
    covxip_list = []; covxim_list = []
    
    nbins=4
    a=[i for i in range(1,nbins+1)]
    #b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2):
    #for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        print(i, j)
        samplesp_i =  samplesp_list[i - 1]; samplesm_i =  samplesm_list[i - 1]
        samplesp_j =  samplesp_list[j - 1]; samplesm_j =  samplesm_list[j - 1]
        binstr = str(i) +'_' + str(j) + '_p_'
        meanr, xip, covxip = getxibiastomo_margin( samplesp_i, samplesp_j, datarhosp, models_combo, nsig=nsig, plots=plots, nameterms=plotspath+binstr + nameterms,
                                                    dxiname=plotspath + binstr+ namedxip )
        meanr, xim, covxim = getxibiastomo_margin( samplesm_i, samplesm_j ,datarhosm, models_combo, nsig=nsig, plots=plots,nameterms=plotspath+binstr + nameterms,
                                                    dxiname=plotspath + binstr+ namedxip )
        ang_list.append(meanr)
        bin1_list.append(np.array( [i]*len(meanr)))
        bin2_list.append(np.array( [j]*len(meanr)))
        angbin_list.append(np.arange(len(meanr)))
        xip_list.append(xip)
        xim_list.append(xim)
        covxip_list.append(np.diag(covxip))
        covxim_list.append(np.diag(covxim))

    covmat = np.diag(np.concatenate([np.concatenate(covxip_list), np.concatenate(covxim_list)]))
    hdul.insert(1, fits.ImageHDU(covmat, name='COVMAT'))
    
    bin1array = np.concatenate(bin1_list)
    bin2array = np.concatenate(bin2_list)
    angbinarray = np.concatenate(angbin_list)
    valuearray = np.concatenate(xip_list)
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
    valuearray = np.concatenate(xim_list)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xim')
    hdul.insert(3, corrhdu)
    hdul.writeto(outpath +'marg_' + filename, overwrite=True)
    print(outpath + 'marg_' + filename,'Written!')
        
def write_tomoxip_overall(parsp_list, parsm_list, rhoscosmo, models_combo,  eq, plots,  outpath , plotspath):
    import itertools
    from src.readfits import read_rhos
    from astropy.io import fits
    import numpy as np
    mflags = getmflags(models_combo)
    namemc, namecont, nameterms, namecovmat, namedxip, filename = getnames(models_combo, eq)
    meanr, rhos,  covrhos =  read_rhos(rhoscosmo)
    datarhosp =  [meanr, [rhos[2*i] for i in range(6)],  [ covrhos[2*i] for i in range(6)]  ]
    datarhosm =  [meanr, [rhos[2*i+1] for i in range(6)],  [ covrhos[2*i + 1] for i in range(6)] ]

    
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    xip_list = []; xim_list = []; bin1_list = []; bin2_list = []; angbin_list = []; ang_list = []
    covxip_list = []; covxim_list = []
    
    nbins=4
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        print(i, j)
        parsp_i =  parsp_list[i - 1]; parsm_i =  parsm_list[i - 1]
        parsp_j =  parsp_list[j - 1]; parsm_j =  parsm_list[j - 1]
        binstr = str(i) +'_' + str(j) + '_p_'
        meanr, xip = getxibiastomo_overall( parsp_i, parsp_j, datarhosp, models_combo,plots=plots, nameterms=plotspath+binstr + nameterms, dxiname=plotspath + binstr+ namedxip)
        meanr, xim = getxibiastomo_overall( parsm_i, parsm_j ,datarhosm, models_combo,plots=plots,nameterms=plotspath+binstr + nameterms, dxiname=plotspath + binstr+ namedxip )
        ang_list.append(meanr)
        bin1_list.append(np.array( [i]*len(meanr)))
        bin2_list.append(np.array( [j]*len(meanr)))
        angbin_list.append(np.arange(len(meanr)))
        xip_list.append(xip)
        xim_list.append(xim)

     
    bin1array = np.concatenate(bin1_list)
    bin2array = np.concatenate(bin2_list)
    angbinarray = np.concatenate(angbin_list)
    valuearray = np.concatenate(xip_list)
    angarray = np.concatenate(ang_list)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f4',  'f4']
    dtype = dict(names = names, formats=forms)
    nrows = len(angarray)
    outdata = np.recarray((nrows, ), dtype=dtype)
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xip')
    hdul.insert(1, corrhdu)
    valuearray = np.concatenate(xim_list)
    array_list = [bin1array, bin2array, angbinarray, valuearray, angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='delta_xim')
    hdul.insert(2, corrhdu)
    hdul.writeto(outpath + filename, overwrite=True)
    print(outpath + filename,'Written!')
    
def RUNTEST_PERTAU(rhofile, taufile, minscale, maxscale, models_combo, eq, nsig,  uwmprior, samecontaminant,  plots,  outpath,  plotspath, zbin=''):
    import numpy as np
    from src.readfits import read_rhos, read_taus
    from src.plot_stats import plotallrhosfits, plotalltausfits, plot_samplesdist

    if (plots):
        xlim = [0.1, 300.]
        taustitle = ''
        plotalltausfits(taufile, outpath=plotspath, title='zbin:' + zbin,  xlim=xlim, zbin=str(zbin))
 
    meanr, rhos,  covrhos =  read_rhos(rhofile, minscale=minscale, maxscale=maxscale)
    meanr, taus,  covtaus =  read_taus(taufile, minscale=minscale, maxscale=maxscale)
    data = {}
    data['rhosp'] = [rhos[2*i] for i in range(6)]; data['rhosm'] = [rhos[2*i+1] for i in range(6)]
    data['covrhosp'] = [ covrhos[2*i] for i in range(6)]; data['covrhosm'] = [ covrhos[2*i + 1] for i in range(6)]
    data['tausp'] = [taus[2*i] for i in range(3)]; data['tausm'] = [taus[2*i + 1] for i in range(3)];
    data['covtausp'] = [covtaus[2*i] for i in range(3)]; data['covtausm'] = [covtaus[2*i + 1] for i in range(3)]; 

    mflags = getmflags(models_combo)
    namemc, namecont, nameterms, namecovmat, namedxip, filename = getnames(models_combo, eq)

    #Finding best alpha beta gamma
    nwalkers,  nsteps = 100,  1000
    moderr = False
    minimize = True
    nsig = nsig
    eq = eq
    print("Using equations: ", eq)
    print('RUNNING tomobin', zbin)
    i_guess0 = [ -0.01, 1, -1 ] #fiducial values
    i_guess = np.array(i_guess0)[np.array(mflags)].tolist()
    if samecontaminant:
        print('Using same contamination parameter for xip and xim')
        fitted_params ,samples, chains = RUNtest(i_guess, data, nwalkers, nsteps,
                                  eq=eq, mflags=mflags, xip=True, xim=True , moderr=moderr,
                                  uwmprior=uwmprior, minimize= minimize)
        
        fitted_paramsp, samplesp, chainsp =  fitted_params, samples, chains
        fitted_paramsm, samplesm, chainsm =  fitted_params, samples, chains
    else:
        print('Using independet contamination parameter for xip and xim')
        print('RUNNING XIP')
        fitted_paramsp, samplesp, chainsp = RUNtest(i_guess, data,
                                                    nwalkers, nsteps, eq=eq,
                                                    mflags=mflags, xip=True ,
                                                    xim=False , moderr=moderr,
                                                    uwmprior=uwmprior, minimize=
                                                    minimize)
        print('RUNNING XIM')
        fitted_paramsm, samplesm, chainsm = RUNtest(i_guess, data,
                                                    nwalkers, nsteps, eq=eq,
                                                    mflags=mflags, xip=False ,
                                                    xim=True , moderr=moderr,
                                                    uwmprior=uwmprior, minimize=
                                                    minimize)
        
    if(plots):
        plot_samplesdist(samplesp, chainsp, mflags, nwalkers, nsteps, plotspath +'p_' + namemc,plotspath + 'p_' +namecont )
        plot_samplesdist(samplesm, chainsm, mflags, nwalkers, nsteps, plotspath +'m_' + namemc,plotspath + 'm_' +namecont )
        plotcovpars(samplesp, namecovmat=plotspath + 'p' + namecovmat)
        plotcovpars(samplesm, namecovmat=plotspath + 'm' + namecovmat)
    return fitted_paramsp, samplesp, chainsp, fitted_paramsm, samplesm, chainsm

                        
def main():
    from src.plot_stats import plotallrhosfits
    
    
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

    models_combo = [args.abe, args.ab, args.ae, args.be, args.a, args.b, args.e]
    eq = args.eq
    nsig = args.nsig

    if args.singletau is not None:
        taufile=args.singletau
        parsp, samplesp, chainsp, parsm, samplesm, chainsm =RUNTEST_PERTAU(args.rhos, taufile, args.minscale,
                                                             args.maxscale, models_combo,eq,  nsig, args.uwmprior,
                                                             args.samecontaminant, args.plots, outpath,
                                                             plotspath)
        mcmcpars = percentiles(samplesp, nsig=nsig) 
        print( ' mcmc parameters xi+',  'nsig=', nsig, ' percentiles: ',  mcmcpars)
        mcmcpars = percentiles(samplesm, nsig=nsig) 
        print( ' mcmc parameters xi-',  'nsig=', nsig, ' percentiles: ',  mcmcpars)
        write_singlexip_margin(samplesp, samplesm, args.rhoscosmo,  models_combo, eq,  nsig, args.plots,  outpath,  plotspath)
        #write_singlexip_overall(parsp, parsm, args.rhoscosmo,  models_combo, eq,  nsig, args.plots,  outpath,  plotspath)
        
    else:
        samplesp_list = []; parsp_list = []
        samplesm_list = []; parsm_list = []
        for i,  taufile in enumerate(args.taus):
            parsp, samplesp, chainsp, parsm, samplesm, chainsm =  RUNTEST_PERTAU(args.rhos, taufile, args.minscale,
                                                                   args.maxscale, models_combo, eq, nsig, args.uwmprior,
                                                                   args.samecontaminant, args.plots, outpath,
                                                                   plotspath,  zbin=str(i + 1))
            
    
            samplesp_list.append(samplesp); parsp_list.append(parsp)
            samplesm_list.append(samplesm); parsm_list.append(parsm)
        write_tomoxip_margin( samplesp_list, samplesm_list, args.rhoscosmo,  models_combo, eq, nsig, args.plots, outpath, plotspath)
        #write_tomoxip_overall( parsp_list, parsm_list, args.rhoscosmo,  models_combo, eq, args.plots, outpath, plotspath)

    

  
    
if __name__ == "__main__":
    main()
