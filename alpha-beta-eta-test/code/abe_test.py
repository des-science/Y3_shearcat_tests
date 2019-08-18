import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta eta test solving the fitting problem of system ofequatiosn, plotting correlations and final correlation function withbias')
    parser.add_argument('--taus',
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                                  '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        #default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1_v2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2_v2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3_v2.fits',
                        #         '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4_v2.fits'],
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    parser.add_argument('--singletau',
                        default=None, 
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3_v2.fits',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/oldbins/TAUS.fits',
                        #default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/oldbins/TAUS_FLASK.fits',
                        #default='/home2/dfa/sobreira/alsina/catalogs/flask/taus/taus_src-cat_s201_z1_ck1.fits',
                        #default=None,
                        help='Fits file containing all tau stats used to estimate abe')
    parser.add_argument('--rhos', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--rhoscosmo', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS_Y3.fits',
                        help='Fits file containing all rho stats used to estimate dxip, the contaminant to be used in cosmosis')
    parser.add_argument('--splitxipxim', default=False,
                        action='store_const', const=True,
                        help='Instead of use only one set of contamination parameter for both xip and xim, treat xip and axim independently')
    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--minscale', default=None, type=float, 
                        help='Limit the analysis to certain minimum scale')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Use Unweighted moments prior')
    parser.add_argument('--overall', default=False, 
                        action='store_const', const=True, help='Use the values of the maximum of the likelihood function, then no covariance matrix will be written in dxi')
    parser.add_argument('--margin', default=False, 
                        action='store_const', const=True, help='Use the values of the best fits of the marginalized likelihood function')
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
def getxibias_margin(samples, datarhos, models_combo, plots=False, nameterms='terms_dxi.png',dxiname='dxi.png'):
    from src.readfits import read_rhos
    from src.maxlikelihood import bestparameters
    from src.plot_stats import pretty_rho
    import numpy as np

    a = b = e = 0; vara =  varb =  vare = 0; covab = covae = covbe = 0
    bestpar = bestparameters(samples)
    #print("Best pars", bestpar)
    par_matcov = np.cov(samples) 
    if (par_matcov.size==1 ): variances = par_matcov
    else: variances = np.diagonal(par_matcov)
    covariances = sum( (par_matcov[i,i+1: ].tolist() for i in range(len(samples) - 1)) , [] )

    eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
    if not (abe_bool or ab_bool or a_bool): abe_bool = True
    if(abe_bool):
        a, b, e = bestpar
        vara, varb, vare =  variances
        covab, covae, covbe =  covariances   
    if(ab_bool):
        a, b = bestpar
        vara, varb =  variances
        covab =  covariances[0]
    if(ae_bool):
        a, e = bestpar
        vara, vare =  variances
        covae =  covariances[0]
    if(be_bool):
        b, e = bestpar
        varb, vare =  variances
        covbe =  covariances[0]
    if(a_bool):
        a =  bestpar[0]
        vara =  variances
    if(b_bool):
        b =  bestpar[0]
        varb =  variances
    if(e_bool):
        e =  bestpar[0]
        vare =  variances

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
        var2 = ((2*e*rho3)**2)*vare +  (n**2)*(sig_rho3**2)
        varab =  vara*(b**2) + varb*(a**2) + 2*covab*(a*b)
        #varab = ((a*b)**2)*( (vara/((a)**2)) + (varb/((b)**2)) + 2*covab/(a*b) )
        var3 = 4*( (rho2**2)*varab + (sig_rho2**2)*((a*b)**2)  )
        #var3 = 4*((a*b*rho2p)**2)*( varab/((a*b)**2) + (sig_rho2/rho2p)**2 )
        varbe =  vare*(b**2) + varb*(e**2) + 2*covbe*(b*e)
        #varbn = ((n*b)**2)*( (varn/((n)**2)) + (varb/((b)**2)) + 2*covbn/(b*n) ) 
        var4 = 4*( (rho4**2)*varbe + (sig_rho4**2)*((e*b)**2)  )
        #var4 = 4*((n*b*rho4p)**2)*(varbn/((b*n)**2) + (sig_rho4/rho4p)**2)
        varae = vare*(a**2) + vara*(e**2) + 2*covan*(a*e)
        #varan = ((n*a)**2)*( (varn/((n)**2)) + (vara/((a)**2)) + 2*covan/(a*n) ) 
        var5 = 4*( (rho5**2)*varae + (sig_rho5**2)*((e*a)**2)  )
        #var5 = 4*((n*a*rho5p)**2)*(varan/((a*n)**2) + (sig_rho5/rho5p)**2) 
        plt.clf()
        lfontsize = 7
        if (abe_bool):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(np.diag(cov_rho0)), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e**2)*rho3, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*e)*rho4, np.sqrt(var4), legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*e*a)*rho5, np.sqrt(var5), legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (ab_bool):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (ae_bool):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(np.diag(cov_rho0)), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e**2)*rho3, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*e*a)*rho5, np.sqrt(var5), legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (be_bool):
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e**2)*rho3, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*e)*rho4, np.sqrt(var4), legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (a_bool):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (b_bool):
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (e_bool):
            pretty_rho(meanr, (e**2)*rho3, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
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

def getxibias_overall(pars, datarhos, models_combo,  plots=False, nameterms='terms_dxi.png',dxiname='dxi.png'):
    from src.plot_stats import pretty_rho
    a = b = e =  0; 
    eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    meanr, rhos, covrhos =  datarhos
    rho0, rho1, rho2, rho3, rho4, rho5 = rhos
    
    if not (abe_bool or ab_bool or a_bool): abe_bool = True
    if(abe_bool):
        a, b, e = pars
    if(ab_bool):
        a, b = pars 
    if(ae_bool):
        a, e = pars
    if(be_bool):
        b, e = pars
    if(a_bool):
        a= pars
    if(b_bool):
        b = pars
    if(e_bool):
        e = pars
    #Ploting each term of the bias
    if(plots):
        xlim = [2., 300.]
        plt.clf()
        lfontsize = 7
        if (abe_bool):
            pretty_rho(meanr, (a**2)*rho0, None, legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, None, legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e**2)*rho3, None, legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, None, legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*e)*rho4, None, legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*e*a)*rho5, None, legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (ab_bool):
            pretty_rho(meanr, (a**2)*rho0, None, legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, None, legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2,None, legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (ae_bool):
            pretty_rho(meanr, (a**2)*rho0, None, legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e**2)*rho3, None, legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*e*a)*rho5, None, legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (be_bool):
            pretty_rho(meanr, (b**2)*rho1, None, legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e**2)*rho3, None, legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*e)*rho4, None, legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (a_bool):
            pretty_rho(meanr, (a**2)*rho0, None, legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (b_bool):
            pretty_rho(meanr, (b**2)*rho1, None, legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (e_bool):
            pretty_rho(meanr, (e**2)*rho3, None, legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)

        
    dxi = (a**2)*rho0 + (b**2)*rho1 + (e**2)*rho3 + (2*a*b)*rho2 + (2*b*e)*rho4 + (2*e*a)*rho5
    if(plots):
        plt.clf()
        pretty_rho(meanr, dxi, None , legend=r"$\delta \xi$",  ylabel=r"$\delta \xi$",  xlim=xlim)
        print('Printing',  dxiname)
        plt.savefig(dxiname, dpi=150)

    return meanr, dxi

def getxibiastomo_margin(samples_i, samples_j, datarhos, models_combo, plots=False, bins=[0, 0], nameterms='terms_dxi.png',dxiname='dxi.png'):
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

    eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
    if not (abe_bool or ab_bool or a_bool): abe_bool = True
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
        xlim = [2., 300.]
        plt.clf()
        i, j = bins
        plt.title(r'$\delta \xi_{+}$ contributions, zbin: (%d,%d)'%(i,j))
        lfontsize = 7
        varterm0 = ((a_j*rho0)**2)*vara_i + ((a_i*rho0)**2)*vara_j + ((a_i*a_j)**2)*np.diag(cov_rho0) +  2*a_i*a_j*(rho0**2)*cova_ia_j
        varterm1 = ((b_j*rho1)**2)*varb_i + ((b_i*rho1)**2)*varb_j + ((b_i*b_j)**2)*np.diag(cov_rho1) +  2*b_i*b_j*(rho1**2)*covb_ib_j
        varterm2 = ((e_j*rho3)**2)*vare_i + ((e_i*rho3)**2)*vare_j + ((e_i*e_j)**2)*np.diag(cov_rho3) +  2*e_i*e_j*(rho3**2)*cove_ie_j    
        varterm3 = ((a_j*rho2)**2)*varb_i + ((b_i*rho2)**2)*vara_j + ((a_i*rho2)**2)*varb_j + ((b_j*rho2)**2)*vara_i + ((b_i*a_j + b_j*a_i)**2)*np.diag(cov_rho2)+ 2*(rho0**2)*(a_j*b_i*cova_jb_i+a_j*a_i*covb_ib_j+a_j*b_j*cova_ib_i + b_i*a_i*cova_jb_j + b_i*b_j*cova_ia_j + a_i*b_j*cova_ib_j)
        varterm4 = ((e_j*rho4)**2)*varb_i + ((b_i*rho4)**2)*vare_j + ((e_i*rho4)**2)*varb_j + ((b_j*rho4)**2)*vare_i + ((b_i*e_j + b_j*e_i)**2)*np.diag(cov_rho4)+ 2*(rho4**2)*(e_j*b_i*covb_ie_j+e_j*e_i*covb_ib_j+e_j*b_j*covb_ie_i + b_i*e_i*covb_je_j + b_i*b_j*cove_ie_j + e_i*b_j*covb_je_i)
        varterm5 = ((a_j*rho5)**2)*vare_i + ((e_i*rho5)**2)*vara_j + ((a_i*rho5)**2)*vare_j + ((e_j*rho5)**2)*vara_i + ((e_i*a_j + e_j*a_i)**2)*np.diag(cov_rho5)+ 2*(rho0**2)*(a_j*e_i*cova_je_i+a_j*a_i*cove_ie_j+a_j*e_j*cova_ie_i + e_i*a_i*cova_je_j + e_i*e_j*cova_ia_j + a_i*e_j*cova_ie_j)
        if (abe_bool):
            pretty_rho(meanr, a_i*a_j*rho0, np.sqrt(varterm0), legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, b_i*b_j*rho1, np.sqrt(varterm1), legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, e_i*e_j*rho3, np.sqrt(varterm2), legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*a_j + b_j*a_i)*rho2 , np.sqrt(varterm3), legend=r'$(\alpha^{(%d)}\beta^{(%d)}+\alpha^{(%d)}\beta^{(%d)})\rho_{2}$'%(i, j, j, i),lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*e_j + b_j*e_i)*rho4, np.sqrt(varterm4), legend=r'$(\beta^{(%d)}\eta^{(%d)}+\beta^{(%d)}\eta^{(%d)})\rho_{4}$'%(i, j, j, i),lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e_i*a_j +e_j*a_i)*rho5, np.sqrt(varterm5), legend=r'$(\eta^{(%d)}\alpha^{(%d)}+\eta^{(%d)}\alpha^{(%d)})\rho_{5}$'%(i, j, j, i), lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig(nameterms, dpi=200)
        if (ab_bool):
            pretty_rho(meanr, a_i*a_j*rho0, np.sqrt(varterm0), legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, b_i*b_j*rho1, np.sqrt(varterm1), legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*a_j + b_j*a_i)*rho2 , np.sqrt(varterm3), legend=r'$(\alpha^{(%d)}\beta^{(%d)}+\alpha^{(%d)}\beta^{(%d)})\rho_{2}$'%(i, j, j, i),lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',nameterms)
            plt.savefig( nameterms, dpi=200)
        if (ae_bool):
            pretty_rho(meanr, a_i*a_j*rho0, np.sqrt(varterm0), legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, e_i*e_j*rho3, np.sqrt(varterm2), legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e_i*a_j +e_j*a_i)*rho5, np.sqrt(varterm5), legend=r'$(\eta^{(%d)}\alpha^{(%d)}+\eta^{(%d)}\alpha^{(%d)})\rho_{5}$'%(i, j, j, i), lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig( nameterms, dpi=200)
        if (be_bool):
            pretty_rho(meanr, b_i*b_j*rho1, np.sqrt(varterm1), legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, e_i*e_j*rho3, np.sqrt(varterm2), legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*e_j + b_j*e_i)*rho4, np.sqrt(varterm4), legend=r'$(\beta^{(%d)}\eta^{(%d)}+\beta^{(%d)}\eta^{(%d)})\rho_{4}$'%(i, j, j, i),lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig(nameterms, dpi=200)
        if (a_bool):
            pretty_rho(meanr, a_i*a_j*rho0, np.sqrt(varterm0), legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig( nameterms, dpi=200)
        if (b_bool):
            pretty_rho(meanr, b_i*b_j*rho1, np.sqrt(varterm1), legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (e_bool):
            pretty_rho(meanr, e_i*e_j*rho3, np.sqrt(varterm2), legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig( nameterms, dpi=200)
        plt.clf()
        pretty_rho(meanr, dxi, np.sqrt(np.diag(covmat_dxi)) , legend=r"$\delta \xi$",  ylabel=r"$\delta \xi$",  xlim=xlim)
        print('Printing',  dxiname)
        plt.savefig(dxiname, dpi=200)
        
    return meanr, dxi, covmat_dxi
 
def getxibiastomo_overall(pars_i, pars_j , datarhos, models_combo, plots=False, bins=[0, 0],  nameterms='terms_dxi.png',dxiname='dxi.png'):
    from src.plot_stats import pretty_rho
    a_i = b_i = e_i = a_j = b_j = e_j =  0; 
    eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    meanr, rhos, covrhos =  datarhos
    rho0, rho1, rho2, rho3, rho4, rho5 = rhos
    
    if not (abe_bool or ab_bool or a_bool): abe_bool = True
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
    if(plots):
        xlim = [2., 300.]
        plt.clf()
        i, j = bins
        plt.title(r'$\delta \xi_{+}$ contributions, zbin: (%d,%d)'%(i,j))
        lfontsize = 7
        if (abe_bool):
            pretty_rho(meanr, a_i*a_j*rho0, None, legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, b_i*b_j*rho1, None, legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, e_i*e_j*rho3, None, legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*a_j + b_j*a_i)*rho2 , None, legend=r'$(\alpha^{(%d)}\beta^{(%d)}+\alpha^{(%d)}\beta^{(%d)})\rho_{2}$'%(i, j, j, i),lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*e_j + b_j*e_i)*rho4, None, legend=r'$(\beta^{(%d)}\eta^{(%d)}+\beta^{(%d)}\eta^{(%d)})\rho_{4}$'%(i, j, j, i),lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e_i*a_j +e_j*a_i)*rho5, None, legend=r'$(\eta^{(%d)}\alpha^{(%d)}+\eta^{(%d)}\alpha^{(%d)})\rho_{5}$'%(i, j, j, i), lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig(nameterms, dpi=200)
        if (ab_bool):
            pretty_rho(meanr, a_i*a_j*rho0, None, legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, b_i*b_j*rho1, None, legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*a_j + b_j*a_i)*rho2 , None, legend=r'$(\alpha^{(%d)}\beta^{(%d)}+\alpha^{(%d)}\beta^{(%d)})\rho_{2}$'%(i, j, j, i),lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',nameterms)
            plt.savefig( nameterms, dpi=200)
        if (ae_bool):
            pretty_rho(meanr, a_i*a_j*rho0, None, legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, e_i*e_j*rho3, None, legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (e_i*a_j +e_j*a_i)*rho5, None, legend=r'$(\eta^{(%d)}\alpha^{(%d)}+\eta^{(%d)}\alpha^{(%d)})\rho_{5}$'%(i, j, j, i), lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig( nameterms, dpi=200)
        if (be_bool):
            pretty_rho(meanr, b_i*b_j*rho1, None, legend=r'$\beta^{%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, e_i*e_j*rho3, None, legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b_i*e_j + b_j*e_i)*rho4, None, legend=r'$(\beta^{(%d)}\eta^{(%d)}+\beta^{(%d)}\eta^{(%d)})\rho_{4}$'%(i, j, j, i),lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig(nameterms, dpi=200)
        if (a_bool):
            pretty_rho(meanr, a_i*a_j*rho0, None, legend=r'$\alpha^{(%d)}\alpha^{(%d)} \rho_{0}$'%(i, j),lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig( nameterms, dpi=200)
        if (b_bool):
            pretty_rho(meanr, b_i*b_j*rho1, None, legend=r'$\beta^{(%d)}\beta^{(%d)} \rho_{1}$'%(i, j),lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (e_bool):
            pretty_rho(meanr, e_i*e_j*rho3, None, legend=r'$\eta^{(%d)}\eta^{(%d)} \rho_{3}$'%(i, j), lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            print('Printing', nameterms)
            plt.savefig( nameterms, dpi=200)
        plt.clf()
        pretty_rho(meanr, dxi, None , legend=r'zbin: (%d,%d)'%(i,j),  ylabel=r'$\delta \xi$',  xlim=xlim)
        print('Printing',  dxiname)
        plt.savefig(dxiname, dpi=200)
        
            
    return meanr, dxi

def getflagsnames(models_combo):
    eq,  abe, ab,  ae, be, a, b, e =  models_combo    
    if not (abe or ab or a): abe = True
    ## ALPHA-BETA-ETA
    if(abe):
        print("### Runing alpha, beta and eta test ### ")
        mflags = [True, True, True]
        namemc = 'mcmc_alpha-beta-eta_eq_' + str(eq) + '.png'
        namecont = 'contours_alpha-beta-eta_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_alpha-beta-eta_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_alpha-beta-eta_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_alpha-beta-eta_eq_' + str(eq) + '.png'
        namedxip = 'xibias_abe_' + str(eq) + '_.png'
        filename =  'abe_dxi_eq'+ str(eq) + '.fits'            
    ## ALPHA-BETA
    if(ab):
        print("### Runing alpha and beta test ### ")
        mflags = [True, True, False] ##alpha,beta,eta
        namemc = 'mcmc_alpha-beta_eq_' + str(eq) + '_png'
        namecont = 'contours_alpha-beta_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_alpha-beta_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_alpha-beta_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_alpha-beta_eq_' + str(eq) + '.png'
        namedxip = 'xibias_ab_' + str(eq) + '.png'
        filename =  'ab_dxi_eq'+ str(eq) + '.fits'
    ## ALPHA-ETA
    if(ae):
        print("### Runing alpha and eta test ### ")
        mflags = [True, False, True]
        namemc = 'mcmc_alpha-eta_eq_' + str(eq) + '.png'
        namecont = 'contours_alpha-eta_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_alpha-eta_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_alpha-eta_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_alpha-eta_eq_' + str(eq) + '.png'
        namedxip = 'xibias_ae_' + str(eq) + '.png'
        filename =  'ae_dxi_eq'+ str(eq) + '.fits'
    ## BETA-ETA
    if(be):
        print("### Runing beta and eta test ### ")
        mflags = [True, False, True]
        namemc = 'mcmc_beta-eta_eq_' + str(eq) + '.png'
        namecont = 'contours_beta-eta_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_beta-eta_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_beta-eta_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_beta-eta_eq_' + str(eq) + '.png'
        namedxip = 'xibias_be_' + str(eq) + '.png'
        filename =  'be_dxi_eq'+ str(eq) + '.fits' 
    ## ALPHA
    if(a):
        print("### Runing alpha test ### ")
        mflags = [True, False, False]
        namemc = 'mcmc_alpha_eq_' + str(eq) + '.png'
        namecont = 'contours_alpha_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_alpha_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_alpha_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_alpha_eq_' + str(eq) + '.png'
        namedxip = 'xibias_a_' + str(eq) + '.png'
        filename =  'a_dxi_eq'+ str(eq) + '.fits'
    ## Beta
    if(b):
        print("### Runing beta test ### ")
        mflags = [False, True, False] 
        namemc = 'mcmc_beta_eq_' + str(eq) + '.png'
        namecont = 'contours_beta_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_beta_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_beta_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_beta_eq_' + str(eq) + '.png'
        namedxip = 'xibias_b_' + str(eq) + '.png'
        filename =  'b_dxi_eq'+ str(eq) + '.fits'
    ## Eta
    if(e):
        print("### Runing eta test ### ")
        mflags = [False, False, True]
        namemc = 'mcmc_eta_eq_' + str(eq) + '.png'
        namecont = 'contours_eta_eq_' + str(eq) + '.png'
        nameterms = 'termsdxi_eta_eq_' + str(eq) + '.png'
        namecovmat = 'covmatrix_eta_eq_' + str(eq) + '.png'
        namebfres = 'Bestfitresiduals_eta_eq_' + str(eq) + '.png'
        namedxip = 'xibias_e_' + str(eq) + '.png'
        filename =  'e_dxi_eq'+ str(eq) + '.fits'
    return [ mflags, namemc, namecont, namecovmat, namebfres, nameterms, namedxip, filename]

def write_singlexip_margin( samplesp, samplesm, rhoscosmo,  models_combo, plots,  outpath,  plotspath):
    import numpy as np
    from src.readfits import read_rhos_plots
    from astropy.io import fits
    nameterms, namedxip, filename = getflagsnames(models_combo)[-3: ] 
    meanr, rhos,  covrhos =  read_rhos_plots(rhoscosmo)
    datarhosp =  [meanr, [rhos[2*i] for i in range(6)],  [ covrhos[2*i] for i in range(6)]  ]
    datarhosm =  [meanr, [rhos[2*i+1] for i in range(6)],  [ covrhos[2*i + 1] for i in range(6)] ]

    meanr, xip, covxip = getxibias_margin(samplesp, datarhosp, plots=plots,
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
def write_singlexip_overall( parsp, pasrsm, rhoscosmo,  models_combo,   plots,  outpath,  plotspath):
    import numpy as np
    from src.readfits import read_rhos_plots
    from astropy.io import fits
    nameterms,  namedxip, filename = getflagsnames(models_combo)[-3: ]
    meanr, rhos,  covrhos =  read_rhos_plots(rhoscosmo)
    datarhosp =  [meanr, [rhos[2*i] for i in range(6)],  [ covrhos[2*i] for i in range(6)]  ]
    datarhosm =  [meanr, [rhos[2*i+1] for i in range(6)],  [ covrhos[2*i + 1] for i in range(6)] ]

    meanr, xip= getxibias_overall(parsp, datarhosp, plots=plots,
                                  nameterms=plotspath + 'p_' +
                                  nameterms, dxiname=plotspath
                                  +'p_' + namedxip)
    meanr, xim = getxibias_overall(parsm, datarhosm, plots=plots,
                                   nameterms=plotspath +'m_' +
                                   nameterms, dxiname=plotspath
                                   +'m_' + namedxip)
    
    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = len(meanr)
    outdata = np.recarray((nrows, ), dtype=dtype)
    
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    
    bin1array = np.array([ -999]*nrows)
    bin2array = np.array([ -999]*nrows)
    angbinarray = np.arange(nrows)
    array_list = [bin1array, bin2array, angbinarray, xip,  meanr ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='xip')
    hdul.insert(1, corrhdu)
    array_list = [bin1array, bin2array, angbinarray, xim,  meanr ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='xim')
    hdul.insert(2, corrhdu)
    hdul.writeto(filename, overwrite=True)
    print(filename,'Written!')
def write_tomoxip_margin(samplesp_list, samplesm_list, rhoscosmo, models_combo, plots, outpath, plotspath ):
    import itertools
    from src.readfits import read_rhos_plots
    from astropy.io import fits
    import numpy as np
    nameterms, namedxip, filename = getflagsnames(models_combo)[-3: ] 
    meanr, rhos,  covrhos =  read_rhos_plots(rhoscosmo)
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
        binstr = 'margin_zbin_%d_%d'%(i, j)
        samplesp_i =  samplesp_list[i - 1]; samplesm_i =  samplesm_list[i - 1]
        samplesp_j =  samplesp_list[j - 1]; samplesm_j =  samplesm_list[j - 1]
        bins=[i, j]
        meanr, xip, covxip = getxibiastomo_margin( samplesp_i, samplesp_j, datarhosp, models_combo,
                                                   plots=plots, bins=bins, nameterms= plotspath + binstr+ nameterms, 
                                                   dxiname=plotspath + binstr+ namedxip )
        meanr, xim, covxim = getxibiastomo_margin( samplesm_i, samplesm_j ,datarhosm, models_combo,
                                                   plots=plots, bins=bins, nameterms= plotspath + binstr+ nameterms, 
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
        
def write_tomoxip_overall(parsp_list, parsm_list, rhoscosmo, models_combo, plots, outpath, plotspath ):
    import itertools
    from src.readfits import read_rhos_plots
    from astropy.io import fits
    import numpy as np
    nameterms, namedxip, filename = getflagsnames(models_combo)[-3: ] 
    meanr, rhos,  covrhos =  read_rhos_plots(rhoscosmo)
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
        parsp_i =  parsp_list[i - 1]; parsm_i =  parsm_list[i - 1]
        parsp_j =  parsp_list[j - 1]; parsm_j =  parsm_list[j - 1]
        bins = [i, j]
        binstr = 'overall_zbin:%d_%d_'%(i, j)
        meanr, xip = getxibiastomo_overall( parsp_i, parsp_j, datarhosp, models_combo, plots=plots, bins=bins, nameterms= plotspath + binstr+ nameterms, 
                                                   dxiname=plotspath + binstr+ namedxip)
        meanr, xim = getxibiastomo_overall( parsm_i, parsm_j ,datarhosm, models_combo, plots=plots, bins=bins, nameterms= plotspath + binstr+ nameterms, 
                                                   dxiname=plotspath + binstr+ namedxip)
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
    hdul.writeto(outpath + 'overall_' + filename, overwrite=True)
    print(outpath + 'overall_' + filename,'Written!')
    
def RUNTEST(i_guess, data, nwalkers, nsteps, eq='All', mflags=[True, True, True], xip=True, xim=False, moderr=False, uwmprior=False, minimize=True, margin=True, overall=False ):
    from src.chi2 import minimizeCHI2
    from src.maxlikelihood import MCMC
    import numpy as np
    if (margin and not overall):
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
    elif (overall and not margin):
        print("Getting best fit of overall likelihood")
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            mflags=mflags, xip=xip,
                                            xim=xim, moderr=moderr)
        
        return fitted_params, chisq
    else:
        print("Not proper configuration, choose only one, or overall or margin")

def RUNTEST_PERTAU(rhofile, taufile, minscale, maxscale, models_combo, nwalkers, nsteps,  uwmprior, splitxipxim, margin, overall,  plots,  plotspath, zbin=0, axs=None):
    import numpy as np
    from src.readfits import  read_rhos, read_taus
    from src.plot_stats import plotallrhosfits, plotalltausfits, plot_samplesdist, plotbestfit, plotbestfitresiduals

    if (plots):
        xlim = [1, 300.]
        taustitle = ''
        plotalltausfits(taufile, outpath=plotspath, title='zbin: %d'%(zbin),  xlim=xlim, zbin=str(zbin))
    
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
    i_guess0 = [ -0.01, 1, -1 ] #fiducial values
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
        
    if (margin and not overall):
        if(plots):
            plot_samplesdist(auxp1, auxp2 , mflags, nwalkers, nsteps, plotspath +'p_zbin_%d_%s'%(zbin,namemc),plotspath + 'p_zbin_%d_%s'%(zbin,namecont),zbin=zbin )
            plot_samplesdist(auxm1 , auxm2, mflags, nwalkers, nsteps, plotspath +'m_zbin_%d_%s'%(zbin,namemc),plotspath + 'm_zbin_%d_%s'%(zbin,namecont),zbin=zbin )
            plotcovpars(auxp1, namecovmat=plotspath + 'p_zbin_%d'%(zbin) + namecovmat)
            plotcovpars(auxm1, namecovmat=plotspath + 'm_zbin_%d'%(zbin) + namecovmat)
            #title='zbin %d %s'%(zbin,namebfres[17:-4])
            #plotbestfitresiduals(auxp1, auxm1, meanr, data, models_combo,  plotspath +'zbin_%d'%(zbin) + namebfres, title=title.replace('_','\_') , margin=margin, overall=overall)
        #samplesp, samplesm
        return auxp1, auxm1 

    if (overall and not margin):
        if (plots):
            if axs is not None:
                plotbestfit(zbin, axs, auxp1, auxm1, meanr, data, models_combo,  plotspath, margin=margin, overall=overall)
            
 
                    
#plt.savefig(plotpath, dpi=150)
            #title='zbin %d %s'%(zbin,namebfres[17:-4])
            #plotbestfitresiduals(auxp1, auxm1, meanr, data, models_combo,  plotspath +'zbin_%d'%(zbin) + namebfres, title=title.replace('_', '\_'), margin=margin, overall=overall)
            
        #parsp,chisqp,parsm,chisqm
        return auxp1, auxp2, auxm1, auxm2

                        
def main():
    from src.plot_stats import plotallrhosfits
    from src.maxlikelihood import percentiles
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
        xlim = [1, 300.]
        ylims = [[1.e-11,5.e-6 ],[1.e-12,1.e-6 ],[3.e-9 ,3.e-4 ],[3.e-12 ,3.e-6 ]]
        rhostitle = ''
        plotallrhosfits(args.rhos, outpath=plotspath, title=rhostitle, xlim=xlim, ylims=ylims)

    models_combo = [args.eq, args.abe, args.ab, args.ae, args.be, args.a, args.b, args.e]
    eq = args.eq
    nsig = args.nsig
    nwalkers,  nsteps = args.nwalkers, args.nsteps

    if args.singletau is not None:
        if args.margin:
            print("STARTING single tau ANALYSIS",  arg.singletau)
            samplesp, samplesm =RUNTEST_PERTAU(args.rhos,args.singletau, args.minscale, args.maxscale,
                                               models_combo ,nwalkers,nsteps, args.uwmprior,
                                               args.splitxipxim, True, False, args.plots,  plotspath)
            mcmcpars = percentiles(samplesp, nsig=nsig) 
            print( ' mcmc parameters xi+',  'nsig=', nsig, ' percentiles: ',  mcmcpars)
            mcmcpars = percentiles(samplesm, nsig=nsig) 
            print( ' mcmc parameters xi-',  'nsig=', nsig, ' percentiles: ',  mcmcpars)
            write_singlexip_margin(samplesp, samplesm, args.rhoscosmo,  models_combo, args.plots,  outpath,  plotspath)
        if args.overall:
            parsp, chi2p_nu, parsm, chi2m_nu = RUNTEST_PERTAU(args.rhos,args.singletau, args.minscale, args.maxscale,
                                               models_combo ,nwalkers,nsteps, args.uwmprior,
                                               args.splitxipxim, False, True, args.plots, plotspath)
            print( ' overall likelihood best fit xi+',  parsp, 'chi2_reduced', chi2p_nu)
            print( ' overall likelihood best fit xi+',  parsm, 'chi2_reduced', chi2m_nu)
            write_singlexip_overall(parsp, parsm, args.rhoscosmo,  models_combo, args.plots,  outpath,  plotspath)
        
    else:
        print("STARTING TOMOGRAPHIC ANALYSIS")
        samplesp_list = []; parsp_list = []
        samplesm_list = []; parsm_list = []
        
        if args.margin:
            for i,  taufile in enumerate(args.taus):
                samplesp, samplesm=RUNTEST_PERTAU(args.rhos,taufile,args.minscale, args.maxscale,
                                                  models_combo ,nwalkers,nsteps, args.uwmprior, args.splitxipxim,
                                                  True, False, args.plots, plotspath,  zbin=(i + 1))
                mcmcparsp = percentiles(samplesp, nsig=nsig) 
                mcmcparsm = percentiles(samplesm, nsig=nsig)
                if args.splitxipxim :
                    print( ' mcmc parameters xi+',  'nsig=', nsig, ' percentiles: ',  mcmcparsp)
                    print( ' mcmc parameters xi-',  'nsig=', nsig, ' percentiles: ',  mcmcparsm)
                else:
                    print( ' mcmc parameters',  'nsig=', nsig, ' percentiles: ',  mcmcparsp)
                    
                samplesp_list.append(samplesp); samplesm_list.append(samplesm)
            write_tomoxip_margin( samplesp_list, samplesm_list, args.rhoscosmo,  models_combo, args.plots,  outpath,  plotspath)

        if args.overall:
            fig1, ax1 = plt.subplots()
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
            fig4, ax4 = plt.subplots()
            fig5, ax5 = plt.subplots()
            fig6, ax6 = plt.subplots()
            fig7, ax7 = plt.subplots()
            axs=[ax1,ax2,ax3,ax4,ax5,ax6,ax7]
            figs=[fig1,fig2,fig3,fig4,fig5,fig6]
            ylabels=[r'$\tau_{0+}$',r'$\tau_{0-}$',r'$\tau_{2+}$',r'$\tau_{2-}$',r'$\tau_{5+}$',r'$\tau_{5-}$']
            for i,  taufile in enumerate(args.taus):
                parsp, chi2p_nu, parsm, chi2m_nu =RUNTEST_PERTAU(args.rhos,taufile,args.minscale, args.maxscale,
                                                                 models_combo ,nwalkers,nsteps, args.uwmprior, args.splitxipxim,
                                                                 False, True, args.plots, plotspath,  zbin=(i + 1), axs=axs)
                if args.splitxipxim :
                    print( 'overall parameters xi+', parsp.tolist()) 
                    print( 'chi2r: ',  chi2p_nu)
                    print( 'overall parameters xi-', parsm.tolist()) 
                    print( 'chi2r: ',  chi2m_nu)
                else:
                    print( 'overall parameters ', parsp.tolist()) 
                    print( 'chi2r: ',  chi2p_nu)
                parsp_list.append(parsp); parsm_list.append(parsm)
            if(args.plots):
                for i, fig in enumerate(figs):
                    print('Printing', plotspath+'%d_tau_bestfit.png'%(i))
                    axs[i].set_ylabel(ylabels[i])
                    axs[i].set_title('Alpha-beta')
                    fig.tight_layout()
                    fig.savefig(plotspath+'%d_tau_bestfit_ab.png'%(i),dpi=200)
                
            write_tomoxip_overall( parsp_list, parsm_list, args.rhoscosmo,  models_combo, args.plots,  outpath,  plotspath )

    

  
    
if __name__ == "__main__":
    main()



        
