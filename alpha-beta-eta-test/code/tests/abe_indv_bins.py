import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import numpy as np
from src.readfits import  read_rhos, read_taus
from src.maxlikelihood import percentiles
from astropy.io import fits
from src.chi2 import minimizeCHI2
from src.maxlikelihood import MCMC

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='run abn test and save files arrays of parameters going out to a maximum bin.  #plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model')
    parser.add_argument('--rhos', default='/data/catalogs/rhos-taus/marco/Y3_03_31_20/rho__JK_Y3_1.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--taus',
                        default='/data/catalogs/rhos-taus/marco/Y3_03_31_20/tau__JK_Y3_1.fits',
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
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
    parser.add_argument('--ab', default=True,
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
    parser.add_argument('--plotspath', default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations',
                        help='location of the output of the files')
    parser.add_argument('--outpath', default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='scaletest.fits', help='Name of the fit file where the info of the files will be saved ')

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

def RUNTEST_PERTAU(rhofile, taufile, minbin, maxbin, models_combo, nwalkers, nsteps,  uwmprior, splitxipxim, margin, overall,  plots,  plotspath, zbin=0, axs=None):
    import numpy as np
    from src.readfits import  read_rhos, read_taus
    from src.plot_stats import plotallrhosfits, plotalltausfits, plot_samplesdist, plotbestfit, plotbestfitresiduals,  plotcovmat

    if (plots):
        xlim = [0.1, 300.]
        if zbin is not None: title='zbin: %d'%(zbin)
        else: title = 'Non-tomographic'
        #plotalltausfits(taufile, outpath=plotspath, title=title,  xlim=xlim, zbin=zbin)
        
    
    meanr, rhos,  covrho =  read_rhos(rhofile, minbin=minbin, maxbin=maxbin)
    meanr, taus,  covtau =  read_taus(taufile, minbin=minbin, maxbin=maxbin)
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
        
    if (margin and not overall):
        if(plots):
            if zbin is not None: title = 'zbin %d'%(zbin)
            else: title =  'Non-tomographic' 
            if (splitxipxim):
                plot_samplesdist(auxp1, auxp2 , mflags, nwalkers, nsteps, os.path.join(plotspath,'%s_%s_p_.png'%(namemc, title)), os.path.join(plotspath, '%s_%s_p_.png'%(namecont, title )), title=title )
                plot_samplesdist(auxm1 , auxm2, mflags, nwalkers, nsteps, os.path.join(plotspath,'%s_%s_m_.png'%(namemc, title)), os.path.join(plotspath, '%s_%s_m_.png'%(namecont, title)), title=title )
                plotcovmat(auxp1, mflags, os.path.join(plotspath, '%s_%s_p_.png'%(namecovmat, title)))
                plotcovmat(auxm1, mflags, os.path.join(plotspath, '%s_%s_m_.png'%(namecovmat, title)))
            else:
                plot_samplesdist(auxp1, auxp2 , mflags, nwalkers, nsteps, os.path.join(plotspath, '%s_%s_.png'%(namemc, title)),  os.path.join(plotspath, '%s_%s_.png'%(namecont, title)), title=title )
                plotcovmat(auxp1, mflags, os.path.join(plotspath, '%s_%s_.png'%(namecovmat, title)))

            if axs is not None:
                plotbestfit(zbin, axs, auxp1, auxm1, meanr, data, models_combo,  plotspath, margin=margin, overall=overall)
                
        
        #samplesp, samplesm
        return auxp1, auxm1 

    if (overall and not margin):
        if (plots):
            if axs is not None:
                plotbestfit(zbin, axs, auxp1, auxm1, meanr, data, models_combo,  plotspath, margin=margin, overall=overall)
                        
        #parsp,chisqp,parsm,chisqm
        return auxp1, auxp2, auxm1, auxm2

def fillparlists(al, ac, ar, bl, bc, br, el, ec, er, mcmcpars,  mflags):
    aflag, bflag, eflag =  mflags
    if(aflag and bflag and eflag):
        al.append(mcmcpars[0][1]); ac.append(mcmcpars[0][0]);ar.append(mcmcpars[0][2])
        bl.append(mcmcpars[1][1]); bc.append(mcmcpars[1][0]);br.append(mcmcpars[1][2])
        el.append(mcmcpars[2][1]); ec.append(mcmcpars[2][0]);er.append(mcmcpars[2][2])
    elif( aflag and bflag and (not eflag)):
        al.append(mcmcpars[0][1]); ac.append(mcmcpars[0][0]);ar.append(mcmcpars[0][2])
        bl.append(mcmcpars[1][1]); bc.append(mcmcpars[1][0]);br.append(mcmcpars[1][2])
        el.append(-999); ec.append(-999);er.append(-999)
    elif(aflag and (not bflag) and eflag ):
        al.append(mcmcpars[0][1]); ac.append(mcmcpars[0][0]);ar.append(mcmcpars[0][2])
        el.append(mcmcpars[1][1]); ec.append(mcmcpars[1][0]);er.append(mcmcpars[1][2])
        bl.append(-999); bc.append(-999);br.append(-999)
    elif( (not aflag) and bflag and eflag ):
        bl.append(mcmcpars[0][1]); bc.append(mcmcpars[0][0]);br.append(mcmcpars[0][2])
        el.append(mcmcpars[1][1]); ec.append(mcmcpars[1][0]);er.append(mcmcpars[1][2])
        al.append(-999); ac.append(-999);ar.append(-999) 
    elif( aflag and (not bflag) and (not eflag)):
        al.append(mcmcpars[0][1]); ac.append(mcmcpars[0][0]);ar.append(mcmcpars[0][2])
        bl.append(-999); bc.append(-999);br.append(-999)
        el.append(-999 ); ec.append(-999);er.append( -999) 
    elif( (not aflag) and bflag and (not eflag)):
        bl.append(mcmcpars[0][1]); bc.append(mcmcpars[0][0]);br.append(mcmcpars[0][2])
        al.append(-999); ac.append(-999);ar.append(-999) 
        el.append(-999); ec.append(-999);er.append(-999)
    elif( (not aflag) and (not bflag) and eflag):
        el.append(mcmcpars[0][1]); ec.append(mcmcpars[0][0]);er.append(mcmcpars[0][2])
        bl.append(-999); bc.append(-999);br.append(-999)
        al.append(-999); ac.append(-999);ar.append(-999)  

def read_pars(filename,  ext, mflags):
    import fitsio
    #print("Reading file",  filename)
    aflag, bflag, eflag =  mflags
    alldata =  fitsio.read(filename, ext=ext)
    meanr = alldata['THETA']
    al = alldata['a_l']; ac = alldata['a_c']; ar = alldata['a_r']
    bl = alldata['b_l']; bc = alldata['b_c']; br = alldata['b_r']
    el = alldata['e_l']; ec = alldata['e_c']; er = alldata['e_r']
    if(aflag and bflag and eflag):
        return meanr, al, ac, ar, bl, bc, br, el, ec,er
    elif( aflag and bflag and (not eflag)):
        return meanr, al, ac, ar, bl, bc, br
    elif(aflag and (not bflag) and eflag ):
        return meanr, al, ac, ar, el, ec,er
    elif( (not aflag) and bflag and eflag ):
        return meanr, bl, bc, br, el, ec,er
    elif( aflag and (not bflag) and (not eflag)):
        return meanr, al, ac, ar
    elif( (not aflag) and bflag and (not eflag)):
        return meanr, bl, bc, br
    elif( (not aflag) and (not bflag) and eflag):
        return meanr, el, ec,er
def plotlineal( args, filename, plotspath):
    import numpy as np
    import itertools
    if not (args.abe or args.ab or args.ae or args.be or args.a or args.b or args.e): args.abe = True
    models_combo = [args.eq, args.abe, args.ab, args.ae, args.be, args.a, args.b, args.e]
    mflags = getflagsnames(models_combo)[0]
    
    ylabels = [r'$\alpha$', r'$\beta$', r'$\eta$']
    #ylims =  [[ - 0.07, 0.07],[ 0 , 5],[ -600, 600] ]
    ylims =  [None, None, None] 
    outputnames = ['alpha_scales.png', 'beta_scales.png', 'eta_scales.png']
    colors = ['black', 'green', 'blue', 'red', 'gray', 'pink']
    ndim =  len(list(itertools.compress(range(len(mflags)),  mflags))) 
    ntomobins = 1
    for t in range(ndim):
        plt.clf()
        for ext in range(ntomobins):
            lists = read_pars(filename, ext+1, mflags)
            meanr = lists[0]
            al = lists[1 + 3*t]; ac = lists[2 + 3*t]; ar = lists[3+3*t]
            #label = 'zbin%d'%(ext + 1)
            label = None
            plt.figure(t)
            plt.scatter(meanr,ac,color=colors[ext],label=label,marker='o',s=10)
            plt.errorbar(meanr,ac,yerr=[-np.array(al),np.array(ar)], color=colors[ext],capsize=2, linestyle='')
            
        plt.figure(t)
        plt.legend(loc='best', fontsize=10)
        plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
        plt.ylabel(ylabels[t], fontsize=24)
        plt.xscale('log')
        #plt.xlim( [ 0.9, 400] )
        plt.ylim(ylims[t] )
        plt.tight_layout()
        filenameout = os.path.join(plotspath, outputnames[t] )
        print("Printing :", filenameout)
        plt.savefig(filenameout,  dpi=200)
def run_allscales(args, outpath, plotspath):
    if not (args.abe or args.ab or args.ae or args.be or args.a or args.b or args.e): args.abe = True
    models_combo = [args.eq, args.abe, args.ab, args.ae, args.be, args.a, args.b, args.e]
    mflags = getflagsnames(models_combo)[0]
    nsig = args.nsig
    nwalkers,  nsteps = args.nwalkers, args.nsteps

    nbins = len(read_rhos(args.rhos)[1][0])
    meanr = read_rhos(args.rhos)[0]

    ##ouputfile format
    names = ['THETA', 'a_l', 'a_c', 'a_r', 'b_l', 'b_c', 'b_r', 'e_l', 'e_c', 'e_r']
    forms = [ 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8','f8', 'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    ac = []; al = []; ar = []
    bc = []; bl = []; br = []
    ec = []; el = []; er = []
    
    for idx in range(nbins):
        minbin = idx
        maxbin = idx + 1
        samplesp, samplesm=RUNTEST_PERTAU(args.rhos,args.taus,minbin,
                                          maxbin, models_combo
                                          ,nwalkers,nsteps,
                                          False, False, True,
                                          False, args.plots,
                                          plotspath)
        mcmcpars = percentiles(samplesp, nsig=nsig) 
        print( ' mcmc parameters xi+',  'nsig=', nsig, ' percentiles: ',  mcmcpars)
        mcmcpars = percentiles(samplesm, nsig=nsig) 
        print( ' mcmc parameters xi-',  'nsig=', nsig, ' percentiles: ',  mcmcpars)
        fillparlists(al, ac, ar, bl, bc, br, el, ec, er, mcmcpars, mflags)
    array_list = [ meanr, np.array(al), np.array(ac), np.array(ar),
                   np.array(bl), np.array(bc), np.array(br), np.array(el),
                   np.array(ec), np.array(er)]
    for array, name in zip(array_list, names): outdata[name] = array
    corrhdu = fits.BinTableHDU(outdata, name='zbin%d'%(0))
    hdul.insert(0, corrhdu)
    
    filename = os.path.join(outpath, args.filename)
    print("Printin file:", filename)
    hdul.writeto(filename, overwrite=True)

def main():    
    args = parse_args()
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
        if not os.path.exists(plotspath): raise

  
    
    run_allscales(args, outpath, plotspath)
    plotlineal(args, os.path.join(outpath, args.filename), plotspath)
    
if __name__ == "__main__":
    main()
