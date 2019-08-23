import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='run abn test and save files arrays of parameters going out to a maximum bin.  #plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model')
    parser.add_argument('--rhos',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/allrhos_4jk.fits',
                        help='location of the file containing all the rhos by patch.')
    parser.add_argument('--taus',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/alltaus_4jk.fits',
                        help='location of the file containing all the taus by patch')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='jkscaletest.fits', help='Name of the fit file where the info of the jkscaletest will be saved ')

    args = parser.parse_args()

    return args

def write_fit(data, filename, count=1,  clobber=False):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    from astropy.io import fits
    #If the first time it enters the file exist does not write
    #anything write only if clobber is set to True
    if (count == 1):
        #If does not exist, create it and write data
        if not os.path.isfile(filename):
            hdu = fits.PrimaryHDU()
            hdul = fits.HDUList([hdu])
            corrhdu = fits.BinTableHDU(data, name='')
            hdul.insert(1, corrhdu)
            hdul.writeto(filename, overwrite=clobber)
            print("Creating file: ", filename)
            #If exist open it and add data
        else:
            if clobber:
                hdu = fits.PrimaryHDU()
                hdul = fits.HDUList([hdu])
                corrhdu = fits.BinTableHDU(data, name='')
                hdul.insert(1, corrhdu)
                hdul.writeto(filename, overwrite=clobber)
                print("Clobering file: ", filename)
            else:
                raise Exception('You tried to overwrite an existing file, without setting clobber')
    else:
        fits = FITS(filename,'rw')
        fits[-1].append(data)
        print("Apending File: ",  filename)
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

def read_pars(filename,  ijk, mflags):
    import fitsio
    #print("Reading file",  filename)
    aflag, bflag, eflag =  mflags
    alldata =  fitsio.read(filename, ext=1)
    alldata = alldata[alldata['JKR']==ijk]
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
                
def run_parspatch(outpath,filename, rhosfile, tausfile, nsig=1,  njk=4, mflags=[True, True, True]  ):
    from readjson import read_rhos, read_taus
    from totalchi2 import minimizeCHI2
    from totalmaxlikelihood import MCMC, percentiles, bestparameters
    import numpy as np
    import fitsio

    ##ouputfile format
    names = ['JKR','THETA', 'a_l', 'a_c', 'a_r', 'b_l', 'b_c', 'b_r', 'e_l', 'e_c', 'e_r']
    forms = ['i4',  'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8','f8', 'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    ## Info previous to the jktest.
    nwalkers,  nsteps = 100,  1000
    eq = 'All'; moderr = False
    uwmprior = False
    i_guess0 = [ -0.01, 1,- 1] #fiducial values
    i_guess = np.array(i_guess0)[np.array(mflags)].tolist()
    data = {}
    data_rhos =  fitsio.read(rhosfile, ext=1)
    data_taus =  fitsio.read(tausfile, ext=1)
    rhosnames =  ['RHO0P','RHO1P','RHO2P','RHO3P','RHO4P','RHO5P']
    varrhosnames =  ['VAR_RHO0','VAR_RHO1','VAR_RHO2','VAR_RHO3','VAR_RHO4','VAR_RHO5']
    tausnames =  ['TAU0P','TAU2P','TAU5P']
    vartausnames =  ['VAR_TAU0','VAR_TAU2','VAR_TAU5']
    meanr = data_rhos['THETA'][data_rhos['JKR']==0]
    
    for i in range(njk):
        jkrhos = (data_rhos['JKR'] == i)
        jktaus = (data_taus['JKR'] == i)

        jkrarr =  np.array([i]*nrows)
        ac = []; al = []; ar = []
        bc = []; bl = []; br = []
        ec = []; el = []; er = []
        for ibin in range(1, 21):
            angrhos = (data_rhos['ANGBIN']< ibin)
            angtaus = (data_taus['ANGBIN']< ibin)
            rhosbool = jkrhos&angrhos
            tausbool = jktaus&angtaus

            data['rhos'] = [ data_rhos[name][rhosbool] for name in rhosnames]
            data['covrhos'] = [ np.diag(data_rhos[name][rhosbool]) for name in varrhosnames]
            data['taus'] = [ data_taus[name][tausbool] for name in tausnames]
            data['covtaus'] = [np.diag(data_taus[name][tausbool]) for name in vartausnames]

            fit_pars, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                           mflags=mflags, moderr=moderr)
            #print("Chi2 Minization", fit_pars, chisq)
            samples, chains = MCMC(fit_pars,data, nwalkers=nwalkers, nsteps=nsteps, eq=eq,
                           mflags=mflags, moderr=moderr, uwmprior=uwmprior)
            #print("Total samples", [len(i) for i in samples] )
            #samples2= np.c_[[par[5000:] for par in samples]].T
            mcmcpars = percentiles(samples, nsig=nsig)
            #bestpars = bestparameters(samples)
            #print("Best pars", bestpars )
            fillparlists(al, ac, ar, bl, bc, br, el, ec, er, mcmcpars, mflags)

        array_list = [jkrarr, meanr, np.array(al), np.array(ac), np.array(ar),
                      np.array(bl), np.array(bc), np.array(br), np.array(el), np.array(ec), np.array(er)]
        for array, name in zip(array_list, names): outdata[name] = array

        write_fit(outdata, outpath + filename, count=i+1, clobber=True)
   
def plot_correlations(outpath,  rhosjk,  tausjk,  njk):
    import fitsio
    from plot_stats import pretty_rho1,  pretty_rho2
    import numpy as np
    data_rhos =  fitsio.read(rhosjk, ext=1)
    data_taus =  fitsio.read(tausjk, ext=1)
    rhosnames =  ['RHO0P','RHO1P','RHO2P','RHO3P','RHO4P','RHO5P']
    varrhosnames =  ['VAR_RHO0','VAR_RHO1','VAR_RHO2','VAR_RHO3','VAR_RHO4','VAR_RHO5']
    tausnames =  ['TAU0P','TAU2P','TAU5P']
    vartausnames =  ['VAR_TAU0','VAR_TAU2','VAR_TAU5']
    meanr = data_rhos['THETA'][data_rhos['JKR']==0]
    
    for i in range(njk):
        title =  'JK Region ' + str(i)
        ylim0, ylim1, ylim2 = [None, None, None]
        xlim =  [ 2, 300] 
        jkrhos = (data_rhos['JKR'] == i)
        jktaus = (data_taus['JKR'] == i)
        rho0p = data_rhos['RHO0P'][jkrhos];rho1p = data_rhos['RHO1P'][jkrhos];rho2p = data_rhos['RHO2P'][jkrhos]
        rho3p = data_rhos['RHO3P'][jkrhos];rho4p = data_rhos['RHO4P'][jkrhos];rho5p = data_rhos['RHO5P'][jkrhos]
        sig_rho0 = np.sqrt(data_rhos['VAR_RHO0'])[jkrhos]; sig_rho1 = np.sqrt(data_rhos['VAR_RHO1'])[jkrhos]
        sig_rho2 = np.sqrt(data_rhos['VAR_RHO2'])[jkrhos]; sig_rho3 = np.sqrt(data_rhos['VAR_RHO3'])[jkrhos]
        sig_rho4 = np.sqrt(data_rhos['VAR_RHO4'])[jkrhos]; sig_rho5 = np.sqrt(data_rhos['VAR_RHO5'])[jkrhos]

        tau0p = data_taus['TAU0P'][jktaus];tau2p = data_taus['TAU2P'][jktaus];tau5p = data_taus['TAU5P'][jktaus]
        sig_tau0 = np.sqrt(data_taus['VAR_TAU0'])[jktaus]; sig_tau2 = np.sqrt(data_taus['VAR_TAU2'])[jktaus];
        sig_tau5 = np.sqrt(data_taus['VAR_TAU5'])[jktaus];
        
        plt.clf()
        pretty_rho1(meanr, rho1p, sig_rho1, rho3p, sig_rho3, rho4p,
                    sig_rho4, title=title, xlim=xlim, ylim=ylim1)
        print("Printing file: ", outpath +'rho1_all_rsrs_' + str(i) +  '_.png')
        plt.savefig(outpath +'rho1_all_rsrs_' + str(i) +  '_.png')
        plt.clf()
        pretty_rho2(meanr, rho0p, sig_rho0, rho2p, sig_rho2, rho5p,
                    sig_rho5, title=title, xlim=xlim, ylim=ylim2)
        print("Printing file: ", outpath +'rho2_all_rsrs_' + str(i) +  '_.png')
        plt.savefig(outpath +'rho2_all_rsrs_' + str(i) +  '_.png')
        plt.clf()
        pretty_rho2(meanr, tau0p, sig_tau0, tau2p, sig_tau2, tau5p,
                    sig_tau5, tauleg=True, title=title, xlim=xlim, ylim=ylim1)
        print("Printing file: ",outpath +'tau_all_rsgal_' + str(i) +  '_.png')
        plt.savefig(outpath +'tau_all_rsgal_' + str(i) +  '_.png')
        
def plotlineal( outpath, filename, mflags, njk):
    import numpy as np
    import itertools
    ylabels = [r'$\alpha$', r'$\beta$', r'$\eta$']
    ylims =  [[ - 0.1, 0.1],[ - 30, 30],[ -600, 600] ] 
    outputnames = ['alpha_quadrants.png', 'beta_quadrants.png', 'eta_quadrants.png']
    colors = ['black', 'green', 'blue', 'red', 'gray', 'pink']
    ndim =  len(list(itertools.compress(xrange(len(mflags)),  mflags))) 
    for t in range(ndim):
        plt.clf()
        for ijk in range(njk):
            lists = read_pars(filename, ijk, mflags)
            meanr = lists[0]
            al = lists[1 + 3*t]; ac = lists[2 + 3*t]; ar = lists[3+3*t]
            label = "P" + str(ijk) 
            plt.figure(t)
            plt.scatter(meanr,ac,color=colors[ijk],label=label,marker='o',s=10)
            plt.errorbar(meanr,ac,yerr=[-np.array(al),np.array(ar)], color=colors[ijk],capsize=0, linestyle='')
            
        plt.figure(t)
        plt.legend(loc='best', fontsize=10)
        plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
        plt.ylabel(ylabels[t], fontsize=24)
        plt.xscale('log')
        #plt.xlim( [ 2, 300] )
        plt.ylim(ylims[t] )
        plt.tight_layout()
        print("Printing :", outpath + outputnames[t])
        plt.savefig(outpath +outputnames[t],  dpi=200)
  
        
def main():
    import os
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    
    
    args = parse_args()
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    mflags = [True, True, False]
    nsig = 1
    njk = 4
    #plot_correlations('', args.rhos, args.taus, njk)
    #run_parspatch(outpath, args.filename, args.rhos, args.taus, nsig=nsig,  njk=njk, mflags=mflags)
    plotlineal('', outpath + args.filename,mflags, njk)
    


    
if __name__ == "__main__":
    main()
