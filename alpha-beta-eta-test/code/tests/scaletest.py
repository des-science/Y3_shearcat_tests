import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='run abn test and save files arrays of parameters going out to a maximum bin.  #plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model')
    parser.add_argument('--rhos', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')
    parser.add_argument('--taus',
                        default=['/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_1.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_2.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_3.fits',
                                 '/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits'],
                        help='Ordered list of fits TAUS, containing all tau stats used to estimate abe')
    parser.add_argument('--nsig', default=1, type=int, 
                        help='How many sigman for the marginalized confidence interval')
    parser.add_argument('--nsteps', default=1000, type=int, 
                        help='nsteps of MCMC')
    parser.add_argument('--nwalkers', default=100, type=int, 
                        help='nwalkers of MCMC')
    parser.add_argument('--eq', default=10, type=int, 
                        help='Select equations to be used for istance --eq=0, 4 represent the whole system of equations')
    parser.add_argument('--abe', default=True,
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
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='scaletest.fits', help='Name of the fit file where the info of the jkscaletest will be saved ')

    args = parser.parse_args()

    return args


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
  


def getmflags(models_combo):
    abe, ab,  ae, be, a, b, e =  models_combo    
    ## ALPHA-BETA-ETA
    if(abe):
        print("### Runing alpha, beta and eta test ### ")
        mflags = [True, True, True]       
    ## ALPHA-BETA
    if(ab):
        print("### Runing alpha and beta test ### ")
        mflags = [True, True, False] ##alpha,beta,eta
    ## ALPHA-ETA
    if(ae):
        print("### Runing alpha and eta test ### ")
        mflags = [True, False, True]
    ## BETA-ETA
    if(be):
        print("### Runing beta and eta test ### ")
        mflags = [False, True, True]
    ## ALPHA
    if(a):
        print("### Runing alpha test ### ")
        mflags = [True, False, False]
    if(b):
        print("### Runing beta test ### ")
        mflags = [False, True, False] 
    ## Eta
    if(e):
        print("### Runing eta test ### ")
        mflags = [False, False, True]
    return mflags
        
def run_allscales(args):
    import numpy as np
    import sys; sys.path.append(".")
    from abe_test import RUNTEST
    from src.readfits import  read_rhos, read_taus
    from src.maxlikelihood import percentiles
    from astropy.io import fits

    ##ouputfile format
    names = ['zbin','THETA', 'a_l', 'a_c', 'a_r', 'b_l', 'b_c', 'b_r', 'e_l', 'e_c', 'e_r']
    forms = ['i4',  'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8','f8', 'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])

    
    rhosfile =  args.rhos
    taus = args.taus
    nwalkers = args.nwalkers
    nsteps = args.nsteps
    eq = args.eq
    mflags = getmflags([args.abe, args.ab, args.ae, args.be, args.a, args.b, args.e])
    nsig = args.nsig
    
    i_guess0 = [ 0, 1, -1 ] #fiducial values
    i_guess = np.array(i_guess0)[np.array(mflags)].tolist()

    for i,  taufile in enumerate(taus):
        zbinarr =  np.array([i]*nrows)
        ac = []; al = []; ar = []
        bc = []; bl = []; br = []
        ec = []; el = []; er = []
        data = {}
        for i in range(20):
            data['rhos'] = read_rhos(rhosfile, maxbin=i)[1]
            data['cov_rhos'] = read_rhos(rhosfile, maxbin = i)[2]
            data['taus'] = read_taus(taufile, maxbin = i)[1]
            data['cov_taus'] = read_taus(taufile, maxbin = i)[2]

            samples, chains = RUNTEST(i_guess, data, nwalkers, nsteps, eq=eq,
                             mflags=mflags, xip=True, xim=True ,
                             moderr=False, uwmprior=False,
                             minimize= True, margin=True,
                                 overall=False)
            mcmcpars = percentiles(samplesp, nsig=nsig)
            fillparlists(al, ac, ar, bl, bc, br, el, ec, er, mcmcpars, mflags)
        array_list = [zbinarr, read_rhos(rhosfile, maxbin=i)[0], np.array(al), np.array(ac), np.array(ar),
                      np.array(bl), np.array(bc), np.array(br), np.array(el), np.array(ec), np.array(er)]
        for array, name in zip(array_list, names): outdata[name] = array
        corrhdu = fits.BinTableHDU(outdata, name='zbin%d'%(zbin))
        hdul.insert(i, corrhdu)

    filename = os.path.join(outpath, args.filename)
    print("Printin file:", filename)
    hdul.writeto(filename, overwrite=True)

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

    run_allscales(args)

    run_parspatch(outpath, args.filename, args.rhos, args.taus, nsig=nsig,  njk=njk, mflags=mflags)
    plotlineal('', outpath + args.filename,mflags, njk)
    


    
if __name__ == "__main__":
    main()
