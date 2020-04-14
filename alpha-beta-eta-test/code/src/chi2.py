def selectipars(params, mflags):
    aflag, bflag, eflag =  mflags
    if(aflag and bflag and eflag):
        alpha, beta, eta = params
    elif( aflag and bflag and (not eflag)):
        alpha, beta = params; eta = 0
    elif(aflag and (not bflag) and eflag ):
        alpha, eta = params; beta = 0
    elif( (not aflag) and bflag and eflag ):
        beta, eta = params; alpha = 0
    elif( aflag and (not bflag) and (not eflag)):
        alpha = params; beta = 0; eta = 0
    elif( (not aflag) and bflag and (not eflag)):
        alpha = 0; beta = params; eta = 0
    elif( (not aflag) and (not bflag) and eflag):
        alpha = 0; beta = 0; eta = params
    return alpha, beta, eta
 
def modelvector(pars, rhos, eq=None, mflags=[True, True, True], xip=True, xim=True):
    import fitsio
    import numpy as np
    alpha, beta, eta = selectipars(pars, mflags)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m =  rhos


    mvec0p = alpha*rho0p + beta*rho2p + eta*rho5p
    mvec1p = alpha*rho2p + beta*rho1p + eta*rho4p 
    mvec2p = alpha*rho5p + beta*rho4p + eta*rho3p

    mvec0m = alpha*rho0m + beta*rho2m + eta*rho5m
    mvec1m = alpha*rho2m + beta*rho1m + eta*rho4m 
    mvec2m = alpha*rho5m + beta*rho4m + eta*rho3m

    if(eq==0):
        if(xip and not xim):
            return mvec0p
        elif(xim and not xip):
            return mvec0m
        else:
            return np.concatenate([mvec0p, mvec0m])
    elif(eq==1):
        if(xip and not xim):
            return mvec1p
        elif(xim and not xip):
            return mvec1m
        else:
            return np.concatenate([mvec1p, mvec1m])
    elif(eq==2):
        if(xip and not xim):
            return mvec2p
        elif(xim and not xip):
            return mvec2m
        else:
            return np.concatenate([mvec2p, mvec2m])
    elif(eq==10):
        if(xip and not xim):
            return np.concatenate([mvec0p, mvec1p])
        elif(xim and not xip):
            return np.concatenate([mvec0m, mvec1m])
        else:
            return np.concatenate([mvec0p, mvec0m, mvec1p, mvec1m])
    elif(eq==20):
        if(xip and not xim):
            return np.concatenate([mvec0p, mvec2p])
        elif(xim and not xip):
            return np.concatenate([mvec0m, mvec2m])
        else:
            return np.concatenate([mvec0p, mvec0m, mvec2p, mvec2m])
    elif(eq==12 or eq==21):
        if(xip and not xim):
            return np.concatenate([mvec1p, mvec2p])
        elif(xim and not xip):
            return np.concatenate([mvec1m, mvec2m])
        else:
            return np.concatenate([mvec1p, mvec1m, mvec2p, mvec2m])   
    else:
        if(xip and not xim):
            return np.concatenate([mvec0p, mvec1p, mvec2p])
        elif(xim and not xip):
            return np.concatenate([mvec0m, mvec1m, mvec2m])
        else:
            return np.concatenate([mvec0p, mvec0m, mvec1p, mvec1m, mvec2p, mvec2m])  
        
def modelcov(pars, cov_rhos, nrows,  eq =None,mflags=[True, True, True], xip=True, xim=True):
    import numpy as np
    alpha, beta, eta = selectipars(pars, mflags)
    covmat = cov_rhos
    #mvar0 = (alpha**2)*covrhos[0]+(beta**2)*covrhos[2]+ (eta**2)*covrhos[5]
    #mvar1 = (alpha**2)*covrhos[2] +(beta**2)*covrhos[1] + (eta**2)*covrhos[4] 
    #mvar2 = (alpha**2)*covrhos[5] +(beta**2)*covrhos[4] + (eta**2)*covrhos[3]   
    if(eq==0):
        if(xip and not xim):
            return covmat[0:nrows,0:nrows]
        elif(xim and not xip):
            return covmat[nrows:2*nrows,nrows:2*nrows] 
        else:
            return covmat[0:2*nrows,0:2*nrows]
    elif(eq==1):
        if(xip and not xim):
            return covmat[2*nrows:3*nrows,2*nrows:3*nrows ] 
        elif(xim and not xip):
            return covmat[3*nrows:4*nrows,3*nrows:4*nrows ] 
        else:
            return covmat[2*nrows:4*nrows,2*nrows:4*nrows ] 
    elif(eq==2):
        if(xip and not xim):
            return covmat[4*nrows:5*nrows,4*nrows:5*nrows ] 
        elif(xim and not xip):
            return covmat[5*nrows:6*nrows,5*nrows:6*nrows ]
        else:
            return covmat[4*nrows:6*nrows,4*nrows:6*nrows ] 
    elif(eq==10):
        if(xip and not xim):
            idx = np.concatenate([np.arange(0,nrows), np.arange(2*nrows, 3*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(nrows, 2*nrows), np.arange(3*nrows, 4*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            return covmat[0:4*nrows,0:4*nrows ] 
    elif(eq==20):
        if(xip and not xim):
            idx = np.concatenate([np.arange(0, nrows), np.arange(4*nrows, 5*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(nrows, 2*nrows), np.arange(5*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            idx = np.concatenate([np.arange(0, nrows), np.arange(4*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
    elif(eq==12 or eq==21):
        if(xip and not xim):
            idx = np.concatenate([np.arange(2*nrows, 3*nrows), np.arange(4*nrows, 5*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(3*nrows, 4*nrows), np.arange(5*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            return covmat[2:6*nrows,2:6*nrows ]   
    else:
        if(xip and not xim):
            idx = np.concatenate([np.arange(0, nrows), np.arange(2*nrows, 3*nrows), np.arange(4*nrows, 5*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(nrows, 2*nrows), np.arange(3*nrows, 4*nrows), np.arange(5*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            return covmat
   
def datavector( taus, eq=None, xip=True, xim=True):
    import fitsio
    import numpy as np
    tau0p, tau0m, tau2p, tau2m, tau5p, tau5m =  taus
    if(eq==0):
        if(xip and not xim):
            return tau0p
        elif(xim and not xip):
            return tau0m
        else:
            return np.concatenate([tau0p, tau0m])
    elif(eq==1):
        if(xip and not xim):
            return tau2p
        elif(xim and not xip):
            return tau2m
        else:
            return np.concatenate([tau2p, tau2m])
    elif(eq==2):
        if(xip and not xim):
            return tau5p
        elif(xim and not xip):
            return tau5m
        else:
            return np.concatenate([tau5p, tau5m])
    elif(eq==10):
        if(xip and not xim):
            return np.concatenate([tau0p, tau2p])
        elif(xim and not xip):
            return np.concatenate([tau0m, tau2m])
        else:
            return np.concatenate([tau0p, tau0m, tau2p, tau2m])
    elif(eq==20):
        if(xip and not xim):
            return np.concatenate([tau0p, tau5p])
        elif(xim and not xip):
            return np.concatenate([tau0m, tau5m])
        else:
            return np.concatenate([tau0p, tau0m, tau5p, tau5m])
    elif(eq==12 or eq == 21):
        if(xip and not xim):
            return np.concatenate([tau2p, tau5p])
        elif(xim and not xip):
            return np.concatenate([tau2m, tau5m])
        else:
            return np.concatenate([tau2p, tau2m, tau5p, tau5m])   
    else:
        if(xip and not xim):
            return np.concatenate([tau0p, tau2p, tau5p])
        elif(xim and not xip):
            return np.concatenate([tau0m, tau2m, tau5m])
        else:
            return np.concatenate([tau0p, tau0m, tau2p, tau2m, tau5p, tau5m])  
def datacov(cov_taus, nrows, eq=None, xip=True, xim=True):
    import numpy as np
    covmat =  cov_taus
    if(eq==0):
        if(xip and not xim):
            return covmat[0:nrows,0:nrows]
        elif(xim and not xip):
            return covmat[nrows:2*nrows,nrows:2*nrows] 
        else:
            return covmat[0:2*nrows,0:2*nrows]
    elif(eq==1):
        if(xip and not xim):
            return covmat[2*nrows:3*nrows,2*nrows:3*nrows ] 
        elif(xim and not xip):
            return covmat[3*nrows:4*nrows,3*nrows:4*nrows ] 
        else:
            return covmat[2*nrows:4*nrows,2*nrows:4*nrows ] 
    elif(eq==2):
        if(xip and not xim):
            return covmat[4*nrows:5*nrows,4*nrows:5*nrows ] 
        elif(xim and not xip):
            return covmat[5*nrows:6*nrows,5*nrows:6*nrows ]
        else:
            return covmat[4*nrows:6*nrows,4*nrows:6*nrows ] 
    elif(eq==10):
        if(xip and not xim):
            idx = np.concatenate([np.arange(0,nrows), np.arange(2*nrows, 3*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(nrows, 2*nrows), np.arange(3*nrows, 4*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            return covmat[0:4*nrows,0:4*nrows ] 
    elif(eq==20):
        if(xip and not xim):
            idx = np.concatenate([np.arange(0, nrows), np.arange(4*nrows, 5*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(nrows, 2*nrows), np.arange(5*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            idx = np.concatenate([np.arange(0, 2*nrows), np.arange(4*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
    elif(eq==12 or eq==21):
        if(xip and not xim):
            idx = np.concatenate([np.arange(2*nrows, 3*nrows), np.arange(4*nrows, 5*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(3*nrows, 4*nrows), np.arange(5*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            return covmat[2*nrows:6*nrows,2*nrows:6*nrows ]   
    else:
        if(xip and not xim):
            idx = np.concatenate([np.arange(0, nrows), np.arange(2*nrows, 3*nrows), np.arange(4*nrows, 5*nrows)])
            return covmat[idx,:][:,idx] 
        elif(xim and not xip):
            idx = np.concatenate([np.arange(nrows, 2*nrows), np.arange(3*nrows, 4*nrows), np.arange(5*nrows, 6*nrows)])
            return covmat[idx,:][:,idx] 
        else:
            return covmat

def CHI2(pars, data, eq=None, mflags=[True, True, True], xip=True, xim=True,  moderr=False):
    taus = data['taus']
    rhos = data['rhos']
    cov_taus = data['cov_taus']
    cov_rhos = data['cov_rhos']
    dvect =  datavector(taus, eq=eq, xip=xip, xim=xim)
    dcov_mat = datacov(cov_taus, len(taus[0]),  eq=eq, xip=xip, xim=xim)
    mvect =  modelvector(pars, rhos,  eq=eq, mflags=mflags, xip=xip, xim=xim)
    mcov_mat = modelcov(pars, cov_rhos,  len(rhos[0]),  eq=eq, mflags=mflags, xip=xip, xim=xim)
    val=chi2(mvect, dvect, mcov_mat, dcov_mat, moderr=moderr )
    return val

def chi2(modelvec, datavec,  covmodel, covdata,  moderr=False ):
    import numpy as np
    d =  np.array([modelvec - datavec])
    if(moderr):
        cov_inv = np.linalg.inv(covdata + covmodel)
        print("ERROR")
    else:
        #cov_inv = np.linalg.pinv(covdata,  rcond=1e-10)
        cov_inv = np.linalg.pinv(covdata)
        
    chisq = np.dot(np.dot(d,cov_inv), d.T)
    #print('chisq:', chisq[0][0])
    #print('np.dot(d,d.T):',  np.dot(d, d.T)[0][0])
    return chisq[0][0]
    
def ndof(data, eq=None, mflags=[True, True, True], xip=True, xim=True):
    import itertools
    npoints = len(data['rhos'][0])
    ndim =  len(list(itertools.compress(range(len(mflags)),  mflags)))
    if (eq ==  0 or eq==1 or eq==2) :
        npoints *=1 
    elif(eq==21 or eq==20 or eq==10 or eq==12):
        npoints *=2
    else:
        npoints *=3
    if(xip and xim):
        npoints *= 2
    dof = npoints- ndim
    return dof
def chi2nu(pars, data, eq=None, mflags=[True, True, True], xip=True, xim=True,  moderr=False):
    dof = ndof(data, eq=eq, mflags=mflags, xip=xip, xim=xim)
    return  CHI2(pars, data, eq=eq, mflags=mflags, xip=xip, xim=xim,  moderr=moderr)/dof
    
def minimizeCHI2(data, initial_guess, eq=None,  mflags=[True, True, True], xip=True, xim=False, moderr=False):
    import scipy.optimize as optimize
    result = optimize.minimize(CHI2, initial_guess,args=(data,eq,mflags, xip, xim, moderr), method='Nelder-Mead', tol=1e-6)

    dof = ndof(data, eq=eq, mflags=mflags, xip=xip, xim=xim)
    print("Degree of freedom", dof)
    print("chi2", result.fun)
    if result.success:
        fitted_params = result.x
        return fitted_params, result.fun/dof
    else:
        raise ValueError(result.message)
