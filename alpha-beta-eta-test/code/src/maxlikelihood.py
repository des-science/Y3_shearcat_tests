#Log natural of the prior
#import matplotlib.pyplot as plt
import numpy as np
#plt.style.use('SVA1StyleSheet.mplstyle')
 
def logprior(pars, mflags=[True, True, True],  uwmprior=False):
    aflag, bflag, eflag =  mflags
    if(uwmprior):
        al = - 2; au =  2
        bl = - 1; bu =  3
        el =- 3; eu = 1
    else:
        l = 1000
        al = -l; au =  l
        bl = -l; bu =  l
        el = -l; eu = l
    alpha_min, beta_min, eta_min  = al,bl,el
    alpha_max, beta_max, eta_max  = au,bu,eu
    if(aflag and bflag and eflag):
        alpha, beta, eta = pars
        if ( alpha_min<alpha< alpha_max)and( beta_min < beta< beta_max)and( eta_min<eta<eta_max): return 0.0
        return -np.inf
    elif( aflag and bflag and (not eflag)):
        alpha, beta = pars
        if (alpha_min<alpha<alpha_max)and( beta_min< beta< beta_max): return 0.0
        return -np.inf
    elif(aflag and (not bflag) and eflag ):
        alpha, eta = pars
        if ( alpha_min<alpha< alpha_max)and( eta_min<eta<eta_max): return 0.0
        return -np.inf
    elif( (not aflag) and bflag and eflag ):
        beta, eta = pars
        if ( beta_min < beta< beta_max)and( eta_min<eta<eta_max): return 0.0
        return -np.inf
    elif( aflag and (not bflag) and (not eflag)):
        alpha = pars
        if ( alpha_min < alpha < alpha_max): return 0.0
        return -np.inf
    elif( (not aflag) and bflag and (not eflag)):
        beta = pars
        if ( beta_min < beta< beta_max): return 0.0
        return -np.inf
    elif( (not aflag) and (not bflag) and eflag):
        eta = pars
        if eta_min<eta<eta_max: return 0.0
        return -np.inf

##Log natural of the likelihood function. Gaussian.
def loglike(chisq):
    return -0.5*chisq
##Log natural of the posterior
def logpost(pars, data, eq=None, mflags=[True, True, True], xip=True, xim=False, moderr=False, uwmprior=False):
    from src.chi2 import CHI2
    chisq = CHI2(pars, data,eq=eq, mflags=mflags, xip=xip, xim=xim,  moderr=moderr )
    lp = logprior(pars, mflags=mflags, uwmprior=uwmprior)
    if not np.isfinite(lp):
        return -np.inf
    return lp + loglike(chisq)
def corner_plot(samples, labels, title):
    import corner
    import numpy as np
    burn = 5000
    samples_burned = np.c_[[par[burn:] for par in samples]]
    fig = corner.corner(samples_burned.T, labels=labels,
                        quantiles=[0.16, 0.5, 0.84],  #-1sigma,0sigma,1sigma
                        levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-9./2)), #1sigma, 2sigma and 3sigma contours
                        show_titles=True, title_kwargs={"fontsize": 12}, title_fmt= '.4f', 
                        smooth1d=None, plot_contours=True,  
                        no_fill_contours=False, plot_density=True, use_math_text=True, )
    print("Printing file:",  title)
    plt.savefig(title)
    plt.close(fig)
    print(title, "Printed")
def MCMC(best_pars,data, nwalkers=50, nsteps=1000, eq=None,
         mflags=[True, True, True], xip=True, xim=False, moderr=False,
         uwmprior=False):
    import emcee
    import itertools
    ndim =  len(list(itertools.compress(range(len(mflags)),  mflags)))
    pos = [best_pars + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost, threads=1,
                                    args=(data, eq, mflags,xip,xim, moderr, uwmprior))
    print("Runing MCMC ...")
    sampler.run_mcmc(pos, nsteps)
    print("Run finished")
    
    if(ndim ==1):
        alpha_chain = sampler.chain[:,:,0]; alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))
        samples = np.c_[alpha_chain_flat].T
        chains = [alpha_chain]
    elif(ndim ==2):
        alpha_chain = sampler.chain[:,:,0]; alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))
        beta_chain = sampler.chain[:,:,1]; beta_chain_flat = np.reshape(beta_chain, (nwalkers*nsteps,))
        samples = np.c_[alpha_chain_flat, beta_chain_flat].T
        chains = [alpha_chain, beta_chain]
    elif(ndim==3):
        alpha_chain = sampler.chain[:,:,0]; alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))
        beta_chain = sampler.chain[:,:,1]; beta_chain_flat = np.reshape(beta_chain, (nwalkers*nsteps,))
        eta_chain = sampler.chain[:,:,2]; eta_chain_flat = np.reshape(eta_chain, (nwalkers*nsteps,))
        samples = np.c_[alpha_chain_flat, beta_chain_flat, eta_chain_flat].T
        chains = [alpha_chain, beta_chain, eta_chain]
    sampler.reset()
    return samples, chains

def bestparameters(samples):
    allpars = []
    for i in range (len(samples)):
        par = np.percentile(samples[i], [50]);
        allpars.append(par[0])
    return allpars
            
def percentiles(samples, nsig=1):
    allpars_percent_list = []
    for i in range (0, len(samples)):
        if (nsig==1):
            a_perc = np.percentile(samples[i], [16, 50, 84]); par_perc_list =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
            allpars_percent_list.append(par_perc_list)
        elif(nsig==2):
            a_perc = np.percentile(samples[i], [2.3, 50, 97.7] ); par_perc_list =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
            allpars_percent_list.append(par_perc_list)
    
    return allpars_percent_list


        
        




   
        
