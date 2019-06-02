import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def pretty_rho(meanr, rho, sig,  legend=None, lfontsize=24, color='black', marker='o', ylabel=r'$\rho(\theta)$',title=None,  xlim=None,  ylim=None):
    '''
    plt.plot(meanr, rho, color=color, label=legend, marker=marker)
    plt.plot(meanr, -rho, color=color, ls=':', marker=marker)
    if sig is not None:
        plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker)
        plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker)
        rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color=color, marker=marker)
    '''
    plt.plot(meanr, rho, color=color, label=legend)
    plt.plot(meanr, -rho, color=color, ls=':')
    #rigth quadrants
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker,  capsize=2)
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker,  capsize=2)
    #leftquadrants
    plt.errorbar( -meanr, rho, yerr=sig, color=color,  marker='^',  capsize=2)
    plt.errorbar( -meanr,-rho, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
    plt.legend(loc='best', fontsize=lfontsize)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    plt.tight_layout()

def pretty_rho0(meanr, rho, sig, title= None,  xlim=None, ylim=None):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    rho0_line = plt.errorbar(-meanr,- rho, yerr=sig, color='blue', marker='o', ls=':')
    plt.legend([rho0_line],[r'$\rho_0(\theta)$'],loc='upper right', fontsize=24)
    
    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    plt.tight_layout()

def pretty_rho1(meanr, rho, sig, rho3=None, sig3=None, rho4=None, sig4=None, title= None,xlim=None, ylim=None):
    import matplotlib.patches as mp
    
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho1_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    #rho1_line = plt.errorbar(-meanr,-rho, yerr=sig, color='blue', marker='o', ls=':')
    if rho3 is not None:
        plt.plot(meanr*1.03, rho3, color='green')
        plt.plot(meanr*1.03, -rho3, color='green', ls=':')
        plt.errorbar(meanr[rho3>0]*1.03, rho3[rho3>0], yerr=sig3[rho3>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho3<0]*1.03, -rho3[rho3<0], yerr=sig3[rho3<0], color='green', ls='', marker='s')
        rho3_line = plt.errorbar(-meanr, rho3, yerr=sig3, color='green', marker='s')
        #rho3_line = plt.errorbar(-meanr,-rho3, yerr=sig3, color='green', marker='s', ls=':')
    if rho4 is not None:
        plt.plot(meanr*1.06, rho4, color='red')
        plt.plot(meanr*1.06, -rho4, color='red', ls=':')
        plt.errorbar(meanr[rho4>0]*1.06, rho4[rho4>0], yerr=sig4[rho4>0], color='red', ls='', marker='^')
        plt.errorbar(meanr[rho4<0]*1.06, -rho4[rho4<0], yerr=sig4[rho4<0], color='red', ls='', marker='^')
        rho4_line = plt.errorbar(-meanr, rho4, yerr=sig4, color='red', marker='^')
        #rho4_line = plt.errorbar(-meanr,-rho4, yerr=sig4, color='red', marker='^', ls=':')

    plt.legend([rho1_line, rho3_line, rho4_line],
               [r'$\rho_1(\theta)$', r'$\rho_3(\theta)$',
                r'$\rho_4(\theta)$'], loc='upper right',
               fontsize=24)

    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig, rho2=None, sig2=None, rho5=None, sig5=None, tauleg=False, title= None, xlim=None, ylim=None):
 
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    #rho2_line = plt.errorbar(-meanr,-rho, yerr=sig, color='blue', marker='o', ls=':')
    if rho2 is not None:
        plt.plot(meanr*1.03, rho2, color='green')
        plt.plot(meanr*1.03, -rho2, color='green', ls=':')
        plt.errorbar(meanr[rho2>0]*1.03, rho2[rho2>0], yerr=sig2[rho2>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho2<0]*1.03, -rho2[rho2<0], yerr=sig2[rho2<0], color='green', ls='', marker='s')
        rho2_line = plt.errorbar(-meanr, rho2, yerr=sig2, color='green', marker='s')
        #rho2_line = plt.errorbar(-meanr,-rho2, yerr=sig2, color='green', marker='s', ls=':')
    if rho5 is not None:
        plt.plot(meanr*1.03, rho5, color='red')
        plt.plot(meanr*1.03, -rho5, color='red', ls=':')
        plt.errorbar(meanr[rho5>0]*1.03, rho5[rho5>0], yerr=sig5[rho5>0], color='red', ls='', marker='^')
        plt.errorbar(meanr[rho5<0]*1.03, -rho5[rho5<0], yerr=sig5[rho5<0], color='red', ls='', marker='^')
        rho5_line = plt.errorbar(-meanr, rho5, yerr=sig5, color='red', marker='^')
        #rho5_line = plt.errorbar(-meanr,-rho5, yerr=sig5, color='red', marker='^', ls=':')
    if tauleg:
        plt.legend([rho0_line, rho2_line, rho5_line],
                   [r'$\tau_0(\theta)$', r'$\tau_2(\theta)$',
                    r'$\tau_5(\theta)$'], loc='upper right',
                   fontsize=24)
        plt.ylabel(r'$\tau(\theta)$', fontsize=24)
    else:
        plt.legend([rho0_line, rho2_line, rho5_line],
                   [r'$\rho_0(\theta)$', r'$\rho_2(\theta)$',
                    r'$\rho_5(\theta)$'], loc='upper right', fontsize=24)
        plt.ylabel(r'$\rho(\theta)$', fontsize=24)

    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    plt.tight_layout()

def plotcorrmat(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    plt.colorbar()
    plt.tight_layout()
    
def plotallrhosfits(stat_file, outpath, title= None, xlim=None, ylims=None):
    import numpy as np
    from readfits import read_rhos
    ylim0p, ylim0m, ylim1p, ylim1m = [None, None,  None, None]
    if ylims is not None: ylim0p, ylim0m, ylim1p, ylim1m = ylims
    meanr, rhos,  cov_rhos =  read_rhos(stat_file)
    rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m, rho4p, rho4m, rho5p, rho5m =  rhos
    cov0p, cov0m, cov1p, cov1m, cov2p, cov2m, cov3p, cov3m, cov4p, cov4m, cov5p, cov5m =  cov_rhos
    sig_rho0p =  np.sqrt(np.diag(cov0p)); sig_rho0m =  np.sqrt(np.diag(cov0m))
    sig_rho1p =  np.sqrt(np.diag(cov1p)); sig_rho1m =  np.sqrt(np.diag(cov1m))
    sig_rho2p =  np.sqrt(np.diag(cov2p)); sig_rho2m =  np.sqrt(np.diag(cov2m))
    sig_rho3p =  np.sqrt(np.diag(cov3p)); sig_rho3m =  np.sqrt(np.diag(cov3m))
    sig_rho4p =  np.sqrt(np.diag(cov4p)); sig_rho4m =  np.sqrt(np.diag(cov4m))
    sig_rho5p =  np.sqrt(np.diag(cov5p)); sig_rho5m =  np.sqrt(np.diag(cov5m))
    plt.clf()
    pretty_rho1(meanr, rho1p, sig_rho1p, rho3p, sig_rho3p, rho4p, sig_rho4p,title=title, xlim=xlim,ylim=ylim0p)
    print("Printing file: ", outpath +'rho1p_all_rsrs.png')
    plt.savefig(outpath +'rho1p_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr,rho0p, sig_rho0p ,rho2p, sig_rho2p,  rho5p, sig_rho5p,title=title, xlim=xlim,ylim=ylim1p)
    print("Printing file: ", outpath +'rho2p_all_rsrs.png')
    plt.savefig(outpath +'rho2p_all_rsrs.png')
    plt.clf()
    pretty_rho1(meanr, rho1m, sig_rho1m, rho3m, sig_rho3m, rho4m, sig_rho4m,title=title, xlim=xlim,ylim=ylim0m)
    print("Printing file: ", outpath +'rho1m_all_rsrs.png')
    plt.savefig(outpath +'rho1m_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr,rho0m, sig_rho0m ,rho2m, sig_rho2m,  rho5m, sig_rho5m,title=title, xlim=xlim,ylim=ylim1m)
    print("Printing file: ", outpath +'rho2m_all_rsrs.png')
    plt.savefig(outpath +'rho2m_all_rsrs.png')

    titles =  [r'$\rho_{0+}(\theta)$', r'$\rho_{0-}(\theta)$',r'$\rho_{1+}(\theta)$',r'$\rho_{1-}(\theta)$',
               r'$\rho_{2+}(\theta)$', r'$\rho_{2-}(\theta)$',r'$\rho_{3+}(\theta)$',r'$\rho_{3-}(\theta)$',
               r'$\rho_{4+}(\theta)$', r'$\rho_{4-}(\theta)$',r'$\rho_{5+}(\theta)$',r'$\rho_{5-}(\theta)$' ]
    names = ['rho0p_covmat.png', 'rho0m_covmat.png', 'rho1p_covmat.png', 'rho1m_covmat.png', 
             'rho2p_covmat.png', 'rho2m_covmat.png', 'rho3p_covmat.png', 'rho3m_covmat.png', 
             'rho4p_covmat.png', 'rho4m_covmat.png', 'rho5p_covmat.png', 'rho5m_covmat.png']
    for i,covmat in enumerate(cov_rhos):
        plt.clf()
        plt.title(titles[i])
        plotcorrmat(covmat)
        plt.savefig(outpath + names[i], dpi=500)
        print(outpath +names[i], 'Printed!')
    
def plotalltausfits(stat_file, outpath, title= None, xlim=None,  ylims=None,  zbin=''):
    import numpy as np
    from readfits import read_taus
    import fitsio
    ylim0p, ylim0m = [None, None]
    if ylims is not None: ylim0p, ylim0m = ylims
    meanr, taus,  cov_taus =  read_taus(stat_file)
    tau0p, tau0m, tau2p, tau2m, tau5p, tau5m =  taus
    cov0p, cov0m, cov2p, cov2m, cov5p, cov5m =  cov_taus
    sig_tau0p =  np.sqrt(np.diag(cov0p)); sig_tau0m =  np.sqrt(np.diag(cov0m))
    sig_tau2p =  np.sqrt(np.diag(cov2p)); sig_tau2m =  np.sqrt(np.diag(cov2m))
    sig_tau5p =  np.sqrt(np.diag(cov5p)); sig_tau5m =  np.sqrt(np.diag(cov5m))

    plt.clf()
    pretty_rho2(meanr, tau0p, sig_tau0p, tau2p, sig_tau2p, tau5p, sig_tau5p, tauleg=True, title=title, xlim=xlim, ylim=ylim0p)
    name = outpath +'taup_all_rsrs' + zbin +  '.png'
    print("Printing file: ", name)
    plt.savefig(name)
    plt.clf()
    pretty_rho2(meanr, tau0m, sig_tau0m, tau2m, sig_tau2m, tau5m, sig_tau5m, tauleg=True, title=title, xlim=xlim, ylim=ylim0m)
    name = outpath +'taum_all_rsrs' + zbin +  '.png'
    print("Printing file: ", name)
    plt.savefig(name)

    '''
    titles =  [r'$\tau_{0+}(\theta)$', r'$\tau_{0-}(\theta)$',
               r'$\tau_{2+}(\theta)$', r'$\tau_{2-}(\theta)$',
               r'$\tau_{5+}(\theta)$',r'$\tau_{5-}(\theta)$' ]
    names = ['tau0p_covmat.png', 'tau0m_covmat.png', 
             'tau2p_covmat.png', 'tau2m_covmat.png', 
             'tau5p_covmat.png', 'tau5m_covmat.png']
    for i,covmat in enumerate(cov_taus):
        plt.clf()
        plt.title(titles[i])
        plotcorrmat(covmat)
        plt.savefig(outpath + names[i], dpi=500)
        print(outpath +names[i], 'Printed!')
    '''

    plt.clf()
    lengths = [len(tau0p),len(tau0m),len(tau2p),len(tau2p), len(tau5p),  len(tau5m)]
    covmatrixfit=fitsio.read(stat_file, ext=1)
    plotcorrmat(covmatrixfit)
    plt.title(r'$\tau_{0+}(\theta) \mid \tau_{0-}(\theta) \mid \tau_{2+}(\theta) \mid \tau_{2-}(\theta) \mid \tau_{5+}(\theta) \mid \tau_{5-}(\theta) $')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = outpath + 'CovariancematrixTaus' + zbin + '.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

def plotalltauscorrmatfits(filenames, outpath):
    import fitsio
    names = ['tau0_covmat.png', 'tau2_covmat.png', 'tau5_covmat.png']
    titles =  [r'$\rho_{0}(\theta)$', r'$\rho_{2}(\theta)$', r'$\rho_{5}(\theta)$' ]
    for i,f in enumerate(filenames):
        covmat = fitsio.read(f, ext=1)
        plt.clf()
        plotcorrmat(covmat)
        plt.title(titles[i])
        plt.savefig(outpath + names[i], dpi=500)
        print(outpath +names[i], 'Printed!')
def corner_plot(samples, labels, title):
    import corner
    import numpy as np
    #burn = 5000
    samples= np.c_[[par[int(0.2 * len(par)):] for par in samples]].T
    fig = corner.corner(samples, labels=labels,
                        quantiles=[0.16, 0.5, 0.84],  #-1sigma,0sigma,1sigma
                        levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-9./2)), #1sigma, 2sigma and 3sigma contours
                        show_titles=True, title_kwargs={"fontsize": 12}, title_fmt= '.4f', 
                        smooth1d=None, plot_contours=True,  
                        no_fill_contours=False, plot_density=True, use_math_text=True, )
    print("Printing file:",  title)
    plt.tight_layout()
    plt.savefig(title)
    plt.close(fig)
    print(title, "Printed")
def plot_samplesdist(samples, chains, mflags, nwalkers, nsteps,  namemc, namecont):
    import numpy as np
    import emcee
    aflag, bflag, eflag =  mflags
    fig = plt.figure(figsize=(16, 12))
    ndim =  len(samples)        
    if(ndim ==1):
         axs = fig.subplots(1, 3)
         if( aflag and (not bflag) and (not eflag)):
             labels = [r"$\alpha$"]
         elif( (not aflag) and bflag and (not eflag)):
             labels = [r"$\beta$"]
         elif( (not aflag) and (not bflag) and eflag):
             labels = [r"$\eta$"]
         axs[0].set_xlabel("Ensemble step")
         axs[1].set_xlabel("Ensemble step")
         axs[2].set_xlabel("Walker Step")
         axs[0].set_title("Ensemble dispersion")
         axs[1].set_title("Ensemble autocorrelation")
         axs[2].set_title("Walker mean and stdev")
         alpha_chain = chains
         alpha_chain_mean = np.mean(alpha_chain, axis=0)
         alpha_chain_err = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
         idx = np.arange(len(alpha_chain_mean))
         axs[2].errorbar(x=idx, y=alpha_chain_mean,
                         yerr=alpha_chain_err, errorevery=50,
                         ecolor='red', lw=0.5, elinewidth=2.,
                         color='k')
    elif(ndim ==2):
        axs = fig.subplots(2, 3)
        if( aflag and bflag and (not eflag)):
            labels = [r"$\alpha$", r"$\beta$"]
        elif(aflag and (not bflag) and eflag ):
            labels = [r"$\alpha$", r"$\eta$"]
        elif( (not aflag) and bflag and eflag ):
            labels = [r"$\beta$", r"$\eta$"]
        axs[1][0].set_xlabel("Ensemble step")
        axs[1][1].set_xlabel("Ensemble step")
        axs[1][2].set_xlabel("Walker Step")
        axs[0][0].set_title("Ensemble dispersion")
        axs[0][1].set_title("Ensemble autocorrelation")
        axs[0][2].set_title("Walker mean and stdev")
        alpha_chain, beta_chain = chains
        alpha_chain_mean = np.mean(alpha_chain, axis=0)
        alpha_chain_err = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
        beta_chain_mean = np.mean(beta_chain, axis=0)
        beta_chain_err = np.std(beta_chain, axis=0) / np.sqrt(nwalkers)
        idx = np.arange(len(alpha_chain_mean))
        axs[0][2].errorbar(x=idx, y=alpha_chain_mean,
                           yerr=alpha_chain_err, errorevery=50,
                           ecolor='red', lw=0.5, elinewidth=2., color='k')
        axs[1][2].errorbar(x=idx, y=beta_chain_mean,
                           yerr=beta_chain_err, errorevery=50,
                           ecolor='red', lw=0.5, elinewidth=2., color='k');
    elif(ndim==3):
        axs = fig.subplots(3, 3)
        labels = [r"$\alpha$", r"$\beta$", r"$\eta$"]
        axs[2][0].set_xlabel("Ensemble step")
        axs[2][1].set_xlabel("Ensemble step")
        axs[2][2].set_xlabel("Walker Step")
        axs[0][0].set_title("Ensemble dispersion")
        axs[0][1].set_title("Ensemble autocorrelation")
        axs[0][2].set_title("Walker mean and stdev")
        alpha_chain, beta_chain, eta_chain = chains
        alpha_chain_mean = np.mean(alpha_chain, axis=0)
        alpha_chain_err = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
        beta_chain_mean = np.mean(beta_chain, axis=0)
        beta_chain_err = np.std(beta_chain, axis=0) / np.sqrt(nwalkers)
        eta_chain_mean = np.mean(eta_chain, axis=0)
        eta_chain_err = np.std(eta_chain, axis=0) / np.sqrt(nwalkers)
        idx = np.arange(len(alpha_chain_mean))
        axs[0][2].errorbar(x=idx, y=alpha_chain_mean,
                           yerr=alpha_chain_err, errorevery=50,
                           ecolor='red', lw=0.5, elinewidth=2.,
                           color='k')
        axs[1][2].errorbar(x=idx, y=beta_chain_mean,
                           yerr=beta_chain_err, errorevery=50,
                           ecolor='red', lw=0.5, elinewidth=2., color='k');
        axs[2][2].errorbar(x=idx, y=eta_chain_mean,
                           yerr=eta_chain_err, errorevery=50,
                           ecolor='red', lw=0.5, elinewidth=2., color='k');

    for i, par in enumerate(samples):
        axs[i][0].set_ylabel(labels[i])
        idx = np.arange(len(par))
        axs[i][0].scatter(idx, par[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
        # Get selfcorrelation using emcee
        ac = emcee.autocorr.function(par)
        idx = np.arange(len(ac),step=1)
        axs[i][1].scatter(idx, ac[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
        axs[i][1].axhline(alpha=1., lw=1., color='red')

    print("Printing file:",  namemc)
    plt.tight_layout()
    plt.savefig(namemc)
    plt.close(fig)
    print(namemc, "Printed")
    corner_plot(samples, labels, namecont)

    
 
  
