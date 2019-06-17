import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def pretty_rho(meanr, rho, sig,  legend=None, lfontsize=24, color='black', marker='o', ylabel=r'$\rho(\theta)$',title=None,  xlim=None,  ylim=None):
    import numpy as np
    '''
    plt.plot(meanr, rho, color=color, label=legend, marker=marker)
    plt.plot(meanr, -rho, color=color, ls=':', marker=marker)
    if sig is not None:
        plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker)
        plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker)
        rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color=color, marker=marker)
    '''
    
    if sig is None: sig =  np.zeros(len(rho))
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
   #plt.tight_layout()

def pretty_rho1(meanr, rho, sig, rho3=None, sig3=None, rho4=None, sig4=None, mlabel=False,  title= None,xlim=None, ylim=None):
    fsize = 24
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

    if mlabel:
        plt.legend([rho1_line, rho3_line, rho4_line],
                   [r'$\rho_{1-}(\theta)$', r'$\rho_{3-}(\theta)$',
                    r'$\rho_{4-}(\theta)$'], loc='upper right',
                   fontsize=fsize)
        plt.ylabel(r'$\rho_{-}(\theta)$', fontsize=fsize)
    else:
        plt.legend([rho1_line, rho3_line, rho4_line],
                   [r'$\rho_{1+}(\theta)$', r'$\rho_{3+}(\theta)$',
                    r'$\rho_{4+}(\theta)$'], loc='upper right',
                   fontsize=fsize)
        plt.ylabel(r'$\rho_{+}(\theta)$', fontsize=fsize)

    

    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=fsize)  
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig, rho2=None, sig2=None, rho5=None, sig5=None, mlabel=False, title= None, xlim=None, ylim=None):
    fsize = 24
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

    if mlabel:
        plt.legend([rho0_line, rho2_line, rho5_line],
                   [r'$\rho_{0-}(\theta)$', r'$\rho_{2-}(\theta)$',
                    r'$\rho_{5-}(\theta)$'], loc='upper right',
                   fontsize=fsize)
        plt.ylabel(r'$\rho_{-}(\theta)$', fontsize=fsize)
    else:
        plt.legend([rho0_line, rho2_line, rho5_line],
                   [r'$\rho_{0+}(\theta)$', r'$\rho_{2+}(\theta)$',
                    r'$\rho_{5+}(\theta)$'], loc='upper right',
                   fontsize=fsize)
        plt.ylabel(r'$\rho_{+}(\theta)$', fontsize=fsize)

    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if title is not None: plt.title(title)
    #plt.tight_layout()

def pretty_tau2(meanr, rho, sig, rho2=None, sig2=None, rho5=None, sig5=None, mlabel=False, title= None, xlim=None, ylim=None):
    fsize = 24
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

    if mlabel:
        plt.legend([rho0_line, rho2_line, rho5_line],
                   [r'$\tau_{0-}(\theta)$', r'$\tau_{2-}(\theta)$',
                    r'$\tau_{5-}(\theta)$'], loc='upper right',
                   fontsize=fsize)
        plt.ylabel(r'$\tau_{-}(\theta)$', fontsize=fsize)
    else:
        plt.legend([rho0_line, rho2_line, rho5_line],
                   [r'$\tau_{0+}(\theta)$', r'$\tau_{2+}(\theta)$',
                    r'$\tau_{5+}(\theta)$'], loc='upper right',
                   fontsize=fsize)
        plt.ylabel(r'$\tau_{+}(\theta)$', fontsize=fsize)

    


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
    from readfits import read_rhos_plots
    ylim0p, ylim0m, ylim1p, ylim1m = [None, None,  None, None]
    if ylims is not None: ylim0p, ylim0m, ylim1p, ylim1m = ylims
    meanr, rhos,  cov_rhos =  read_rhos_plots(stat_file)
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
    pretty_rho1(meanr, rho1m, sig_rho1m, rho3m, sig_rho3m, rho4m, sig_rho4m,mlabel=True, title=title, xlim=xlim,ylim=ylim0m)
    print("Printing file: ", outpath +'rho1m_all_rsrs.png')
    plt.savefig(outpath +'rho1m_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr,rho0m, sig_rho0m ,rho2m, sig_rho2m,  rho5m, sig_rho5m, mlabel=True, title=title, xlim=xlim,ylim=ylim1m)
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
    from readfits import read_taus_plots
    import fitsio
    ylim0p, ylim0m = [None, None]
    if ylims is not None: ylim0p, ylim0m = ylims
    meanr, taus,  cov_taus =  read_taus_plots(stat_file)
    tau0p, tau0m, tau2p, tau2m, tau5p, tau5m =  taus
    cov0p, cov0m, cov2p, cov2m, cov5p, cov5m =  cov_taus
    sig_tau0p =  np.sqrt(np.diag(cov0p)); sig_tau0m =  np.sqrt(np.diag(cov0m))
    sig_tau2p =  np.sqrt(np.diag(cov2p)); sig_tau2m =  np.sqrt(np.diag(cov2m))
    sig_tau5p =  np.sqrt(np.diag(cov5p)); sig_tau5m =  np.sqrt(np.diag(cov5m))

    plt.clf()
    pretty_tau2(meanr, tau0p, sig_tau0p, tau2p, sig_tau2p, tau5p, sig_tau5p, mlabel=True, title=title, xlim=xlim, ylim=ylim0p)
    name = outpath +'taup_all_rsrs' + zbin +  '.png'
    print("Printing file: ", name)
    plt.savefig(name)
    plt.clf()
    pretty_tau2(meanr, tau0m, sig_tau0m, tau2m, sig_tau2m, tau5m, sig_tau5m, mlabel=False, title=title, xlim=xlim, ylim=ylim0m)
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
    plt.clf()
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

def pretty_residuals(axs, idx, meanr,res,sig,label='Corr',color='black'):
    axs[idx].plot(meanr, res, color=color, label=label)
    axs[idx].plot(meanr, -res, color=color, ls=':')
    axs[idx].errorbar(meanr[res>0], res[res>0], yerr=sig[res>0], color=color, ls='', marker='.',  capsize=2)
    axs[idx].errorbar(meanr[res<0], -res[res<0], yerr=sig[res<0], color=color, ls='', marker='.',  capsize=2)   
    axs[idx].errorbar( -meanr, res, yerr=sig, color=color,  marker='^',  capsize=2)
    axs[idx].errorbar( -meanr,-res, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
    axs[idx].legend(loc='best', fontsize=5) 
    axs[idx].tick_params(axis='both', which='major', labelsize=10)
    axs[idx].set_xlabel(r'$\theta$ (arcmin)', fontsize=7)
    axs[idx].set_ylim(1.e-12,1.e-5 )
    axs[idx].set_xscale('log')
    axs[idx].set_yscale('log', nonposy='clip')

def pretty_residuals_tomo(axs,  meanr,res,sig,label='Corr',color='black'):
    #axs.plot(meanr, res, color=color, label=label)
    #axs.plot(meanr, -res, color=color, ls=':')
    axs.errorbar(meanr, res, yerr=sig, color=color, ls='', marker='.',  capsize=2, label=label)
    #axs.errorbar(meanr[res>0], res[res>0], yerr=sig[res>0], color=color, ls='', marker='.',  capsize=2)
    #axs.errorbar(meanr[res<0], -res[res<0], yerr=sig[res<0], color=color, ls='', marker='.',  capsize=2)   
    #axs.errorbar( -meanr, res, yerr=sig, color=color,  marker='^',  capsize=2)
    #axs.errorbar( -meanr,-res, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
    axs.legend(loc='best', fontsize=5) 
    axs.tick_params(axis='both', which='major', labelsize=15)
    axs.set_xlabel(r'$\theta$ (arcmin)', fontsize=17)
    #axs.set_ylim(-3.e-6,6.e-6 )
    axs.set_xscale('log')
    #axs.tight_layout()
    #axs.set_yscale('log', nonposy='clip')

def plotbestfitresiduals(samplesp, samplesm,  meanr, data, models_combo, plotname,  margin=False, overall=False):
    from maxlikelihood import bestparameters
    import numpy as np
    import matplotlib.gridspec as gridspec

    rhosp =[data['rhos'][2*i] for i in range(6)]; nrows = len(rhosp[0])
    covrhosp =[data['cov_rhos'][2*i*nrows:(2*i + 1)*nrows , 2*i*nrows:(2*i + 1)*nrows] for i in range(6)]
    rhosm =[data['rhos'][2*i + 1] for i in range(6)]
    covrhosm =[data['cov_rhos'][(2*i + 1)*nrows:2*(i + 1)*nrows , (2*i + 1)*nrows:2*(i + 1)*nrows] for i in range(6)]
    tausp =[data['taus'][2*i] for i in range(3)]; nrows = len(rhosp[0])
    covtausp =[data['cov_taus'][2*i*nrows:(2*i + 1)*nrows , 2*i*nrows:(2*i + 1)*nrows] for i in range(3)]
    tausm =[data['taus'][2*i + 1] for i in range(3)]
    covtausm =[data['cov_taus'][(2*i + 1)*nrows:2*(i + 1)*nrows , (2*i + 1)*nrows:2*(i + 1)*nrows] for i in range(3)]
    
    if(overall and not margin):
        a = b = e = 0; 
        am = bm = em = 0; 
        eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
        if not (abe_bool or ab_bool or a_bool): abe_bool = True
        if(abe_bool):
            a, b, e = samplesp
            am, bm, em = samplesm
        if(ab_bool):
            a, b = samplesp
            am, bm = samplesm
        if(ae_bool):
            a, e = samplesp
            am, em = samplesm
        if(be_bool):
            b, e = samplesp
            bm, em = samplesm
        if(a_bool):
            a = samplesp
            am= samplesm
        if(b_bool):
            b = samplesp
            bm= samplesm
        if(e_bool):
            e = samplesp
            em= samplesm

        res0p = tausp[0] - a*rhosp[0] - b*rhosp[2] - e*rhosp[5]
        res0m = tausm[0] - am*rhosm[0] - b*rhosm[2] - e*rhosm[5]
        v0p= np.diag(covtausp[0])+(a**2)*np.diag(covrhosp[0])+(b**2)*np.diag(covrhosp[2])+(e**2)*np.diag(covrhosp[5])
        v0m= np.diag(covtausm[0])+(am**2)*np.diag(covrhosm[0])+(bm**2)*np.diag(covrhosm[2])+(em**2)*np.diag(covrhosm[5])
        res1p = tausp[1] - a*rhosp[2] - b*rhosp[1] - e*rhosp[4]
        res1m = tausm[1] - am*rhosm[2] - b*rhosm[1] - e*rhosm[4]
        v1p= np.diag(covtausp[1]) + (a**2)*np.diag(covrhosp[2]) + (b**2)*np.diag(covrhosp[1])+(e**2)*np.diag(covrhosp[4])  
        v1m=  np.diag(covtausm[1]) + (am**2)*np.diag(covrhosm[2]) + (bm**2)*np.diag(covrhosm[1])+(em**2)*np.diag(covrhosm[4]) 
        res2p = tausp[2] - a*rhosp[5] - b*rhosp[4] - e*rhosp[3]
        res2m = tausm[2] - am*rhosm[5] - b*rhosm[4] - e*rhosm[3]
        v2p= np.diag(covtausp[2]) + (a**2)*np.diag(covrhosp[5]) + (b**2)*np.diag(covrhosp[4]) + (e**2)*np.diag(covrhosp[3])
        v2m= np.diag(covtausm[2]) + (am**2)*np.diag(covrhosm[5]) + (bm**2)*np.diag(covrhosm[4]) + (e**2)*np.diag(covrhosm[3])

        plt.clf()
        #fig, axs = plt.subplots(1, 3, figsize=(1.6*3,1.6), sharey=True)
        fig, axs = plt.subplots(1, 6, sharey=True)
        fig.subplots_adjust(hspace=0, wspace=0)
        pretty_residuals(axs, 0, meanr,res0p, np.sqrt(v0p),label=r'$\delta \tau_{0+}$', color='red')
        pretty_residuals(axs, 1, meanr,res0m, np.sqrt(v0m),label=r'$\delta \tau_{0-}$', color='green')
        pretty_residuals(axs, 2, meanr,res1p, np.sqrt(v1p),label=r'$\delta \tau_{2+}$', color='blue')
        pretty_residuals(axs, 3, meanr,res1m, np.sqrt(v1m),label=r'$\delta \tau_{2-}$', color='black')
        pretty_residuals(axs, 4, meanr,res2p, np.sqrt(v2p),label=r'$\delta \tau_{5+}$', color='gray')
        pretty_residuals(axs, 5, meanr,res2m, np.sqrt(v2m),label=r'$\delta \tau_{5-}$', color='purple')
        #plt.tight_layout()

        '''
        pretty_rho(meanr, res0p, np.sqrt(v0p),  legend=r'$\delta \tau_{0+}$', lfontsize=14, color='red', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res0m, np.sqrt(v0m),  legend=r'$\delta \tau_{0-}$', lfontsize=14, color='blue', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res1p, np.sqrt(v1p),  legend=r'$\delta \tau_{2+}$', lfontsize=14, color='green', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res1m, np.sqrt(v1m),  legend=r'$\delta \tau_{2-}$', lfontsize=14, color='pink', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res2p, np.sqrt(v2p),  legend=r'$\delta \tau_{5+}$', lfontsize=14, color='gray', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res2m, np.sqrt(v2m),  legend=r'$\delta \tau_{5-}$', lfontsize=14, color='black', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        '''
 
    elif(margin and not overall):
        a = b = e = 0; vara =  varb =  vare = 0; covab = covae = covbe = 0
        am = bm = em = 0; varam =  varbm =  varem = 0; covabm = covaem = covbem = 0
        bestpar = bestparameters(samplesp)
        bestparm = bestparameters(samplesm)
        #print("Best pars", bestpar)
        par_matcov = np.cov(samplesp) 
        if (par_matcov.size==1 ): variances = par_matcov
        else: variances = np.diagonal(par_matcov)
        covariances = sum( (par_matcov[i,i+1: ].tolist() for i in range(len(samplesp) - 1)) , [] )

        par_matcovm = np.cov(samplesm) 
        if (par_matcov.size==1 ): variancesm = par_matcovm
        else: variancesm = np.diagonal(par_matcovm)
        covariancesm = sum( (par_matcovm[i,i+1: ].tolist() for i in range(len(samplesm) - 1)) , [] )

        eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
        if not (abe_bool or ab_bool or a_bool): abe_bool = True
        if(abe_bool):
            a, b, e = bestpar
            vara, varb, vare =  variances
            covab, covae, covbe =  covariances
            am, bm, em = bestparm
            varam, varbm, varem =  variancesm
            covabm, covaem, covbem =  covariancesm 
        if(ab_bool):
            a, b = bestpar
            vara, varb =  variances
            covab =  covariances[0]
            am, bm = bestparm
            varam, varbm =  variancesm
            covabm =  covariancesm[0]
        if(ae_bool):
            a, e = bestpar
            vara, vare =  variances
            covae =  covariances[0]
            am, em = bestparm
            varam, varem =  variancesm
            covaem =  covariancesm[0]
        if(be_bool):
            b, e = bestpar
            varb, vare =  variances
            covbe =  covariances[0]
            bm, em = bestparm
            varbm, varem =  variancesm
            covbem =  covariancesm[0]
        if(a_bool):
            a =  bestpar[0]
            vara =  variances
            am =  bestparm[0]
            varam =  variancesm
        if(b_bool):
            b =  bestpar[0]
            varb =  variances
            bm =  bestparm[0]
            varbm =  variancesm
        if(e_bool):
            e =  bestpar[0]
            vare =  variances
            em =  bestparm[0]
            varem =  variancesm

        res0p = tausp[0] - a*rhosp[0] - b*rhosp[2] - e*rhosp[5]
        res0m = tausm[0] - am*rhosm[0] - b*rhosm[2] - e*rhosm[5]
        v0p= np.diag(covtausp[0])+(rhosp[0]**2)*vara +(rhosp[2]**2)*varb+(rhosp[5]**2)*vare+2*(rhosp[0]*rhosp[2]*covab+rhosp[2]*rhosp[5]*covbe+rhosp[0]*rhosp[5]*covae)+(a**2)*np.diag(covrhosp[0])+(b**2)*np.diag(covrhosp[2])+(e**2)*np.diag(covrhosp[5])
        v0m= np.diag(covtausm[0])+(rhosm[0]**2)*varam +(rhosm[2]**2)*varbm+(rhosm[5]**2)*varem+2*(rhosm[0]*rhosm[2]*covabm+rhosm[2]*rhosm[5]*covbem+rhosm[0]*rhosm[5]*covaem)+(am**2)*np.diag(covrhosm[0])+(bm**2)*np.diag(covrhosm[2])+(em**2)*np.diag(covrhosm[5])
        res1p = tausp[1] - a*rhosp[2] - b*rhosp[1] - e*rhosp[4]
        res1m = tausm[1] - am*rhosm[2] - b*rhosm[1] - e*rhosm[4]
        v1p= np.diag(covtausp[1]) + (rhosp[2]**2)*vara +(rhosp[1]**2)*varb + (rhosp[4]**2)*vare + 2*(rhosp[1]*rhosp[2]*covab+rhosp[1]*rhosp[4]*covbe+rhosp[2]*rhosp[4]*covae) + (a**2)*np.diag(covrhosp[2]) + (b**2)*np.diag(covrhosp[1])+(e**2)*np.diag(covrhosp[4])  
        v1m=  np.diag(covtausm[1]) + (rhosm[2]**2)*varam +(rhosm[1]**2)*varbm + (rhosm[4]**2)*varem + 2*(rhosm[1]*rhosm[2]*covabm+rhosm[1]*rhosm[4]*covbem+rhosm[2]*rhosm[4]*covaem) + (am**2)*np.diag(covrhosm[2]) + (bm**2)*np.diag(covrhosm[1])+(em**2)*np.diag(covrhosm[4]) 
        res2p = tausp[2] - a*rhosp[5] - b*rhosp[4] - e*rhosp[3]
        res2m = tausm[2] - am*rhosm[5] - b*rhosm[4] - e*rhosm[3]
        v2p= np.diag(covtausp[2]) + (rhosp[5]**2)*vara +(rhosp[4]**2)*varb + (rhosp[3]**2)*vare + 2*(rhosp[5]*rhosp[4]*covab+rhosp[3]*rhosp[4]*covbe+rhosp[5]*rhosp[3]*covae) + (a**2)*np.diag(covrhosp[5]) + (b**2)*np.diag(covrhosp[4]) + (e**2)*np.diag(covrhosp[3])
        v2m= np.diag(covtausm[2]) + (rhosm[5]**2)*varam +(rhosm[4]**2)*varbm + (rhosm[3]**2)*varem + 2*(rhosm[5]*rhosm[4]*covabm+rhosm[3]*rhosm[4]*covbem+rhosm[5]*rhosm[3]*covaem) + (am**2)*np.diag(covrhosm[5]) + (bm**2)*np.diag(covrhosm[4]) + (e**2)*np.diag(covrhosm[3])

        plt.clf()
        pretty_rho(meanr, res0p, np.sqrt(v0p),  legend=r'$\delta \tau_{0+}$', lfontsize=14, color='red', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res0m, np.sqrt(v0m),  legend=r'$\delta \tau_{0-}$', lfontsize=14, color='blue', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res1p, np.sqrt(v1p),  legend=r'$\delta \tau_{2+}$', lfontsize=14, color='green', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res1m, np.sqrt(v1m),  legend=r'$\delta \tau_{2-}$', lfontsize=14, color='pink', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res2p, np.sqrt(v2p),  legend=r'$\delta \tau_{5+}$', lfontsize=14, color='gray', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
        pretty_rho(meanr, res2m, np.sqrt(v2m),  legend=r'$\delta \tau_{5-}$', lfontsize=14, color='black', marker='o', ylabel=r'Residuals',title=None,  xlim=None,  ylim=None)
            
    else:
        print("Not valid configuration for overall and margin flags")
    
    
    print('Printing',  plotname)
    plt.savefig(plotname, dpi=150)       
 
  
def plotbestfit(zbin,axs,samplesp, samplesm,  meanr, data, models_combo, plotpath, title=None,  margin=False, overall=False):
    from maxlikelihood import bestparameters
    import numpy as np
    import matplotlib.gridspec as gridspec

    rhosp =[data['rhos'][2*i] for i in range(6)]; nrows = len(rhosp[0])
    covrhosp =[data['cov_rhos'][2*i*nrows:(2*i + 1)*nrows , 2*i*nrows:(2*i + 1)*nrows] for i in range(6)]
    rhosm =[data['rhos'][2*i + 1] for i in range(6)]
    covrhosm =[data['cov_rhos'][(2*i + 1)*nrows:2*(i + 1)*nrows , (2*i + 1)*nrows:2*(i + 1)*nrows] for i in range(6)]
    tausp =[data['taus'][2*i] for i in range(3)]; nrows = len(tausp[0])
    covtausp =[data['cov_taus'][2*i*nrows:(2*i + 1)*nrows , 2*i*nrows:(2*i + 1)*nrows] for i in range(3)]
    tausm =[data['taus'][2*i + 1] for i in range(3)]
    covtausm =[data['cov_taus'][(2*i + 1)*nrows:2*(i + 1)*nrows , (2*i + 1)*nrows:2*(i + 1)*nrows] for i in range(3)]
    
    if(overall and not margin):
        a = b = e = 0; 
        am = bm = em = 0; 
        eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
        if not (abe_bool or ab_bool or a_bool): abe_bool = True
        if(abe_bool):
            a, b, e = samplesp
            am, bm, em = samplesm
        if(ab_bool):
            a, b = samplesp
            am, bm = samplesm
        if(ae_bool):
            a, e = samplesp
            am, em = samplesm
        if(be_bool):
            b, e = samplesp
            bm, em = samplesm
        if(a_bool):
            a = samplesp
            am= samplesm
        if(b_bool):
            b = samplesp
            bm= samplesm
        if(e_bool):
            e = samplesp
            em= samplesm

        res0p = a*rhosp[0] + b*rhosp[2] + e*rhosp[5]
        res0m = am*rhosm[0] + b*rhosm[2] + e*rhosm[5]
        
        res1p = a*rhosp[2] + b*rhosp[1] + e*rhosp[4]
        res1m = am*rhosm[2] + bm*rhosm[1] + em*rhosm[4]
        
        res2p = a*rhosp[5]+ b*rhosp[4] + e*rhosp[3]
        res2m = am*rhosm[5] + bm*rhosm[4] + em*rhosm[3]
        

        
        colors=['red','green','blue','black']
        pretty_residuals_tomo(axs[0], meanr,tausp[0], np.sqrt(np.diag(covtausp[0])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[0].plot( meanr,res0p, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[1], meanr,tausm[0], np.sqrt(np.diag(covtausm[0])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[1].plot( meanr,res0m, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[2],  meanr,tausp[1], np.sqrt(np.diag(covtausp[1])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[2].plot( meanr,res1p, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[3], meanr,tausm[1], np.sqrt(np.diag(covtausm[1])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[3].plot( meanr,res1m, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[4],  meanr,tausp[2], np.sqrt(np.diag(covtausp[2])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[4].plot( meanr,res2p, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[5], meanr,tausm[2], np.sqrt(np.diag(covtausm[2])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[5].plot( meanr,res2m, 'k', color=colors[zbin-1])
                

    elif(margin and not overall):
        a = b = e = 0; vara =  varb =  vare = 0; covab = covae = covbe = 0
        am = bm = em = 0; varam =  varbm =  varem = 0; covabm = covaem = covbem = 0
        bestpar = bestparameters(samplesp)
        bestparm = bestparameters(samplesm)
        #print("Best pars", bestpar)
        par_matcov = np.cov(samplesp) 
        if (par_matcov.size==1 ): variances = par_matcov
        else: variances = np.diagonal(par_matcov)
        covariances = sum( (par_matcov[i,i+1: ].tolist() for i in range(len(samplesp) - 1)) , [] )

        par_matcovm = np.cov(samplesm) 
        if (par_matcov.size==1 ): variancesm = par_matcovm
        else: variancesm = np.diagonal(par_matcovm)
        covariancesm = sum( (par_matcovm[i,i+1: ].tolist() for i in range(len(samplesm) - 1)) , [] )

        eq, abe_bool, ab_bool,  ae_bool, be_bool, a_bool, b_bool, e_bool =  models_combo
    
    
        if not (abe_bool or ab_bool or a_bool): abe_bool = True
        if(abe_bool):
            a, b, e = bestpar
            vara, varb, vare =  variances
            covab, covae, covbe =  covariances
            am, bm, em = bestparm
            varam, varbm, varem =  variancesm
            covabm, covaem, covbem =  covariancesm 
        if(ab_bool):
            a, b = bestpar
            vara, varb =  variances
            covab =  covariances[0]
            am, bm = bestparm
            varam, varbm =  variancesm
            covabm =  covariancesm[0]
        if(ae_bool):
            a, e = bestpar
            vara, vare =  variances
            covae =  covariances[0]
            am, em = bestparm
            varam, varem =  variancesm
            covaem =  covariancesm[0]
        if(be_bool):
            b, e = bestpar
            varb, vare =  variances
            covbe =  covariances[0]
            bm, em = bestparm
            varbm, varem =  variancesm
            covbem =  covariancesm[0]
        if(a_bool):
            a =  bestpar[0]
            vara =  variances
            am =  bestparm[0]
            varam =  variancesm
        if(b_bool):
            b =  bestpar[0]
            varb =  variances
            bm =  bestparm[0]
            varbm =  variancesm
        if(e_bool):
            e =  bestpar[0]
            vare =  variances
            em =  bestparm[0]
            varem =  variancesm
            
    
        res0p = a*rhosp[0] + b*rhosp[2] + e*rhosp[5]
        res0m = am*rhosm[0] + b*rhosm[2] + e*rhosm[5]
        
        res1p = a*rhosp[2] + b*rhosp[1] + e*rhosp[4]
        res1m = am*rhosm[2] + bm*rhosm[1] + em*rhosm[4]
        
        res2p = a*rhosp[5]+ b*rhosp[4] + e*rhosp[3]
        res2m = am*rhosm[5] + bm*rhosm[4] + em*rhosm[3]
        
        colors=['red','green','blue','black']
        pretty_residuals_tomo(axs[0], meanr,tausp[0], np.sqrt(np.diag(covtausp[0])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[0].plot( meanr,res0p, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[1], meanr,tausm[0], np.sqrt(np.diag(covtausm[0])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[1].plot( meanr,res0m, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[2],  meanr,tausp[1], np.sqrt(np.diag(covtausp[1])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[2].plot( meanr,res1p, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[3], meanr,tausm[1], np.sqrt(np.diag(covtausm[1])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[3].plot( meanr,res1m, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[4],  meanr,tausp[2], np.sqrt(np.diag(covtausp[2])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[4].plot( meanr,res2p, 'k', color=colors[zbin-1])
        
        pretty_residuals_tomo(axs[5], meanr,tausm[2], np.sqrt(np.diag(covtausm[2])), label=r'Bin %d'%(zbin), color=colors[zbin-1])
        axs[5].plot( meanr,res2m, 'k', color=colors[zbin-1])
        
            
    else:
        print("Not valid configuration for overall and margin flags")
    
    
