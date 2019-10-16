import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    parser.add_argument('--tausmetacal',
                        default='/home/dfa/sobreira/alsina/catalogs/JK/taus_jk1000_01-250/taus_src-cat_jk1_z4.fits',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--tausflask',
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/taus_flask_01-250/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--tausjk',
                        default='/home/dfa/sobreira/alsina/catalogs/JK/taus_jk1000_01-250/',
                        help='Full Path to the taus measurement of flask catalogs')
    parser.add_argument('--zbin', default=4 , type=int,
                        help='seed used, useful to run parallel')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()

    return args



def covariance_jck(TOTAL_PHI,jk_r,type_cov):
    import numpy as np
    if type_cov=='jackknife':
        fact=(jk_r-1.)/(jk_r)
    elif type_cov=='bootstrap':
        fact=1./(jk_r)
    #  Covariance estimation
    average=np.zeros(TOTAL_PHI.shape[0])
    cov_jck=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
    err_jck=np.zeros(TOTAL_PHI.shape[0])
    for kk in range(jk_r):
        average+=TOTAL_PHI[:,kk]
    average=average/(jk_r)
    # print average
    for ii in range(TOTAL_PHI.shape[0]):
        for jj in range(ii+1):
            for kk in range(jk_r):
                cov_jck[jj,ii]+=TOTAL_PHI[ii,kk]*TOTAL_PHI[jj,kk]
            cov_jck[jj,ii]=(-average[ii]*average[jj]*jk_r+cov_jck[jj,ii])*fact
            cov_jck[ii,jj]=cov_jck[jj,ii]

    for ii in range(TOTAL_PHI.shape[0]):
        err_jck[ii]=np.sqrt(cov_jck[ii,ii])

    #compute correlation
    corr=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
    for i in range(TOTAL_PHI.shape[0]):
        for j in range(TOTAL_PHI.shape[0]):
            corr[i,j]=cov_jck[i,j]/(np.sqrt(cov_jck[i,i]*cov_jck[j,j]))
    average=average*fact
    return {'cov' : cov_jck, 'err' : err_jck, 'corr':corr, 'mean':average}


def plotmetacal(axs, tausname,  color, label):
    from src.readfits import read_taus
    import numpy as np

    ax1, ax2, ax3, ax4, ax5, ax6 = axs
    
    name = os.path.join(tausname)
      
    meanr, taus, covmat = read_taus(name)
           
    nrows = len(meanr)
   
    print('matrix covariance shape', covmat.shape)
    sig_tau0p = np.sqrt(np.diag(covmat[0:nrows, 0:nrows]))
    sig_tau0m = np.sqrt(np.diag(covmat[nrows:2*nrows, nrows:2*nrows]))
    sig_tau2p = np.sqrt(np.diag(covmat[2*nrows:3*nrows, 2*nrows:3*nrows]))
    sig_tau2m = np.sqrt(np.diag(covmat[3*nrows:4*nrows, 3*nrows:4*nrows]))
    sig_tau5p = np.sqrt(np.diag(covmat[4*nrows:5*nrows, 4*nrows:5*nrows]))
    sig_tau5m = np.sqrt(np.diag(covmat[5*nrows:6*nrows, 5*nrows:6*nrows]))

    taumeans = taus
    sig_taus = [sig_tau0p,sig_tau0m,sig_tau2p,sig_tau2m,sig_tau5p,sig_tau5m  ]

    #ylabels = [r'$\tau_{0+}$', r'$\tau_{0-}$', r'$\tau_{2+}$', r'$\tau_{2-}$', r'$\tau_{5+}$', r'$\tau_{5-}$']
    ylabels = [r'$\sigma \tau_{0+}$', r'$\sigma \tau_{0-}$', r'$\sigma \tau_{2+}$', r'$\sigma \tau_{2-}$', r'$\sigma \tau_{5+}$', r'$\sigma \tau_{5-}$']
    for i, ax in enumerate(axs):
        ax.errorbar(meanr,sig_taus[i],yerr=None,color=color, ls='', marker='.', capsize=2, label=label)
        #ax.errorbar(meanr,taumeans[i],yerr=sig_taus[i],color=color, ls='', marker='.', capsize=2, label=label)
        ax.legend(loc='best', fontsize=10)
        ax.set_ylabel(ylabels[i]); ax1.set_xlabel(r'$\theta$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        #ax.set_ylim([ -2.e-6,2.e-6 ])
def plotflask(axs, zbin, tausflask,color, label, ndraws):
    from src.readfits import read_taus
    import numpy as np

    ax1, ax2, ax3, ax4, ax5, ax6 = axs
    
    veclist = []
    count = 0; ck = 1; seed = 1
    while(count <ndraws):
        name = os.path.join(tausflask, 'taus_src-cat_s%d_z%d_ck%d.fits'%(seed,zbin, ck  ))
        exist =  os.path.isfile(name)
        if exist:
            meanr, taus, covtaus = read_taus(name)
            if (np.count_nonzero(taus) == 0):
                print("Warning, weird measurement, skipping", name)
            else:
                veclist.append(np.concatenate(np.c_[taus]))
                count +=1
        else:
            print(name, 'Does not exist')          
        seed += 1

    print(count, "FLASK catalogs were read")
    meanvec = np.mean(veclist, axis=0)
    
    nrows = len(meanr)
    tau0pmean =  meanvec[0:nrows]
    tau0mmean =  meanvec[nrows:2*nrows]
    tau2pmean =  meanvec[2*nrows:3*nrows]
    tau2mmean =  meanvec[3*nrows:4*nrows]
    tau5pmean =  meanvec[4*nrows:5*nrows]
    tau5mmean =  meanvec[5*nrows:6*nrows]

    ranveclist = np.c_[veclist].T
    covmat = np.cov(ranveclist)
    print("COVMAT FLASK")
    print(covmat)
    
    print('matrix covariance shape', covmat.shape)
    sig_tau0p = np.sqrt(np.diag(covmat[0:nrows, 0:nrows]))
    sig_tau0m = np.sqrt(np.diag(covmat[nrows:2*nrows, nrows:2*nrows]))
    sig_tau2p = np.sqrt(np.diag(covmat[2*nrows:3*nrows, 2*nrows:3*nrows]))
    sig_tau2m = np.sqrt(np.diag(covmat[3*nrows:4*nrows, 3*nrows:4*nrows]))
    sig_tau5p = np.sqrt(np.diag(covmat[4*nrows:5*nrows, 4*nrows:5*nrows]))
    sig_tau5m = np.sqrt(np.diag(covmat[5*nrows:6*nrows, 5*nrows:6*nrows]))

    taumeans = [tau0pmean,tau0mmean,tau2pmean,tau2mmean,tau5pmean,tau5mmean ]
    sig_taus = [sig_tau0p,sig_tau0m,sig_tau2p,sig_tau2m,sig_tau5p,sig_tau5m  ]

    #ylabels = [r'$\tau_{0+}$', r'$\tau_{0-}$', r'$\tau_{2+}$', r'$\tau_{2-}$', r'$\tau_{5+}$', r'$\tau_{5-}$']
    ylabels = [r'$\sigma \tau_{0+}$', r'$\sigma \tau_{0-}$', r'$\sigma \tau_{2+}$', r'$\sigma \tau_{2-}$', r'$\sigma \tau_{5+}$', r'$\sigma \tau_{5-}$']
    for i, ax in enumerate(axs):
        ax.errorbar(meanr,sig_taus[i],yerr=None,color=color, ls='', marker='.', capsize=2, label=label)
        #ax.errorbar(meanr,taumeans[i],yerr=sig_taus[i],color=color, ls='', marker='.', capsize=2, label=label)
        ax.legend(loc='best', fontsize=10)
        ax.set_ylabel(ylabels[i]); ax1.set_xlabel(r'$\theta$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        #ax.set_ylim([ -2.e-6,2.e-6 ])
def plotjk(axs, zbin, tausflask,  color, label, ndraws, numpycov=True):
    from src.readfits import read_taus
    import numpy as np

    ax1, ax2, ax3, ax4, ax5, ax6 = axs
    
    veclist = []
    count = 0; ck = 1; jkid = 1
    while(count <ndraws):
            name = os.path.join(tausflask, 'taus_src-cat_jk%d_z%d.fits'%(jkid,zbin))
            exist =  os.path.isfile(name)
            if exist:
                meanr, taus, covtaus = read_taus(name)
                if (np.count_nonzero(taus) == 0):
                    print("Warning, weird measurement, skipping", name)
                else:
                    veclist.append(np.concatenate(np.c_[taus]))
                    count +=1
            else:
                print(name, 'Does not exist')          
            jkid += 1
    print(count, "JK catalogs were read")
    meanvec = np.mean(veclist, axis=0)
    
    nrows = len(meanr)
    tau0pmean =  meanvec[0:nrows]
    tau0mmean =  meanvec[nrows:2*nrows]
    tau2pmean =  meanvec[2*nrows:3*nrows]
    tau2mmean =  meanvec[3*nrows:4*nrows]
    tau5pmean =  meanvec[4*nrows:5*nrows]
    tau5mmean =  meanvec[5*nrows:6*nrows]

    ranveclist = np.c_[veclist].T

    if numpycov == True:
        print("Using numpy.cov")
        covmat = np.cov(ranveclist)
        covmat *= (ndraws - 1)
    else:
        print("Calculating manually cov")
        covmat = covariance_jck(ranveclist,ndraws,'jackknife')['cov']
    print("COVMAT JK")
    print(covmat)
    print('matrix covariance shape', covmat.shape)
    sig_tau0p = np.sqrt(np.diag(covmat[0:nrows, 0:nrows]))
    sig_tau0m = np.sqrt(np.diag(covmat[nrows:2*nrows, nrows:2*nrows]))
    sig_tau2p = np.sqrt(np.diag(covmat[2*nrows:3*nrows, 2*nrows:3*nrows]))
    sig_tau2m = np.sqrt(np.diag(covmat[3*nrows:4*nrows, 3*nrows:4*nrows]))
    sig_tau5p = np.sqrt(np.diag(covmat[4*nrows:5*nrows, 4*nrows:5*nrows]))
    sig_tau5m = np.sqrt(np.diag(covmat[5*nrows:6*nrows, 5*nrows:6*nrows]))

    taumeans = [tau0pmean,tau0mmean,tau2pmean,tau2mmean,tau5pmean,tau5mmean ]
    sig_taus = [sig_tau0p,sig_tau0m,sig_tau2p,sig_tau2m,sig_tau5p,sig_tau5m  ]

    #ylabels = [r'$\tau_{0+}$', r'$\tau_{0-}$', r'$\tau_{2+}$', r'$\tau_{2-}$', r'$\tau_{5+}$', r'$\tau_{5-}$']
    ylabels = [r'$\sigma \tau_{0+}$', r'$\sigma \tau_{0-}$', r'$\sigma \tau_{2+}$', r'$\sigma \tau_{2-}$', r'$\sigma \tau_{5+}$', r'$\sigma \tau_{5-}$']
    for i, ax in enumerate(axs):
        ax.errorbar(meanr,sig_taus[i],yerr=None,color=color, ls='', marker='.', capsize=2, label=label)
        #ax.errorbar(meanr,taumeans[i],yerr=sig_taus[i],color=color, ls='', marker='.', capsize=2, label=label)
        ax.legend(loc='best', fontsize=10)
        ax.set_ylabel(ylabels[i]); ax1.set_xlabel(r'$\theta$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        #ax.set_ylim([ -2.e-6,2.e-6 ])

def main():
    import sys; sys.path.append(".")
    from src.readfits import read_taus
    import numpy as np
    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise
        

    figs = [];  axs = []; filenames = []
    names = ['taus0p', 'taus0m', 'taus2p', 'taus2m' , 'taus5p' , 'taus5m']
    for i in range(6):
        figaux, axaux = plt.subplots()
        figs.append(figaux); axs.append(axaux)
        filenames.append(os.path.join(plotspath,'%s_flask_JK_zbin%d%s'%(names[i], args.zbin, '.png') ))

    '''
    measurespath = os.path.expanduser('/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/' )
    plotmetacal(axs, os.path.join(measurespath,'TAUS_FLASK_zbin_4.fits'), 'green', 'MY SHAPE NOISE')
    plotmetacal(axs, os.path.join(measurespath,'TAUS_FLASK_zbin_old_4.fits'), 'red', 'OLD SHAPE NOISE')
    plotmetacal(axs, os.path.join(measurespath,'TAUS_zbin_4_1-250.fits'), 'blue', 'CURRENT METACAL')
    plotmetacal(axs, os.path.join(measurespath,'tau_4.fits'), 'gray', 'MARCO JK')
    '''
    
    '''
    ndraws = 290
    args.zbin = 3
    #plotmetacal(axs, args.tausmetacal, 'green', 'Taus treecorr')
    plotflask(axs, args.zbin, args.tausflask, 'red', 'Taus flask', ndraws)
    plotjk(axs, args.zbin, args.tausjk, 'blue', 'Taus JK ', 300, numpycov=True)
    '''
    
    ndraws = 250
    plotflask(axs, 1, args.tausflask, 'blue', 'Taus flask zbin1', ndraws)
    plotjk(axs, 1, args.tausjk, 'blue', 'Taus JK zbin1 ', 250, numpycov=True)
    plotflask(axs, 2, args.tausflask, 'red', 'Taus flask  zbin2', ndraws)
    plotjk(axs, 2, args.tausjk, 'red', 'Taus JK  zbin2', 250, numpycov=True)
    plotflask(axs, 3, args.tausflask, 'green', 'Taus flask  zbin3', ndraws)
    plotjk(axs, 3, args.tausjk, 'green', 'Taus JK  zbin3', 250, numpycov=True)
    plotflask(axs, 4, args.tausflask, 'yellow', 'Taus flask  zbin4', ndraws)
    plotjk(axs, 4, args.tausjk, 'yellow', 'Taus JK  zbin4', 250, numpycov=True)
    
    
    
    

    '''
    plotjk(axs, args.zbin, args.tausjk, 'blue', 'Taus JK 100', 100)
    plotjk(axs, args.zbin, args.tausjk, 'green', 'Taus JK 200', 200)
    plotjk(axs, args.zbin, args.tausjk, 'yellow', 'Taus JK 300', 300)
    plotjk(axs, args.zbin, args.tausjk, 'red', 'Taus JK 400', 400)
    '''


    '''
    measurespath = os.path.expanduser('/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/' )
    plotmetacal(axs, os.path.join(measurespath,'TAUS_FLASK_zbin_4.fits'),  'green', 'Flask 1-250')
    plotmetacal(axs, os.path.join(measurespath,'tau_4.fits'), 'gray', 'JK (Y1 zbin) Marco 1-250')
    plotflask(axs, args.zbin, args.tausflask, 'red', 'Flask 0.1-250', 250)
    plotjk(axs, args.zbin, args.tausjk, 'yellow', 'JK 0.1-250', 690)
    '''
    
    '''
    measurespath = os.path.expanduser('/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/' )
    plotmetacal(axs, os.path.join(measurespath,'TAUS_zbin_4.fits'),  'green', 'Sigma Taus default bin-slop')
    plotmetacal(axs, os.path.join(measurespath,'TAUS_zbin_4_binslope003.fits'), 'gray', 'Sigma Taus 0.03 bin-slop')
    '''

    for i, fig in enumerate(figs):
        fig.tight_layout()
        fig.savefig(filenames[i], dpi=200)
        plt.close(fig)
        print(filenames[i], 'Printed!')
    
if __name__ == "__main__":
    main()
