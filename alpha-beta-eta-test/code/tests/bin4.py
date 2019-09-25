import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='plotting bestfit  manually varying alpha and beta')
    parser.add_argument('--taus',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/TAUS_FLASK_zbin_4.fits',
                        help='Fits file containing all tau stats used to estimate abe')
    parser.add_argument('--rhos', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/RHOS.fits',
                        help='Fits file containing all rho stats used to estimate abe')

    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    
    args = parser.parse_args()

    return args
                        
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

def plotbestfit(axs, meanr, data,  plotpath):
    import numpy as np
  
    rhosp =[data['rhos'][2*i] for i in range(6)]; nrows = len(rhosp[0])
    covrhosp =[data['cov_rhos'][2*i*nrows:(2*i + 1)*nrows , 2*i*nrows:(2*i + 1)*nrows] for i in range(6)]
    rhosm =[data['rhos'][2*i + 1] for i in range(6)]
    covrhosm =[data['cov_rhos'][(2*i + 1)*nrows:2*(i + 1)*nrows , (2*i + 1)*nrows:2*(i + 1)*nrows] for i in range(6)]
    tausp =[data['taus'][2*i] for i in range(3)]; nrows = len(tausp[0])
    covtausp =[data['cov_taus'][2*i*nrows:(2*i + 1)*nrows , 2*i*nrows:(2*i + 1)*nrows] for i in range(3)]
    tausm =[data['taus'][2*i + 1] for i in range(3)]
    covtausm =[data['cov_taus'][(2*i + 1)*nrows:2*(i + 1)*nrows , (2*i + 1)*nrows:2*(i + 1)*nrows] for i in range(3)]
    

 
    a, b, e = [0, 0, 0]
    am, bm, em = [0, 0, 0]
    
    for i, b in enumerate([0.5, 1.0, 1.5, 2.0]):
        bm = b
    
        res0p = a*rhosp[0] + b*rhosp[2] + e*rhosp[5]
        res0m = am*rhosm[0] + bm*rhosm[2] + em*rhosm[5]
        
        res1p = a*rhosp[2] + b*rhosp[1] + e*rhosp[4]
        res1m = am*rhosm[2] + bm*rhosm[1] + em*rhosm[4]
    
        res2p = a*rhosp[5]+ b*rhosp[4] + e*rhosp[3]
        res2m = am*rhosm[5] + bm*rhosm[4] + em*rhosm[3]
        
        
        colors=['red','green','blue','black']
        pretty_residuals_tomo(axs[0], meanr,tausp[0],  np.sqrt(np.diag(covtausp[0])), label=r'$\alpha$=0 , $\beta$=  %.1f, $\chi_{\nu}^{rel}$ = %.2f'%(b, np.power((tausp[0] - res0p)/np.sqrt(np.diag(covtausp[0]))  , 2).sum()/19. ), color=colors[i-1])
        axs[0].plot( meanr,res0p, 'k', color=colors[i-1])
                              
        pretty_residuals_tomo(axs[1], meanr,tausm[0], np.sqrt(np.diag(covtausm[0])), label=r'$\alpha$=0 , $\beta$=%.1f, $\chi_{\nu}^{rel}$ = %.2f'%(b, np.power((tausm[0] - res0m)/np.sqrt(np.diag(covtausm[0])) , 2).sum()/19.  ),  color=colors[i-1])
        axs[1].plot( meanr,res0m, 'k', color=colors[i-1])
                              
        pretty_residuals_tomo(axs[2],  meanr,tausp[1], np.sqrt(np.diag(covtausp[1])), label=r'$\alpha$=0 ,$\beta$=%.1f, $\chi_{\nu}^{rel}$ = %.2f'%(b, np.power((tausp[1] - res1p)/np.sqrt(np.diag(covtausp[1])), 2).sum()/19.  ),   color=colors[i-1])
        axs[2].plot( meanr,res1p, 'k', color=colors[i-1])
                              
        pretty_residuals_tomo(axs[3], meanr,tausm[1], np.sqrt(np.diag(covtausm[1])), label=r'$\alpha$=0 , $\beta$=%.1f, $\chi_{\nu}^{rel}$ = %.2f'%(b, np.power((tausm[1] - res1m)/np.sqrt(np.diag(covtausm[1])), 2).sum()/19.  ) ,  color=colors[i-1])
        axs[3].plot( meanr,res1m, 'k', color=colors[i-1])
                              
        pretty_residuals_tomo(axs[4],  meanr,tausp[2], np.sqrt(np.diag(covtausp[2])), label=r'$\alpha$=0 , $\beta$=%.1f, $\chi_{\nu}^{rel}$ = %.2f'%(b, np.power((tausp[2] - res2p)/np.sqrt(np.diag(covtausp[2])), 2).sum()/19.  ) ,  color=colors[i-1])
        axs[4].plot( meanr,res2p, 'k', color=colors[i-1])
                              
        pretty_residuals_tomo(axs[5], meanr,tausm[2], np.sqrt(np.diag(covtausm[2])), label=r'$\alpha$=0 , $\beta$=%.1f, $\chi_{\nu}^{rel}$ = %.2f'%(b, np.power((tausm[2] - res2m)/np.sqrt(np.diag(covtausm[2])), 2).sum()/19.  ) ,  color=colors[i-1])
        axs[5].plot( meanr,res2m, 'k', color=colors[i-1])
                              
                              
def main():
    import sys; sys.path.append(".")
    import numpy as np
    from src.readfits import  read_rhos, read_taus
    
    args = parse_args()
    
    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise
        
        
        
    meanr, rhos,  covrho =  read_rhos(args.rhos )
    meanr, taus,  covtau =  read_taus(args.taus )
    data = {}
    data['rhos'] = rhos
    data['cov_rhos'] = covrho
    data['taus'] = taus
    data['cov_taus'] = covtau
    
    
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()
    fig7, ax7 = plt.subplots()
    figs = [fig1, fig2, fig3, fig4, fig5, fig6]
    axs=[ax1,ax2,ax3,ax4,ax5,ax6]
    ylabels=[r'$\tau_{0+}$',r'$\tau_{0-}$',r'$\tau_{2+}$',r'$\tau_{2-}$',r'$\tau_{5+}$',r'$\tau_{5-}$']
    
    plotbestfit(axs, meanr, data,  plotspath)
    for i, fig in enumerate(figs):
        print('Printing', plotspath+'tau%d_bestfit.png'%(i))
        axs[i].set_ylabel(ylabels[i])
        axs[i].set_title('Alpha-beta')
        fig.tight_layout()
        fig.savefig(plotspath+'tau%d_bestfit_ab.png'%(i),dpi=200)
                              
                              
                              
                              
if __name__ == "__main__":
    main()
                              
                              
                              
                              
                              
