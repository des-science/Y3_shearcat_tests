import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')
import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='shape noise flask vs metacal')
    
    parser.add_argument('--metacal_cat',
                        #default='/home/dfa/sobreira/alsina/catalogs/y3_master/Y3fullmaster/Y3_mastercat_v2_6_20_18.h5',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/Y3_mastercat_7_24_19.h5',
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--flask_cats',
                        #default='/home/dfa/sobreira/alsina/catalogs/FLASK/desy3_6x2pt_lognormal-maskedcats_v3/', 
                        default='/home/dfa/sobreira/alsina/catalogs/FLASK/desy3_v3_mysp/',
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--seed', default=1, type=int, 
                        help='Selecting a particular seed')
    parser.add_argument('--nz_source',
                        #default='/home/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        default='/home/dfa/sobreira/alsina/catalogs/Y3_mastercat_7_24/nz_source_zbin.h5',
                        help='Indexes catalog to select galaxies in a particular redshift bin in Metacal')
    parser.add_argument('--plotspath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/plots/',
                        help='location of the plots.')
    args = parser.parse_args()
    return args


def main():
    import sys; sys.path.append(".")
    import numpy as np
    from src.read_cats import read_metacal,  read_flask
    from matplotlib.pyplot import contour
    import getdist
    from getdist import plots, MCSamples,  loadMCSamples

    '''
    sns.jointplot(e1_met, e2_met, kind='kde')
    if zbin is not None: filename = os.path.join(plotspath,'ee_zbin%d.png'%(zbin))
    else: filename =  os.path.join(plotspath,'ee.png')
    print("Printing File", filename)
    sns.savefig(filename)
    '''

    args = parse_args()

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(plotspath): raise

    galkeys = ['ra','dec','e_1','e_2','R11','R22','T', 'psf_T']
    
    for zbin in [1, 2, 3, 4]:
    #for zbin in [None]:
        print('Plotting histogram for zbin=', zbin)
        data_galaxies = read_metacal(args.metacal_cat, galkeys, zbin=zbin,nz_source_file=args.nz_source)
        data_gal_flask = read_flask(args.flask_cats,  args.seed, zbin=zbin, ck=1)
        
        e1_met = data_galaxies['e_1']; e1_flask = data_gal_flask['e_1']
        e2_met = data_galaxies['e_2']; e2_flask = data_gal_flask['e_2']
        e_met = np.sqrt(e1_met**2 +  e2_met**2); e_flask = np.sqrt(e1_flask**2 +  e2_flask**2)
        ee_met = e1_met**2 +  e2_met**2; ee_flask = e1_flask**2 +  e2_flask**2
    
        nbins = 10000; xmin = -2; xmax = 2  
        if zbin is not None:
            label = 'zbin:%d metacal'%(zbin)
            label2 = 'zbin:%d flask'%(zbin)
        else: label = None

        plt.clf()
        samps_met = np.column_stack((e1_met, e2_met))
        names_metcal = ['e1', 'e2']
        labels = [r'$e_{1}$', r'$e_{2}$']
        samples_metacal = MCSamples(samples=samps_met,names = names_metcal, labels = labels,  label=r'M%d $\langle e^{2} \rangle$=%.3f'%(zbin, np.mean(ee_met)))
        samps_flask = np.column_stack((e1_flask, e2_flask))
        samples_flask = MCSamples(samples=samps_flask,names = names_metcal, labels = labels,  label=r'F%d $\langle e^{2} \rangle$=%.3f'%(zbin, np.mean(ee_flask)))
        #g = plots.getSinglePlotter(width_inch=4, ratio=1)
        #g.plot_2d([samples_metacal, samples_metacal], 'e1', 'e2', filled=True)
        #g.add_legend(['metacal 1', 'metacal 2'], colored_text=True)
        g = plots.getSubplotPlotter()
        g.triangle_plot([samples_metacal, samples_flask], filled_compare=True, contour_colors=['green','darkblue'])
        if zbin is not None: filename = os.path.join(plotspath,'ee_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'ee.png')
        plt.savefig(filename, dpi=200)
        print('printed', filename)


        ##WITHOUT SHAPENOISE
        '''
        flaskeys = ['GAMMA1_TRUE', 'GAMMA2_TRUE']
        data_gal_flask_t = read_flask(args.flask_cats,  args.seed, zbin=zbin, ck=1,  keys=flaskeys)
        e1_flask_t = data_gal_flask_t['GAMMA1_TRUE']
        e2_flask_t = data_gal_flask_t['GAMMA2_TRUE']
        ee_flask_t = e1_flask_t**2 +  e2_flask_t**2
        
        plt.clf()
        samps_met = np.column_stack((e1_met, e2_met))
        names_metcal = ['e1', 'e2']
        labels = [r'$e_{1}$', r'$e_{2}$']
        samples_metacal = MCSamples(samples=samps_met,names = names_metcal, labels = labels,  label=r'M%d $\langle e^{2} \rangle$=%.3f'%(zbin, np.mean(ee_met)))
        samps_flask = np.column_stack((e1_flask_t, e2_flask_t))
        samples_flask = MCSamples(samples=samps_flask,names = names_metcal, labels = labels,  label=r'F%d $\langle e^{2} \rangle$=%.3f'%(zbin, np.mean(ee_flask_t)))
        g = plots.getSubplotPlotter()
        g.triangle_plot([samples_metacal, samples_flask], filled_compare=True, contour_colors=['green','darkblue'])
        if zbin is not None: filename = os.path.join(plotspath,'ee_zbin%d_true.png'%(zbin))
        else: filename =  os.path.join(plotspath,'ee_true.png')
        plt.savefig(filename, dpi=200)
        '''

        
        '''
        print('Len e1_met', len(e1_met), 'Min e1_met:',  min(e1_met), 'Max e1_met', max(e1_met), 'Mean e1_met', np.mean(e1_met) )
        plt.hist(e1_met, bins=np.linspace(xmin, xmax, nbins), label=label,  weights=np.ones(len(e1_met))/len(e1_met) )
        print('Len e1_flask', len(e1_flask), 'Min e1_flask:',  min(e1_flask), 'Max e1_flask', max(e1_flask), 'Mean e1_flask', np.mean(e1_flask)  )
        plt.hist(e1_flask, bins=np.linspace(xmin, xmax, nbins), label=label2,  weights=np.ones(len(e1_flask))/len(e1_flask), alpha=0.5, color='red' )
        plt.ylabel('Counts'); plt.xlabel(r'$e_{1}$')
        #plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        #plt.ylim(0, 1.2)
        plt.xlim( -2.2, 2.2)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'e1_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'e1.png')
        print("Printing File", filename)
        plt.savefig(filename)

        plt.clf()
        print('Len e2_met', len(e2_met), 'Min e2_met:',  min(e2_met), 'Max e2_met', max(e2_met), 'Mean e2_met', np.mean(e2_met) )
        plt.hist(e2_met, bins=np.linspace(xmin, xmax, nbins), label=label,  weights=np.ones(len(e2_met))/len(e2_met) )
        print('Len e2_flask', len(e2_flask), 'Min e2_flask:',  min(e2_flask), 'Max e2_flask', max(e2_flask), 'Mean e2_flask', np.mean(e2_flask) )
        plt.hist(e2_flask, bins=np.linspace(xmin, xmax, nbins), label=label2,  weights=np.ones(len(e2_flask))/len(e2_flask), alpha=0.5, color='red' )
        plt.ylabel('Counts'); plt.xlabel(r'$e_{2}$')
        #plt.ylim(0, 1.2)
        plt.xlim( -2.2, 2.2)
        #plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'e2_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'e2.png')
        print("Printing File", filename)
        plt.savefig(filename)

        nbins = 10000; xmin = 0; xmax = 2  
        plt.clf()
        print('Len e_met', len(e_met), 'Min e_met:',  min(e_met), 'Max e_met', max(e_met), 'Mean e_met', np.mean(e_met) )
        plt.hist(e_met, bins=np.linspace(xmin, xmax, nbins), label=label,  weights=np.ones(len(e_met))/len(e_met) )
        print('Len e_flask', len(e_flask), 'Min e_flask:',  min(e_flask), 'Max e_flask', max(e_flask), 'Mean e_met', np.mean(e_met) )
        plt.hist(e_flask, bins=np.linspace(xmin, xmax, nbins), label=label2, weights=np.ones(len(e_flask))/len(e_flask), alpha=0.5, color='red' )
        plt.ylabel('Counts'); plt.xlabel('e')
        #plt.ylim(0, 1.2)
        plt.xlim(0, 2.2)
        #plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'e_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'e.png')
        print("Printing File", filename)
        plt.savefig(filename)

        plt.clf()
        print('Len ee_met', len(ee_met), 'Min ee_met:',  min(ee_met), 'Max ee_met', max(ee_met), 'Mean ee_met', np.mean(ee_met) )
        plt.hist(ee_met, bins=np.linspace(xmin, xmax, nbins), label=label,  weights=np.ones(len(ee_met))/len(ee_met) )
        print('Len ee_flask', len(ee_flask), 'Min ee_flask:',  min(ee_flask), 'Max ee_flask', max(ee_flask), 'Mean ee_flask', np.mean(ee_flask) )
        plt.hist(ee_flask, bins=np.linspace(xmin, xmax, nbins), label=label2, weights=np.ones(len(ee_flask))/len(ee_flask), alpha=0.5, color='red' )
        plt.ylabel('Counts'); plt.xlabel(r'$e^{2}$')
        #plt.ylim(0, 1.2)
        plt.xlim(0, 2.2)
        #plt.yscale('log')
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'ee_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'ee.png')
        print("Printing File", filename)
        plt.savefig(filename)
        '''
        
        
        '''
        fig = plt.figure(figsize = (6, 6))
        grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
        main_ax = fig.add_subplot(grid[:-1, 1: ])
        y_hist = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
        x_hist = fig.add_subplot(grid[-1, 1: ], yticklabels=[], sharex=main_ax) 
        
        
        nbins = 100; xmin =-2; xmax = 2 
        bins = [np.linspace(xmin, xmax, nbins), np.linspace(xmin, xmax, nbins)]
        weights_met = None
        weights_flask = None
        hist = main_ax.hist2d(e1_met, e2_met, bins=bins, label=label,  weights=weights_met,  cmap='gnuplot')
        main_ax.contour(hist[0],extent=[hist[1].min(),hist[1].max(),hist[2].min(),hist[2].max()],linewidths=3)
        #main_ax.hist2d(e1_flask, e2_flask, bins=bins, label=label2, weights=weights_flask, alpha=0.5, color='red' )
        #main_ax.plot(e1_met, e2_met,  'ok',  markersize=3, alpha =0.2)
        #main_ax.legend(loc='best', fontsize=10)
        #fig.colorbar(hist[3], ax=main_ax)

        y_hist.set_ylabel(r'$e_{2}$')
        y_hist.hist(e2_met, bins=np.linspace(xmin, xmax, nbins), label=label,  weights=np.ones(len(e2_met))/len(e2_met), orientation='horizontal' )
        y_hist.hist(e2_flask, bins=np.linspace(xmin, xmax, nbins), label=label2, weights=np.ones(len(e2_flask))/len(e2_flask), orientation='horizontal',  alpha=0.5, color='red')
        #y_hist.legend(loc='best', fontsize=10)
        y_hist.invert_xaxis()
        
        x_hist.set_xlabel(r'$e_{1}$')
        x_hist.hist(e1_met, bins=np.linspace(xmin, xmax, nbins), label=label,  weights=np.ones(len(e1_met))/len(e1_met) )
        x_hist.hist(e1_flask, bins=np.linspace(xmin, xmax, nbins), label=label2, weights=np.ones(len(e1_flask))/len(e1_flask), alpha=0.5, color='red' )
        x_hist.legend(loc='best', fontsize=10)
        x_hist.invert_yaxis()

        fig.tight_layout()
        if zbin is not None: filename = os.path.join(plotspath,'ee_zbin%d.png'%(zbin))
        else: filename =  os.path.join(plotspath,'ee.png')
        print("Printing File", filename)
        fig.savefig(filename)
        '''
        

        
        
if __name__ == "__main__":
    main()
