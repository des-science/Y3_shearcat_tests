import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--samplesfile_continf', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/runs/contamination/2sig-inf/2pt_sim_1110_baseline_Y3cov_contaminated_inf_2sig.fits_d_l_chain.txt', help='txt file with the samples contaminated after running cosmosis')
    parser.add_argument('--samplesfile_contsup', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/runs/contamination/2sig-sup/2pt_sim_1110_baseline_Y3cov_contaminated_sup_2sig.fits_d_l_chain.txt', help='txt file with the samples contaminated after running cosmosis')
    parser.add_argument('--samplesfile_forecast', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/runs/nocontamination/5/2pt_sim_1110_baseline_Y3cov.fits_d_l_chain.txt', help='txt file with the samples after running cosmosis')
    parser.add_argument('--outpath',
                        default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/',
                        help='path where output will be send') 
    
    args = parser.parse_args()
    return args

def get_nsample(filename):
    with open(filename,"r") as fi:
        for ln in fi:
            if ln.startswith("#nsample="):
                nsamples = int(ln[9:])
        return nsamples
def main():
    import getdist
    from getdist import plots, MCSamples,  loadMCSamples
    import numpy as np
    args = parse_args()
    out = os.path.expanduser(args.outpath)
    out = os.path.join(out,'plots')
    if not os.path.isdir(out):
        os.makedirs(out)

    allnames = np.array(['Om','h0','Ob','ns','a_s','Onuh2','b1','b2','b3','b4','b5','m1','m2','m3','m4','ia_a','ia_alpha', 'wpz_b1','wpz_b2','wpz_b3','wpz_b4','lpz_b1','lpz_bin2','lpz_bin3','lpz_bin4','lpz_bin5','s8','like','post','weight'])
    alllabels = np.array(['\Omega_m', 'h_0', '\Omega_b', 'n_s','a_s', r'\Omega_{\nu}h^{2}','b1','b2','b3','b4','b5','m1','m2','m3','m4','ia_a','ia_alpha', 'wpz_b1','wpz_b2','wpz_b3','wpz_b4','lpz_b1','lpz_bin2','lpz_bin3','lpz_bin4','lpz_bin5',r'\sigma_{8}','like','post','weight'])
    useindex = [0, 1, 2, 3, 4, 5,- 4]
    #Only change this indixes to selec particular parameters
    #useindex = [0, - 4]
    usednames = allnames[useindex]
    usedlabels = alllabels[useindex]

    ##FIDUCIAL FORECAST
    nsample =  get_nsample(args.samplesfile_forecast)
    allsamplestable = np.loadtxt(args.samplesfile_forecast)
    allsamplestable =  allsamplestable[ -nsample:, : ]                     
    usedsamples = allsamplestable[:, useindex]
    usedweights = allsamplestable[: , -1]
    usedpost =  allsamplestable[:, -2]
    samples = MCSamples(samples=usedsamples, names=usednames,
                        labels=usedlabels, weights=usedweights , loglikes=usedpost,
                        label='Forecast' )
    samples.removeBurn(remove=0.1)
    filename = os.path.join(out,'fiducial.tex')
    print('Printing', filename)
    print(samples.getTable().tableTex(), file=open(filename, "w"))

    #INF (MIN) contamination
    nsample_continf =  get_nsample(args.samplesfile_continf)
    allsamplestable_continf = np.loadtxt(args.samplesfile_continf)
    allsamplestable_continf =  allsamplestable_continf[-nsample_continf:, : ]
    usedsamples_continf = allsamplestable_continf[:, useindex]
    usedweights_continf = allsamplestable_continf[: , -1]
    usedpost_continf =  allsamplestable_continf[:, -2]
    samples_continf = MCSamples(samples=usedsamples_continf,
                             names=usednames, labels=usedlabels, weights=usedweights_continf ,loglikes=usedpost_continf,
                             label='Min contamination' )
    samples_continf.removeBurn(remove=0.1)
    filename = os.path.join(out,'fiducial_contaminatedinf.tex')
    print('Printing', filename)
    print(samples_continf.getTable().tableTex(), file=open(filename, "w"))

    #sup (Max) contamination
    nsample_contsup =  get_nsample(args.samplesfile_contsup)
    allsamplestable_contsup = np.loadtxt(args.samplesfile_contsup)
    allsamplestable_contsup =  allsamplestable_contsup[-nsample_contsup:, : ]
    usedsamples_contsup = allsamplestable_contsup[:, useindex]
    usedweights_contsup = allsamplestable_contsup[: , -1]
    usedpost_contsup =  allsamplestable_contsup[:, -2]
    samples_contsup = MCSamples(samples=usedsamples_contsup,
                             names=usednames, labels=usedlabels, weights=usedweights_contsup ,loglikes=usedpost_contsup,
                             label='Max contamination' )
    samples_contsup.removeBurn(remove=0.1)
    filename = os.path.join(out,'fiducial_contaminatedsup.tex')
    print('Printing', filename)
    print(samples_contsup.getTable().tableTex(), file=open(filename, "w"))




    
    g = plots.getSubplotPlotter()
    g.settings.plot_meanlikes = False
    g.settings.alpha_factor_contour_lines = True
    #g.settings.axis_marker_lw = 5
    g.settings.figure_legend_frame = True
    g.settings.alpha_filled_add=0.1
    g.settings.title_limit_fontsize = 16
    g.settings.figure_legend_loc = 'best'
    g.settings.rcSizes(axes_fontsize = 12, lab_fontsize=20, legend_fontsize =40)
    g.triangle_plot([samples, samples_continf, samples_contsup], filled_compare=[True, False, False],  shaded=False, Filled=False,
                    line_args=[{'ls':'solid', 'lw':2, 'color':'green'},{'ls':'--','lw':2, 'color':'blue'},{'ls':'dotted','lw':2, 'color':'red'}],
                    contour_colors=['green','blue', 'red'], title_limit=1)
    #g.add_legend(legend_labels=[legend_name], fontsize=36, legend_loc=(-3.5,7))
    filename = os.path.join(out,'getdistplot.png')
    print('Printing', filename)
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    #g.export(filename)
    

    


if __name__ == "__main__":
    main()
