import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--samplesfile_contaminated', default='/home/dfa/sobreira/alsina/Y3_shearcat_tests/alpha-beta-eta-test/forecast/runs/nocontamination/4/2pt_sim_1110_baseline_Y3cov.fits_d_l_chain.txt', help='txt file with the samples contaminated after running cosmosis')
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

    
    nsample =  get_nsample(args.samplesfile_forecast)
    allsamplestable = np.loadtxt(args.samplesfile_forecast)
    allsamplestable =  allsamplestable[ -nsample:, : ]                     
    usedsamples = allsamplestable[:, useindex]
    usedweights = allsamplestable[: , -1]
    usedpost =  allsamplestable[:, -2]
    samples = MCSamples(samples=usedsamples, names=usednames,
                        labels=usedlabels, weights=usedweights , loglikes=usedpost,
                        label='Forecast Run1' )
    samples.removeBurn(remove=0.1)
    filename = os.path.join(out,'fiducial.tex')
    print('Printing', filename)
    print(samples.getTable().tableTex(), file=open(filename, "w"))

    nsample_cont =  get_nsample(args.samplesfile_contaminated)
    allsamplestable_cont = np.loadtxt(args.samplesfile_contaminated)
    allsamplestable_cont =  allsamplestable_cont[-nsample_cont:, : ]
    usedsamples_cont = allsamplestable_cont[:, useindex]
    usedweights_cont = allsamplestable_cont[: , -1]
    usedpost_cont =  allsamplestable_cont[:, -2]
    samples_cont = MCSamples(samples=usedsamples_cont,
                             names=usednames, labels=usedlabels, weights=usedweights_cont ,loglikes=usedpost_cont,
                             label='Forecast Run2' )
    samples_cont.removeBurn(remove=0.1)
    filename = os.path.join(out,'fiducial_contaminated.tex')
    print('Printing', filename)
    print(samples_cont.getTable().tableTex(), file=open(filename, "w"))




    
    g = plots.getSubplotPlotter()
    g.settings.plot_meanlikes = False
    g.settings.alpha_factor_contour_lines = True
    #g.settings.axis_marker_lw = 5
    g.settings.figure_legend_frame = True
    g.settings.alpha_filled_add=0.35
    g.settings.title_limit_fontsize = 16
    g.settings.figure_legend_loc = 'best'
    g.settings.rcSizes(axes_fontsize = 12, lab_fontsize=20, legend_fontsize =40)
    g.triangle_plot([samples, samples_cont], filled_compare=[True,  True], 
                    line_args=[{'ls':'solid', 'lw':2, 'color':'green'},{'ls':'--','lw':2, 'color':'blue'}],
                    contour_colors=['green','blue'], title_limit=1)
    #g.add_legend(legend_labels=[legend_name], fontsize=36, legend_loc=(-3.5,7))
    filename = os.path.join(out,'getdistplot.png')
    print('Printing', filename)
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    #g.export(filename)
    

    


if __name__ == "__main__":
    main()
