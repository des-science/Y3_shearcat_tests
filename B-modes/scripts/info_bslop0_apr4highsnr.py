import numpy as np
import os
from colormaps import plasma,viridis

'''
Information about paths and parameters used in pipeline. 
'''

"""
BASIC
---------------------------------
Define the basic parameters here.
Define running on data (mode = 'data') or running or sims,
for now implemented on mice simulations (mode = 'mice').
Config and corresponding paths will be automatically
defined depending on the mode used. 

BLINDING Instructions: set blind = True and plot_blinded = False 
below and run the scripts normally (i.e. run call_ggl_y3kp). 
This will save a TwoPointFile with the unblinded measurements. 
Then, source the setup file of cosmosis 
(i.e. source cosmosis/config/setup-cosmosis-nersc, 
using your own path to the file). This will enable the 
cosmosis environment. Once this is done, run the blind.py script 
(you might need to pull from 2pt_pipeline repository first, 
and edit the path to that directory in blind.py). 
The blind.py script will blind the measurements using 
cosmosis, will save a twopoint file with the blinded measurements 
and will delete the unblinded ones. Then, to plot the blinded 
measurements, set plot_blinded = True and run call_ggl_y3kp.py again. 
Then you are done. All this process only affects the gammat 
measurements, not the systematics tests. 
"""

basic = {
    'mode':'data',
    'blind': False,
    'plot_blinded': False
}


"""
CONFIGURATION
---------------------------------
Define the parameters of your run here. 
If you are running on data edit config_data, 
and if you are using MiCE simulations use
config_mice.
"""

config_data = {
    'njk': 2,
    'bslop': 0.0,
    'nthbins': 120000,
    #'nthbins': 1000,
    'thlims': np.array([2.5,250.]),
    #'nthbins': 10,
    #'thlims': np.array([2.5,250.]),
    #'filename_mastercat': '/fs/scratch/cond0080/y3_catalogs/Y3_mastercat_6_15_19_subsampled.h5',
    #'filename_mastercat': '/fs/scratch/cond0080/y3_catalogs/Y3_mastercat_7_24_19.h5',
    'filename_mastercat': '/fs/scratch/cond0080/y3_catalogs/March16_2020/Y3_mastercat_03_16_20_highsnr.h5',
    'redmagic_v': 'combined_sample_fid',
    'zslim_v': 'y1',
    'zs_v': 'bpz',
    #'zs_v': 'sompz',  # This doesn't seem to work
    'zllim_v': 'y3'
    }

# Really should update this to accommodate different path lengths...
config_data['mastercat_v'] = config_data['filename_mastercat'][53:-3]

config_mice = {
    'njk': 300,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'version': '2',
    'redmagic_v': 'y1',
    'zslim_v': 'y1',
    'zs_v': 'bpz',
    'zllim_v': 'y1'
    }

if basic['mode'] == 'data':
    config = config_data
if basic['mode'] == 'mice':
    config = config_mice

#print '\nChosen configuration:\n--------------------------\n', config

"""
PATHS
---------------------------------
Define the paths dictionary. This dictionary is imported
in the other scripts. Add more paths as necessary.
"""

paths = {}

paths['y1'] = '/Volumes/Data/y1_shear_tests/cats/jk/test_mcal_bpzmof_unblind/'
paths['runs'] =  'runs/' 
paths['plots'] = 'plots/'
paths['y3'] = 'ggl_results/'
paths['y3_exp'] = 'ggl_data/'
paths['redmagic', config['redmagic_v']] = '../../src/xcorr/lens_cats/redmagic/%s/%s/njk_%d/'%(config_data['mastercat_v'], config['redmagic_v'], config['njk'])
paths['redmagic', 'y1'] = '../../src/xcorr/lens_cats/redmagic/y1/'
paths['lens'] = paths['redmagic', '%s'%config['redmagic_v']] + 'lens.fits'
paths['randoms'] = paths['redmagic', '%s'%config['redmagic_v']] + 'random.fits'
paths['lens_nz'] = 'y1_2pt_NG_mcal_1110.fits'
paths['source_nz'] = 'y1_2pt_NG_mcal_1110.fits'
paths['mice'] = '/global/project/projectdirs/des/y3-bias/mice2/' 
paths['lens_mice'] = paths['mice'] + 'lens.fits'
paths['source_mice'] = paths['mice'] + 'source.fits'
paths['randoms_mice'] = paths['mice'] + 'random.fits'
#paths['yaml'] = 'destest/' 
paths['yaml'] = ''
paths['config_data'] = os.path.join('mastercat_%s'%config_data['mastercat_v'], 'zslim_%s'%config_data['zslim_v'], 'zs_%s'%config_data['zs_v'],
                        'redmagic_%s'%config_data['redmagic_v'], 'zllim_%s'%config_data['zllim_v'], 'njk_%d'%config_data['njk'],
                        'thbin_%0.1f_%d_%d'%(config_data['thlims'][0], config_data['thlims'][1], config_data['nthbins']),
                        'bslop_%0.1g'%config_data['bslop']) 

paths['config_mice'] = os.path.join('mice', 'v_%s'%config_mice['version'], 'zslim_%s'%config_mice['zslim_v'], 'zs_%s'%config_mice['zs_v'],
                        'redmagic_%s'%config_mice['redmagic_v'], 'zllim_%s'%config_mice['zllim_v'], 'njk_%d'%config_mice['njk'],
                        'thbin_%0.1f_%d_%d'%(config_mice['thlims'][0], config_mice['thlims'][1], config_mice['nthbins']),
                        'bslop_%0.1g'%config_mice['bslop']) 

# Where we save the runs and plots for one particular configuration:
paths['runs_config'] = os.path.join(paths['runs'], paths['config_%s'%basic['mode']]) + '/'
paths['plots_config'] = os.path.join(paths['plots'], paths['config_%s'%basic['mode']]) + '/'

"""
ZBINS
---------------------------------
Define the zbins dictionary. This dictionary is imported
in the other scripts and it defines the number of lens and source
redshift bins and their limits. 
"""

zbins = {}
zbins['lbins'] = ['l1', 'l2', 'l3', 'l4', 'l5']
zbins['sbins'] = ['s1', 's2', 's3', 's4']
zbins['lsbins'] = [l + '_' + s for l in zbins['lbins'] for s in zbins['sbins']]
# Updated bins for Y3 of the redmagic sample
zbins['l1'] = [0.15, 0.35]
zbins['l2'] = [0.35, 0.5] 
zbins['l3'] = [0.5, 0.65] 
zbins['l4'] = [0.65, 0.85] 
zbins['l5'] = [0.85, 0.95] 
zbins['lims'] = [zbins['l1'][0], zbins['l2'][0], zbins['l3'][0], zbins['l4'][0], zbins['l5'][0], zbins['l5'][1]]
#zbins['s1'] = [0.90, 1.30]
#zbins['s1'] = [0.20, 0.90]
#zbins['s1'] = [0.20, 0.43]
zbins['s1'] = [0.20, 1.30]
zbins['s2'] = [0.43, 0.63] 
zbins['s3'] = [0.63, 0.90] 
zbins['s4'] = [0.90, 1.30] 
zbins['source_lims'] = [zbins['s1'][0], zbins['s2'][0], zbins['s3'][0], zbins['s4'][0], zbins['s4'][1]]

"""
PLOTTING
---------------------------------
Define the plotting dictionary. This dictionary is imported
in the other scripts and it defines useful quantites that are used
accross several plots, to ensure consistency. 
"""

plotting = {}
if basic['mode'] == 'data':
    plotting['catname'] = r'Metacalibration PSF ' + config['mastercat_v'][0:2]
if basic['mode'] == 'mice':
    plotting['catname'] = r'\textsc{MICE}'

plotting['latex'] = False
plotting['cmap'] = viridis
plotting['redshift_l'] = [r'$%0.2f < z_l < %0.2f $'%(zbins['lims'][i], zbins['lims'][i+1]) for i in range(len(zbins['lims'])-1)]
plotting['redshift_s'] = [r'$%0.2f < z_s < %0.2f $'%(zbins['source_lims'][i], zbins['source_lims'][i+1]) for i in range(len(zbins['source_lims'])-1)]
plotting['titles_redmagic'] = ['redMaGiC HiDens', 'redMaGiC HiDens', 'redMaGiC HiDens', 'redMaGiC HiLum', 'redMaGiC HiLum']
plotting['th_limit'] = [64.,40.,30., 24., 21.] 


"""
SIZE AND S/N NOISE TESTS
---------------------------------
Defines a dictionary used for the size and S/N. 
Currently values from y1. Needs to be updated for Y3.
"""

source_nofz_pars = {}
source_nofz_pars['dzs','size'] = [-0.00071855,0.003097] # true - bpz for [low, high] size
source_nofz_pars['dzs','snr'] = [0.032083,-0.003159] # true - bpz for [low, high] snr
source_nofz_pars['dzs_sigma'] = 1.6 #sqrt(2) times the non-tomographic uncertainty as estimated in Hoyle et al. using COSMOS
source_nofz_pars['thetamin'] = 64.

"""
SYSTEMATICS MAPS TESTS
---------------------------------
Defines a dictionary used for the systematics maps tests.
Currently names for y1. 
"""

sysmaps = {}
sysmaps['nside'] = 4096
sysmaps['nested_bool'] = False # means that it is ring
sysmaps['airmass'] = 'AIRMASS_coaddweights3_mean'
sysmaps['count'] = 'count__fracdet'
sysmaps['exptime'] = 'EXPTIME__total'
sysmaps['fwhm'] = 'FWHM_MEAN_coaddweights3_mean'
sysmaps['maglimit'] = 'maglimit3__'
sysmaps['skybrite'] = 'SKYBRITE_coaddweights3_mean'
sysmaps['skysigma'] = 'SKYSIGMA_coaddweights3_mean'
