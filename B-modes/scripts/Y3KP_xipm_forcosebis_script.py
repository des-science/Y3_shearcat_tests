import yaml
import destest
import treecorr
import numpy as np

# basic dict props
destest_dict_ = {
    'output_exists' : True,
    'use_mpi'       : False,
    'source'        : 'hdf5',
    'dg'            : 0.01
    }

# Populates a full destest yaml dict for each catalog selection based on the limited catalog 
# input info provided in the common cats.yaml file
def create_destest_yaml( params, name, cal_type, group, table, select_path ):
    """
    Creates the input dictionary structure from a passed dictionary rather than\
 reading froma yaml file.
    """
    destest_dict = destest_dict_.copy()
    destest_dict['load_cache'] = params['load_cache']
    destest_dict['output'] = params['output']
    destest_dict['name'] = name
    destest_dict['filename'] = params['datafile']
    destest_dict['param_file'] = params['param_file']
    destest_dict['cal_type'] = cal_type
    destest_dict['group'] = group
    destest_dict['table'] = table
    destest_dict['select_path'] = select_path
    destest_dict['e'] = ['e_1','e_2']
    destest_dict['Rg'] = ['R11','R22']
    destest_dict['w'] = 'weight'

    return destest_dict

# Build selector (and calibrator) classes from destest for the catalog.
def load_catalog(
        pipe_params, name, cal_type, group, table, select_path, inherit=None, return_calibrator=None):
    """
    Loads data access and calibration classes from destest for a given yaml set\
up file.
    """

    # Input yaml file defining catalog
    params = create_destest_yaml(pipe_params, name, cal_type, group, table, select_path)
    # Load destest source class to manage access to file
    source = destest.H5Source(params)

    # Load destest selector class to manage access to data in a structured way
    if inherit is None:
        sel = destest.Selector(params,source)
    else:
        sel = destest.Selector(params,source,inherit=inherit)

    # Load destest calibrator class to manage calibration of the catalog
    if return_calibrator is not None:
        cal = return_calibrator(params,sel)
        return sel, cal
    else:
        return sel

# Read yaml file that defines all the catalog selections used
params = yaml.load(open('cats_cosebis.yaml'))
params['param_file'] = 'cats_cosebis.yaml'

# Source catalog
source_selector, source_calibrator = load_catalog(
    params, 'mcal', 'mcal', params['source_group'], params['source_table'], 
    params['source_path'], return_calibrator=destest.MetaCalib)

# Gold catalog
gold_selector = load_catalog(
    params, 'gold', 'mcal', params['gold_group'], params['gold_table'], 
    params['gold_path'], inherit=source_selector)

# BPZ (or DNF) catalog, depending on paths in cats.yaml file (exchange bpz and dnf)
pz_selector = load_catalog(
    params, 'pz', 'mcal', params['pz_group'], params['pz_table'], 
    params['pz_path'], inherit=source_selector)

# Load ra,dec from gold catalog
ra  = gold_selector.get_col('ra')[0]
dec = gold_selector.get_col('dec')[0]

# Get e1,e2 (all five version)
g1=source_selector.get_col('e_1')
g2=source_selector.get_col('e_2')
R1,c,w = source_calibrator.calibrate('e_1')
R2,c,w = source_calibrator.calibrate('e_2')

corr = treecorr.GGCorrelation(bin_type='Linear',nbins=120000, min_sep=2.5,
                                          max_sep=250.0, sep_units='arcmin',
                                          bin_slop=0.0)

cat_s = treecorr.Catalog(
    ra=ra, dec=dec, g1=(g1[0]-np.average(g1[0],weights=w))/R1, 
    g2=(g2[0]-np.average(g2[0],weights=w))/R2,
    w=w,ra_units='deg', dec_units='deg')

corr.process(cat_s, num_threads=40)  # <--num_threads could be removed or set

theta = corr.meanr
gts = corr.xip
gxs = corr.xim
errsp = np.sqrt(np.abs(corr.varxip))
errsm = np.sqrt(np.abs(corr.varxim))
wts = corr.weight
npairs = corr.npairs
R=np.average((R1,R2))
#print(repr(np.asarray(theta)))
###print(repr(np.asarray(gts/R**2)))
###print(repr(np.asarray(gxs/R**2)))
#print(repr(np.asarray(gts)))
#print(repr(np.asarray(gxs)))
#print(repr(np.asarray(errsp)))
#print(repr(np.asarray(errsm)))

out_ascii='xi_weighted_120klinth2.5_250.0_bslop0.0_final'
fileo = open(out_ascii,'w')
np.savetxt(fileo, np.vstack((theta,gts,gxs, errsp, errsm, wts, npairs)).T)

# End of script
