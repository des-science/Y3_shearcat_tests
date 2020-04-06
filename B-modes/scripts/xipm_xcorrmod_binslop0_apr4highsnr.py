""" Adapted from xcorr's ggl.py
"""
import sys
sys.path.append('/users/PCON0003/cond0080/src/xcorr/src')
sys.path.append('/users/PCON0003/cond0080/src/xcorr/destest')
import psutil
import os
import numpy as np
import astropy.io.fits as pf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.stats
import healpy as hp
import treecorr
import signal
import twopoint
from scipy import interpolate
import functions
import yaml
import destest
#import ipdb

def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

class GGL(object):
    """
    Basic class that has all the functions shared for several tests.
    """

    def __init__(self, basic, config, paths):
        self.basic = basic
        self.config = config
        self.paths = paths

    def load_metacal(self):
        """
        Loads metacal data for Y3 catalog using h5 interface.
        Combines with Gold and BPZ.
        Returns: two dictionaries for the sources, with all the relevant columns and the calibrator
                 in order to be used for other functions too. 
                 These are the two dictionaries it returns, which are used for different things:

        - source_5sels: It is a nested dictionary with 'sheared' and 'unsheared' quantities. 
                        *'unsheared': Dictionary with the unsheared version of each 
                                      quantity with the selections from: unsheared, 1p, 1m, 2p, 2m.
                                      The unsheared version of a quantity means that has been measured on
                                      images which are unsheared, or using fluxes obtained from unsheared images.
                                      This dictionary is useful to compute the selection response manually.

                        *'sheared': Dictionary with the 5 selections (1p, 1m, 2p, 2m), in which the 
                                    quantities are obtained from sheared images or using fluxes measured
                                    on sheared images. 
                                    This dictionary is useful to make further selections besides the baseline selection,
                                    for instance to select of redshift bins, size or S/N, which you need to do also 
                                    using the sheard quanties to compute the selection response. 
        
        - source: Simplest one. Dictionary with the baseline unsheared selection and quantites.
        
        """

        mcal_file = self.paths['yaml'] + 'destest_mcal_highsnr.yaml'
        params_mcal = yaml.load(open(mcal_file))
        params_mcal['param_file'] = mcal_file
        params_mcal['filename'] = self.config['filename_mastercat']
        source_mcal = destest.H5Source(params_mcal)
        source_selector = destest.Selector(params_mcal,source_mcal)
        source_calibrator = destest.MetaCalib(params_mcal,source_selector)

        gold_file = self.paths['yaml'] + 'destest_gold.yaml'
        params_gold = yaml.load(open(gold_file))
        params_gold['param_file'] = gold_file
        params_gold['filename'] = self.config['filename_mastercat']
        source_gold = destest.H5Source(params_gold)
        gold_selector = destest.Selector(params_gold,source_gold,inherit=source_selector)

        param_file = self.paths['yaml'] + './destest_pz.yaml'
        params_pz = yaml.load(open(param_file))
        params_pz['filename'] = self.config['filename_mastercat']
        source_pz = destest.H5Source(params_pz)
        pz_selector = destest.Selector(params_pz, source_pz, inherit=source_selector)


        # Dictionary with the unsheared version of each quantity with the selections from: unsheared, 1p, 1m, 2p, 2m. 
        source_5sels = {}
        source_5sels['unsheared'] = {}
        source_5sels['unsheared']['ra'] = [gold_selector.get_col('ra', uncut=True)[0][gold_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['dec'] = [gold_selector.get_col('dec', uncut=True)[0][gold_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['e1'] = [source_selector.get_col('e_1', uncut=True)[0][source_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['e2'] = [source_selector.get_col('e_2', uncut=True)[0][source_selector.get_mask()[i]] for i in range(5)]

        # Dictionary with the 5 selections (1p, 1m, 2p, 2m), in which the quantities are obtained from sheared images or using fluxes measured on sheared images.   
        source_5sels['sheared'] = {}
        source_5sels['sheared']['e1'] = source_selector.get_col('e_1')
        source_5sels['sheared']['e2'] = source_selector.get_col('e_2')
        source_5sels['sheared']['snr'] = source_selector.get_col('snr')
        source_5sels['sheared']['size'] = source_selector.get_col('T')
        #source_5sels['sheared']['bpz_mean'] = pz_selector.get_col('zmean_sof')
        #source_5sels['sheared']['bpz_zmc'] = pz_selector.get_col('zmc_sof')
        # This change updates to SOM pz selection using bhat
        source_5sels['sheared']['somz_bhat'] = pz_selector.get_col('bhat')

        # Dictionary with the unsheared version and selection only:
        source = {}
        source['ra'] = source_5sels['unsheared']['ra'][0]
        source['dec'] = source_5sels['unsheared']['dec'][0]
        source['e1'] = source_5sels['unsheared']['e1'][0]
        source['e2'] = source_5sels['unsheared']['e2'][0]
        source['psf_e1'] = source_selector.get_col('psf_e1')[0]
        source['psf_e2'] = source_selector.get_col('psf_e2')[0]
        source['snr'] = source_5sels['sheared']['snr'][0]
        source['size'] = source_5sels['sheared']['size'][0]
        #source['bpz_mean'] = source_5sels['sheared']['bpz_mean'][0]
        #source['bpz_zmc'] = source_5sels['sheared']['bpz_zmc'][0]
        source['somz_bhat'] = source_5sels['sheared']['somz_bhat'][0]

        calibrator = destest.MetaCalib(params_mcal, source_selector)
        if 'v1' in self.config['mastercat_v']:
            R11, _, _ = calibrator.calibrate('e1')
            R22, _, _ = calibrator.calibrate('e2')
            source['Rgamma'] = calibrator.calibrate('e1', return_wRg=True)

        #if ('v2' in self.config['mastercat_v'])|('_subsampled' in self.config['mastercat_v']):
        #if '_subsampled' in self.config['mastercat_v']:
        else:
            R11, _, _ = calibrator.calibrate('e_1')
            R22, _, _ = calibrator.calibrate('e_2')
            # Note: This returns the Rgamma per galaxy already averaged over e1, e2.
            # So source['Rmean'] is equal to source['Rgamma'].mean()
            source['Rgamma'] = calibrator.calibrate('e_1', return_wRg=True)
            source['R11'] = calibrator.calibrate('e_1', return_full=True)[0]
            source['R22'] = calibrator.calibrate('e_2', return_full=True)[0]

        source['Rmean'] = np.mean([R11, R22])
        #print np.mean(np.mean(source['R11']), 
        #              np.mean(source['R22']))
        print np.mean(source['R11']), np.mean(source['R22'])
        print 'Response full sample', source['Rmean']

        return source, source_5sels, calibrator

    def load_metacal_bin(self, source, source_5sels, calibrator, zlim_low, zlim_high):
        """
        source: dictionary containing relevant columns for the sources, with the baseline selection applied already.
        source_5sels: dictionary containing relevant columns for the sources, with the baseline selection applied already,
                     for each of the 5 selections 1p, 1m, 2p, 2m. 
        calibrator: class to compute the response. Taken from baseline selection.
        zlim_low, zlim_high: limits to select the tomographic bin.
        Obtains 5 masks (unsheared, sheared 1p, 1m, 2p, 2m) to obtain the new selection response.
        Returns: Source dictionary masked with the unsheared mask and with the mean response updated.
        """
        #photoz_masks = [(source_5sels['sheared']['bpz_mean'][i] > zlim_low) & (source_5sels['sheared']['bpz_mean'][i] < zlim_high) for i in range(5)]
        #photoz_masks = [(source_5sels['sheared']['somz_bhat'][i] >= 0) & (
        #        source_5sels['sheared']['somz_bhat'][i] <= 3) for i in range(5)]
        photoz_masks = [(source_5sels['sheared']['snr'][i]>=0) for i in range(5)]

        source_bin = {}
        source_bin['ra'] = source['ra'][photoz_masks[0]]
        source_bin['dec'] = source['dec'][photoz_masks[0]]
        source_bin['e1'] = source['e1'][photoz_masks[0]]
        source_bin['e2'] = source['e2'][photoz_masks[0]]
        source_bin['psf_e1'] = source['psf_e1'][photoz_masks[0]]
        source_bin['psf_e2'] = source['psf_e2'][photoz_masks[0]]
        source_bin['snr'] = source['snr'][photoz_masks[0]]
        source_bin['size'] = source['size'][photoz_masks[0]]
        #source_bin['bpz_mean'] = source['bpz_mean'][photoz_masks[0]]
        #source_bin['bpz_zmc'] = source['bpz_zmc'][photoz_masks[0]]
        source_bin['somz_bhat'] = source['somz_bhat'][photoz_masks[0]]
        source_bin['Rgamma'] = source['Rgamma'][photoz_masks[0]]

        if 'v1' in self.config['mastercat_v']:
            R11, _, _ = calibrator.calibrate('e1', mask=photoz_masks)
            R22, _, _ = calibrator.calibrate('e2', mask=photoz_masks)
        #if 'subsampled' in self.config['mastercat_v']:
        else:
            R11, _, _ = calibrator.calibrate('e_1', mask=photoz_masks)
            R22, _, _ = calibrator.calibrate('e_2', mask=photoz_masks)
            source_bin['R11'] = calibrator.calibrate('e_1', return_full=True, mask=photoz_masks)[0]
            source_bin['R22'] = calibrator.calibrate('e_2', return_full=True, mask=photoz_masks)[0]
            
        source_bin['Rmean'] = np.mean([R11, R22])
        print 'Mean response redshift bin (%0.2f, %0.2f):'%(zlim_low, zlim_high), source_bin['Rmean'], np.mean(source_bin['Rgamma']), np.mean(source_bin['R11']), np.mean(source_bin['R22'])
        return source_bin


    def load_metacal_bin_sels_responses(self, source_5sels, zlim_low, zlim_high):
        """
        source_5sels: Nested dictionary containing relevant columns for the sources, with the baseline selection applied already,
                   for each of the unsheared and sheared versions and for each selection.
        zlim_low, zlim_high: limits to select the tomographic bin.
        Returns: Source dictionary masked each time with one of the photo-z masks, to compute the selection response manually and
                 test its scale dependence. 
        """
        photoz_masks = [(source_5sels['sheared']['bpz_mean'][i] > zlim_low) & (source_5sels['sheared']['bpz_mean'][i] < zlim_high) for i in range(5)]
        source_bin_sels = {}
        source_bin_sels['ra'] = [source_5sels['unsheared']['ra'][i][photoz_masks[i]] for i in range(1, 5)]
        source_bin_sels['dec'] = [source_5sels['unsheared']['dec'][i][photoz_masks[i]] for i in range(1, 5)]
        source_bin_sels['e1'] = [source_5sels['unsheared']['e1'][i][photoz_masks[i]] for i in range(1, 5)]
        source_bin_sels['e2'] = [source_5sels['unsheared']['e2'][i][photoz_masks[i]] for i in range(1, 5)]

        return source_bin_sels


    def get_source(self, source):
        """
        Given a source sample, returns ra, dec, and weight, in case it exists.
        """

        ra_s = source['ra']
        dec_s = source['dec']
        try:
            w = source['w']
            print 'Weights found in source catalog.'
        except:
            print 'There are no identified weights for the sources.'
            w = np.ones(len(ra_s))

        return ra_s, dec_s, w

    def run_treecorr_jackknife(self, source, meanR, type_corr):
        """
        Function that runs treecorr for a given lens and source sample,
        and a given configuration, and paralalizes the measurements
        in different jackknife patches for the lenses and randoms.
        Returns the measurements for each jackknife patch.
        type_corr: string, type of correlation, i.e. NG, NN, NK_Rgamma
        NG for gammat, NN for wtheta, NK for scalar quantities, like
        the responses. The string after NK_ is used to load the column.
        """
        assert type_corr == 'NG' or type_corr == 'NN' or 'NK' in type_corr, 'This type_corr of correlation is not accepted by this function.'

        if type_corr == 'NG':
            corr = treecorr.GGCorrelation(bin_type='Linear',nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
            #corr = treecorr.GGCorrelation(bin_type='Log',nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                          max_sep=self.config['thlims'][1], sep_units='arcmin',
                                          bin_slop=self.config['bslop'])

        if type_corr == 'NN':
            if jk == 0: print 'Doing NN correlation.'
            corr = treecorr.NNCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                          max_sep=self.config['thlims'][1], sep_units='arcmin',
                                          bin_slop=self.config['bslop'])
                
        if 'NK' in type_corr:
            if jk == 0: print 'Doing NK correlation with variable %s.'%type_corr[3:]
            corr = treecorr.NKCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                          max_sep=self.config['thlims'][1], sep_units='arcmin',
                                          bin_slop=self.config['bslop'])
            
        ra_s, dec_s, w = self.get_source(source)
        if type_corr == 'NG' or type_corr == 'NN':
            e1 = source['e1']
            e2 = source['e2']
            # We still need to define these variables because otherwise
            # multiprocessing gives an error (cell is empty), because one has t\o define
            # all variables even when they are not needed for correct execution\.
            scalar = np.zeros(len(e1))
            #print 'Doing correlations with responses!!!! e1=R11 and e2=R22!!'
        if 'NK' in type_corr:
            print type_corr, type_corr[3:]
            scalar = source['%s'%type_corr[3:]]
            # We still need to define these variables because otherwise
            # multiprocessing gives an error (cell is empty), because one has t\o define
            # all variables even when they are not needed for correct execution\.
            e1 = np.zeros(len(scalar))
            e2 = np.zeros(len(scalar))

        if type_corr == 'NG' or type_corr == 'NN':
            if self.basic['mode'] == 'data' or self.basic['mode'] =='data_y1sources':
                cat_s = treecorr.Catalog(
                    ra=ra_s, dec=dec_s, g1=(e1-e1.mean()), g2=(e2-e2.mean()),
                    w=w,ra_units='deg', dec_units='deg')

        if 'NK' in type_corr:
            cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, k=scalar, w=w, 
                                     ra_units='deg', dec_units='deg')
 
        corr.process(cat_s)#, num_threads=4)
        #theta = corr.logr
        theta = corr.meanr
        if type_corr == 'NG':
            gts = corr.xip
            gxs = corr.xim
            errsp = np.sqrt(np.abs(corr.varxip))
            errsm = np.sqrt(np.abs(corr.varxim))
            wts = corr.weight
            npairs = corr.npairs
        else:
            print "Dying"
            df

        return theta, gts, gxs, errsp, errsm, wts, npairs


    def process_run(self, all, theta, path_test, end):
        """
        From the jackknife measurements in all jackknife regions but all, constructs covariance, mean and stats.
        Saves them into file.
        all: gt_all or gx_all.
        theta: in arcmin.
        path_test: where to save the files.
        end: string to save the files: gt, gx, randoms, etc.
        """
        mean = np.mean(all, axis=0)
        cov = functions.covariance(all, mean)
        err = np.sqrt(np.diag(cov))
        chi2 = np.dot(mean.T, np.dot(np.linalg.inv(cov), mean))
        ndf = len(mean)

        std = np.sqrt((len(all) - 1.)) * np.std(all, axis=0)

        # Hartlap factor
        N = len(all)  # number of jackknife regions
        p = len(mean)  # number of angular bins
        factor = (N - p - 2) / float(N - 1)
        chi2_hartlap = chi2 * factor

        stats = np.array([chi2_hartlap, ndf])

        np.savetxt(path_test + 'mean_%s' % end, zip(theta, mean, err), header='th, %s, err_%s' % (end, end))
        np.savetxt(path_test + 'cov_%s' % end, cov)
        np.savetxt(path_test + 'all_%s' % end, all,
                   header='%s (sum of %s from all patches except for one, different each time)' % (end, end))
        np.savetxt(path_test + 'null_chi2_%s' % end, stats.reshape(1, stats.shape[0]),
                   fmt='%0.1f  %d', header='chi2_hartlap  ndf')


    def load_data_or_sims(self):
        '''
        Loads and returns lens, randoms and sources, either from data or from simulations.
        '''

        if self.basic['mode'] == 'data':
            #lens_all = pf.getdata(self.paths['lens'])
            #random_all = pf.getdata(self.paths['randoms'])
            source_all, source_all_5sels, calibrator = self.load_metacal()

        if self.basic['mode'] == 'data_y1sources':
            lens_all = pf.getdata(self.paths['lens'])
            random_all = pf.getdata(self.paths['randoms'])

        if self.basic['mode'] == 'mice':
            lens_all = pf.getdata(self.paths['lens_mice'])
            random_all = pf.getdata(self.paths['randoms_mice'])
            source_all = pf.getdata(self.paths['source_mice'])

        return source_all, source_all_5sels, calibrator


class Measurement(GGL):
    """
    Subclass that runs the gglensing measurement for all the lens-source bin combinations.
    Includes:
    - Mean response calculation.
    - Random points subtraction.
    - Jackknife covariance calculation.
    - Boost factors calculation.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self, sbin):
        return os.path.join(self.paths['runs_config'], 'measurement', sbin) + '/'

    def get_path_test_allzbins(self):
        return os.path.join(self.paths['runs_config'], 'measurement') + '/'

    def get_twopointfile_name(self):
        return os.path.join(self.get_path_test_allzbins() + 'gammat_twopointfile.fits')

    def run(self):

        source_all, source_all_5sels, calibrator = self.load_data_or_sims()
        for sbin in self.zbins['sbins'][:1]:
        #for sbin in self.zbins['sbins']:
    		print 'Running measurement for source %s.' % sbin

		if self.basic['mode'] == 'data':
		    source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, zlim_low=self.zbins[sbin][0], zlim_high=self.zbins[sbin][1])
		    R = source['Rmean']
                    R11 = source['R11']
                    R22 = source['R22']

                path_test = self.get_path_test(sbin)
                make_directory(path_test)

                theta, gts, gxs, errsp, errsm, wts, npairs = self.run_treecorr_jackknife(source, R, 'NG')

                #gt_all = (gtnum / wnum) / R - (gtnum_r / wnum_r) / R
                #gx_all = (gxnum / wnum) / R - (gxnum_r / wnum_r) / R
                print repr(np.asarray(theta))
                print repr(np.asarray(gts/R**2))
                print repr(np.asarray(gxs/R**2))
                print repr(np.asarray(errsp))
                print repr(np.asarray(errsm))
                # theta, xi_p, xi_m, random 5 cols "xi_zall_nBins_1_Bin1_Bin1"
                #fileo = open('xi_z0.2_1.3_60kth4.0_250.0_bslop0','w')
                #fileo = open('xi_z0.2_1.3_120kth4.0_250.0_bslop0','w')
                #fileo = open('xi_z0.2_1.3_1kth2.5_250.0_bslop0','w')
                fileo = open('xi_z0.2_1.3_120klinth2.5_250.0_bslop0_highsnr_apr4','w')
                np.savetxt(fileo, np.vstack((theta, (gts/(R**2)), 
                                             (gxs/(R**2)), errsp,
                                             errsm, wts, npairs)).T)

                #self.process_run(gt_all, theta, path_test, 'gt')
                #self.process_run(gx_all, theta, path_test, 'gx')
                #self.process_run((gtnum_r / wnum_r) / R, theta, path_test, 'randoms')


    def save_2pointfile(self, string):
        """
        Save the correlation function measurements (i.e. gammat, boost factors, etc),
        N(z)'s and jackknife covariance into the 2point format file.
        Creates a:
        - SpectrumMeasurement obejct: In which the 2PCF measurements are saved. 
        - Kernel object: The N(z)'s are here. 
        - CovarianceMatrixInfo object: In which the jackknife covariance is saved.

        Then, it builds the TwoPointFile objet from the above objects,
        and saves it to a file.
        """
        
        # Preparing spectrum
        length = self.config['nthbins'] * len(self.zbins['lbins']) * len(self.zbins['sbins'])
        values = np.zeros(length, dtype=float)
        bin1 = np.zeros(length, dtype=int)
        bin2 = np.zeros_like(bin1)
        angular_bin = np.zeros_like(bin1)
        angle = np.zeros_like(values)
        dv_start = 0
        cov = np.zeros((length, length))

        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, xi, xi_err = np.loadtxt(path_test + 'mean_%s'%string, unpack=True)
                cov_ls = np.loadtxt(path_test + 'cov_%s'%string)

                bin_pair_inds = np.arange(dv_start, dv_start + self.config['nthbins'])
                values[bin_pair_inds] = xi
                bin1[bin_pair_inds] = l+1
                bin2[bin_pair_inds] = s+1
                angular_bin[bin_pair_inds] = np.arange(self.config['nthbins'])
                angle[bin_pair_inds] = theta
                dv_start += self.config['nthbins']
                cov[bin_pair_inds[0]:bin_pair_inds[-1]+1, bin_pair_inds] = cov_ls

        # Preparing N(z) for the blinding script
        if 'y1_2pt_NG_mcal_1110' in self.paths['lens_nz']:
            file_lens_nz = twopoint.TwoPointFile.from_fits(self.paths['lens_nz'])
            lens_nz = file_lens_nz.get_kernel('nz_lens')

        if 'stellar_mass' in self.paths['lens_nz']:
            import astropy
            ext = astropy.io.fits.open(self.paths['lens_nz'])['nz_lens']
            lens_nz = twopoint.NumberDensity.from_fits(ext)

        if 'y1_2pt_NG_mcal_1110' in self.paths['source_nz']:
            file_source_nz = twopoint.TwoPointFile.from_fits(self.paths['source_nz'])
            source_nz = file_source_nz.get_kernel('nz_source')

        if string == 'gt':
            gammat = twopoint.SpectrumMeasurement('gammat', (bin1, bin2),
                                                  (twopoint.Types.galaxy_position_real,
                                                   twopoint.Types.galaxy_shear_plus_real),
                                                  ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, values,
                                                  angle=angle, angle_unit='arcmin')

            cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['gammat'], [length], cov)
            print 'Saving TwoPointFile with lens N(z) from %s'%(self.paths['lens_nz'])
            print 'Saving TwoPointFile with source N(z) from %s'%(self.paths['source_nz'])
            gammat_twopoint = twopoint.TwoPointFile([gammat], [lens_nz, source_nz], windows=None, covmat_info=cov_mat_info)
            twopointfile_unblind = self.get_twopointfile_name()

            # Remove file if it exists already because to_fits function doesn't overwrite
            if os.path.isfile(twopointfile_unblind):
                os.system('rm %s'%(twopointfile_unblind))

            gammat_twopoint.to_fits(twopointfile_unblind)


        if string == 'boost_factor':
            boost_factor = twopoint.SpectrumMeasurement('boost_factor', (bin1, bin2),
                                                  (twopoint.Types.galaxy_position_real,
                                                   twopoint.Types.galaxy_shear_plus_real),
                                                  ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, values,
                                                  angle=angle, angle_unit='arcmin')

            cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['boost_factor'], [length], cov)
            print 'Saving TwoPointFile with Y1 N(z)s'
            boost_factor_twopoint = twopoint.TwoPointFile([boost_factor], [y1_lensnz, y1_sourcenz], windows=None, covmat_info=cov_mat_info)
            save_path = os.path.join(self.get_path_test_allzbins() + '%s_twopointfile.fits'%string)
            boost_factor_twopoint.to_fits(save_path)


    def save_spectrum_measurement_file(self):
        """
        Save the boost factors into a SpectrumMeasurment file, from the twopoint code.
        Caveat: cannot save the errors, use TwoPointFile instead (see save_2pointfile function above).
        """

        bf_length = self.config['nthbins'] * len(self.zbins['lbins']) * len(self.zbins['sbins'])
        bf_values = np.zeros(bf_length, dtype=float)
        bin1 = np.zeros(bf_length, dtype=int)
        bin2 = np.zeros_like(bin1)
        angular_bin = np.zeros_like(bin1)
        angle = np.zeros_like(bf_values)
        dv_start = 0

        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, bf, bf_err = np.loadtxt(path_test + 'mean_boost_factor', unpack=True)

                bin_pair_inds = np.arange(dv_start, dv_start + self.config['nthbins'])
                bf_values[bin_pair_inds] = bf
                bin1[bin_pair_inds] = l
                bin2[bin_pair_inds] = s
                angular_bin[bin_pair_inds] = np.arange(self.config['nthbins'])
                angle[bin_pair_inds] = theta
                dv_start += self.config['nthbins']

        boost_factors = twopoint.SpectrumMeasurement('boost_factors', (bin1, bin2),
                                                     (twopoint.Types.galaxy_position_real,
                                                      twopoint.Types.galaxy_shear_plus_real),
                                                     ['no_nz', 'no_nz'], 'SAMPLE', angular_bin, bf_values,
                                                     angle=angle, angle_unit='arcmin', extra_cols=None)

        path_save = self.get_path_test_allzbins()
        filename = 'boost_factors.fits'
        (boost_factors.to_fits()).writeto(path_save + filename)


    def plot_NO_BLINDED(self):
        """"
        Makes plot of the fiducial measurement for all redshift bins.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        title_source = self.plotting['catname']

        # Figure
        fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharey=True, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for l in range(0, len(self.zbins['lbins'])):

            # To iterate between the three columns and two lines
            j = 0 if l < 3 else 1
            ax[j][l % 3].axvspan(self.config['thlims'][0]*0.8, self.plotting['th_limit'][l], color='gray', alpha=0.2)
            for s in range(len(self.zbins['sbins'])):

                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                if os.path.isfile(path_test + 'mean_gt'):
                    th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)

                    mask_neg = gt < 0
                    mask_pos = gt > 0

                    #chi2, ndf = self.get_chi2(path_test, 'gt')
                    ax[j][l % 3].errorbar(th[mask_neg] * (1 + 0.05 * s), -gt[mask_neg], err[mask_neg], fmt='.', mfc='None',
                                          mec=plt.get_cmap(cmap)(cmap_step * s), ecolor=plt.get_cmap(cmap)(cmap_step * s), capsize=2)
                    ax[j][l % 3].errorbar(th[mask_pos] * (1 + 0.05 * s), gt[mask_pos], err[mask_pos], fmt='.',
                                          color=plt.get_cmap(cmap)(cmap_step * s),
                                          mec=plt.get_cmap(cmap)(cmap_step * s), label=self.plotting['redshift_s'][s], capsize=2)

                    ax[j][l % 3].set_xlim(self.config['thlims'][0]*0.8, self.config['thlims'][1]*1.2)
                    ax[j][l % 3].set_xscale('log')
                    ax[j][l % 3].set_yscale('log')

                    ax[j][l % 3].text(0.5, 0.9, self.plotting['redshift_l'][l], horizontalalignment='center',
                                      verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)
                    #ax[j][l % 3].text(0.5, 0.93, self.plotting['titles_redmagic'][l], horizontalalignment='center',
                    #                  verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)

                    #if l % 3 > 0:  # In case we want to keep labels on the left y-axis
                    ax[j][l % 3].yaxis.set_ticklabels([])  # to remove the ticks labels
                    if l < 2:
                        ax[0][l].xaxis.set_ticklabels([])  # to remove the ticks labels

                    ax[j][l % 3].set_xlabel(r'$\theta$ [arcmin]', size='large')
                    ax[j][0].set_ylabel(r'$\gamma_t (\theta)$', size='large')
                    #ax[j][l % 3].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

                    """
                    # Chi2
                    ax[j][l%3].text(0.25,0.3,r'Null $\chi^2$/ndf',
                                    horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 10)
                    ax[j][l%3].text(0.25,0.23 -0.06*s,r'$%0.1f/%d$'%(chi2, ndf),
                             horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 12, color = plt.get_cmap(cmap)(cmap_step*s))
                    """

        #ax[1][0].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        #ax[1][1].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        # handles, labels = ax[0][0].get_legend_handles_labels()
        fig.delaxes(ax[1, 2])
        # ax[1][1].legend(handles[::-1], labels[::-1], frameon=True, fancybox = True,prop={'size':12}, numpoints = 1, loc='center left', bbox_to_anchor=(1, 0.5))
        ax[0][0].legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='center',
                        bbox_to_anchor=(2.45, -0.52))
        fig.suptitle(title_source, fontsize=16)
        fig.subplots_adjust(top=0.93)
        self.save_plot('plot_measurement')

    def load_twopointfile(self):
        '''
        Loads TwoPointFile and returns it.
        '''
        filename = self.get_twopointfile_name()
        if self.basic['blind']: gammat_file = twopoint.TwoPointFile.from_fits('%s_BLINDED.fits'%filename[:-5])
        else: gammat_file = twopoint.TwoPointFile.from_fits('%s.fits'%filename[:-5])
        return gammat_file

    def compute_sn_ratio(self):
        '''
        Compute S/N ratio using null chi2. 
        S/N = sqrt(null chi2 - ndof)
        Uses full jackknife covariance.
        '''

        gammat_file = self.load_twopointfile()
        gammat = gammat_file.spectra[0].value
        cov = gammat_file.covmat
        
        COV = np.array(cov)
        INVCOV = np.linalg.inv(COV)
        null_chi2 = np.mat(gammat)*INVCOV*np.transpose(np.mat(gammat))
        dof = gammat.shape[0]
        sn = np.sqrt(null_chi2 - dof)
        
        path_save = self.get_path_test_allzbins()
        print 'S/N of the full measurements:', float(sn)
        np.savetxt(path_save + 'sn_ratio_full_measurements', sn, fmt = '%0.5g', header = 'S/N ratio computed with JK covariance, S/N = sqrt(null chi2 - ndof)')


    def plot_from_twopointfile(self):
        
        """"
        Makes plot of the fiducial measurement for all redshift bins, from a twopoint file, like Y1 style.
        It also uses the twopoint plotting functions to plot the measurements in a different style,
        the covariance and the N(z)'s.
        Useful to plot the blinded measurements (now the default). 
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        gammat_file = self.load_twopointfile()

        gammat = gammat_file.spectra[0]
        pairs = gammat.bin_pairs
        npairs = len(pairs)
        # It starts with 1, not with 0
        bins_l = np.transpose(pairs)[0]
        bins_s = np.transpose(pairs)[1]
        nbins_l = np.max(bins_l)
        nbins_s = np.max(bins_s)
        assert len(self.zbins['lbins']) == nbins_l, 'Number of lens bins in info does not match with the one in the two-point file.'
        assert len(self.zbins['sbins']) == nbins_s, 'Number of source bins in info does not match with the one in the two-point file.'

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        title_source = self.plotting['catname']

        # Figure
        fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharey=False, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for l in range(0, len(self.zbins['lbins'])):

            # To iterate between the three columns and two lines
            j = 0 if l < 3 else 1
            ax[j][l % 3].axvspan(self.config['thlims'][0]*0.8, self.plotting['th_limit'][l], color='gray', alpha=0.2)

            for s in range(len(self.zbins['sbins'])):

                    path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                    th, gt = gammat.get_pair(l+1, s+1)
                    err = gammat.get_error(l+1, s+1)

                    mask_neg = gt < 0
                    mask_pos = gt > 0

                    ax[j][l % 3].errorbar(th[mask_neg] * (1 + 0.05 * s), -gt[mask_neg], err[mask_neg], fmt='.', mfc='None',
                                          mec=plt.get_cmap(cmap)(cmap_step * s), ecolor=plt.get_cmap(cmap)(cmap_step * s), capsize=2)
                    ax[j][l % 3].errorbar(th[mask_pos] * (1 + 0.05 * s), gt[mask_pos], err[mask_pos], fmt='.',
                                          color=plt.get_cmap(cmap)(cmap_step * s),
                                          mec=plt.get_cmap(cmap)(cmap_step * s), label=self.plotting['redshift_s'][s], capsize=2)

                    ax[j][l % 3].set_xlim(self.config['thlims'][0]*0.8, self.config['thlims'][1]*1.2)
                    ax[j][l % 3].set_xscale('log')
                    ax[j][l % 3].set_yscale('log')

                    ax[j][l % 3].text(0.5, 0.9, self.plotting['redshift_l'][l], horizontalalignment='center',
                                      verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)
                    #ax[j][l % 3].text(0.5, 0.93, self.plotting['titles_redmagic'][l], horizontalalignment='center',
                    #                  verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)

                    #if l % 3 > 0:  # In case we want to keep labels on the left y-axis
                    ax[j][l % 3].yaxis.set_ticklabels([])  # to remove the ticks labels
                    if l < 2:
                        ax[0][l].xaxis.set_ticklabels([])  # to remove the ticks labels

                    ax[j][l % 3].set_xlabel(r'$\theta$ [arcmin]', size='large')
                    ax[j][0].set_ylabel(r'$\gamma_t (\theta)$', size='large')
                    ax[j][l % 3].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

                    """
                    # Chi2
                    ax[j][l%3].text(0.25,0.3,r'Null $\chi^2$/ndf',
                                    horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 10)
                    ax[j][l%3].text(0.25,0.23 -0.06*s,r'$%0.1f/%d$'%(chi2, ndf),
                             horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 12, color = plt.get_cmap(cmap)(cmap_step*s))
                    """

        fig.delaxes(ax[1, 2])
        # ax[1][1].legend(handles[::-1], labels[::-1], frameon=True, fancybox = True,prop={'size':12}, numpoints = 1, loc='center left', bbox_to_anchor=(1, 0.5))
        ax[0][0].legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='center',
                        bbox_to_anchor=(2.45, -0.52))
        fig.suptitle(title_source, fontsize=16)
        fig.subplots_adjust(top=0.93)
        self.save_plot('plot_measurement_BLINDED')

        # Use twopoint library to make the rest of the plots
        gammat_file.plots(self.paths['plots_config'] + 'gammat_twopointfile_BLINDED', blind_yaxis=self.basic['blind'], latex = self.plotting['latex'])


    def plot_boostfactors(self):

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        fig, ax = plt.subplots(4, 5, figsize=(12.5, 9.375), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        c1 = plt.get_cmap(cmap)(0.)

        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):

                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, bf, bf_err = np.loadtxt(path_test + 'mean_boost_factor', unpack=True)

                ax[s][l].axhline(y=1, ls=':', color='k', alpha=1)
                ax[s][l].errorbar(theta, bf, bf_err, fmt='.', color=c1, mec=c1, capsize=2)
                ax[s][l].set_xscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'Boost factors', size='larger',
                                        linespacing=3)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], size='larger')

                ax[s][l].axvspan(2.5, self.plotting['th_limit'][l], color='gray', alpha=0.2)

        ax[0][4].legend(frameon=False, fontsize=16, loc='lower right')
        self.save_plot('boost_factors')

    def plot_randoms(self):
        """
        Makes plot of the tangential shear around random points.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        labels = [self.plotting['catname']]
        c = 0  # If adding im3shape, s=1
        markers = ['o', '^']
        fs = 18  # fontsize
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        fig, ax = plt.subplots(4, 5, figsize=(17.25, 13.8), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.0, wspace=0.0)

        for l in range(0, len(self.zbins['lbins'])):

            for s in range(len(self.zbins['sbins'])):

                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                th, gt, err = np.loadtxt(path_test + 'mean_randoms', unpack=True)

                ax[s][l].axhline(y=0, ls=':', color='k')
                ax[s][l].errorbar(th * (1 + 0.07 * c), gt, err, fmt=markers[c], color=plt.get_cmap(cmap)(cmap_step * s),
                                  mec=plt.get_cmap(cmap)(cmap_step * s), label=labels[c], capsize=1.3)

                ax[s][l].set_xlim(2.5, 300)
                ax[s][l].set_ylim(-2.3 * 10 ** (-4), 2.3 * 10 ** (-4))
                ax[s][l].set_xscale('log')

                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize=fs)
                ax[s][l].tick_params(axis='both', which='minor', labelsize=fs)

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', fontsize=fs)
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'$\gamma_t (\theta)$', fontsize=fs,
                                        linespacing=3)
                    ax[s][l].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                if s > 0:
                    ax[s][l].yaxis.get_offset_text().set_visible(False)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], fontsize=fs)

        ax[0][4].legend(frameon=False, fancybox=True, prop={'size': 13}, numpoints=1, loc='upper right')
        self.save_plot('plot_randoms')

    def plot_gammax(self):
        """
        Makes plot of the cross-component gammax.
        Top panel: Plot of the gammax measurement for a single lens-source combination.
        Bottom panel: chi2 distribution from all lens-source combinations.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        labels = [self.plotting['catname']]
        c = 0  # for metacal
        markers = ['o', '^']
        fs = 12  # fontsize
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)
        colors = [c1, c2]
        fig, ax = plt.subplots(2, 1, figsize=(5, 8), sharey=False, sharex=False)
        fig.subplots_adjust(hspace=0.2, wspace=0.)

        # TOP panel
        l = 0
        s = 0
        path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
        th, gx, err = np.loadtxt(path_test + 'mean_gx', unpack=True)
        ax[0].axhline(y=0, ls=':', color='k')
        ax[0].errorbar(th * (1 + 0.07 * c), gx * th, err * th, fmt=markers[c], color=colors[c],
                       mec=colors[c], label=labels[c], capsize=1.3)
        ax[0].set_xlim(2.5, 300)
        ax[0].set_ylim(-8 * 10 ** (-3), 8 * 10 ** (-3))
        ax[0].set_xscale('log')

        ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
        ax[0].set_xlabel(r'$\theta$ [arcmin]', fontsize=fs)
        ax[0].set_ylabel(r'$\gamma_{\times} \times \theta$', fontsize=fs)
        ax[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[0].text(0.5, 0.84, self.plotting['redshift_l'][l], horizontalalignment='center',
                   verticalalignment='center', transform=ax[0].transAxes, fontsize=12)
        ax[0].text(0.5, 0.92, self.plotting['redshift_s'][s], horizontalalignment='center',
                   verticalalignment='center', transform=ax[0].transAxes, fontsize=12)
        ax[0].legend(frameon=False, fancybox=True, prop={'size': 10}, numpoints=1, loc='lower left')

        # BOTTOM panel
        chi2s = []
        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                chi2, ndf = self.get_chi2(path_test, 'gx')
                chi2s.append(chi2)
        chi2s = np.array(chi2s)
        print len(chi2s)
        ax[1].hist(chi2s, bins=7, color=colors[c], ec=colors[c], lw=2, normed=True, histtype='step', alpha=1,
                   label=labels[c])

        # Chi2 distribution
        ndf = len(gx)
        t = np.linspace(0., 60, 300)
        ax[1].fill_between(t, 0, scipy.stats.chi2.pdf(t, ndf), color='gray', alpha=0.5,
                           label='$\chi^2$ pdf (ndf$ = %d$)' % ndf)
        ax[1].set_ylim(0, 0.1)
        ax[1].legend(frameon=False, fancybox=True, prop={'size': 10}, numpoints=1, loc='best')
        ax[1].set_xlabel(r'$\chi^2_\mathrm{Null}$')
        self.save_plot('plot_gammax')


class ResponsesScale(GGL):
    """
    Subclass that obtains the scale dependence responses (NK correlations) for all the lens-source bin combinations. Both for randoms and lenses.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test_allzbins(self):
        return os.path.join(self.paths['runs_config'], 'responses_nk') + '/'

    def get_path_test(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'responses_nk', lbin + '_' + sbin) + '/'

    def save_responses_nk(self, path_test, responses, end):
        np.savetxt(path_test + 'responses_nk_no_jackknife_%s' % end, responses, header='theta(arcmin), R_nk, Rgamma_nk, Rs_nk')

    def load_responses_nk(self, path_test, end):
        theta, R_nk, Rgamma_nk, Rs_nk = np.loadtxt(path_test + 'responses_nk_%s' % end, unpack=True)
        return theta, R_nk, Rgamma_nk, Rs_nk

    def load_responses_nk_errors(self, path_test, end):
        theta, R_nk, err = np.loadtxt(path_test + 'mean_R_nk_JK_%s' % end, unpack=True)
        return theta, R_nk, err

    def save_responses_mean(self, responses_mean, end):
        responses_mean = np.array([responses_mean[sbin] for sbin in self.zbins['sbins']])
        np.savetxt(self.get_path_test_allzbins() + 'responses_mean_%s'%end, responses_mean,
                   header='R_mean, Rgamma_mean, Rs_mean')

    def load_responses_mean(self, end):
        R_mean, Rgamma, Rs = np.loadtxt(self.get_path_test_allzbins() + 'responses_mean_%s'%end, unpack=True)
        return R_mean 

    def compute_Rs(self, e_ix, delta_gamma):
        """
        Computes R_s.
        e_ix: Dictionary per each component 1p, 1m, 2p, 2m.
        delta_gamma: value of the artificially applied shear to the images. 
        It can be averaged over all angular scales, or averaged in angular bins using NK correlation.
        """
        Rs11_mean = (e_ix['1p'] - e_ix['1m'])/delta_gamma
        Rs22_mean = (e_ix['2p'] - e_ix['2m'])/delta_gamma
        Rs_mean = 0.5 * (Rs11_mean + Rs22_mean)
        return Rs_mean

    def compute_chi2_consistency(self, datavector, constant, cov):
        diff = datavector-constant
        chi2 = np.dot(diff.T, np.dot(np.linalg.inv(cov), diff))

        # Hartlap factor
        N = self.config['njk']  # number of jackknife regions
        p = len(diff)  # number of angular bins
        factor = (N - p - 2) / float(N - 1)
        chi2_hartlap = chi2 * factor
        return chi2_hartlap

    def fit_constant_minuit(self, datavector, cov):
        """
        Function to fit a constant to some data.
        Minimizes chi2 using Minuit
        """
        from iminuit import Minuit
        data_mat = np.mat(datavector)
        def f(c):
            return np.abs(np.array((data_mat - c)*np.linalg.inv(cov)*(data_mat.T - c))[0][0])

        m = Minuit(f, print_level=0, errordef=1, pedantic = False)
        m.migrad()
        fit, err_fit = m.values['c'], m.errors['c']
        hartlap_factor = (self.config['njk'] - len(datavector) - 2) / float(self.config['njk'] - 1)
        chi2_fit = f(fit)*hartlap_factor
        
        print 'fit, err_fit = %0.3e +- %0.3e (Minuit)'%(fit, err_fit)
        print 'chi2_fit/ndf: %0.2f/%d (Minuit)'%(chi2_fit, (len(datavector)-1))
        return fit, err_fit, chi2_fit, len(datavector)-1


    def fit_constant_least_squares(self, datavector, cov, R0):
        """
        Function to fit a constant to some data.
        Minimizes chi2 using least squares from scipy
        R0: initial guess for the constant
        """
        from scipy import optimize
        data_mat = np.mat(datavector)
        def f(c):
            return np.abs(np.array((data_mat - c)*np.linalg.inv(cov)*(data_mat.T - c))[0][0])
        #ipdb.set_trace()
        opt = optimize.least_squares(f, R0)
        fit = opt.x
        print 'Fit with least squares:', fit

        hartlap_factor = (self.config['njk'] - len(datavector) - 2) / float(self.config['njk'] - 1)
        chi2_fit = f(fit[0])*hartlap_factor
        
        print 'fit = %0.3e (least squares)'%(fit[0])
        print 'chi2_fit/ndf: %0.2f/%d (least squares)'%(chi2_fit, (len(datavector)-1))
        return fit 


    def build_dictionary_e_ix(self, lens, source_sels, average_type):
        """
        Function to build the dictionary e_ix, which is the mean ellipticity for 
        the component i=(1,2) for a given selection x=(p,m). Therefore, it will have
        the mean of e_1p['1p'], e_1m['1m'], e_2p['2p'], e_2m['2m'], which will be obtained
        from the selections 1p, 1m, 2p, 2m respectively. This dictionary is used 
        to compute the selection response Rs.
        """
    
        e_ix = {}
        components = ['1p', '1m', '2p', '2m']
        for i, comp in enumerate(components):
            source_component_ix = {
                'ra': source_sels['ra'][i],
                'dec': source_sels['dec'][i],
                'e_ix': source_sels['e%s'%comp[0]][i]} # Choose e1 for 1p, 1m selections, and e2 for 2p, 2m selections. 
            if average_type == 'mean':
                e_ix[comp] = np.mean(source_component_ix['e_ix'])
            if average_type == 'NK_no_jackknife':
                theta, e_ix[comp] = self.run_nk_no_jackknife(lens, source_component_ix, scalar = source_component_ix['e_ix'])

            if average_type == 'NK_jackknife':
                theta, xi_nk, weights, npairs = self.run_treecorr_jackknife(lens, source_component_ix, type_corr = 'NK_e_ix')
                e_ixnum, _, wnum = self.numerators_jackknife(xi_nk, xi_nk, weights)
                e_ix[comp] = e_ixnum / wnum # contains all the jackknife regions measurements

        return e_ix

    def run_jackkniferesponses_tomo(self, lens, source, source_sels, delta_gamma, average_type, path_test, lens_or_random):
        """
        Function that computes mean response for each source bin or scale dependent responses for each lens-source combination.
        Uses NK TreeCorr correlation for scale dependence.
        Source: Source catalog for a certain redshift bin for the unsheared selection. 
        Source_sels: List of source catalogs for each of the four sheared selections (sheared 1p, 1m, 2p, 2m)
                     but with the unsheared quantities ra, dec, e1, e2, to obtain the new selection response.  
        delta_gamma: value of the artificially applied shear to the images. 
        average_type: way to average the ellipticities, can be:
                     - 'mean': np.mean
                     - 'NK_no_jackknife': using NK correlations without jackknife
                     - 'NK_jackknife': using NK correlations with jackknife
        lens_or_random: string for saving.
        Rgamma: (R11 + R22)/2 for each galaxy.
        """

        if average_type == 'mean':
            Rgamma = np.mean(source['Rgamma'])

        if average_type == 'NK_no_jackknife':
            theta, Rgamma = self.run_nk_no_jackknife(lens, source, scalar = source['Rgamma'])

        if average_type == 'NK_jackknife':
            theta, xi_nk, weights, npairs = self.run_treecorr_jackknife(lens, source, type_corr = 'NK_Rgamma')
            Rgammanum, _, wnum = self.numerators_jackknife(xi_nk, xi_nk, weights)
            Rgamma = Rgammanum / wnum # contains all the jackknife regions measurements (i.e. Rgamma_all)
            
        e_ix = self.build_dictionary_e_ix(lens, source_sels, average_type)
        Rs = self.compute_Rs(e_ix, delta_gamma)

        R = Rgamma + Rs

        print 'average_type, R, Rgamma, Rs', average_type, R, Rgamma, Rs

        if average_type == 'mean':
            responses = [R, Rgamma, Rs]
            return responses
        
        if average_type == 'NK_no_jackknife':
            responses = zip(theta, R, Rgamma, Rs)
            self.save_responses_nk(path_test, responses, lens_or_random)

        if average_type == 'NK_jackknife':
            self.process_run(R, theta, path_test, 'R_nk_JK_%s'%lens_or_random)
            self.process_run(Rgamma, theta, path_test, 'Rgamma_nk_JK_%s'%lens_or_random)
            self.process_run(Rs, theta, path_test, 'Rs_nk_JK_%s'%lens_or_random)

    def run(self):
        """
        Runs the NK responses between lenses and sources.
        Runs for lenses and randoms too.
        """
        lens_all, random_all, source_all, source_all_5sels, calibrator = self.load_data_or_sims()
        resp = {}
        responses_mean = {}
        for sbin in self.zbins['sbins']:

            print 'Running responses test for source %s.' % sbin
            source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, zlim_low=self.zbins[sbin][0], zlim_high=self.zbins[sbin][1])
            resp[sbin] = source['Rmean']
            print 'R = source[Rmean]', resp
            source_sels = self.load_metacal_bin_sels_responses(source_all_5sels, zlim_low=self.zbins[sbin][0], zlim_high=self.zbins[sbin][1])
            delta_gamma = 2*0.01

            for lbin in self.zbins['lbins']:
                print 'Running responses test for lens %s.' % lbin
                path_test = self.get_path_test(lbin, sbin)
                make_directory(path_test)

                lens = lens_all[(lens_all['z'] > self.zbins[lbin][0]) & (lens_all['z'] < self.zbins[lbin][1])]
                responses_mean[sbin] = self.run_responses_tomo(lens, source, source_sels, delta_gamma, average_type='mean', path_test=path_test, lens_or_random='lens')
                #comment:responses_nk_no_jackknife = self.run_responses_tomo(lens, source, source_sels, delta_gamma, average_type='NK_no_jackknife', path_test=path_test, lens_or_random='lens') #much slower
                self.run_responses_tomo(lens, source, source_sels, delta_gamma, average_type='NK_jackknife',path_test=path_test, lens_or_random='lens')

                #comment:random = random_all[(random_all['z'] > self.zbins[lbin][0]) & (random_all['z'] < self.zbins[lbin][1])]
                #comment:responses_nk_no_jackknife = self.run_responses_tomo(random, source, source_sels, delta_gamma, average_type='NK_no_jackknife', path_test=path_test, lens_or_random='random') #super slow
                #comment:self.run_responses_tomo(random, source, source_sels, delta_gamma, average_type='NK_jackknife', path_test=path_test, lens_or_random='random')


        print resp
        print responses_mean
        self.save_responses_mean(resp, 'destest')
        self.save_responses_mean(responses_mean, 'xcorr')
            
    def plot(self, lens_random, mask_scales):
        """
        Makes plot comparing the NK responses to the mean ones.
        lens_random: string, can be lens or random.
        mask_scales: boolean, scales to use for the chi2
        Indicates which is the foreground sample when computing the NK correlations.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(0.)
        c2 = plt.get_cmap(cmap)(0.25)
        c3 = plt.get_cmap(cmap)(0.5)
        fig, ax = plt.subplots(len(self.zbins['sbins']), len(self.zbins['lbins']), figsize=(16.5, 13.2), sharey='row', sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)

        R_mean_all = self.load_responses_mean('xcorr')

        err_Rmean = np.array([0.0003528, 0.0004943, 0.0004696, 0.0005272]) #values from Marco from JK, see slack on 14th of March, propagate errors
        print R_mean_all
        for l in range(len(self.zbins['lbins'])):

            for s in range(len(self.zbins['sbins'])):
                R_mean = R_mean_all[s]
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                #theta, R_nk, _, _ = self.load_responses_nk(path_test, lens_random)
                theta, R_nk_jk, err = self.load_responses_nk_errors(path_test, lens_random)
                cov = np.loadtxt(path_test + 'cov_R_nk_JK_%s'%lens_random)

                if mask_scales:
                    ax[s][l].axvspan(self.config['thlims'][0]*0.8, self.plotting['th_limit'][l], color='gray', alpha=0.2)                
                    mask = theta>self.plotting['th_limit'][l]
                    chi2 = self.compute_chi2_consistency(datavector=R_nk_jk[mask], constant=R_mean, cov=(cov[mask].T)[mask].T)
                    save = 'mask_scales'
                    ndf = len(R_nk_jk[mask])
                    c, err_c, chi2_c, ndf_c = self.fit_constant_minuit(datavector=R_nk_jk[mask], cov=(cov[mask].T)[mask].T)
                    
                else:
                    chi2 = self.compute_chi2_consistency(datavector=R_nk_jk, constant=R_mean, cov=cov)
                    save = 'all_scales'
                    ndf = len(R_nk_jk)
                    c, err_c, chi2_c, ndf_c = self.fit_constant_minuit(datavector=R_nk_jk, cov=cov)
                    c_ls = self.fit_constant_least_squares(datavector=R_nk_jk, cov=cov, R0=R_mean)

                ax[s][l].margins(x=0, y=0.3)
                t = np.linspace(0,max(theta)*2, 10)
                ax[s][l].plot(t, [R_mean] * len(t), '-', lw=2, color=c1, mec=c1, label=r'$R_{\mathrm{mean}}$')
                ax[s][l].fill_between(t, np.array([R_mean - err_Rmean[s] for i in t]), 
                              np.array([R_mean + err_Rmean[s] for i in t]), alpha=0.4, edgecolor=c1, facecolor=c1)
                #ax[s][l].plot(theta, R_nk, '-', lw=2, color=c2, mec=c2, label=r'$R_{\mathrm{nk}}$')
                #ax[s][l].plot(theta, [np.mean(R_nk)] * len(theta), '--', lw=2, color=c2, mec=c2,
                #              label=r'$\overline{R_{\mathrm{nk}}}$')
                ax[s][l].errorbar(theta, R_nk_jk, err, fmt = 'o', color=c3, mec=c3, markersize=3., label=r'$R_{\mathrm{nk, jk}}$')
                #ax[s][l].plot(theta, [np.mean(R_nk_jk)] * len(theta), ':', lw=2, color=c3, mec=c3,
                #              label=r'$\overline{R_{\mathrm{nk, jk}}}$')


                #ax[s][l].plot(theta, [c_ls] * len(theta), '--', lw=2, color=c3, mec=c3,
                #              label=r'Fit with least squares')
                ax[s][l].fill_between(t, np.array([c - err_c for i in t]), 
                              np.array([c + err_c for i in t]), alpha=0.4, edgecolor=c3, facecolor=c3, 
                              label='Fit to a constant')

                ax[s][l].set_xscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'Responses', size='larger', linespacing=3)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], size='larger')

                diff = R_mean/R_nk_jk - 1

                #ax[s][l].text(0.5, 0.88,
                #              r'Mean $R_{\mathrm{mean}}/R_{\mathrm{nk}}-1 = %0.2f \%%$' % (100 * np.mean(diff)),
                #              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                #              fontsize='medium')

                #ax[s][l].text(0.5, 0.79,
                #              r'Max $R_{\mathrm{mean}}/R_{\mathrm{nk}}-1 = %0.2f \%%$' % (100 * np.max(np.absolute(diff))),
                #              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                #              fontsize='medium')

                ax[s][l].text(0.5, 0.85,
                              r'$\chi^2$/ndf$ = %0.1f/%d$' % (chi2, ndf),
                              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                              fontsize='medium')

        ax[0][4].legend(frameon=False, fontsize=10, loc='lower right')
        self.save_plot('plot_responses_scale_dependence_%s_%s'%(save, lens_random))


    def plot_sigmas(self, lens_random, mask_scales):
        """
        Makes plot comparing the NK responses to the mean ones, divided by the uncertainty on the measurement. 
        lens_random: string, can be lens or random.
        Indicates which is the foreground sample when computing the NK correlations.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(0.)
        c2 = plt.get_cmap(cmap)(0.6)
        fig, ax = plt.subplots(len(self.zbins['sbins']), len(self.zbins['lbins']), figsize=(16.5, 13.2), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)

        R_mean_all = self.load_responses_mean('xcorr')
        print R_mean_all

        measurement = Measurement(self.basic, self.config, self.paths, self.zbins, self.plotting)
        gammat_file = measurement.load_twopointfile()
        gammat = gammat_file.spectra[0]


        for l in range(len(self.zbins['lbins'])):

            for s in range(len(self.zbins['sbins'])):
                R_mean = R_mean_all[s]
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, R_nk, _ = self.load_responses_nk_errors(path_test, lens_random)
                path_test_measurement = measurement.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                cov = np.loadtxt(path_test_measurement + 'cov_gt')
                _, gt = gammat.get_pair(l+1, s+1)
                err = gammat.get_error(l+1, s+1)
                diff_err = (R_mean/R_nk-1)*gt/err
                diff = (R_mean/R_nk-1)*gt

                if mask_scales:
                    ax[s][l].axvspan(self.config['thlims'][0]*0.8, self.plotting['th_limit'][l], color='gray', alpha=0.2)                
                    mask = theta>self.plotting['th_limit'][l]
                    chi2 = self.compute_chi2_consistency(datavector=R_mean/R_nk[mask]*gt[mask], constant=gt[mask], cov=(cov[mask].T)[mask].T)
                    save = 'mask_scales'
                    ndf = len(R_nk[mask])
                else:
                    chi2 = self.compute_chi2_consistency(datavector=R_mean/R_nk*gt, constant=gt, cov=cov)
                    save = 'all_scales'
                    ndf = len(R_nk)

                print 'l, s, chi2:', l, s, chi2
                ax[s][l].plot(theta, diff_err, lw=2, color=c1, mec=c1)
                ax[s][l].axhline(y=0, color= 'k', ls=':')
                #ax[s][l].plot(theta, err, lw=2, color=c2, mec=c2, label=r'$\sigma_{\gamma_t, \mathrm{JK}}$')
                #ax[s][l].plot(theta, diff/err, lw=2, color=c2, mec=c2, label=r'$(R_{\mathrm{nk}} - R_{\mathrm{mean}})/\sigma_{\gamma_t}$')
                ax[s][l].set_xscale('log')
                #ax[s][l].set_yscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'$(R_{\mathrm{mean}}/R_{\mathrm{nk}}-1)\gamma_t/\sigma_{\gamma_t, \mathrm{JK}}$', size='larger', linespacing=3)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], size='larger')

                ax[s][l].text(0.5, 0.87,
                              r'$\chi^2$/ndf$ = %0.3f/%d$' % (chi2, ndf),
                              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                              fontsize='medium')

        ax[0][4].legend(frameon=False, fontsize=16, loc='lower right')
        self.save_plot('plot_responses_%s_diff' % lens_random)



class ResponsesProjection(GGL):
    """
    Subclass that obtains the tangential compomenent projection of the responses (NG correlations) for all the lens-source bin combinations. Both for randoms and lenses.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test_allzbins(self):
        return os.path.join(self.paths['runs_config'], 'responses_ng') + '/'

    def get_path_test(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'responses_ng', lbin + '_' + sbin) + '/'

    def run(self):
        """
        Runs the projection of the responeses (in the same way as in the shears)
        """

        lens_all, random_all, source_all, source_all_5sels, calibrator = self.load_data_or_sims()

        for sbin in self.zbins['sbins']:
    		print 'Running measurement for source %s.' % sbin

		if self.basic['mode'] == 'data':
		    source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, zlim_low=self.zbins[sbin][0], zlim_high=self.zbins[sbin][1])
		    R = source['Rmean']
                    print 'R = source[Rmean]', R, sbin

		if self.basic['mode'] == 'data_y1sources':
		    source = pf.getdata(self.paths['y1'] + 'metacal_sel_sa%s.fits'%sbin[1])

    		for l, lbin in enumerate(self.zbins['lbins']):
    		    print 'Running measurement for lens %s.' % lbin
    		    path_test = self.get_path_test(lbin, sbin)
    		    make_directory(path_test)

    		    lens = lens_all[(lens_all['z'] > self.zbins[lbin][0]) & (lens_all['z'] < self.zbins[lbin][1])]

    		    theta, R, Rx, errs, weights, npairs = self.run_treecorr_jackknife(lens, source, 'NG')
    		    #self.save_runs(path_test, theta, R, Rx, errs, weights, npairs, False)
    		    Rnum, Rxnum, wnum = self.numerators_jackknife(R, Rx, weights)

                    random = random_all[(random_all['z'] > self.zbins[lbin][0]) & (random_all['z'] < self.zbins[lbin][1])]

    		    theta, R, Rx, errs, weights, npairs = self.run_treecorr_jackknife(random, source, 'NG')
    		    #self.save_runs(path_test, theta, R, Rx, errs, weights, npairs, True)
    		    Rnum_r, Rxnum_r, wnum_r = self.numerators_jackknife(R, Rx, weights)

    		    R_all = Rnum/wnum
                    R_r_all = Rnum_r/wnum_r
    		    Rx_all = Rxnum/wnum
                    Rx_r_all = Rxnum_r/wnum_r

    		    self.process_run(R_all, theta, path_test, 'R')
    		    self.process_run(Rx_all, theta, path_test, 'Rx')
    		    self.process_run(R_r_all, theta, path_test, 'R_randoms')
    		    self.process_run(Rx_r_all, theta, path_test, 'Rx_randoms')

        

class TestStars(GGL):
    """
    SubClass to test if the tangential shear around stars is consistent with zero.
    Uses no tomography for the source sample.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self, typestars):
        return os.path.join(self.paths['runs_config'], 'stars_%s' % typestars) + '/'

    def run(self, typestars):

        """
        Runs the gglensing measurement with stars as lenses and the full source sample.
        typestars: string that indicates which stars to load. Either 'bright' or 'faint'.
        Includes:
        - Mean response calculation for the full sample.
        - Random points subtraction.
        - Jackknife covariance calculation.
        """

        stars = pf.getdata(self.paths['y1'] + 'star_%s.fits' % typestars)
        randoms = pf.getdata(self.paths['y1'] + 'random.fits')
        mask_randoms = np.random.randint(len(randoms), size=len(stars) * 10)
        sources = pf.getdata(self.paths['y1'] + 'metacal_sel_allbins.fits')

        path_test = self.get_path_test(typestars)
        make_directory(path_test)

        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(stars, sources, 'NG')
        gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(randoms[mask_randoms], sources, 'NG')
        gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

        R = self.run_responses_mean_notomo(sources['Rgamma'])

        gt_all = (gtnum / wnum) / R - (gtnum_r / wnum_r) / R
        gx_all = (gxnum / wnum) / R - (gxnum_r / wnum_r) / R

        self.process_run(gt_all, theta, path_test, 'gt')
        self.process_run(gx_all, theta, path_test, 'gx')

    def plot(self):
        """
        Make plot for all the stars test.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        types = ['bright', 'faint']
        styles = ['o', '^']
        shifts = [0, 1, 2, 3]
        colors = [plt.get_cmap(self.plotting['cmap'])(0), plt.get_cmap(self.plotting['cmap'])(0.7)]
        titles_c = [self.plotting['catname']]
        s = 0  # If adding im3shape, s=1

        fig, ax = plt.subplots(1, 1, figsize=(5., 5.), sharey=False, sharex=False)
        ax.axhline(y=0, ls=':', color='k')

        for t, type in enumerate(types):

            path_test = self.get_path_test(type)
            th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)

            if s == 0:
                ax.errorbar(th * (1 + 0.05 * (shifts[t])), gt, err, fmt=styles[t], color=colors[s], mec=colors[s],
                            label=titles_c[s] + ' ' + type, markersize=4.5)
            if s == 1:
                ax.errorbar(th * (1 + 0.05 * (shifts[t + 2])), gt, err, fmt=styles[t], color=colors[s], mec=colors[s],
                            label=titles_c[s] + ' ' + type, markersize=4.5)

            chi2, ndf = self.get_chi2(path_test, 'gt')

        ax.set_xlim(2.5, 300)
        ax.set_ylim(-2 * 10 ** (-4), 2 * 10 ** (-4))
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))
        ax.set_xlabel(r'$\theta$ [arcmin]', size='large')
        ax.set_ylabel(r'$\gamma_{t} (\theta)$', size='large')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='best')
        self.save_plot('plot_stars')


class TestPSF(GGL):
    """
    SubClass to test if the psf residuals are compatible with zero.
    Uses no tomography for the lens sample.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self):
        return os.path.join(self.paths['runs_config'], 'psfresiduals') + '/'

    def save_psf_residuals(self, ra_lims, dec_lims):
        """
        Computes the psf residuals for the r band and saves them to a file.
        """

        info = pf.open(paths['y1base'] + 'cats/y1a1-v13/exposure_info_y1a1-v13.fits')[1].data
        exp = info['exp']
        print 'Original info:', len(exp)
        fil = info['filter']
        flag_info = info['flag']

        min_ra, max_ra = ra_lims
        min_dec, max_dec = dec_lims

        # Use only r-band
        mask_info = (fil == 'r')
        exp = exp[mask_info]
        print 'Only r band:', len(exp)

        # To not have repeated exposure names due to different ccds (in info file every ccd has a line)
        exp_unique = np.unique(exp)
        print 'All ccds:', len(exp_unique)
        print exp_unique

        ra_all, dec_all, psf1_all, psf2_all, res1_all, res2_all, mag_all = [], [], [], [], [], [], []

        for i in range(len(exp_unique)):
            # Load every exposure
            data = pf.open(paths['y1base'] + 'cats/y1a1-v13/psf_cats/%s_exppsf.fits' % exp_unique[i])[1].data
            ra = data['ra']
            dec = data['dec']
            e1 = data['e1']
            e2 = data['e2']
            psf1 = data['psf_e1']
            psf2 = data['psf_e2']
            # Resiual psf is the difference between the measurement of the psf (e1) and the model of the psf(psfex)
            # at the position of the stars, that is the only place you can measure the psf
            res1 = e1 - psf1
            res2 = e2 - psf2
            flag = data['flag']
            mag = data['mag']
            # Use only reserved stars and spt
            mask = ((dec > min_dec) & (dec < max_dec) & (flag == 64))
            ra_all.extend(ra[mask])
            dec_all.extend(dec[mask])
            psf1_all.extend(psf1[mask])
            psf2_all.extend(psf2[mask])
            res1_all.extend(res1[mask])
            res2_all.extend(res2[mask])
            mag_all.extend(mag[mask])
        # Convert to numpy array and simplify name
        # plt.hist(mag_all, bins=25)
        ra = np.array(ra_all)
        dec = np.array(dec_all)
        psf1 = np.array(psf1_all)
        psf2 = np.array(psf2_all)
        res1 = np.array(res1_all)
        res2 = np.array(res2_all)
        # w = np.ones(len(ra))
        print 'Number of exposures:', len(ra)

        c1 = pf.Column(name='RA', format='E', array=ra)
        c2 = pf.Column(name='DEC', format='E', array=dec)
        c3 = pf.Column(name='E1', format='E', array=res1)
        c4 = pf.Column(name='E2', format='E', array=res2)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra))
        hdu.writeto(self.paths['y1'] + 'psfresiduals.fits', clobber=True)

    def save_psf_residuals_y3(self, ra_lims, dec_lims):
        """
        Computes the psf residuals for the r band and saves them to a file.
        """

        # info = pf.open(paths['y1base'] + 'cats/y1a1-v13/exposure_info_y1a1-v13.fits')[1].data
        # exp = info['exp']
        # print 'Original info:', len(exp)
        # fil = info['filter']
        # flag_info = info['flag']
        #
        min_ra, max_ra = ra_lims
        min_dec, max_dec = dec_lims
        #
        # # Use only r-band
        # mask_info = (fil == 'r')
        # exp = exp[mask_info]
        # print 'Only r band:', len(exp)
        #
        # # To not have repeated exposure names due to different ccds (in info file every ccd has a line)
        # exp_unique = np.unique(exp)
        # print 'All ccds:', len(exp_unique)
        # print exp_unique

        ra_all_bandr, dec_all_bandr, psf1_all_bandr, psf2_all_bandr, res1_all_bandr, res2_all_bandr, mag_all_bandr, exp_all_bandr, band_all_bandr = [], [], [], [], [], [], [], [], []

        ra_all_bandi, dec_all_bandi, psf1_all_bandi, psf2_all_bandi, res1_all_bandi, res2_all_bandi, mag_all_bandi, exp_all_bandi, band_all_bandi = [], [], [], [], [], [], [], [], []

        ra_all_bandz, dec_all_bandz, psf1_all_bandz, psf2_all_bandz, res1_all_bandz, res2_all_bandz, mag_all_bandz, exp_all_bandz, band_all_bandz = [], [], [], [], [], [], [], [], []

        exposures = []

        y3_exp_dir = self.paths['y3_exp']

        for root, directories, filenames in os.walk(y3_exp_dir):
            for directory in directories:
                exposures.append(directory)

        exposures_int = np.array([int(exp) for exp in exposures])
        exp_unique = np.unique(exposures_int)

        print len(exp_unique)

        for i in range(len(exp_unique)):

            if np.mod(i,100) == 0:
                print i

            # Load every exposure
            data = pf.open(y3_exp_dir + str(exp_unique[i]) + '/exp_psf_cat_%s.fits' % exp_unique[i])[1].data
            info = pf.open(y3_exp_dir + str(exp_unique[i]) + '/exp_psf_cat_%s.fits' % exp_unique[i])[2].data

            flag = np.array(info['flag'])

            if np.sum(flag) == 0:
                band = info['band'][0][0]
                ra = data['ra']
                dec = data['dec']
                obs_e1 = data['obs_e1']
                obs_e2 = data['obs_e2']
                piff_e1 = data['piff_e1']
                piff_e2 = data['piff_e2']
                # Resiual psf is the difference between the measurement of the psf (e1) and the model of the psf(psfex)
                # at the position of the stars, that is the only place you can measure the psf
                res1 = obs_e1 - piff_e1
                res2 = obs_e2 - piff_e2
                obs_flag = data['obs_flag']
                mag = data['mag']
                # Use only reserved stars and spt
                mask = ((dec > min_dec) & (dec < max_dec) & (obs_flag == 64))

                if band == 'r':

                    ra_all_bandr.extend(ra[mask])
                    dec_all_bandr.extend(dec[mask])
                    psf1_all_bandr.extend(piff_e1[mask])
                    psf2_all_bandr.extend(piff_e2[mask])
                    res1_all_bandr.extend(res1[mask])
                    res2_all_bandr.extend(res2[mask])
                    mag_all_bandr.extend(mag[mask])

                elif band == 'i':
                    ra_all_bandi.extend(ra[mask])
                    dec_all_bandi.extend(dec[mask])
                    psf1_all_bandi.extend(piff_e1[mask])
                    psf2_all_bandi.extend(piff_e2[mask])
                    res1_all_bandi.extend(res1[mask])
                    res2_all_bandi.extend(res2[mask])
                    mag_all_bandi.extend(mag[mask])

                elif band == 'z':
                    ra_all_bandz.extend(ra[mask])
                    dec_all_bandz.extend(dec[mask])
                    psf1_all_bandz.extend(piff_e1[mask])
                    psf2_all_bandz.extend(piff_e2[mask])
                    res1_all_bandz.extend(res1[mask])
                    res2_all_bandz.extend(res2[mask])
                    mag_all_bandz.extend(mag[mask])

                else:
                    print 'no correct band alloted'

        # Convert to numpy array and simplify name
        # plt.hist(mag_all, bins=25)
        ra_bandr = np.array(ra_all_bandr)
        dec_bandr = np.array(dec_all_bandr)
        psf1_bandr = np.array(psf1_all_bandr)
        psf2_bandr = np.array(psf2_all_bandr)
        res1_bandr = np.array(res1_all_bandr)
        res2_bandr = np.array(res2_all_bandr)

        print 'Number of exposures:', len(ra_bandr)

        c1 = pf.Column(name='RA', format='E', array=ra_bandr)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandr)
        c3 = pf.Column(name='E1', format='E', array=res1_bandr)
        c4 = pf.Column(name='E2', format='E', array=res2_bandr)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra_bandr))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandr.fits', clobber=True)

        ra_bandi = np.array(ra_all_bandi)
        dec_bandi = np.array(dec_all_bandi)
        psf1_bandi = np.array(psf1_all_bandi)
        psf2_bandi = np.array(psf2_all_bandi)
        res1_bandi = np.array(res1_all_bandi)
        res2_bandi = np.array(res2_all_bandi)

        print 'Number of exposures bandi:', len(ra_bandi)

        c1 = pf.Column(name='RA', format='E', array=ra_bandi)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandi)
        c3 = pf.Column(name='E1', format='E', array=res1_bandi)
        c4 = pf.Column(name='E2', format='E', array=res2_bandi)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra_bandi))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandi.fits', clobber=True)

        ra_bandz = np.array(ra_all_bandz)
        dec_bandz = np.array(dec_all_bandz)
        psf1_bandz = np.array(psf1_all_bandz)
        psf2_bandz = np.array(psf2_all_bandz)
        res1_bandz = np.array(res1_all_bandz)
        res2_bandz = np.array(res2_all_bandz)

        # w = np.ones(len(ra))
        print 'Number of exposures bandz:', len(ra_bandz)

        c1 = pf.Column(name='RA', format='E', array=ra_bandz)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandz)
        c3 = pf.Column(name='E1', format='E', array=res1_bandz)
        c4 = pf.Column(name='E2', format='E', array=res2_bandz)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra_bandz))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandz.fits', clobber=True)

    def run(self):
        """
        Obtains the tangential component of the psf residuals around lenses, with random point subtraction.
        Obtains the corresponding jackknife covariance.
        """

        lens = pf.getdata(self.paths['y1'] + 'lens.fits')
        masklens = ((lens['z'] > self.zbins['l1'][0]) & (lens['z'] < self.zbins['l5'][1]))
        lens = lens[masklens]

        random = pf.getdata(self.paths['y1'] + 'random.fits')
        maskrandom = ((random['z'] > self.zbins['l1'][0]) & (random['z'] < self.zbins['l5'][1]))
        print self.zbins['l1'][0]
        print self.zbins['l5'][1]
        random = random[maskrandom]

        psfres = pf.getdata(self.paths['y1'] + 'psfresiduals.fits')

        path_test = self.get_path_test()
        make_directory(path_test)

        print 'PSF residuals around lenses...'
        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(lens, psfres, 'NG')
        gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

        print 'PSF residuals around randoms...'
        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(random, psfres, 'NG')
        gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

        gt_all = gtnum / wnum - gtnum_r / wnum_r

        self.process_run(gt_all, theta, path_test, 'gt')

    def run_y3(self,bands):
        """
        Obtains the tangential component of the psf residuals around lenses, with random point subtraction.
        Obtains the corresponding jackknife covariance.
        """

        lens = pf.getdata(self.paths['y3'] + 'lens.fits')
        masklens = ((lens['z'] > self.zbins['l1'][0]) & (lens['z'] < self.zbins['l5'][1]))
        lens = lens[masklens]

        random = pf.getdata(self.paths['y3'] + 'random.fits')
        maskrandom = ((random['z'] > self.zbins['l1'][0]) & (random['z'] < self.zbins['l5'][1]))
        print self.zbins['l1'][0]
        print self.zbins['l5'][1]
        random = random[maskrandom]

        for band in bands:

            psfres = pf.getdata(self.paths['y3'] + 'psfresiduals_band' + band + '.fits')

            path_test = self.get_path_test()
            make_directory(path_test)

            print 'PSF residuals around lenses...'
            theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(lens, psfres, 'NG')
            gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

            print 'PSF residuals around randoms...'
            theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(random, psfres, 'NG')
            gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

            gt_all = gtnum / wnum - gtnum_r / wnum_r

            self.process_run(gt_all, theta, path_test, 'gt')

    def plot(self):
        """
        Makes plot of the psf resdiuals.
        """
        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        c1 = plt.get_cmap(cmap)(0)
        titles_l = r'$0.15 < z < 0.90 $'
        title_redmagic = 'redMaGiC'
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        path_test = self.get_path_test()
        th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)
        ax.errorbar(th, gt, err, fmt='.', color=c1, mec=c1, markersize=5.7, capsize=1.4)
        ax.set_xlim(2.5, 250)
        ax.set_ylim(-1.8 * 10 ** (-5), 1.8 * 10 ** (-5))
        ax.set_xscale('log')
        ax.axhline(y=0, ls=':', color='k')
        ax.text(0.5, 0.85, titles_l, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.5, 0.92, title_redmagic, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

        ax.set_xlabel(r'$\theta$ [arcmin]', size='large')
        ax.set_ylabel(r'$\gamma_{t,\mathrm{PSF\, residuals}}$', size='large')

        chi2, ndf = self.get_chi2(path_test, 'gt')

        ax.text(0.7, 0.18, r'Null $\chi^2$/ndf ',
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
        ax.text(0.7, 0.11, r'$%0.1f/%d$' % (chi2, ndf),
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=12, color=c1)

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.save_plot('plot_psfresiduals')


class TestSizeSNR(GGL):
    """
    SubClass to test if the tangential shear signal has no dependence on source size or S/N. Using the first lens bin and all sources, as in Y1.
    The variable size_snr on run will be 'size' or 'snr' and will specify the test.
    """

    def __init__(self, basic, config, paths, zbins, plotting, source_nofz_pars):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting
        self.source_nofz_pars = source_nofz_pars

    def get_path_test(self, size_snr):
        return os.path.join(self.paths['runs_config'], '%s' % size_snr) + '/'

    def save_size_snr(self, path_test, size_snr, result_data, result_theory):
        np.savetxt(path_test + '%s_data' % size_snr, result_data, header='ratio err_ratio')
        np.savetxt(path_test + '%s_theory' % size_snr, result_theory, header='ratio err_ratio')

    def load_size_snr(self, path_test, size_snr):
        ratio_data, err_ratio_data = np.loadtxt(path_test + '%s_data' % size_snr, unpack=True)
        ratio_theory, err_ratio_theory = np.loadtxt(path_test + '%s_theory' % size_snr, unpack=True)
        return ratio_data, err_ratio_data, ratio_theory, err_ratio_theory

    def run_responses_mean_notomo_size_snr(self, Rgamma, size_snr, cut, high_low, delta_gamma):
        """
        Computes responses when there is an extra selection on size or snr.
        For all the sources (no tomography).
        - Rgamma: Rgamma for the high or low selection in size or snr.
        - size_snr: string, either size or snr.
        - cut: median of size or snr for the whole source sample.
        - high_low: string, either high or low.
        - delta_gamma: value of the artificially applied shear to the images.      
        """

        e_ix = {}  # ellipticities for component i for a given selection s, divided by Delta gamma.
        # x: p, m
        components = ['1p', '1m', '2p', '2m']
        for comp in components:
            print comp
            e_ix_allbins = np.zeros(0)
            par_ix_allbins = np.zeros(0)  # size or snr in the selection ix, for instance, size_1p, size_2p, etc.
            for sbin in self.zbins['sbins']:
                # Appending responses source bin sbin.
                source_selection = pf.getdata(
                    self.paths['y1'] + 'metacal_sel_responses/metacal_sel_responses_sa%s_%s.fits' % (sbin[1], comp))
                print 'dgamma is now dividing in the function that loads metacal and writes the file on it!'
                e_ix_allbins = np.append(e_ix_allbins, source_selection['Riisx'])
                par_ix_allbins = np.append(par_ix_allbins, source_selection['%s_ix' % size_snr])

            # For each component, mask the high or low part of e_ix
            if high_low == 'high':
                mask = np.where(par_ix_allbins > cut)
            if high_low == 'low':
                mask = np.where(par_ix_allbins <= cut)

            e_ix[comp] = np.mean(e_ix_allbins[mask])

        Rgamma_mean = np.mean(Rgamma)
        Rs_mean = self.compute_Rs(e_ix, delta_gamma)
        R_mean = self.save_responses_mean(Rgamma_mean, Rs_mean, 'notomo_%s_%s' % (size_snr, high_low))

        return R_mean

    def load_nzs(self, size_snr):
        """
	Loads redshift distributions for lenses, sources and source splits, low(l) and high(h).
	"""
        zl, nzl = np.loadtxt(self.paths['nz_lens'] + 'lens', unpack=True, usecols=(0, 1))
        zs, nzsl, nzsh = np.loadtxt(self.paths['nz_source_notomo_%s' % size_snr], unpack=True)
        nzsl = interpolate.interp1d(zs + self.source_nofz_pars['dzs', size_snr][0], nzsl, bounds_error=False,
                                    fill_value=0)(zs)
        nzsh = interpolate.interp1d(zs + self.source_nofz_pars['dzs', size_snr][1], nzsh, bounds_error=False,
                                    fill_value=0)(zs)

        return zl, nzl, zs, nzsl, nzsh

    def run(self, size_snr):
        lens_all = pf.getdata(self.paths['y1'] + 'lens.fits')
        lens = lens_all[(lens_all['z'] > self.zbins['l1'][0]) & (lens_all['z'] < self.zbins['l1'][1])]
        random_all = pf.getdata(self.paths['y1'] + 'random.fits')
        random = random_all[(random_all['z'] > self.zbins['l1'][0]) & (random_all['z'] < self.zbins['l1'][1])]

        sources = pf.getdata(self.paths['y1'] + 'metacal_sel_allbins.fits')

        # Source size or snr split, selecting the two halves of the split
        par = sources[size_snr]
        cut = np.median(par)
        print 'len(sources)', len(sources)
        maskl = par <= cut
        maskh = par > cut

        path_test = self.get_path_test(size_snr)
        make_directory(path_test)
        print size_snr
        print 'NEW len(sources[maskl])', len(sources[maskl])
        print 'NEW len(sources[maskh])', len(sources[maskh])

        # Computing the measurements for the split halves, both around lenses and randoms
        theta, gtsl, gxsl, errsl, weightsl, npairsl = self.run_treecorr_jackknife(lens, sources[maskl], 'NG')
        self.save_runs(path_test, theta, gtsl, gxsl, errsl, weightsl, npairsl, False)
        gtlnum, gxlnum, wlnum = self.numerators_jackknife(gtsl, gxsl, weightsl)

        theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r = self.run_treecorr_jackknife(random, sources[maskl],
                                                                                            'NG')
        self.save_runs(path_test, theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r, True)
        gtlnum_r, gxlnum_r, wlnum_r = self.numerators_jackknife(gtsl_r, gxsl_r, weightsl_r)

        theta, gtsh, gxsh, errsh, weightsh, npairsh = self.run_treecorr_jackknife(lens, sources[maskh], 'NG')
        gthnum, gxhnum, whnum = self.numerators_jackknife(gtsh, gxsh, weightsh)

        theta, gtsh_r, gxsh_r, errsh_r, weightsh_r, npairsh_r = self.run_treecorr_jackknife(random, sources[maskh],
                                                                                            'NG')
        gthnum_r, gxhnum_r, whnum_r = self.numerators_jackknife(gtsh_r, gxsh_r, weightsh_r)

        # Computing the responses for the split halves
        Rl = self.run_responses_mean_notomo_size_snr(sources['Rgamma'][maskl], size_snr, cut, 'low')
        Rh = self.run_responses_mean_notomo_size_snr(sources['Rgamma'][maskh], size_snr, cut, 'high')

        # Combining measurements and responses to get gammat
        gtl_all = (gtlnum / wlnum) / Rl - (gtlnum_r / wlnum_r) / Rl
        gth_all = (gthnum / whnum) / Rh - (gthnum_r / whnum_r) / Rh
        np.savetxt(path_test + 'gtl_all', gtl_all)
        np.savetxt(path_test + 'gth_all', gth_all)

        # Getting the data ratio using the simulations
        sims, cov_sims = self.load_sims()
        ratio, err_ratio = self.ratio_from_sims(theta, gtl_all, gth_all, sims, cov_sims)

        # Load N(z)'s and corrects the mean using Cosmos calibration.
        zl, nzl, zs, nzsl, nzsh = self.load_nzs(size_snr)

        # Computing inverse sigma_crit for the splits
        isch = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsh)
        iscl = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsl)

        # Computing the error on the theory prediction for the ratios, based on photo-z uncertainties
        # Shifting the N(z)'s for the splits up and down by the dzs_sigma parameter, and computing the theory ratio for each of these four cases
        shifts = [[-0.01, 0.], [0.01, 0.], [0., -0.01], [0., +0.01]]
        shifts = np.array(shifts) * self.source_nofz_pars['dzs_sigma']
        ratios = np.zeros(5)
        ratios[0] = isch / iscl
        sigmas = np.zeros(5)
        sigmas[0] = 0
        i = 1
        for s in shifts:
            nzsli = interpolate.interp1d(zs + s[0], nzsl, bounds_error=False, fill_value=0)
            nzsls = nzsli(zs)

            nzshi = interpolate.interp1d(zs + s[1], nzsh, bounds_error=False, fill_value=0)
            nzshs = nzshi(zs)

            ischs = functions.inv_sigma_crit_eff(zl, zs, nzl, nzshs)
            iscls = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsls)
            ratios[i] = ischs / iscls
            sigmas[i] = abs(ratios[0] - ratios[i])
            i = i + 1

        # Taking the mean of the sigma up and down for each of the splits, and computing the total error by adding them in quadrature
        sigma_tot = np.sqrt(np.mean([sigmas[1:3]]) ** 2 + np.mean(sigmas[3:5]) ** 2)
        print ratio, err_ratio
        print '%0.2f +- %0.2f' % (ratios[0], sigma_tot)

        # Saving the data ratio
        result_data = [ratio, err_ratio]
        result_theory = [ratios[0], sigma_tot]
        self.save_size_snr(path_test, size_snr, result_data, result_theory)

    def plot(self):
        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        lss = ['-', ':']
        labels_c = [self.plotting['catname']]
        tests = ['snr', 'size']
        labels_t = [r'S/N', r'Size']

        # Text
        titles_l = r'$(0.15 < z_l < 0.30)$'
        titles_s = r'$0.20 < z_s < 1.30$'

        fig, ax = plt.subplots(2, 2, figsize=(8, 6), sharey='row', sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.32, wspace=0.0)

        for t, size_snr in enumerate(tests):
            color = plt.get_cmap(cmap)(0)
            # Load N(z)'s and corrects the mean using Cosmos calibration.
            zl, nzl, zs, nzsl, nzsh = self.load_nzs(size_snr)
            ax[0][t].fill_between(zl, 0, nzl / np.max(nzl) * 2.04696244339, color='gray', alpha=0.2)
            ax[0][t].set_xlim(0, 1.8)
            ax[0][t].set_xlabel('Redshift')
            ax[0][t].plot(zs + self.source_nofz_pars['dzs', size_snr][0], nzsl, ls='-', lw=1.5, c=color,
                          label='Low %s' % (labels_c[0]))
            ax[0][t].plot(zs + self.source_nofz_pars['dzs', size_snr][1], nzsh, ls='-.', lw=1.5, c=color,
                          label='High %s' % (labels_c[0]))
            ax[0][t].set_title('%s split' % labels_t[t])

        ax[0][0].plot([], [], color='gray', alpha=0.2, linewidth=10, label='redMaGiC ' + titles_l)
        ax[0][0].set_ylabel('$n(z)$', fontsize=14)
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)
        c3 = 'k'

        handles, labels = ax[0][0].get_legend_handles_labels()

        # Sort the legend by the labels names
        # ax[0][1].legend(handles, labels, frameon=False,loc='best', prop={'size':9})
        # ax[0][0].set_xlim(0,1.79)
        ax[0][0].set_xlim(0, 1.79)
        ax[0][0].set_ylim(0, 3.)
        ax[0][1].legend(frameon=False, loc='best', prop={'size': 9})

        def compute_sigmas(q1, err1, q2, err2):
            # q1, err1: quantity and corresponding error
            print 'q1,err1,q2,err2:', q1, err1, q2, err2
            return np.abs(q1 - q2) / np.sqrt(err1 ** 2 + err2 ** 2)

        for t, size_snr in enumerate(tests):
            path_test = self.get_path_test(size_snr)
            ratio, err_ratio, invsc, err_invsc = self.load_size_snr(path_test, size_snr)
            sigma = compute_sigmas(ratio, err_ratio, invsc, err_invsc)

            color, x = c1, 1
            ax[1][t].errorbar(x, ratio, err_ratio, fmt='o', color=color, mec=color, label='%s' % (labels_c[0]))
            x_lin = [x - 0.2, x, x + 0.2]
            ax[1][t].fill_between(x_lin, np.array([invsc - err_invsc for i in x_lin]),
                                  np.array([invsc + err_invsc for i in x_lin]), alpha=0.4, edgecolor=color,
                                  facecolor=color,
                                  label='$\Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{high}}/ \Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{low}}$ %s' % (
                                      labels_c[0]))

            ax[1][t].set_xlim(0.5, 1.5)
            ax[1][t].set_ylim(0.4, 1.6)
            ax[1][t].set_xticklabels([])
            ax[1][t].set_xticks([])

        ax[1][0].legend(frameon=False, prop={'size': 10.5}, loc='best', numpoints=1)
        ax[1][0].set_ylabel(r'$\gamma^{\mathrm{high}}_t /\ \gamma^{\mathrm{low}}_t $', fontsize=14)

        self.save_plot('snr_size')


class TestSysMaps(GGL):
    """
    SubClass to test if the tangential shear signal has no dependence on observational conditions such as seeing, airmass etc, in each band griz.
    Using the first lens bin and all sources, as in Y1.
    The variable map can be: 'airmass', 'fwhm', 'maglimit', 'skybrite'. We iterate over them.
    The variable band can be: 'g', 'r', 'i', 'z'. We iterate over them. In Y1 we only used r band in the end, because it was used by im3shape.
    """

    def __init__(self, basic, config, paths, zbins, plotting, source_nofz_pars, sysmaps):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting
        self.source_nofz_pars = source_nofz_pars
        self.sysmaps = sysmaps

    def get_path_test(self, map, band):
        return os.path.join(self.paths['runs_config'], 'systematics_maps', map, band) + '/'

    def save_systematics_maps_ratios(self, path_test, result_data, result_theory):
        np.savetxt(path_test + 'data', result_data, header='ratio err_ratio')
        np.savetxt(path_test + 'theory', result_theory, header='ratio err_ratio')

    def load_systematics_maps_ratios(self, path_test):
        data, data_err = np.loadtxt(path_test + 'data', unpack=True)
        theory, theory_err = np.loadtxt(path_test + 'theory', unpack=True)
        return data, data_err, theory, theory_err

    def visualize_map(self, pix, signal, map, band, nside, nested_bool, name):
        """
        Plots the map and saves it.
        """
        ma = np.array([-1.6375 * 10 ** 30] * hp.nside2npix(nside))
        naturals = np.arange(0, len(ma))
        pix = np.in1d(naturals, pix)
        ma[pix] = signal
        # hp.gnomview(ma, min= mini, max=maxi, rot =  (73, -52, 0), xsize = 750 , title = '%s i'%map ,  notext = True, cmap = self.plotting['cmap'], nest = nested_bool)
        hp.mollview(ma, title='%s %s' % (map, band), notext=True, cmap=self.plotting['cmap'], nest=nested_bool)
        path = self.paths['plots'] + 'systematics_maps/'
        make_directory(path)
        plt.savefig(path + '%s_%s%s.pdf' % (map, band, name))

    def load_systematics_map(self, map, band, nside, nested_bool):
        '''
        Loads the systematics map, splits into high and low parts, and plots all maps.
        nested_bool: True if nest, False if ring.
        Returns: pixels corresponding the low half and the high half, for each map.
        '''
        path = os.path.join(self.paths['y1base'], 'cats', 'systematics_maps') + '/'
        sys_map = pf.getdata(
            path + 'Y1A1NEW_COADD_SPT_band_%s/' % band + 'Y1A1NEW_COADD_SPT_band_%s_nside%s_oversamp4_' % (
                band, nside) + self.sysmaps[map] + '.fits')
        pix = sys_map['PIXEL']
        sig = sys_map['SIGNAL']
        self.visualize_map(pix, sig, map, band, nside, nested_bool, '')

        mask_low = sig <= np.median(sig)
        mask_hi = sig > np.median(sig)

        pix_low = pix[mask_low]
        sig_low = sig[mask_low]

        pix_hi = pix[mask_hi]
        sig_hi = sig[mask_hi]

        self.visualize_map(pix_low, sig_low, map, band, nside, nested_bool, '_low')
        self.visualize_map(pix_hi, sig_hi, map, band, nside, nested_bool, '_high')

        return pix_low, pix_hi

    def load_nzs(self, map, band):
        """
	Loads redshift distributions for lenses, sources and source splits, low(l) and high(h).
	"""
        zl, nzl = np.loadtxt(self.paths['nz_lens'] + 'lens', unpack=True, usecols=(0, 1))
        zs, nzsl, nzsh = np.loadtxt(
            self.paths['y1base'] + 'runs/test_mcal_bpzmof_unblind/nofzs/source_%s_%s_notomo' % (map, band), unpack=True)
        return zl, nzl, zs, nzsl, nzsh

    def radec_to_thetaphi(self, ra, dec):
        """
        Converts ra and dec in degrees to theta and phi.
        Returns theta and phi in radians.
        """
        theta = (90. - dec) * np.pi / 180.
        phi = ra * np.pi / 180.
        return theta, phi

    def run(self, maps, bands):
        """
        Runs gglensing measurment for all maps and bands.
        """

        lens_all = pf.getdata(self.paths['y1'] + 'lens.fits')
        lens = lens_all[(lens_all['z'] > self.zbins['l1'][0]) & (lens_all['z'] < self.zbins['l1'][1])]
        random_all = pf.getdata(self.paths['y1'] + 'random.fits')
        random = random_all[(random_all['z'] > self.zbins['l1'][0]) & (random_all['z'] < self.zbins['l1'][1])]
        sources = pf.getdata(self.paths['y1'] + 'metacal_sel_allbins.fits')
        R = self.run_responses_mean_notomo(sources['Rgamma'])

        for map in maps:
            print 'Running map %s...' % map
            for band in bands:
                print 'Band %s' % band

                path_test = self.get_path_test(map, band)
                make_directory(path_test)

                # Load and split the systematics map
                pix_low, pix_hi = self.load_systematics_map(map, band, self.sysmaps['nside'],
                                                            self.sysmaps['nested_bool'])

                # Building lenses masks, low and high
                theta_l, phi_l = self.radec_to_thetaphi(lens['ra'], lens['dec'])
                pix_all_l = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_l, phi_l, self.sysmaps['nested_bool'])
                maskl_low = np.in1d(pix_all_l, pix_low)
                maskl_hi = np.in1d(pix_all_l, pix_hi)

                # Building randoms masks, low and high
                theta_r, phi_r = self.radec_to_thetaphi(random['ra'], random['dec'])
                pix_all_r = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_r, phi_r, self.sysmaps['nested_bool'])
                maskr_low = np.in1d(pix_all_r, pix_low)
                maskr_hi = np.in1d(pix_all_r, pix_hi)

                # Building sources masks, low and high
                theta_s, phi_s = self.radec_to_thetaphi(sources['ra'], sources['dec'])
                pix_all_s = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_s, phi_s, self.sysmaps['nested_bool'])
                masks_low = np.in1d(pix_all_s, pix_low)
                masks_hi = np.in1d(pix_all_s, pix_hi)

                # Computing the measurements for the split halves, both around lenses and randoms
                print 'Lenses, low.'
                theta, gtsl, gxsl, errsl, weightsl, npairsl = self.run_treecorr_jackknife(lens[maskl_low],
                                                                                          sources[masks_low], 'NG')
                gtlnum, gxlnum, wlnum = self.numerators_jackknife(gtsl, gxsl, weightsl)

                print 'Randoms, low.'
                theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r = self.run_treecorr_jackknife(random[maskr_low],
                                                                                                    sources[masks_low],
                                                                                                    'NG')
                gtlnum_r, gxlnum_r, wlnum_r = self.numerators_jackknife(gtsl_r, gxsl_r, weightsl_r)

                print 'Lenses, high.'
                theta, gtsh, gxsh, errsh, weightsh, npairsh = self.run_treecorr_jackknife(lens[maskl_hi],
                                                                                          sources[masks_hi], 'NG')
                gthnum, gxhnum, whnum = self.numerators_jackknife(gtsh, gxsh, weightsh)

                print 'Randoms, high.'
                theta, gtsh_r, gxsh_r, errsh_r, weightsh_r, npairsh_r = self.run_treecorr_jackknife(random[maskr_hi],
                                                                                                    sources[masks_hi],
                                                                                                    'NG')
                gthnum_r, gxhnum_r, whnum_r = self.numerators_jackknife(gtsh_r, gxsh_r, weightsh_r)

                # Combining measurements and responses to get gammat
                gtl_all = (gtlnum / wlnum) / R - (gtlnum_r / wlnum_r) / R
                gth_all = (gthnum / whnum) / R - (gthnum_r / whnum_r) / R

                # Getting the data ratio using the simulations
                sims, cov_sims = self.load_sims()
                ratio, err_ratio = self.ratio_from_sims(theta, gtl_all, gth_all, sims, cov_sims)

                # Load N(z)'s
                zl, nzl, zs, nzsl, nzsh = self.load_nzs(map, band)

                # Computing inverse sigma_crit for the splits
                isch = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsh)
                iscl = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsl)

                # Computing the error on the theory prediction for the ratios, based on photo-z uncertainties
                # Shifting the N(z)'s for the splits up and down by the dzs_sigma parameter, and computing the theory ratio for each of these four cases
                shifts = [[-0.01, 0.], [0.01, 0.], [0., -0.01], [0., +0.01]]
                shifts = np.array(shifts) * self.source_nofz_pars['dzs_sigma']
                ratios = np.zeros(5)
                ratios[0] = isch / iscl
                sigmas = np.zeros(5)
                sigmas[0] = 0
                i = 1
                for s in shifts:
                    nzsli = interpolate.interp1d(zs + s[0], nzsl, bounds_error=False, fill_value=0)
                    nzsls = nzsli(zs)

                    nzshi = interpolate.interp1d(zs + s[1], nzsh, bounds_error=False, fill_value=0)
                    nzshs = nzshi(zs)

                    ischs = functions.inv_sigma_crit_eff(zl, zs, nzl, nzshs)
                    iscls = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsls)
                    ratios[i] = ischs / iscls
                    sigmas[i] = abs(ratios[0] - ratios[i])
                    i = i + 1

                # Taking the mean of the sigma up and down for each of the splits, and computing the total error by adding them in quadrature
                sigma_tot = np.sqrt(np.mean([sigmas[1:3]]) ** 2 + np.mean(sigmas[3:5]) ** 2)
                print ratio, err_ratio
                print '%0.2f +- %0.2f' % (ratios[0], sigma_tot)

                # Saving the data ratio
                result_data = [ratio, err_ratio]
                result_theory = [ratios[0], sigma_tot]
                self.save_systematics_maps_ratios(path_test, result_data, result_theory)

    def plot(self):
        """
        Plots measurments compared to expected theory.
        Set up for four maps, and a single band.
        """
        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        #labels_cshort = r'\textsc{Metacal}'
        labels_cshort = 'Metacal'
        fontsize = 16
        maps = ['airmass', 'fwhm', 'maglimit', 'skybrite']
        band = 'r'
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)

        fig, ax = plt.subplots(2, 2, figsize=(6, 6), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for m, map in enumerate(maps):

            path_test = self.get_path_test(map, band)
            data, data_err, theory, theory_err = self.load_systematics_maps_ratios(path_test)
            if m % 2 == 0: c = 0
            if m % 2 == 1: c = 1
            if m == 0 or m == 1: f = 0
            if m == 2 or m == 3: f = 1

            color, x = c1, 1

            ax[f][c].errorbar(x, data, data_err, fmt='o', color=color, mec=color,
                              label='%s' % (self.plotting['catname']))
            x_lin = [x - 0.2, x, x + 0.2]
            ax[f][c].fill_between(x_lin, np.array([theory - theory_err for i in x_lin]),
                                  np.array([theory + theory_err for i in x_lin]), alpha=0.4, edgecolor=color,
                                  facecolor=color,
                                  label=r'$\Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{high}}/ \Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{low}} %s$' % (
                                      labels_cshort))

            ax[f][c].set_xlim(0.2, 2.8)
            ax[f][c].set_xticklabels([])
            ax[f][c].set_xticks([])
            ax[f][c].set_ylim(0.5, 1.5)
            ax[f][0].set_ylabel(r'$\gamma^{\mathrm{high}}_t /\ \gamma^{\mathrm{low}}_t $', fontsize=fontsize)
            ax[f][c].text(0.5, 0.9, r'\textsc{%s}' % map,
                          horizontalalignment='center', verticalalignment='center', transform=ax[f][c].transAxes,
                          fontsize=12, color='k')

        ax[1][0].legend(frameon=False, loc='best', numpoints=1, fontsize=9.8)
        self.save_plot('systematics_maps')

