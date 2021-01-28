import numpy as np
from astropy.table import Table, vstack, hstack
from astropy.stats import jackknife_stats
from matplotlib import pyplot as plt
from matplotlib import rc
import galsim
import pickle
import pdb

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=11)

plt.close('all') # tidy up any unshown plots

def make_response_plots(mastercat_version, mcal_bpz_dir, plot_dir='./plots/', pickle_dir='./', rseed=1234, n_realisations=int(1e5)):

  print('Making response plots...')

  cat = Table.read(mcal_bpz_dir + 'total_mcal_NearestNeighbors_class_match_Y3_mastercat_{0}.fits'.format(mastercat_version))

  bands = ['u', 'g', 'r', 'i', 'z', 'J', 'H', 'K']
  mag = {}
  mag_arr =np.zeros([len(cat), len(bands)])
  for ii in range(0, len(bands)):
    mag[bands[ii]] = cat['bdf_mag_dered'][:,ii]
    mag_arr[:,ii] = cat['bdf_mag_dered'][:,ii]

  cut_multiband = (mag['u'] > 0.) & (mag['u'] < 35.) & \
                  (mag['r'] > 0.) & (mag['r'] < 35.) & \
                  (mag['z'] > 0.) & (mag['z'] < 35.) & \
                  (mag['K'] > 0.) & (mag['K'] < 35.)

  ug = mag['u'] - mag['g']
  gr = mag['g'] - mag['r']
  ri = mag['r'] - mag['i']
  iz = mag['i'] - mag['z']
  zJ = mag['z'] - mag['J']
  JH = mag['J'] - mag['H']
  HK = mag['H'] - mag['K']

  ur = mag['u'] - mag['r']
  zK = mag['z'] - mag['K']

  binary_cut = (np.log10(cat['T']) < (22.25 - mag['r'])/2.5)
  highe_cut = np.sqrt(cat['e_1']**2. + cat['e_2']**2.) > 0.8
  binary_cut = binary_cut*highe_cut

  #shape_cut = make_shape_cuts(cat, snr=[10,1000])
  shape_cut = cat['select_bool']
  shape_cut = shape_cut*(~binary_cut)

  # kNN star/gals
  stars_kNN = cat['NearestNeighbors_class']==2
  gals_kNN = cat['NearestNeighbors_class']==1

  #Make histograms of responses for all classes of objects.
  star_R11 = cat['R11'][stars_kNN]
  gal_R11 = cat['R11'][gals_kNN]

  shape_star_R11 = cat['R11'][stars_kNN*shape_cut]
  shape_gal_R11 = cat['R11'][gals_kNN*shape_cut]

  star_R_66 = jackknife_stats(star_R11, np.mean, 0.66)[2]
  shape_star_R_66 = jackknife_stats(shape_star_R11, np.mean, 0.66)[2]

  gal_R_66 = 0.
  shape_gal_R_66 = 0.

  n_sub_r = 100

  for i_r in np.arange(n_sub_r):
    print('Jackknife of subsample {0}/{1}'.format(i_r, n_sub_r))

    gal_R11_sub = np.random.choice(gal_R11, star_R11.shape[0], replace=False)
    shape_gal_R11_sub = np.random.choice(shape_gal_R11, shape_star_R11.shape[0], replace=False)

    gal_R_66 += jackknife_stats(gal_R11_sub, np.mean, 0.66)[2] / n_sub_r
    shape_gal_R_66 += jackknife_stats(shape_gal_R11_sub, np.mean, 0.66)[2] / n_sub_r

  plt.figure(1, figsize=(2*4.5, 3.75))
  plt.subplot(121)
  plt.hist(star_R11, histtype='step', bins=np.linspace(-3,3,50), density=True, label='kNN Stars $\\langle R \\rangle  = {0:.4f} \pm {1:.4f}\,(66\%)$'.format(star_R11.mean(), star_R_66))
  plt.hist(gal_R11, histtype='step', bins=np.linspace(-3,3,50), density=True, label='kNN Galaxies $\\langle R \\rangle  = {0:.4f} \pm {1:.4f}\,(66\%)$'.format(gal_R11.mean(), gal_R_66))
  plt.title('All objects')
  plt.axvline(0, alpha=0.4, color='k', linestyle='dashed', zorder=-1)
  plt.legend(loc='upper left', fontsize='small')
  plt.xlabel('$R_{11}$')
  plt.subplot(122)
  star_R_counts, star_R_bins, _ = plt.hist(shape_star_R11, histtype='step', bins=np.linspace(-3,3,50), density=True, label='kNN Stars $\\langle R \\rangle  = {0:.4f} \pm {1:.4f}\,(66\%)$'.format(shape_star_R11.mean(), shape_star_R_66))
  galaxy_R_counts, galaxy_R_bins, _ = plt.hist(shape_gal_R11, histtype='step', bins=np.linspace(-3,3,50), density=True, label='kNN Galaxies $\\langle R \\rangle  = {0:.4f} \pm {1:.4f}\,(66\%)$'.format(shape_gal_R11.mean(), shape_gal_R_66))
  star_R_bins = (star_R_bins[:-1]+ star_R_bins[1:])/2
  galaxy_R_bins = (galaxy_R_bins[:-1]+ galaxy_R_bins[1:])/2
  plt.axvline(0, alpha=0.4, color='k', linestyle='dashed', zorder=-1)
  plt.title('Shear cat. objects. $N_{\\rm gal} = $'+'{0}'.format(np.sum(gals_kNN*shape_cut))+', $N_{\\rm star} = $'+'{0}'.format(np.sum(stars_kNN*shape_cut)))
  plt.legend(loc='upper left', fontsize='small')
  plt.xlabel('$R_{11}$')
  # this one
  plt.savefig(plot_dir + 'stellar-class-responses_{0}.png'.format(mastercat_version), dpi=300, bbox_inches='tight')
  plt.savefig(plot_dir + 'stellar-class-responses_{0}.pdf'.format(mastercat_version), bbox_inches='tight')

  star_mean_R = np.average(star_R11)
  gal_mean_R = np.average(gal_R11)
  shape_star_mean_R = np.average(shape_star_R11)
  shape_gal_mean_R = np.average(shape_gal_R11)

  N_star = np.sum(stars_kNN)
  N_shape_star = np.sum(stars_kNN*shape_cut)
  N_gal = np.sum(gals_kNN)
  N_shape_gal = np.sum(gals_kNN*shape_cut)

  fstar = N_shape_star / (N_shape_star + N_shape_gal)
  fgal = N_shape_gal / (N_shape_star + N_shape_gal)

  mean_m = (fstar/fgal)*(shape_star_mean_R/shape_gal_mean_R)

  # randomly sample an m
  # N have posson distribution
  # R have observed distribution

  star_R_rng = galsim.DistDeviate(rseed, function=galsim.LookupTable(star_R_bins,star_R_counts))
  galaxy_R_rng = galsim.DistDeviate(rseed, function=galsim.LookupTable(galaxy_R_bins,galaxy_R_counts))

  N_shape_star_rng = galsim.PoissonDeviate(rseed, mean=N_shape_star)

  star_R_arr = np.zeros(n_realisations)
  galaxy_R_arr = np.zeros(n_realisations)
  N_shape_star_arr = np.zeros(n_realisations)

  print('MC-ing errors...')
  for ir in np.arange(n_realisations):
    star_R_arr[ir] = star_R_rng()
    galaxy_R_arr[ir] = galaxy_R_rng()
    N_shape_star_arr[ir] = N_shape_star_rng()

  fstar_arr = N_shape_star_arr / (N_shape_star + N_shape_gal)
  fgal_arr = 1. - fstar_arr

  m = (fstar_arr/fgal_arr)*(star_R_arr/galaxy_R_arr)

  mask = (m>-0.05)*(m<0.05)

  plt.figure(2, figsize=(4.5, 3.75))
  plt.hist(m[mask], histtype='step', bins=np.linspace(-0.02, 0.04, 50), density=True)
  plt.axvline(np.median(m[mask]), linestyle='dashed', alpha=0.4)
  plt.axhline(0, color='k', linewidth=1, zorder=-10)
  plt.xlabel('Stellar contamination shear bias $m$')
  plt.savefig(plot_dir + 'stellar-contamination-biases_{0}.png'.format(mastercat_version), dpi=300, bbox_inches='tight')
  plt.savefig(plot_dir + 'stellar-contamination-biases_{0}.pdf'.format(mastercat_version), bbox_inches='tight')

  sav_obj = {
           'star_R11' : star_R11,
           'gal_R11' : gal_R11,
           'star_R_66' : star_R_66,
           'gal_R_66' : gal_R_66,
           'shape_star_R11' : shape_star_R11,
           'shape_gal_R11' : shape_gal_R11,
           'shape_star_R_66' : shape_star_R_66,
           'shape_gal_R_66' : shape_gal_R_66,
           'gals_kNN' : gals_kNN,
           'stars_kNN' : stars_kNN,
           'shape_cut' : shape_cut,
           'mastercat_version' : mastercat_version,
           'm' : m,
           'mask' : mask,
           }
    
  with open(pickle_dir + 'everything_you_need_for_stellar_contamination.pkl', 'wb') as f:
    pickle.dump(sav_obj, f, pickle.HIGHEST_PROTOCOL)

if __name__=='__main__':

  make_response_plots(mastercat_version='03_31_20',
                      mcal_bpz_dir='./data/')