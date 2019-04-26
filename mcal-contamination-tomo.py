from astropy.table import Table, vstack
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import os
import pdb

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=11)

plt.close('all') # tidy up any unshown plots

def make_shape_cuts(cat, snr=[10,100], flags=0, size_ratio=0.5):
  flag_cut = cat['flags']==flags
  snr_cut = np.less(cat['snr'], snr[1])*np.greater(cat['snr'], snr[0])
  size_cut = np.greater(cat['T']/cat['psf_T'], size_ratio)
  shape_cut = flag_cut*snr_cut*size_cut

  return shape_cut

extended_class_coadd_list = [0, 1, 2, 3]
cmap_list = ['Blues', 'Oranges', 'Greens', 'Reds']
extended_class_coadd_names = {'0': 'High confidence stars',
                              '1': 'Candidate stars',
                              '2': 'Mostly galaxies',
                              '3': 'High confidence galaxies'}

field_list = ['sxds',
              'vvds',
              'deep23']

m_per_sub = {'sxds' : [],
              'vvds' : [],
              'deep23' : []}

sigma_m_per_sub = {'sxds' : [],
              'vvds' : [],
              'deep23' : []}

N_per_sub = {'sxds' : [],
              'vvds' : [],
              'deep23' : []}

for field_name in field_list:

  plt.close('all')

  hsc_fname = 'data/{0}.hsc.fits'.format(field_name)
  gold_fname = 'data/{0}.y3_gold_2_2.fits'.format(field_name)
  mcal_fname = 'data/{0}.metacal_unsheared_Y3_mastercat_v2_6_20_18.fits'.format(field_name)

  hsc_gold_match_fname = 'data/{0}.hsc_gold_matches.fits'.format(field_name)
  mcal_gold_match_fname = 'data/{0}.mcal_gold_matches.fits'.format(field_name)
  mcal_hsc_gold_match_fname = 'data/{0}.hsc_mcal_gold_matches.fits'.format(field_name)

  #hsc = Table.read(hsc_fname)
  #gold = Table.read(gold_fname)
  #mcal = Table.read(mcal_fname)

  stilts_cmd = 'java -jar /Applications/TOPCAT.app/Contents/Resources/Java/topcat-full.jar -stilts'
  stilts_sky_join_cmd = ' tmatch2 in1={0} in2={1} out={2} matcher=sky values1="ra dec" values2="ra dec" params="0.5"'
  stilts_id_join_cmd = ' tmatch2 in1={0} in2={1} out={2} matcher=exact values1="coadd_object_id" values2="COADD_OBJECT_ID"'


  # cross match all the catalogues
  if not os.path.exists(hsc_gold_match_fname):
    os.system(stilts_cmd+stilts_sky_join_cmd.format(hsc_fname, gold_fname, hsc_gold_match_fname))

  hsc_gold_matches = Table.read(hsc_gold_match_fname)

  if not os.path.exists(mcal_gold_match_fname):
    os.system(stilts_cmd+stilts_id_join_cmd.format(mcal_fname, gold_fname, mcal_gold_match_fname))

  mcal_gold_matches = Table.read(mcal_gold_match_fname)

  if not os.path.exists(mcal_hsc_gold_match_fname):
    os.system(stilts_cmd+stilts_id_join_cmd.format(mcal_fname, hsc_gold_match_fname, mcal_hsc_gold_match_fname))

  mcal_hsc_gold_matches = Table.read(mcal_hsc_gold_match_fname)

  z_bins = np.array([0.2,0.43,0.63,0.90,1.30])
  hsc_gold_sample_bins = np.digitize(hsc_gold_matches['BPZ_ZMEAN_MOF'], z_bins)
  mcal_gold_sample_bins = np.digitize(mcal_gold_matches['BPZ_ZMEAN_MOF'], z_bins)
  mcal_hsc_gold_sample_bins = np.digitize(mcal_hsc_gold_matches['BPZ_ZMEAN_MOF'], z_bins)

  for ibin in np.arange(1,len(z_bins)):
    plt.close('all')
    hsc_gold_matches_bin = hsc_gold_matches[hsc_gold_sample_bins==ibin]
    mcal_gold_matches_bin = mcal_gold_matches[mcal_gold_sample_bins==ibin]
    mcal_hsc_gold_matches_bin = mcal_hsc_gold_matches[mcal_hsc_gold_sample_bins==ibin]

    # shape cuts to catalogues
    shape_mcal_gold_cut = make_shape_cuts(mcal_gold_matches_bin)
    shape_mcal_hsc_gold_cut = make_shape_cuts(mcal_hsc_gold_matches_bin)

    # y1 spreadmodel cuts to mcal-gold
    star_spread_cut = np.less(np.abs(mcal_gold_matches_bin['SPREAD_MODEL_I'] + (5./3.)*mcal_gold_matches_bin['SPREADERR_MODEL_I']), 0.002)
    gal_spread_cut = ~star_spread_cut

    # hsc star cuts to mcal-hsc-gold
    star_hsc_cut = mcal_hsc_gold_matches_bin['iclassification_extendedness']==0
    gal_hsc_cut = ~star_hsc_cut

    R11_gold_cut = np.less(mcal_gold_matches_bin['R11'], 75.)*np.greater(mcal_gold_matches_bin['R11'], -75)*(mcal_gold_matches_bin['flags']==0)
    R11_hsc_cut = np.less(mcal_hsc_gold_matches_bin['R11'], 75.)*np.greater(mcal_hsc_gold_matches_bin['R11'], -75)*(mcal_hsc_gold_matches_bin['flags']==0)

    print('############### Field: {0}, z bin: {1} ###############'.format(field_name, ibin))
    # need to think about completeness across this
    #N_mcal = len(mcal)
    #N_gold = len(gold)
    #N_hsc = len(hsc)

    N_mcal_gold = len(mcal_gold_matches_bin)
    N_hsc_gold = len(hsc_gold_matches_bin)
    N_mcal_hsc_gold = len(mcal_hsc_gold_matches_bin)

    N_shape = np.sum(shape_mcal_gold_cut)
    N_shape_hsc = np.sum(shape_mcal_hsc_gold_cut)
    N_hsc_shape_stars = np.sum(star_hsc_cut*shape_mcal_hsc_gold_cut)
    N_hsc_shape_gal = np.sum(gal_hsc_cut*shape_mcal_hsc_gold_cut)
    N_spread_shape_stars = np.sum(star_spread_cut*shape_mcal_gold_cut)

    N_high_shape_stars = np.sum(shape_mcal_gold_cut*(mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==0))
    N_cand_shape_stars = np.sum(shape_mcal_gold_cut*(mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==1))

    #print(N_mcal, N_gold, N_hsc)
    #print(N_mcal_gold, N_hsc_gold, N_mcal_hsc_gold)

    #print(N_shape, N_hsc_shape_stars, N_spread_shape_stars)
    #print(N_shape, N_hsc_shape_stars/N_shape, N_spread_shape_stars/N_shape)

    print('HSC contamination: {0:.2f}%'.format(100*N_hsc_shape_stars/N_shape))
    print('Y1Spreadmodel contamination: {0:.2f}%'.format(100*N_spread_shape_stars/N_shape))

    print('HSC Star <R>: {0}'.format(np.mean(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*star_hsc_cut])))
    print('HSC Star std(R): {0}'.format(np.std(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*star_hsc_cut])))

    print('HSC Galaxy <R>: {0}'.format(np.mean(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*gal_hsc_cut])))
    print('HSC Galaxy std(R): {0}'.format(np.std(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*gal_hsc_cut])))

    hsc_star_mean_R = np.mean(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*star_hsc_cut])
    hsc_galaxy_mean_R = np.mean(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*gal_hsc_cut])

    hsc_star_var_R = np.var(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*star_hsc_cut])
    hsc_galaxy_var_R = np.var(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*gal_hsc_cut])

    m_hsc = (N_hsc_shape_stars / N_hsc_shape_gal)*(hsc_star_mean_R / hsc_galaxy_mean_R)
    #sigma_m_hsc = np.sqrt((hsc_star_var_R/hsc_star_mean_R)**2. + (hsc_galaxy_var_R/hsc_galaxy_mean_R)**2.)*m_hsc
    sigma_m_hsc = np.sqrt((hsc_star_var_R/(hsc_star_mean_R*np.sqrt(N_hsc_shape_stars)))**2. + (hsc_galaxy_var_R/(hsc_galaxy_mean_R*np.sqrt(N_hsc_shape_gal)))**2.)*m_hsc
  
    m_per_sub[field_name].append(m_hsc)
    sigma_m_per_sub[field_name].append(sigma_m_hsc)

    print('m from HSC stars = {0:.5f} +- {1:.5f}'.format(m_hsc, sigma_m_hsc))

    plt.figure(1, figsize=(2*4.5, 3.75))
    plt.subplot(121)
    plt.title('Galaxies')
    plt.hist2d(mcal_hsc_gold_matches_bin['MAG_AUTO_I'][gal_hsc_cut], mcal_hsc_gold_matches_bin['SPREAD_MODEL_I'][gal_hsc_cut], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap='Blues')
    plt.xlabel('MAG\_AUTO\_I')
    plt.ylabel('SPREAD\_MODEL\_I')
    plt.subplot(122)
    plt.title('Stars')
    plt.hist2d(mcal_hsc_gold_matches_bin['MAG_AUTO_I'][star_hsc_cut], mcal_hsc_gold_matches_bin['SPREAD_MODEL_I'][star_hsc_cut], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap='Oranges')
    plt.xlabel('MAG\_AUTO\_I')
    #plt.ylabel('SPREAD MODEL I')
    plt.suptitle('Field {0}, '.format(field_name)+'Redshift bin {0}'.format(ibin)+' HSC \\textsc{iclassification\_extendedness}')
    plt.savefig('plots/{1}-stargal-hsc-bin{0}.png'.format(ibin, field_name), dpi=300, bbox_inches='tight')

    plt.figure(2, figsize=(4*4.5, 3.75))
    for iext, extended_class_coadd in enumerate(extended_class_coadd_list):
      plt.subplot(1,len(extended_class_coadd_list),iext+1)
      plt.title(extended_class_coadd_names[str(extended_class_coadd)])
      plt.hist2d(mcal_gold_matches_bin['MAG_AUTO_I'][mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==extended_class_coadd], mcal_gold_matches_bin['SPREAD_MODEL_I'][mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==extended_class_coadd], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap=cmap_list[iext])
      plt.xlabel('MAG\_AUTO\_I')
      if iext==0:
        plt.ylabel('SPREAD\_MODEL\_I')
    plt.suptitle('Field {0}, '.format(field_name)+'Redshift bin {0}'.format(ibin)+' Gold \\textsc{EXTENDED\_CLASS\_COADD}')
    plt.savefig('plots/{1}-stargal-y3extendedness-bin{0}.png'.format(ibin, field_name), dpi=300, bbox_inches='tight')

    plt.figure(3, figsize=(3*4.5, 3.75))
    for iext, extended_class_coadd in enumerate(extended_class_coadd_list):
      class_name = extended_class_coadd_names[str(extended_class_coadd)]
      if class_name.endswith('stars'):
        plt.subplot(131)
      else:
        plt.subplot(132)
      ext_class_cut = mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==extended_class_coadd
      plt.hist(mcal_gold_matches_bin['R11'][R11_gold_cut*ext_class_cut],  bins=50, histtype='step', normed=True, label=class_name)
      plt.xlabel('Shear Response $R_{11}$')

    plt.subplot(131)
    plt.title('Stars')
    plt.hist(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*star_hsc_cut], bins=50, histtype='step', normed=True, label='HSC Stars')
    plt.legend(fontsize='xx-small', ncol=1, loc='upper left')
    plt.ylabel('Normalised counts')
    plt.xlim([-3,3])
    plt.subplot(132)
    plt.title('Galaxies')
    plt.hist(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*gal_hsc_cut], bins=50, histtype='step', normed=True, label='HSC Galaxies')
    plt.hist(mcal_gold_matches_bin['R11'][R11_gold_cut*shape_mcal_gold_cut], bins=50, histtype='step', normed=True, label='Mcal Shape objects')
    plt.xlabel('Shear Response $R_{11}$')
    plt.legend(fontsize='xx-small', ncol=1, loc='upper left')
    plt.xlim([-3,3])
    #plt.yscale('log')
    plt.subplot(133)
    plt.title('Stars in Shape Catalogue')
    plt.hist(mcal_gold_matches_bin['R11'][R11_gold_cut*shape_mcal_gold_cut*(mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==0)], bins=30, histtype='step', normed=True, label='{0:.2f}\% High confidence stars'.format(100*N_high_shape_stars/N_shape))
    plt.hist(mcal_hsc_gold_matches_bin['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*star_hsc_cut], bins=30, histtype='step', normed=True, label='{0:.2f}\% HSC Stars, $m\\approx${1:.4f}'.format(100*N_hsc_shape_stars/N_shape, m_hsc))
    #plt.hist(mcal_gold_matches_bin['R11'][R11_gold_cut*shape_mcal_gold_cut*mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==0], bins=30, histtype='step', normed=True, label='high confidence stars')
    plt.hist(mcal_gold_matches_bin['R11'][R11_gold_cut*shape_mcal_gold_cut*(mcal_gold_matches_bin['EXTENDED_CLASS_COADD']==1)], bins=30, histtype='step', normed=True, label='{0:.2f}\% Candidate stars'.format(100*N_cand_shape_stars/N_shape))
    plt.xlabel('Shear Response $R_{11}$')
    plt.legend(fontsize='xx-small', ncol=1, loc='upper left')
    plt.suptitle('Field {1}, Redshift bin {0}'.format(ibin, field_name))
    plt.xlim([-3,3])
    plt.savefig('plots/{1}-R11-objecttype-bin{0}.png'.format(ibin, field_name), dpi=300, bbox_inches='tight')
  
  del(hsc_gold_matches)
  del(mcal_gold_matches)
  del(mcal_hsc_gold_matches)

colors = ['C0', 'C1', 'C2']
plt.close('all')
plt.figure(4, figsize=(4.5, 3.75))
x = np.array([1,2,3,4])
x = x - 0.1
for ifld,field_name in enumerate(field_list):
  plt.errorbar(x, np.array(m_per_sub[field_name]), yerr=np.asarray(sigma_m_per_sub[field_name], dtype=float), fmt='o', color=colors[ifld])#, label=field_name)
  plt.plot(x, np.array(m_per_sub[field_name]), '-', color=colors[ifld])#, label=field_name)
  x = x+0.1
plt.ylim([-0.008,0.008])
plt.legend(['SXDS', 'VVDS', 'DEEP2\_3'])
plt.xlabel('Redshift bin')
plt.ylabel('Stellar contamination shear bias $m$')
plt.savefig('./plots/stellar-contamination-biases.png', dpi=300, bbox_inches='tight')