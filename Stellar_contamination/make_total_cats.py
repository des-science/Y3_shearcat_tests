import numpy as np
from astropy.table import Table, vstack, hstack
import os
import pdb


def make_total_cats(mastercat_version, mcal_bpz_dir, deep_dir, output_dir, field_list):
  
  total_mcal_cat = Table()
  total_class_cat = Table()
  total_cat = Table()

  for field_name in field_list:

    print('Adding {0} to total catalogue...'.format(field_name))

    mcal_cat = Table.read(mcal_bpz_dir + '{0}.mcal_bpz_Y3_mastercat_{1}.h5'.format(field_name, mastercat_version))
    class_cat = Table.read(deep_dir + '{0}.ugriz-mof02-JHK-ff03_extcorr_NearestNeighbors_class.fits'.format(field_name))
    match_cat = Table.read(output_dir + '{0}_mcal_NearestNeighbors_class_match_Y3_mastercat_{1}.fits'.format(field_name, mastercat_version))
    
    total_mcal_cat = vstack([total_mcal_cat, mcal_cat])
    total_class_cat = vstack([total_class_cat, class_cat])

    idx_mcal = match_cat['idx_mcal']
    idx_class = match_cat['idx_class']

    mcal_cat = mcal_cat[idx_mcal]
    class_cat = class_cat[idx_class]
    total_cat = vstack([total_cat, hstack([mcal_cat, class_cat])])

  total_mcal_cat.write(output_dir + 'total.mcal_bpz_Y3_mastercat_{0}.fits'.format(mastercat_version))
  class_cat.write(output_dir + 'total.ugriz-mof02-JHK-ff03_extcorr_NearestNeighbors_class.fits')
  total_cat.write(output_dir + 'total_mcal_NearestNeighbors_class_match_Y3_mastercat_{0}.fits'.format(mastercat_version))
  print('...done')
