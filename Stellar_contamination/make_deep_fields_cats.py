import numpy as np
import h5py as h5
from astropy.table import Table

from region_defs import deep_fields

def make_deep_cats(mastercat_version, mastercat_dir, output_dir, field_list):

  master = h5.File(mastercat_dir + 'Y3_mastercat_{}.h5'.format(mastercat_version), 'r')

  index_keys = ['select']
  gold_keys = ['extended_class_mash_sof']
  mcal_keys = ['R11', 'flags', 'snr', 'size_ratio', 'ra', 'dec', 'coadd_object_id', 'e_1', 'e_2', 'T']
  bpz_keys = ['zmean_sof']

  select_bool = np.zeros(len(master['catalog/metacal/unsheared']['R11']), dtype=bool)
  select_bool[master['index']['select']] = True

  for field_name in field_list:

    field = deep_fields[field_name]

    ra_cut = (np.greater(master['catalog/gold']['ra'], field['ra_min']))*(np.less(master['catalog/gold']['ra'], field['ra_max']))
    dec_cut = (np.greater(master['catalog/gold']['dec'], field['dec_min']))*(np.less(master['catalog/gold']['dec'], field['dec_max']))

    field_cut = ra_cut*dec_cut

    output_cat = Table()

    output_cat['select_bool'] = select_bool[field_cut]

    for key in mcal_keys:
      print(field['name']+' '+key+'...')
      output_cat[key] = master['catalog/metacal/unsheared'][key][field_cut]

    for key in bpz_keys:
      print(field['name']+' '+key+'...')
      output_cat[key] = master['catalog/bpz/unsheared'][key][field_cut]

    for key in gold_keys:
      print(field['name']+' '+key+'...')
      output_cat[key] = master['catalog/gold'][key][field_cut]

    output_cat.write(output_dir + '{0}.mcal_bpz_Y3_mastercat_{1}.h5'.format(field['name'], mastercat_version), format='hdf5', path='catalog/mcal_bpz', overwrite=True)
