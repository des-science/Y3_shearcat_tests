import numpy as np
import healpy as hp
import h5py as h5
import pdb
from astropy.table import Table

def flux2mag(flux, zero_pt=30):
    return zero_pt - 2.5 * np.log10(flux)

def make_shape_cuts(cat, snr=[10,1000], flags=0, size_ratio=0.5, size=10, large_faint=[2, 30]):
  flag_cut = np.array(cat['catalog']['metacal']['unsheared']['flags'])==flags
  snr_cut = np.less(np.array(cat['catalog']['metacal']['unsheared']['snr']), snr[1])*np.greater(np.array(cat['catalog']['metacal']['unsheared']['snr']), snr[0])
  size_ratio_cut = np.greater(np.array(cat['catalog']['metacal']['unsheared']['size_ratio']), size_ratio)
  
  size_cut = np.less(np.array(cat['catalog']['metacal']['unsheared']['T']), size)
  large_faint_cut = ~(np.greater(np.array(cat['catalog']['metacal']['unsheared']['T']), large_faint[0])*np.less(np.array(cat['catalog']['metacal']['unsheared']['snr']), large_faint[1]))

  shape_cut = flag_cut*snr_cut*size_ratio_cut*size_cut#*large_faint_cut

  return shape_cut

def make_gold_cuts(cat, flags_foreground=0, flags_badregion=2, flags_gold=8, flags_footprint=1):
  foreground_cut = np.array(cat['catalog']['gold']['flags_foreground'])==flags_foreground
  badregion_cut = np.less(np.array(cat['catalog']['gold']['flags_badregions']), flags_badregion)
  gold_cut = np.less(np.array(cat['catalog']['gold']['flags_gold']), flags_gold)
  footprint_cut = np.array(cat['catalog']['gold']['flags_footprint'])==flags_footprint

  gold_cut = foreground_cut*badregion_cut*gold_cut*footprint_cut
  return gold_cut

def make_binary_cut(logT, magr, c=22.25, m=2.5, return_plotline=False):
  binaries = logT < (c - magr)/m

  if return_plotline:

    y = np.linspace(5,35,128)
    x = np.power(10.,(c - y)/m)

    return binaries, x, y
  else:
    return binaries

mastercat_version = '03_31_20'

#mastercat_dir = '/global/cscratch1/sd/troxel/cats_des_y3/'
mastercat_dir = '/project/projectdirs/des/www/y3_cats/'

master = h5.File(mastercat_dir + 'Y3_mastercat_{}.h5'.format(mastercat_version), 'r')

gpix  = master['catalog/gold/hpix_16384'][:]//(hp.nside2npix(16384)//hp.nside2npix(4096))
mask_select = ~np.in1d(gpix//(hp.nside2npix(16384)//hp.nside2npix(4096)),master['index/mask/hpix'][:],assume_unique=False)

highe_cut = np.greater(np.sqrt(np.power(master['catalog/metacal/unsheared']['e_1'],2.) + np.power(master['catalog/metacal/unsheared']['e_2'],2)), 0.8)

highe_shape_cut = highe_cut

magr = flux2mag(master['catalog/metacal/unsheared']['flux_r'][mask_select*highe_shape_cut])
magi = flux2mag(master['catalog/metacal/unsheared']['flux_i'][mask_select*highe_shape_cut])
magz = flux2mag(master['catalog/metacal/unsheared']['flux_z'][mask_select*highe_shape_cut])
coadd_id = master['catalog/metacal/unsheared']['coadd_object_id'][mask_select*highe_shape_cut]
highe_shape_weights = master['catalog/metacal/unsheared']['weight'][mask_select*highe_shape_cut]

logT = np.log10(master['catalog']['metacal']['unsheared']['T'][mask_select*highe_shape_cut])
snr = master['catalog']['metacal']['unsheared']['snr'][mask_select*highe_shape_cut]

R_11 = master['catalog']['metacal']['unsheared']['R11'][mask_select*highe_shape_cut]

binaries, x, y = make_binary_cut(logT, magr, c=22.25, m=3.5, return_plotline=True)

output_cat = Table()

output_cat['coadd_id'] = coadd_id
output_cat['logT'] = logT
output_cat['magr'] = magr
output_cat['magi'] = magi
output_cat['magz'] = magz
output_cat['snr'] = snr
output_cat['R11'] = R_11
output_cat['binaries'] = binaries

output_cat.write('./data/binary_responses_{0}.fits'.format(mastercat_version))
#np.savetxt('./data/binary_cut_line.txt', np.column_stack([x, y]))
