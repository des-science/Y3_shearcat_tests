import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as uns


def match_deep_master(mastercat_version, mcal_bpz_dir, deep_dir, output_dir, field_list, match_distance=1.*uns.arcsec):

  for field_name in field_list:

      mcal_cat = Table.read(mcal_bpz_dir + '{0}.mcal_bpz_Y3_mastercat_{1}.h5'.format(field_name, mastercat_version))
      class_cat = Table.read(deep_dir + '{0}.ugriz-mof02-JHK-ff03_extcorr_NearestNeighbors_class.fits'.format(field_name))

      match_fname = output_dir + '{0}_mcal_NearestNeighbors_class_match_Y3_mastercat_{1}.fits'.format(field_name, mastercat_version)

      print('Matching {0} mcal to class catalogue...'.format(field_name))
      mcal_coord = SkyCoord(mcal_cat['ra']*uns.deg, mcal_cat['dec']*uns.deg)
      class_coord = SkyCoord(class_cat['ra']*uns.deg, class_cat['dec']*uns.deg)
      idx_class, idx_mcal, d2d, d3d = SkyCoord.search_around_sky(mcal_coord, class_coord, match_distance)
      match_flag = Table()
      match_flag['idx_mcal']= idx_mcal
      match_flag['idx_class']= idx_class
      match_flag.write(match_fname)
      print('...done')