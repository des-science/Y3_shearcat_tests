import numpy as np

from make_deep_fields_cats import *
from match_deep_master import *
from make_total_cats import *
from make_response_plots import *

field_list = [
              'c3',
              'x3',
              'e2'
              ]

mastercat_version = '03_31_20'

mastercat_dir = '/global/cscratch1/sd/troxel/cats_des_y3/'
mcal_bpz_dir = '/global/cscratch1/sd/itrharri/cats_des_y3/'
#mcal_bpz_dir = './data/'
#deep_dir = './data/'
deep_dir = '/project/projectdirs/des/www/y3_deep_cats/'

make_deep_cats(mastercat_version=mastercat_version,
               mastercat_dir=mastercat_dir,
               output_dir=mcal_bpz_dir,
               field_list=field_list)

match_deep_master(mastercat_version=mastercat_version,
                  mcal_bpz_dir=mcal_bpz_dir,
                  deep_dir=deep_dir,
                  output_dir=mcal_bpz_dir,
                  field_list=field_list)

make_total_cats(mastercat_version=mastercat_version,
                mcal_bpz_dir=mcal_bpz_dir,
                deep_dir=deep_dir,
                output_dir=mcal_bpz_dir,
                field_list=field_list)

make_response_plots(mastercat_version=mastercat_version,
                    mcal_bpz_dir=mcal_bpz_dir)