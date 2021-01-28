import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
from matplotlib import rc

import healpy as hp
import pickle
import pdb

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=11)

plt.close('all') # tidy up any unshown plots

seeing_map = hp.read_map('./data/y3a2_i_o.4096_t.32768_FWHM.WMEAN_EQU.fits')
nside = hp.npix2nside(seeing_map.shape[0])

c3_cen = SkyCoord(52.6484, -28.1000, frame='icrs', unit='deg')
x3_cen = SkyCoord(36.4500, -4.6000, frame='icrs', unit='deg')
e2_cen = SkyCoord(9.5000, -43.9980, frame='icrs', unit='deg')

c3_pix = hp.ang2pix(nside, c3_cen.ra.deg, c3_cen.dec.deg, lonlat=True)
x3_pix = hp.ang2pix(nside, x3_cen.ra.deg, x3_cen.dec.deg, lonlat=True)
e2_pix = hp.ang2pix(nside, e2_cen.ra.deg, e2_cen.dec.deg, lonlat=True)

plt.close('all') # tidy up any unshown plots
plt.figure(1, figsize=(2*4.5, 3.75))
#plt.subplot(121)
hp.mollview(seeing_map, fig=1, sub=(121), title='Seeing FWHM', unit='arcsec')
hp.projscatter(c3_cen.ra.deg, c3_cen.dec.deg, lonlat=True, c='c')
hp.projtext(c3_cen.ra.deg, c3_cen.dec.deg, 'C3', lonlat=True)

hp.projscatter(x3_cen.ra.deg, x3_cen.dec.deg, lonlat=True, c='m')
hp.projtext(x3_cen.ra.deg, x3_cen.dec.deg, 'X3', lonlat=True)

hp.projscatter(e2_cen.ra.deg, e2_cen.dec.deg, lonlat=True, c='y')
hp.projtext(e2_cen.ra.deg, e2_cen.dec.deg, 'E2', lonlat=True)

plt.subplot(122)
plt.hist(seeing_map[seeing_map>0], bins=30, density=True, histtype='step')
plt.axvline(seeing_map[c3_pix], c='c', label='C3')
plt.axvline(seeing_map[x3_pix], c='m', label='X3')
plt.axvline(seeing_map[e2_pix], c='y', label='E2')
plt.legend()
plt.xlabel('Seeing FWHM [arcsec]')
plt.xticks([])
plt.savefig('./plots/seeing.png', dpi=300, bbox_inches='tight')