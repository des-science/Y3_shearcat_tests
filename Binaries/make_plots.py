import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib import rc
import pickle

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=11)

def frexp10(x):
    # convenience function to print standard form numbers in latex
    exp = int(np.floor(np.log10(abs(x))))
    return x / 10**exp, exp

mastercat_version = '03_31_20'

cat = Table.read('./data/binary_cut_{0}.fits'.format(mastercat_version))
stellar_loci = Table.read('data/covey_medianlocus.tbl', format='cds') # Covey et al (2007) stellar loci
cut_x, cut_y = np.loadtxt('./data/binary_cut_line.txt', unpack=True)

h, xedges, yedges, _ = plt.hist2d(cat['logT'], cat['magr'], bins=[np.linspace(-1,1,100), np.linspace(20,25,100)], cmap='Blues')

plt.close('all')
plt.figure(1, figsize=(3*4.5, 3.75))

plt.subplot(131)
plt.title('All high $\\left|e\\right|$')
plt.hist2d(cat['logT'][~cat['binaries']], cat['magr'][~cat['binaries']], bins=(xedges, yedges), cmap='Blues', cmin=1.e-9, label='Pass')
plt.hist2d(cat['logT'][cat['binaries']], cat['magr'][cat['binaries']], bins=(xedges, yedges), cmap='Reds', cmin=1.e-9, label='Fail')
#plt.hist2d(cat['logT'], cat['magr'], bins=(xedges, yedges), cmap='Blues')#, cmin=1.e-9)
plt.plot(np.log10(cut_x), cut_y, 'k-', label='Binary cut')
plt.legend(loc='upper right')
plt.xlabel('$\log_{10}(T)$')
plt.ylabel('$r \,$[mag]')

plt.subplot(132)
plt.title('Pass binary cut')
h, xedges, yedges, _ = plt.hist2d(cat['magi'][~cat['binaries']]-cat['magz'][~cat['binaries']], cat['magr'][~cat['binaries']]-cat['magi'][~cat['binaries']], bins=[np.linspace(-0.5,2.5,50),np.linspace(-0.5,2.5,50)], cmap='Blues')
plt.plot(stellar_loci['i-z'], stellar_loci['r-i'], 'ko', ms=1, alpha=0.4, label='Stellar locus')
plt.legend(loc='upper right')
plt.xlabel('$i - z\,$[mag]')
plt.ylabel('$r - i\,$[mag]')

plt.subplot(133)
plt.title('Fail binary cut $N_{\\rm binaries}$'+'$={0:.2f}\\times 10^{{{1}}}$'.format(*frexp10(np.sum(cat['binaries']))))
plt.hist2d(cat['magi'][cat['binaries']]-cat['magz'][cat['binaries']], cat['magr'][cat['binaries']]-cat['magi'][cat['binaries']], bins=(xedges, yedges), cmap='Reds')
plt.plot(stellar_loci['i-z'], stellar_loci['r-i'], 'ko', ms=1, alpha=0.4)
plt.xlabel('$i - z\,$[mag]')

plt.subplots_adjust(wspace=0.25)
plt.savefig('./plots/binaries.png', dpi=300, bbox_inches='tight')
plt.savefig('./plots/binaries.pdf', bbox_inches='tight')

sav_obj = {
           'cat' : cat,
           'stellar_loci' : stellar_loci,
           'cut_x' : cut_x,
           'cut_y' : cut_y
           }

with open('everything_you_need_for_binaries.pkl', 'wb') as f:
  pickle.dump(sav_obj, f, pickle.HIGHEST_PROTOCOL)