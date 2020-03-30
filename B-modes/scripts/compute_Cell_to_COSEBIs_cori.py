# Script to plot WnLog(ell) and convert PCells to COSEBIs
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import pickle

col_arr = ['#d55e00','#cc79a7','#0072b2','#009e73','#f0e442','#000000']

#nmode = 1
# Function to read the WnLog tables
def get_wnlog(nmode):
    wdir = '/project/projectdirs/des/y3-sheartests//WnLog' #location on cori
    file = open('%s/WnLog%s-2.50-250.00.table'%(wdir, str(nmode)))
    ell, wn_ell = np.loadtxt(file, unpack=True)
    file.close()
    return np.exp(ell), wn_ell

def get_cbb_ell():
    # Set up to read PCl data vector from Cyrille in repo
    cdir = '../data/'
    cfile = '%s/cls_Y3_mastercat_12_3_19_v3_Bmode_JK_non-tomo_highell_nopureE_countmode_C1_apo0.0_nside1024.pickle'%cdir
    cdata = pickle.load(open(cfile, 'rb'))
    ell = np.linspace(0,3071,num=3072)
    # cdata['ell'] and cdata[('galaxy_shear', 'galaxy_shear')][(0,0)]['data'][3,:] have 32 values, which could be taken as stepwise or interpolated
    # Interpolation is extremely noisy, try a stepwise function instead
    #P_BB_model = ius(cdata['ell'], cdata[('galaxy_shear', 'galaxy_shear')][(0,0)]['data'][3,:])
    P_BB = cdata[('galaxy_shear', 'galaxy_shear')][(0,0)]['data'][3,:]
    mid = cdata['ell'][:-1]+np.diff(cdata['ell'])/2. # midpoint
    ell_bins = np.append(np.append(0, mid), 3071) # bin edges from 0 to 3071
    P_BB_model = np.zeros_like(ell)
    # Probs there is a better way of writing this for loop
    for idx in range(ell_bins.size-1):
        P_BB_model[(ell>ell_bins[idx])&(ell<=ell_bins[idx+1])]=P_BB[idx]
    return ell, P_BB_model #P_BB_model(ell)
    
def old_get_cbb_ell():
    # This function read in a previous data vector from Cyrille
    cdir = '/users/PCON0003/cond0080/scripts/y3_bmode' # in Ohio
    cfile = open('%s/Y3_mastercat_12_3_19_v2_Bmode_highell_nopureE_countmode_C1_apo0.0_nside1024_cls.txt'%cdir)
    data = np.loadtxt(cfile)
    cfile.close()
    ell = np.linspace(0,3071,num=3072)
    #P_BB = data[0,:].copy() #select (stricter cut)
    P_BB = data[5,:].copy() #select_shape (less strict cut)
    return ell, P_BB

cd_cosebis = np.zeros(20)

# Loop over each COSEBI mode
for ndx in range(20):
    n_ell, n_wn_ell = get_wnlog(ndx+1)
    cd_ell, cd_p_bb = get_cbb_ell()
    # Interpolate
    wnlog_model = ius(n_ell, n_wn_ell)
    pbb_model = ius(cd_ell, cd_p_bb)
    # Get COSEBI for this mode
    ell_eval = np.linspace(0,3071,num=3072)
    B_integrand = ius(ell_eval, (ell_eval/(2*np.pi))*wnlog_model(ell_eval)*pbb_model(ell_eval))
    cd_cosebis[ndx] = B_integrand.integral(0,3071)

plt.plot(np.arange(1,21), cd_cosebis/1.E-10, 'k-')
print(cd_cosebis)
np.savetxt('../data/cd_cosebis_from_Pbb_12_3_19_v3.txt',cd_cosebis, fmt='%.18e')
plt.show()

"""
plt.title("Filters for 2.5'-250.0'",fontsize=14)
plt.xlabel('ell',fontsize=14)
plt.ylabel('Wn(ell)',fontsize=14)
plt.legend(loc='best')

#plt.xlim(1., 1.E5)
#plt.xlim(1., 10.)
plt.xscale('log')

plt.show()

"""
