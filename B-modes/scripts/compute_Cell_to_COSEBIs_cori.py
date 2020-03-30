# Script to plot WnLog(ell) and convert PCells to COSEBIs
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius

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
    # This function needs to be completely updated to Cyrille's latest
    cdir = '/users/PCON0003/cond0080/scripts/y3_bmode'
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

#plt.plot(np.arange(1,21), cd_cosebis/1.E-10, 'k-')
#print cd_cosebis
#np.savetxt('cd_cosebis_from_Pbb_select_shape.txt',cd_cosebis, fmt='%.18e')
#plt.show()

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
