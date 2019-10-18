import numpy as np
import scipy
from scipy import spatial

def covariance_scalar_jck(TOTAL_PHI,jk_r, type_c = 'jackknife'):

  #  Covariance estimation
  if type_c == 'jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_c=='bootstrap':
      fact=1./(jk_r)
        
  average=0.
  cov_jck=0.
  err_jck=0.


  for kk in range(jk_r):
    average+=TOTAL_PHI[kk]
  average=average/(jk_r)

  for kk in range(jk_r):
    #cov_jck+=TOTAL_PHI[kk]#*TOTAL_PHI[kk]

    cov_jck+=(-average+TOTAL_PHI[kk])*(-average+TOTAL_PHI[kk])


  err_jck=np.sqrt(cov_jck*fact)


  #average=average*(jk_r)/(jk_r-1)
  return {'cov' : cov_jck*fact,
          'err' : err_jck,
          'mean': average}


def analysis(save_output_folder,kE_label,sys_map,mask,mapa, hpix, fract_limits,len_hist=10,mapa_weight = None, cov_ext = False):
    from scipy import linalg
    print ('\n****************\n')
    print ('TEST [method1] {1}: -> {0}'.format(sys_map['title'],kE_label))
    
    if cov_ext:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[cov_ext['indexes']])
                                      
        mask = cov_ext['indexes'] # it can be slightly different
    else:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[mask])     
    fractional = sys_map['fractional']
    
    # plot typical ranges ***********************
    area_tot = compute_area(mask, nside)
    y,x = np.histogram(sys_map['map'][mask], bins = 50)
    sys_map1 =  sys_map['map'][mask]
    mapa1 = mapa[mask]-mapa_weight[mask]

    #print len(mapa1),len(sys_map1)
    #mapa1 = mapa1 - np.mean(mapa1)
    
    ave = sys_map1.mean()
    std = sys_map1.std()
    if fractional:
    
        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
    
    
    min_s, max_s = find_limits(sys_map1, fract_limits)
    if not fractional:
        if 'min_f' in sys_map.keys():
            min_s = max([min_s,sys_map['min_f']-ave])
        if 'max_f' in sys_map.keys():
            max_s = min([max_s,sys_map['max_f']-ave])   
    if fractional:
        min_ss = min_s*ave + ave
        max_ss = max_s*ave + ave    
    else:
        min_ss = min_s + ave
        max_ss = max_s + ave

    
    def histedges_equalN(x, nbin):
        npt = len(x)
        return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

    mh = (sys_map1>min_s) & (sys_map1<max_s)
    bins = histedges_equalN(sys_map1[mh], len_hist)

    
    #print np.histogram(sys_map1, bins= bins)
    
    #define bin centers as the median point.
    x = np.zeros(len(bins)-1)
    for i in range(len(x)):

        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        x[i] = np.median(sys_map1[mask_bin])
        
        
        #print bins[i],x[i],bins[i+1]
    #x = bins[:-1]+0.5*(bins[1]-bins[0])
     
    mask_asd = (sys_map1 > min_s ) & (sys_map1 < max_s )
    mask_pixels = (sys_map1 < min_s) ^ ( sys_map1 > max_s)
    area_rest = compute_area(mask_pixels, nside)
    print ('Area excluded: {0:.1f}/{1:.1f}'.format(area_rest,area_tot))
    

    # plot maps proerties *************************
    value_jck = np.zeros((len_hist,n_jck+1))



    for i in range(len_hist):
        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        value_jck[i,0] = (mapa1[mask_bin]).sum()
        for j in range(n_jck):
            #print ((bins[i],bins[i+1]))
            
            mask_bin1 = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) & (hpix[mask_hpix] == j)
            value_jck[i,j+1] = (value_jck[i,0] - (mapa1[mask_bin1]).sum())/(len(mapa1[mask_bin])-len(mapa1[mask_bin1]))
        value_jck[i,0] = value_jck[i,0]/ len(mapa1[mask_bin])
    dict_jck = covariance_jck(value_jck[:,1:], n_jck, 'jackknife')


    

    params_guess= np.array([0.0,0])
    params_guess2= np.array([.0,0.0,0])
    params_guess3= np.array([.0,.0,0.0,0])
    if not cov_ext:

        ccov = dict_jck#['cov']
        
    else:
        vector_cov = np.zeros((len_hist,cov_ext['cov'].shape[1]))
        for i in range(len_hist):
            mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1])# & cov_ext['indexes']
            vector_cov[i,:] = np.mean(cov_ext['cov'][mask_bin,:],axis = 0)
        mask_bin = (sys_map1 > bins[0]) & (sys_map1 < bins[-1])

        ccov = covariance_jck(vector_cov, cov_ext['cov'].shape[1], 'bootstrap')
        '''
        mute= ccov['cov'].diagonal()
        cov_mute = np.zeros((len(mute),len(mute)))
        for ff in range(len(mute)):
            cov_mute[ff,ff] = mute[ff]
        ccov['cov'] = cov_mute
        '''
        
    dict_output=dict()
    if not cov_ext:

        poptbts, pcovbts = curve_fit(fitting_linear, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess)
        poptbts2, pcovbts2 = curve_fit(fitting_2nd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess2)
        poptbts3, pcovbts3 = curve_fit(fitting_3rd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess3)

        paramsbts=poptbts
        paramsbts2=poptbts2
        paramsbts3=poptbts3
        #print ('slope = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))
        wbt = value_jck[:,0]
        cns = np.zeros(len(wbt))* np.mean(mapa1[mask_asd])
        inv_covbt = linalg.inv(ccov['cov'])
        
        # jackknife values ********************************************************
        cns = fitting_2nd(x,poptbts2[0],poptbts2[1],poptbts2[2]) 
        chi2redbt2 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        cns = fitting_3rd(x,poptbts3[0],poptbts3[1],poptbts3[2],poptbts3[3]) 
        chi2redbt3 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = fitting_linear(x,poptbts[0],poptbts[1]) 
        chi2redbt =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = np.ones(len(wbt))* np.mean(mapa1[mask_asd])
        chi2redbt_null = (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        #print 'pull:'
        ##print (wbt-cns)/ccov['err'],np.sum(((wbt-cns)**2/ccov['err']**2.))
        chi2redbt_diff = chi2redbt_null - chi2redbt 


  
        #********************************************************
        dict_output.update({'chi2red_null': chi2redbt_null})
        dict_output.update({'chi2red_diff': chi2redbt_diff})
        dict_output.update({'chi2red': chi2redbt})
        dict_output.update({'chi2red2': chi2redbt2})
        dict_output.update({'chi2red3': chi2redbt3})

        popt, pcov = curve_fit(fitting_linear, x, value_jck[:,0],sigma = dict_jck['cov'],p0=params_guess)

        params=popt
        print ('slope jck = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))


        
    else:
        
        # all the mocks *****************
        chi2_vect = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect1 = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect2 = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect3 = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect4 = np.zeros(cov_ext['cov'].shape[1]) 
        
        #print 'mean mocks',np.mean(cov_ext['cov'][:,0])
        for jj in range(cov_ext['cov'].shape[1]):

            popt, pcov = curve_fit(fitting_linear, x, vector_cov[:,jj],sigma = ccov['cov'],p0=params_guess)

            params=popt
            #print ('slope jck = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))
            w = vector_cov[:,jj]
            #print np.mean(cov_ext['cov'][:,j])
            cns = np.ones(len(w))*  np.mean(cov_ext['cov'][:,jj])
            inv_cov = linalg.inv(ccov['cov'])
            chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            #plt.plot(x,cns)

            cns = fitting_linear(x,popt[0],popt[1]) 

            chi2_vect[jj] = chi2red - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            chi2_vect1[jj] = cimhi2red# - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            chi2_vect2[jj] = (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            

            mutw = chi2red/(len(x)-2) - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))//(len(x)-1) 
            df = 1
            p = 1 - stats.t.cdf(mutw,df=df)
            chi2_vect3[jj] = p2s(p)
            
            mutw = chi2red - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))
            df = 1
            p = 1 - stats.t.cdf(mutw,df=df)
            chi2_vect4[jj] = p2s(p)


    
        m1 = np.sort(chi2_vect1)
        m2 = np.sort(chi2_vect2)
        m3 = np.sort(chi2_vect)
        m4 = np.sort(chi2_vect3)
        m5 = np.sort(chi2_vect4)
        chi681 = m1[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        chi682 = m2[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        chi683 = m4[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        chi684 = m5[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        # chi68 out of 1000 sims
        chi68 = m3[int(0.68*np.float(cov_ext['cov'].shape[1]))]

        poptbts, pcovbts = curve_fit(fitting_linear, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess)
        poptbts2, pcovbts2 = curve_fit(fitting_2nd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess2)
        poptbts3, pcovbts3 = curve_fit(fitting_3rd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess3)

        paramsbts=poptbts
        paramsbts2=poptbts2
        paramsbts3=poptbts3
        #print ('slope resampling = {0:.4f} +- {1:.4f}').format(paramsbts[0],np.sqrt(pcovbts[0,0]))
        wbt = value_jck[:,0]
        #print 'mean data',np.mean(mapa1[mask_asd])
        
        inv_covbt = linalg.inv(ccov['cov'])
        
        # bootstrap values ********************************************************
        cns = fitting_2nd(x,poptbts2[0],poptbts2[1],poptbts2[2]) 
        chi2redbt2 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        cns = fitting_3rd(x,poptbts3[0],poptbts3[1],poptbts3[2],poptbts3[3]) 
        chi2redbt3 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = fitting_linear(x,poptbts[0],poptbts[1]) 
        chi2redbt =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = np.ones(len(wbt))* np.mean(mapa1[mask_asd])
        chi2redbt_null = (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        #print 'pull:'
        ##print (wbt-cns)/ccov['err'],np.sum(((wbt-cns)**2/ccov['err']**2.))
        chi2redbt_diff = chi2redbt_null - chi2redbt 
        
        #print "                        chi2/dof (d chi2/dof) [d chi2/dof _68]"
        #print ('chi2/dof resampling :  {0:.2f} ({2:.3f}) [{3:.3f}]').format(chi2redbt,sys_map['title'],chi2redbt_diff,chi2redbt_diff/chi68)
         
            
        #********************************************************
        dict_output.update({'chi2red_boot_null': chi2redbt_null})
        dict_output.update({'chi2red_diff_boot': chi2redbt_diff})
        dict_output.update({'chi2red_boot': chi2redbt})
        dict_output.update({'chi2red_boot2': chi2redbt2})
        dict_output.update({'chi2red_boot3': chi2redbt3})
        dict_output.update({'chi2red_diff_boot68': chi2redbt_diff/chi68}) 
        dict_output.update({'chi68': chi68}) 
        popt, pcov = curve_fit(fitting_linear, x, value_jck[:,0],sigma = dict_jck['cov'],p0=params_guess)

        params=popt
        print ('slope jck = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))
        w = value_jck[:,0]
        cns = np.ones(len(w))* np.mean(mapa1[mask_asd])
        inv_cov = linalg.inv(dict_jck['cov'])
        chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
        
        cns = fitting_linear(x,popt[0],popt[1]) 
        chi2red_diff = chi2red - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
        
        print ('chi2/dof jck  <{1}>0 =<{1}>:  {0:.2f} ({2:.3f})  [{3:.3f}]').format(chi2red,sys_map['title'],chi2red_diff,chi2red_diff/chi68)
        
        dict_output.update({'chi2red_diff': chi2red_diff})
        dict_output.update({'chi2red_diff68': chi2red_diff/chi68})
    

  
    
    if not cov_ext:
        dict_output.update({'dict_cov' : ccov})
        dict_output.update({'params_bts' :paramsbts})
        dict_output.update({'params_bts2' :paramsbts2})
        dict_output.update({'params_bts3' :paramsbts3})
        dict_output.update({'params_cov_bts' : pcovbts})  
    else:
        dict_output.update({'dict_cov' : dict_jck})
        dict_output.update({'dict_cov_bts' : ccov})
        dict_output.update({'params_bts' :paramsbts})
        dict_output.update({'params_bts2' :paramsbts2})
        dict_output.update({'params_bts3' :paramsbts3})
        dict_output.update({'params_cov_bts' : pcovbts})        
        
        dict_output.update({'v':vector_cov})
        dict_output.update({'chi68_new':chi683})
        dict_output.update({'chi68_new1':chi684})
        dict_output.update({'chi68_array':m4})
    dict_output.update({'w' :value_jck[:,0]})
    dict_output.update({'x' :x})   
    dict_output.update({'bins' :bins})       

    dict_output.update({'params' :params})
    dict_output.update({'params_cov' : pcov})
    dict_output.update({'min_s': min_s})
    dict_output.update({'max_s': max_s})
    dict_output.update({'min_ss': min_ss})
    dict_output.update({'max_ss': max_ss})
    dict_output.update({'ave_value': ave})
    dict_output.update({'mask_pixel': mask_pixels})
    dict_output.update({'value_jck': value_jck})
    dict_output.update({'Area': area_rest/area_tot})
    dict_output.update({'mean shear': np.mean(mapa1[mask_asd])})


    if not cov_ext:
        label_boot =''
    else:
        label_boot ='_boot'
    if 1==1:
        chi = dict_output['chi2red'+label_boot+'_null'] 
        df1 = len(dict_output['x']) - 1
        chi=chi/df1
        p = 1 - stats.t.cdf(chi,df=df1)
        #print p2s(p)
        dict_output.update({'s_null': p2s(p)})
        chi = dict_output['chi2red'+label_boot] 
        df2 = len(dict_output['x']) - 2
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
    
        dict_output.update({'s_fit': p2s(p)})
        #print p2s(p)
    
        chi = dict_output['chi2red'+label_boot+'2'] 
        df2 = len(dict_output['x']) - 3
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
    
        dict_output.update({'s_fit2': p2s(p)})
        #print p2s(p)
    
        chi = dict_output['chi2red'+label_boot+'3'] 
        df2 = len(dict_output['x']) - 4
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
    
        dict_output.update({'s_fit3': p2s(p)})
        #print p2s(p)
        chi = dict_output['chi2red'+label_boot+'_null']/df1 - dict_output['chi2red'+label_boot+'']/df2 
        df = 1
        p = 1 - stats.t.cdf(chi,df=df)
        

        
        #print p2s(p)
        dict_output.update({'s_H1': p2s(p)})
        
        
        
        
        #******************************************
        chi = dict_output['chi2red'+label_boot] 
        df2 = len(dict_output['x']) - 2
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'_null'] ,dict_output['chi2red'+label_boot] ,len(dict_output['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)
        dict_output.update({'F01': p2s(pf)})
    
        #******************************************
        chi = dict_output['chi2red'+label_boot+'2'] 
        df2 = len(dict_output['x']) - 3
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'_null'] ,dict_output['chi2red'+label_boot+'2'] ,len(dict_output['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot] ,dict_output['chi2red'+label_boot+'2'] ,len(dict_output['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        dict_output.update({'F02': p2s(pf)})
        dict_output.update({'F12': p2s(pf2)})
    
        #******************************************
        chi = dict_output['chi2red'+label_boot+'3'] 
        df2 = len(dict_output['x']) - 4
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'_null'] ,dict_output['chi2red'+label_boot+'3'] ,len(dict_output['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot] ,dict_output['chi2red'+label_boot+'3'] ,len(dict_output['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
 
        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'2'] ,dict_output['chi2red'+label_boot+'3'] ,len(dict_output['x']),1,2)  
        pf3 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        print ('ok')
        dict_output.update({'F03': p2s(pf)})
        dict_output.update({'F13': p2s(pf2)})
        dict_output.update({'F23': p2s(pf3)})  
        

    return dict_output

def do_analysis(out_file,add_label,kE,kE_label,weight,systematic_maps,info,fract_limits,len_hist=10,hpix_type = 'normal',cov_ext = False, fit ='linear'):
    systematics_dict = dict() 
    systematics_dict_2 = dict() 
    for key in systematic_maps.keys():
      
      if key != kE_label:
       if ((key == 'snr' and kE_label =='kE') or (key == 'size_ratio' and kE_label =='kE') or (key == 'E1' and kE_label =='kE') or (key == 'E2' and kE_label =='kE')):
        pass
       else:
        if hpix_type == 'normal':
            dict_output = analysis(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix'], fract_limits,len_hist=len_hist,mapa_weight = weight, cov_ext = cov_ext)
            dict_output_2 = analysis_2(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix'], fract_limits,len_hist=len_hist,mapa_weight = weight, cov_ext = cov_ext)
           # print  dict_output_2['b_arr'][0], ' +- ',  dict_output_2['b_err']
      
        if hpix_type == 'marc':
            dict_output = analysis(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix_f'], fract_limits,len_hist=len_hist, mapa_weight = weight, cov_ext = cov_ext)
            dict_output_2 = analysis_2(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix_f'], fract_limits,len_hist=len_hist,mapa_weight = weight, cov_ext = cov_ext)

        make_plot(out_file,add_label,dict_output,systematic_maps[key]['title'],kE_label,systematic_maps[key]['map'],systematic_maps[key]['fractional'],info, fract_limits,len_hist = len_hist,cov_ext = cov_ext,fit=fit,dict_output_2=dict_output_2)
        systematics_dict.update({systematic_maps[key]['title'] : dict_output})
        systematics_dict_2.update({systematic_maps[key]['title'] : dict_output_2})

    return systematics_dict,systematics_dict_2


def derivative(x,f):
            xmi = int(len(x)/2.)
            dx_i = x[xmi+1]-x[xmi]
            dx_mi = x[xmi]-x[xmi-1]
            f_i = f[xmi+1]
            f_0 = f[xmi]
            f_mi = f[xmi-1]
            b = -dx_i*f_mi / (dx_mi*(dx_i+dx_mi)) + (dx_i-dx_mi)*f_0/(dx_i*dx_mi)+dx_mi*f_i/(dx_i*(dx_i+dx_mi))
            return b,xmi

        
def analysis_2(save_output_folder,kE_label,sys_map,mask,mapa, hpix, fract_limits,len_hist=10,mapa_weight = None, cov_ext = False):
    from scipy import linalg
    print ('\n****************\n')
    print ('TEST [method2] {1}: -> {0}'.format(sys_map['title'],kE_label))
    
    if cov_ext:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[cov_ext['indexes']])                        
        mask = cov_ext['indexes'] # it can be slightly different
    else:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[mask])   
    fractional = sys_map['fractional']
    
    # plot typical ranges ***********************
    area_tot = compute_area(mask, nside)
    y,x = np.histogram(sys_map['map'][mask], bins = 50)
    sys_map1 =  sys_map['map'][mask]
    mapa1 = mapa[mask]-mapa_weight[mask]

    
    ave = sys_map1.mean()
    std = sys_map1.std()
    if fractional:
    
        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
    
    
    min_s, max_s = find_limits(sys_map1, fract_limits)
    if not fractional:
        if 'min_f' in sys_map.keys():
            min_s = max([min_s,sys_map['min_f']-ave])
        if 'max_f' in sys_map.keys():
            max_s = min([max_s,sys_map['max_f']-ave])   
    if fractional:
        min_ss = min_s*ave + ave
        max_ss = max_s*ave + ave    
    else:
        min_ss = min_s + ave
        max_ss = max_s + ave

    
    def histedges_equalN(x, nbin):
        npt = len(x)
        return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

    mh = (sys_map1>min_s) & (sys_map1<max_s)
    bins = histedges_equalN(sys_map1[mh], len_hist)

    
    #print np.histogram(sys_map1, bins= bins)
    
    #define bin centers as the median point.
    x = np.zeros(len(bins)-1)
    for i in range(len(x)):

        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        x[i] = np.median(sys_map1[mask_bin])
        

    mask_asd = (sys_map1 > min_s ) & (sys_map1 < max_s )
    mask_pixels = (sys_map1 < min_s) ^ ( sys_map1 > max_s)
    area_rest = compute_area(mask_pixels, nside)
    print ('Area excluded: {0:.1f}/{1:.1f}'.format(area_rest,area_tot))
    

    # plot maps proerties *************************
    value_jck = np.zeros((len_hist,n_jck+1))



    for i in range(len_hist):
        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        value_jck[i,0] = (mapa1[mask_bin]).sum()
        for j in range(n_jck):
            #print ((bins[i],bins[i+1]))
            
            mask_bin1 = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) & (hpix[mask_hpix] == j)
            value_jck[i,j+1] = (value_jck[i,0] - (mapa1[mask_bin1]).sum())/(len(mapa1[mask_bin])-len(mapa1[mask_bin1]))
        value_jck[i,0] = value_jck[i,0]/ len(mapa1[mask_bin])
    dict_jck = covariance_jck(value_jck[:,1:], n_jck, 'jackknife')



    if not cov_ext:
        ccov = dict_jck#['cov']
    else:
        vector_cov = np.zeros((len_hist,cov_ext['cov'].shape[1]))
        for i in range(len_hist):
            mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1])# & cov_ext['indexes']
            vector_cov[i,:] = np.mean(cov_ext['cov'][mask_bin,:],axis = 0)
        mask_bin = (sys_map1 > bins[0]) & (sys_map1 < bins[-1])

        ccov = covariance_jck(vector_cov, cov_ext['cov'].shape[1], 'bootstrap')

        
    
    
    # ***********************************************************************************************
    b_arr = np.zeros(n_jck+1)
    if not cov_ext:

        # compute b with non linear
        b_arr[0],xmi = derivative(x,value_jck[:,0])
        for jk in range(n_jck):
            b_arr[1+jk],xmi = derivative(x,value_jck[:,1+jk])
        b_err = covariance_scalar_jck(b_arr[1:],n_jck)['err']
    
    else:
        # compute b with non linear
        b_arr[0],xmi = derivative(x,value_jck[:,0])
        for jk in range(n_jck):
            b_arr[1+jk],xmi = derivative(x,vector_cov[:,jk])
        b_err = covariance_scalar_jck(b_arr[1:],n_jck, type_c = 'bootstrap')['err']
           
    dict_output=dict()
    dict_output.update({'b_arr' : b_arr})
    dict_output.update({'xmi' : xmi})
    dict_output.update({'x' : x})
    dict_output.update({'b_err' : b_err})   
    dict_output.update({'bins' :bins})       


    dict_output.update({'min_s': min_s})
    dict_output.update({'max_s': max_s})
    dict_output.update({'min_ss': min_ss})
    dict_output.update({'max_ss': max_ss})
    dict_output.update({'ave_value': ave})
    dict_output.update({'mask_pixel': mask_pixels})
    dict_output.update({'value_jck': value_jck})
    if cov_ext:
        dict_output.update({'vector_cov': vector_cov})
    dict_output.update({'Area': area_rest/area_tot})
    dict_output.update({'mean shear': np.mean(mapa1[mask_asd])})
        

    return dict_output





def make_plot_final_method1(dict_systematics_tot_1,label,z_minz,z_maxz,fit ='linear'):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(10,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    if fit =='linear':
        fit_f ='F01'
    elif fit =='quadratic':
        fit_f ='F02'
    elif fit =='cubic':
        fit_f ='F03'
        
    hpix_type = 'normal'

    
    xxx = []
    for kk in dict_systematics_tot_1['{0}_{1}'.format(z_minz[0],z_maxz[0])][label]['_w0'].keys():
        xxx.append(kk)
        
    for i in range(len(z_minz)):

        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        #for count in range(20):
        add_label = '_w0'#.format(count)

        if 1==1:
                mute = dict_systematics_tot_1[binx][label][add_label]
                
                ax[i].set_ylabel('bin : {0} \n F-test'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                yy =[]
                for kk in xxx:
                    yy.append(dict_systematics_tot_1[binx][label][add_label][kk][fit_f])
                    
                ax[i].scatter(np.arange(len(xxx)), yy,label = '{0} fit'.format(fit))
                #ax[i].scatter(mute['x'], mute['original'],label = 'linear_fit (no weight)')
                ll = []
        
                '''
                for xx in (xxx):
                    try:
                        ll.append( '{0} {1}'.format(str(xx.split('_')[2]),str(xx.split('_')[3])))
                    except:
                        ll.append( '{0}'.format(str(xx.split('_')[2])))#,str(xx.split('_')[3])))
                '''
                ax[i].plot([-1.,len(xxx)+1],[2,2], linestyle = 'dashed',color = 'black', alpha = 0.7)
                ax[i].plot([-1.,len(xxx)+1],[1,1], linestyle = 'dashed',color = 'black', alpha = 0.7)
                plt.xticks(np.arange(len(xxx)),xxx, rotation='90')
                ax[0].legend()
                plt.suptitle('mean {0} systematic test'.format(label))
                plt.xlim([-1,len(xxx)])
                #print count, mute.keys()
        #except:
        #        pass
        fig.subplots_adjust(wspace=0.,hspace=0.)
        
def make_plot_final(dict_final1, label,z_minz,z_maxz,type_c = 'jackknife',fact = False,Except=False,only=False):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(10,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    bd = dict()
    for i in range(len(z_minz)):
        #i = 4-i
        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        fold = './output_new1_{0}_{1}_{2}/'.format(z_min,z_max,nside)
        #for count in range(20):
        if 1==1:
                if 1==1:
                    ll = []
                    ll1 = []
                    key = binx
                    m1 = 0
                    count = 0
                    for key1 in dict_final1[key].keys():
                        if only: 
                            new_key = only[key1]
                            m1+=len(only[key1])
                        elif Except:
                            new_key = []
                            for ii,jj in enumerate(~np.in1d(dict_final1[key][key1]['_w0'].keys(),Except[key1])):
                                if jj:
                                    new_key.append(dict_final1[key][key1]['_w0'].keys()[ii])
                            m1 += len(new_key)
                        else:
                            m1 += len(dict_final1[key][key1]['_w0'].keys())
                            new_key = dict_final1[key][key1]['_w0'].keys()
                    b1 = np.zeros((m1,n_jck+1))
                    for key1 in dict_final1[key].keys():
                        
                        if only: 
                            new_key = only[key1]
             
                        elif Except:
                            new_key = []
                            for ii,jj in enumerate(~np.in1d(dict_final1[key][key1]['_w0'].keys(),Except[key1])):
                                if jj:
                                    new_key.append(dict_final1[key][key1]['_w0'].keys()[ii])
                        else:
                            new_key = dict_final1[key][key1]['_w0'].keys()
                            
                        for key2 in new_key:
                            b1[count,:]= dict_final1[key][key1]['_w0'][key2]['b_arr']
                            count+=1
                            ll.append('{0}_{1}'.format(key1,key2))


                if fact:
                    b1 = b1
                    b_dict1 = covariance_jck(b1[:,1:]*fact[binx], n_jck, type_cov = type_c)
                else:
                    b_dict1 = covariance_jck(b1[:,1:], n_jck, type_cov = type_c)
                
                
                ax[i].set_ylabel('bin : {0} \n b/b_err'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                
                ax[i].scatter(np.arange(len( b1[:,0]))+1, b1[:,0]/b_dict1['err'],label = 'b/b_err')
                
                
               
                #ax[i].scatter(np.arange(len( b1[:,0])),  p2s(1 - stats.t.cdf((b1[:,0]/b_dict1['err'])**2.,df=1.)) ,label = 'sigma')

                #print b1[:,0]/b_dict1['err']
                    #print ll
                   # #print '{0} {1}'.format(str(xx.split('_')[2]),str(xx.split('_')[3])) ,mute['label'][ii]
                #ax[i].plot([0.,len( b1[:,0])+1],[2,2], linestyle = 'dashed',color = 'black', alpha = 0.7)
                #ax[i].plot([0.,len( b1[:,0])+1],[1,1], linestyle = 'dashed',color = 'black', alpha = 0.7)
                plt.xticks(np.arange(len( b1[:,0]))+1,ll, rotation='90')
                #ax[0].legend()
                #ax[0].suptitle('{0}'.format(label))
                ax[0].text(0.4,1.05 , '{0}'.format(label),transform=ax[0].transAxes)
                
                w = b1[:,0]
                cns = np.zeros(len(b1[:,0]))
                from scipy import linalg
                inv_cov = linalg.inv(b_dict1['cov'])
                N_p = n_jck
                p_p = len(cns)
                f_hartlap = (N_p-1.)/(N_p-p_p-2.)
                #print f_hartlap
                chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/(len(cns)*f_hartlap)

                ax[i].text(1.1 , 0.5 , 'chi/dof : {0:.2f}; sigma : {1:.2f}'.format(chi2red,p2s(1 - stats.t.cdf(chi2red,df=len(cns)))),transform=ax[i].transAxes)
     #plt.text(0, 5.6, ' Average mean shear = {0:.5f}'.format(np.mean(points_2)))
                plt.xlim([0.,len( b1[:,0])+1])
                #print count, mute.keys()

        fig.subplots_adjust(wspace=0.,hspace=0.)
        bd.update({binx:[b_dict1,b1]})
    return bd, ll


def final_plot_Pearson(dict_final1, label,z_minz,z_maxz,type_c = 'jackknife',fact = False):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(10,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    bd = dict()
    for i in range(len(z_minz)):
        #i = 4-i
        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        fold = './output_new1_{0}_{1}_{2}/'.format(z_min,z_max,nside)
        #for count in range(20):
        if 1==1:
                if 1==1:
                    ll = []
                    
                    m1 = 0
                    for key1 in dict_final1[binx].keys():
                        m1 += len(dict_final1[binx][key1].keys())
                    b1 = np.zeros((m1,n_jck+1))

                    count = 0
                    for key1 in dict_final1[binx].keys():
                        for key2 in dict_final1[binx][key1].keys():
                            b1[count,1:]= dict_final1[binx][key1][key2]['Pj']
                            b1[count,0]= dict_final1[binx][key1][key2]['P']
                            count+=1
                            ll.append('{0}_{1}'.format(key1,key2))
                if fact:
                    b1 = b1
                    b_dict1 = covariance_jck(b1[:,1:]*fact[binx], n_jck, type_cov = type_c)
                else:
                    b_dict1 = covariance_jck(b1[:,1:], n_jck, type_cov = type_c)
                
                
                ax[i].set_ylabel('bin : {0} \n P/P_err'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                
                ax[i].scatter(np.arange(len( b1[:,0])), b1[:,0]/b_dict1['err'],label = 'b/b_err')
                import copy 
                points2 = copy.copy(b1[:,0])
                
                xxt = np.arange(len( b1[:,0]))
                plt.xticks(xxt,ll, rotation='90')
                ax[0].text(0.4,1.05 , '{0}'.format(label),transform=ax[0].transAxes)
                
                w = b1[:,0]
                cns = np.zeros(len(b1[:,0]))
                from scipy import linalg
                inv_cov = linalg.inv(b_dict1['cov'])
                N_p = n_jck
                p_p = len(cns)
                f_hartlap = (N_p-1.)/(N_p-p_p-2.)
                #print f_hartlap
                chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/(len(cns)*f_hartlap)

                ax[i].text(1.1 , 0.5 , 'chi/dof : {0:.2f}; sigma : {1:.2f}'.format(chi2red,p2s(1 - stats.t.cdf(chi2red,df=len(cns)))),transform=ax[i].transAxes)
                plt.xlim([0.,len( b1[:,0])+1])
                plt.ylim([0.,3.5])
                
                # plot above lims
                
                mm = b1[:,0]>3.5

                for hh in range(len(mm[mm])):
                    ax.annotate( '{0:.2f}'.format(ll[mm][hh]), ((xxt[mm][hh]-0.5),3.2),transform=ax.transAxes)
            
            
                if len(mm[mm])>0:
                    points2[mm] = 3.4
                    plt.scatter(xxt[mm],points2[mm], color = 'b', marker = '^')

            
        fig.subplots_adjust(wspace=0.,hspace=0.)
        bd.update({binx:b_dict1})
    return bd, ll

import healpy as hp


def convert_to_pix_coord(ra, dec, nside=1024):
    """
    Converts RA,DEC to hpix coordinates
    """

    theta = (90.0 - dec) * np.pi / 180.
    phi = ra * np.pi / 180.
    pix = hp.ang2pix(nside, theta, phi, nest=False)

    return pix
print ('done') 


# compile also this !


import matplotlib.pyplot as plt
import pyfits as pf
import healpy as hp
import copy
import scipy
from scipy import stats
from scipy.optimize import curve_fit
from scipy import linalg

def covariance_scalar_jck(TOTAL_PHI,jk_r, type_c = 'jackknife'):

  #  Covariance estimation
  if type_c == 'jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_c=='bootstrap':
      fact=1./(jk_r)
        
  average=0.
  cov_jck=0.
  err_jck=0.


  for kk in range(jk_r):
    average+=TOTAL_PHI[kk]
  average=average/(jk_r)

  for kk in range(jk_r):
    #cov_jck+=TOTAL_PHI[kk]#*TOTAL_PHI[kk]

    cov_jck+=(-average+TOTAL_PHI[kk])*(-average+TOTAL_PHI[kk])


  err_jck=np.sqrt(cov_jck*fact)


  #average=average*(jk_r)/(jk_r-1)
  return {'cov' : cov_jck*fact,
          'err' : err_jck,
          'mean': average}



def covariance_jck(TOTAL_PHI,jk_r,type_cov):
  if type_cov=='jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_cov=='bootstrap':
      fact=1./(jk_r)
  #  Covariance estimation

  average=np.zeros(TOTAL_PHI.shape[0])
  cov_jck=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  err_jck=np.zeros(TOTAL_PHI.shape[0])

  for kk in range(jk_r):
    average+=TOTAL_PHI[:,kk]
  average=average/(jk_r)

 # print average
  for ii in range(TOTAL_PHI.shape[0]):
     for jj in range(ii+1):
          for kk in range(jk_r):
            cov_jck[jj,ii]+=TOTAL_PHI[ii,kk]*TOTAL_PHI[jj,kk]

          cov_jck[jj,ii]=(-average[ii]*average[jj]*jk_r+cov_jck[jj,ii])*fact
          cov_jck[ii,jj]=cov_jck[jj,ii]

  for ii in range(TOTAL_PHI.shape[0]):
   err_jck[ii]=np.sqrt(cov_jck[ii,ii])
 # print err_jck

  #compute correlation
  corr=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  for i in range(TOTAL_PHI.shape[0]):
      for j in range(TOTAL_PHI.shape[0]):
        corr[i,j]=cov_jck[i,j]/(np.sqrt(cov_jck[i,i]*cov_jck[j,j]))

  average=average*fact
  return {'cov' : cov_jck,
          'err' : err_jck,
          'corr':corr,
          'mean':average}


def compute_area(map, nside):
    map1 = copy.copy(map)
    #map1[map1 > 0] = 1
    area = np.sum(map1) * 1.0 * hp.pixelfunc.nside2pixarea(nside, degrees=True)
    return area

def find_limits(mapp, fract):
    min_ms = min(mapp)
    max_ms = max(mapp)
    his,arr=np.histogram(mapp, bins = np.linspace(min_ms,max_ms,1000))

    tot = np.sum(his)
    low = tot*(fract)/100.
    mute = 0.
    ll =0
    while ll <(len(his)):
        mute += his[ll]
        
        if mute> low:
            min_s = arr[ll]
            ll = (len(his))
        ll += 1    
    ll = 0
    mute = 0.
    while ll <(len(his)):
        mute += his[len(his)-ll-1]

        if mute> low:
            max_s = arr[len(his)-ll]
            ll = (len(his))
        ll += 1    
        
    return min_s,max_s

# linear fit
def fitting_linear(x, para0,para1):
    return para0*(x)+para1

def fitting_2nd(x, para0,para1,para2):
    return para2*(x*x)+para0*(x)+para1
def fitting_3rd(x, para0,para1,para2,para3):
    return para3*(x*x*x)+para2*(x*x)+para0*(x)+para1

from scipy import stats
from scipy import special


def p2s(p):
    return np.sqrt(2.)*special.erfinv((1-p))
def s2p(sigma):
    return 1-special.erf(sigma/np.sqrt(2.))



def show_only_plots(name_folder,add_label,kE,kE_label,weight,systematic_maps,systematics_dict,info,fract_limits, len_hist = 10, cov_ext = False,fit ='linear'):
    
    for key in systematic_maps.keys():
      if key != kE_label:
        
        out_file = name_folder
        make_plot(out_file,add_label,systematics_dict[key],systematic_maps[key]['title'],kE_label,systematic_maps[key]['map'],systematic_maps[key]['fractional'],info,fract_limits, len_hist = len_hist,cov_ext = cov_ext,fit=fit)
        #make_plot(out_file,add_label,dict_output,systematic_maps[key]['title'],kE_label,systematic_maps[key]['map'],systematic_maps[key]['fractional'],info, fract_limits,len_hist = len_hist,cov_ext = cov_ext,fit=fit,dict_output_2=dict_output_2)
        


from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})



def def_max_val(dictu,lim = 2.0, fit = 'linear'):
    F01 = dictu['F01']
       
    if fit == 'linear':
        h = 0.
        m = 0.
    elif fit == 'quadratic':
        h = 1.
        m = 0.
    elif fit == 'cubic':
        h = 1.
        m = 1.
    
    F02 = dictu['F02']*h
    F03 = dictu['F03']*m
    
    F13 = dictu['F13']*m
    F12 = dictu['F12']*h
    F23 = dictu['F23']*m
    

    
    v = 0
    value = 0
    
    if F01>lim:
        value = F01
        v = 1
        if F02>lim and F12>lim:
            value = F02
            v = 2
            if F03>lim and F23>lim:
                value = F03
                v = 3
        else:
            if F03>lim and F13>lim:
                value = F03
                v = 3
            
    else:
        if F02>lim:   
            value = F02
            v = 2    
            if F03>lim and F23>lim:
                value = F03
                v = 3            
        else:
            if F03>lim:
                value = F03
                v = 3            
    return v,value

def correct_weight(count,systematics_dict_total,systematic_maps, info,weights , lim = 2.0, keys = False, fit = 'linear'):
    #order weights!
    from scipy import linalg     
    if 1==1:
        key = '_w{0}'.format(count)
        if keys:
            syst_labb = keys
        else:
            syst_labb = systematics_dict_total[key].keys()

        points = np.zeros(len(syst_labb))
        vv_arr = np.zeros(len(syst_labb))
        label = []
        string = []
     
        for num, key2 in enumerate(syst_labb):
            if key2 != 'mask_tot':
                points[num] = def_max_val( systematics_dict_total[key][key2],lim = lim,fit=fit)[1]
                vv_arr[num] = def_max_val( systematics_dict_total[key][key2],lim = lim,fit=fit)[0]
                label.append(key2)
                vv=''
                if vv_arr[num] == 1:
                    vv = 'linear'
                if vv_arr[num] == 2:
                    vv = 'quadratic'
                if vv_arr[num] == 3:
                    vv = 'cubic'
                
                string.append('{0} ({1:.2f}, {2}), '.format(key2,points[num],vv))
        
        ind = np.argsort(points)

        mask = points[ind] > lim
        weights_add = np.zeros(len(systematic_maps[key2]['map']))
        if len(points[ind][mask])>0:
            

            syst_l = [str(x) for x in np.array(label)[ind][mask]][::-1][0]
            nn = int(vv_arr[ind][mask][::-1][0])
            
            print (' CORRECTION: {0} {1}'.format(syst_l,nn))
            string = ''.join([str(x) for x in np.array(string)[ind][mask]][::-1])
        
            print ("systematics > lim : ",string)#[mask].flatten()
        
            # define weights:
            import copy
            #weights_add = np.zeros(len(systematic_maps[key2]['map']))
        
            
            
            fractional = systematic_maps[syst_l]['fractional']
            x = systematic_maps[syst_l]['map'][info['mask_sims']]
            ave  = x.mean()                      
            if fractional:
                x -= ave
                x /= ave
            else:
                x -= ave
                  
            #mask_sys = (systematic_maps[syst_l]['map'] > systematics_dict_total[key][key2]['min_ss']) & (systematic_maps[syst_l]['map'] < systematics_dict_total[key][key2]['max_ss'])
            
            mask_nn = (systematic_maps[syst_l]['map'][info['mask_sims']] > systematics_dict_total[key][syst_l]['min_ss']) & (systematic_maps[syst_l]['map'][info['mask_sims']] < systematics_dict_total[key][syst_l]['max_ss'])
            
            mutenn = np.zeros(len(systematic_maps[syst_l]['map'][info['mask_sims']]))
            print ("efficiency: ", (np.float(len(x[mask_nn])))/(np.float(len(x))),len(x[mask_nn]))
            #print len(x[mask_nn]), len(weights_add[info['mask_sims']&mask_sys])
            if nn ==3:
                #cubic correction
                print ('cubic')
                dict_out = systematics_dict_total[key][syst_l]#['params_bts3']
                #the corrections should not correct for the shift
                weights_add[info['mask_sims']][mask_nn] = fitting_3rd(x[mask_nn],dict_out['params_bts3'][0],dict_out['params_bts3'][1],dict_out['params_bts3'][2],dict_out['params_bts3'][3])  

            if nn ==2:
                #cubic correction
                print ('quadratic')
                dict_out = systematics_dict_total[key][syst_l]#['params_bts3']

                nw = fitting_2nd(x[mask_nn],dict_out['params_bts2'][0],dict_out['params_bts2'][1],dict_out['params_bts2'][2])  
                mutenn[mask_nn] = nw
                weights_add[info['mask_sims']] = mutenn
             
            if nn ==1:
                #cubic correction
                print ('linear')
                dict_out = systematics_dict_total[key][syst_l]#['params_bts3']
                #the corrections should not correct for the shift
                #print max(fitting_linear(x[mask_nn],dict_out['params_bts'][0],dict_out['params_bts'][1])  ), min(fitting_linear(x[mask_nn],dict_out['params_bts'][0],dict_out['params_bts'][1])  )
                nw= fitting_linear(x[mask_nn],dict_out['params_bts'][0],dict_out['params_bts'][1])  
                mutenn[mask_nn] = nw
                weights_add[info['mask_sims']] = mutenn
                print (dict_out['params_bts'][0],dict_out['params_bts'][1],len(weights_add[info['mask_sims']]),len(weights_add[info['mask_sims']][mask_nn]), max(weights_add[info['mask_sims']][mask_nn]))
            mute_l = False
        else:
            mute_l = True
        return weights_add,mute_l

    
def compute_fisher(chi_n,chi_fit,N,p_null,p_fit):
    #F = -((N-p_fit)/(p_null-p_fit))*(chi_n/chi_fit-1.),
    F = ((chi_n-chi_fit)/(p_fit-p_null))/((chi_fit)/(N-p_fit))
    v1 = -(p_null-p_fit)
    v2 = (N-p_fit)
    return F, v1,v2


def make_plot(save_output_folder,add_label,dict_out, title, kE_label, sys_map1,fractional,info,fract_limits,len_hist=10,cov_ext = False,fit='linear', dict_output_2 = False):
    
    print ("Jackknife: {0}".format(not cov_ext))
    plt.figure(1,figsize=(6,5))
    ax = plt.subplot(211)
    if fractional:
        plt.xlabel('{0} [{1:.2f},{2:.2f}] (fractional)'.format(title, dict_out['min_s'],dict_out['max_s']))
    else:
        plt.xlabel('{0}[{1:.2f},{2:.2f}] (mean sub)'.format(title, dict_out['min_s'],dict_out['max_s']))
    plt.title('systematic map : {0}'.format(title))
    plt.ylabel('counts')
    
    plt.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        
    sys_map1 = sys_map1[info['mask_sims']]
    ave = sys_map1.mean()
    if fractional:

        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
    plt.hist(sys_map1, bins = dict_out['bins'])
    
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width* 1.4 , box.height])   
    
    plt.figure(1)
    ax = plt.subplot(212)
    plt.title('systematic map : {0}'.format(title))
    if fractional:
        plt.xlabel('{0} (mean subtracted, fractional)'.format(title))
    else:
        plt.xlabel('{0} (mean subtracted)'.format(title))        
    plt.ylabel('<{0}>'.format(kE_label))
    
    bins = dict_out['bins']
    x = dict_out['x']
    
    plt.errorbar(x, dict_out['w'],dict_out['dict_cov']['err'])
    plt.tight_layout()
    
    w = dict_out['w']
    cns = np.ones(len(w))* dict_out['mean shear']
    inv_cov = linalg.inv(dict_out['dict_cov']['cov'])
    chi2red = (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/len(w)
    

    #plt.plot(x, fitting_linear(x,dict_out['params'][0],dict_out['params'][1]),label='fit')
              
    if cov_ext:         
        label_bst ="_boot"
        inv_cov = linalg.inv(dict_out['dict_cov_bts']['cov'])
    else:
        label_bst =""
        inv_cov = linalg.inv(dict_out['dict_cov']['cov'])
    if 1==1:
        
        w = dict_out['w']
        cns = np.ones(len(w))* dict_out['mean shear']
        
        chi2red = (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/len(w)
    
    
            
        scipy.stats.f.cdf
        xx = 1.05
        yy = 1.10

        chi = dict_out['chi2red'+label_bst+'_null'] 
        df1 = len(dict_out['x']) - 1
        chi=chi/df1
        p = 1 - stats.t.cdf(chi,df=df1)
        #print chi,df,p
        #print p2s(p)
        plt.text(xx, yy+1.35, 'summary stat', transform=ax.transAxes)
        plt.text(xx, yy+1.2, '[method 1]', transform=ax.transAxes)

        plt.text(xx, yy+0.95, 'sigma_null: {0:.3f}'.format(p2s(p) ), transform=ax.transAxes)

        chi = dict_out['chi2red'+label_bst] 
        df2 = len(dict_out['x']) - 2
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        #print chi,df,p
        #print p2s(p)
        
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'_null'] ,dict_out['chi2red'+label_bst] ,len(dict_out['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)


        plt.text(xx, yy+0.8, 's_1: {0:.2f}'.format(p2s(p),p2s(pf)), transform=ax.transAxes)
        plt.text(xx, yy+0.65, 'F: [{1:.2f}]'.format(p2s(p),p2s(pf)), transform=ax.transAxes)
  


        chi = dict_out['chi2red'+label_bst+'2'] 
        df2 = len(dict_out['x']) - 3
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        
        
        
        #print chi,df,p
        #print p2s(p)
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'_null'] ,dict_out['chi2red'+label_bst+'2'] ,len(dict_out['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)

        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst] ,dict_out['chi2red'+label_bst+'2'] ,len(dict_out['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        plt.text(xx, yy+0.5, 's_2: {0:.3f}'.format(p2s(p),p2s(pf),p2s(pf2)), transform=ax.transAxes)
        plt.text(xx, yy+0.35, 'F: [{1:.2f},{2:.2f}]'.format(p2s(p),p2s(pf),p2s(pf2)), transform=ax.transAxes)
   
        
        chi = dict_out['chi2red'+label_bst+'3'] 
        df2 = len(dict_out['x']) - 4
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        #print chi,df,p
        #print p2s(p)
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'_null'] ,dict_out['chi2red'+label_bst+'3'] ,len(dict_out['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)


        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst] ,dict_out['chi2red'+label_bst+'3'] ,len(dict_out['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
 
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'2'] ,dict_out['chi2red'+label_bst+'3'] ,len(dict_out['x']),1,2)  
        pf3 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        plt.text(xx, yy+0.2, 's_3: {0:.3f}'.format(p2s(p),p2s(pf),p2s(pf2),p2s(pf3)), transform=ax.transAxes)
        plt.text(xx, yy+0.05, 'F: [{1:.2f},{2:.2f},{3:.2f}]'.format(p2s(p),p2s(pf),p2s(pf2),p2s(pf3)), transform=ax.transAxes)
                
                 
        chi = dict_out['chi2red'+label_bst+'_null']/df1 - dict_out['chi2red'+label_bst]/df2
        df = 1
        p = 1 - stats.t.cdf(chi,df=df)
        
        '''
        plt.text(xx, yy+0.35, 'sigma_H1 (fit): {0:.3f}'.format(p2s(p)), transform=ax.transAxes)

        plt.text(xx, yy+0.2, 'sigma_68 : {0:.3f}'.format(dict_out['chi68_new']), transform=ax.transAxes)
        plt.text(xx, yy+0.05, 'sigma_68_ratio : {0:.3f}'.format(p2s(p)/dict_out['chi68_new']), color='r',transform=ax.transAxes)
        plt.text(xx, yy-0.1, 'sigma_68 (nd) : {0:.3f}'.format(dict_out['chi68_new1']),transform=ax.transAxes)

        chi = dict_out['chi2red_boot_null'] - dict_out['chi2red_boot']
        df = 1
        p = 1 - stats.t.cdf(chi,df=df)
        plt.text(xx, yy-0.25, 'sigma_68_ratio (nd: {0:.3f}'.format(p2s(p)/dict_out['chi68_new1']), color='r',transform=ax.transAxes)


        '''
        plt.text(xx, 0.65, 'chi2 null: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst+'_null'],len(dict_out['x']) - 1), transform=ax.transAxes)
        plt.text(xx, 0.5, 'chi2 fi1t: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst],len(dict_out['x']) - 2), transform=ax.transAxes)
        plt.text(xx, 0.35, 'chi2 fit2: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst+'2'],len(dict_out['x']) - 2), transform=ax.transAxes)
        plt.text(xx, 0.2, 'chi2 fit3: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst+'3'],len(dict_out['x']) - 2), transform=ax.transAxes)
      
        '''
        plt.text(xx, 0.35, 'diff chi2: {0:.4f}'.format(dict_out['chi2red_diff_boot'],dict_out['chi68']), transform=ax.transAxes)
        #plt.text(xx, 0.2, 'diff chi2/dof 68:  {0:.4f}'.format(dict_out['chi68']), transform=ax.transAxes)
        #plt.text(xx, 0.05, 'diff chi2/ diff chi2 68: {0:.4f}'.format(dict_out['chi2red_diff_boot68']),color='r', transform=ax.transAxes)
        '''
        plt.text(0.05, 0.05, 'slope [method_1]= {0:.4f} +- {1:.4f}'.format(dict_out['params_bts'][0],np.sqrt(dict_out['params_cov_bts'][0,0])), transform=ax.transAxes)
        
        if dict_output_2:
            plt.text(0.05, 0.2, 'slope [method_2]= {0:.4f} +- {1:.4f}'.format(dict_output_2['b_arr'][0],dict_output_2['b_err']), transform=ax.transAxes)
            
        plt.plot(x, fitting_2nd(x,dict_out['params_bts2'][0],dict_out['params_bts2'][1],dict_out['params_bts3'][2]),label='fit 2nd')
        plt.plot(x, fitting_3rd(x,dict_out['params_bts3'][0],dict_out['params_bts3'][1],dict_out['params_bts3'][2],dict_out['params_bts3'][3]),label='fit 3rd')

        plt.plot(x, fitting_linear(x,dict_out['params_bts'][0],dict_out['params_bts'][1]),label='fit')
        plt.plot(x, cns,label='null')
        



        
        
        
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width* 1.4, box.height])




    plt.ylim([-0.004,0.004])
    plt.savefig(save_output_folder+kE_label+'_'+title+' '+str(fract_limits)+'_'+str(len_hist)+'_'+add_label+'_'+fit+'.png')   

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95),
          ncol=4, fancybox=True)
    #plt.legend()
    plt.show()


def plot_resume(systematics_dict_total,z_min,z_max,nside,kE_label,d68 = False,type_c='jackknife') :
 from scipy import linalg     
 for key in systematics_dict_total.keys():
    fig = plt.figure(figsize=(8, 10)) 
    
    
    
    if type_c == 'jackknife':
        d68 = False
        
    syst_labb = systematics_dict_total[key].keys()
    if d68:
        
        ax = plt.subplot(313)
        
        plt.title('sys map impact wrt null hypotesys - {0}   bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    
        points = np.zeros(len(syst_labb))
        points_2 = np.zeros(len(syst_labb))  
    
        points_err = np.zeros(len(syst_labb))
        x = np.array(range(len(syst_labb)))+1
        label = []

        plt.show
        for num, key2 in enumerate(syst_labb):
          if key2 != 'mask_tot':
        
            points_2[num] = systematics_dict_total[key][key2]['chi2red_diff_boot68']
            key2 = '_'+syst_labb[num]
            key2 =  '   '+key+key2
            label.append(key2)
    
        plt.xticks(x, label, rotation='90')
        plt.xlim([0.,len(syst_labb)+1])
        plt.ylim([0,6.])
        mm = points_2>6

        for hh in range(len(mm[mm])):
            ax.annotate( '{0:.2f}'.format(points_2[mm][hh]), ((x[mm][hh]-0.5),5.4),transform=ax.transAxes)
            
            #    plt.text(x[mm[hh]], points_2[mm[hh]], 'cc')#,transform=ax.transAxes)
        plt.scatter(x,points_2,label = '<e>_0 = <e>', color = 'b')
        if len(mm[mm])>0:
            points_2[mm] = 5.9
            plt.scatter(x[mm],points_2[mm], label = '<e>_0 = 0.', color = 'b', marker = '^')

            
        plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
        plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
        plt.ylabel('delta chi2/delta chi2 (68)')
    

    
    
        
        
        
        
    plt.subplot(312)#,figsize=(10,20))
    #plt.figure(1,figsize=(10,100))
    

    points = np.zeros(len(syst_labb))
    points_2 = np.zeros(len(syst_labb))
    x = np.array(range(len(syst_labb)))+1
    label = []


  
    for num, key2 in enumerate(syst_labb):
     if key2 != 'mask_tot':
        points[num] = np.sqrt((systematics_dict_total[key][key2]['params_bts'][0]/np.sqrt(systematics_dict_total[key][key2]['params_cov'][0,0]))**2.)
        points_2[num] = systematics_dict_total[key][key2]['mean shear'] 
        if num == 0:
            mask_tot = systematics_dict_total[key][key2]['mask_pixel']
        else:
            mask_tot = mask_tot ^ systematics_dict_total[key][key2]['mask_pixel']
        key2 = '_'+syst_labb[num]
        key2 = '   '+ key+key2
        label.append(key2)
    

    plt.xticks(x, label,  rotation='90')
    plt.xlim([0.,len(syst_labb)+1])
    plt.ylim([0.,8])
    plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
    plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
    plt.title('Linear fit slope: {0}   bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    plt.scatter(x,points)
    plt.ylabel('|b|/b_err')
    #fig.tight_layout()
    p_value = 1 - stats.chi2.cdf(np.sum(points**2), len(points))
    #plt.text(0, 7, '  chi2/dof = {0:.2f}'.format(np.sum(points)/len(points)))
    #plt.text(0, 5.6, ' Average mean shear = {0:.5f}'.format(np.mean(points_2)))
    #plt.text(0, 4, '  p-value = {0:.2f}'.format(p_value))
    tot_mask = compute_area(mask_tot,nside)
    #plt.text(10, 7, 'area masked = {0:.1f} deg2'.format(tot_mask))
    #systematics_dict_total[key].update({'mask_tot' : tot_mask})
    #plt.tight_layout()
    
    

    ax = plt.subplot(311)
    #plt.figure(1,figsize=(10,20))
    plt.title('Chi2 null hypothesis {0} - bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    
    points = np.zeros(len(syst_labb))
    points_2 = np.zeros(len(syst_labb))  
    
    points_err = np.zeros(len(syst_labb))
    x = np.array(range(len(syst_labb)))+1
    label = []

    plt.show
    for num, key2 in enumerate(syst_labb):
      if key2 != 'mask_tot':
        
        w = systematics_dict_total[key][key2]['w']
        cns = np.ones(len(w))* systematics_dict_total[key][key2]['mean shear']
        inv_cov = linalg.inv(systematics_dict_total[key][key2]['dict_cov']['cov'])
        points[num] = (np.matmul(w,np.matmul(inv_cov,w)))/len(w)
        if type_c =='jackknife':
            points_2[num] = systematics_dict_total[key][key2]['chi2red']/(len(systematics_dict_total[key][key2]['x']) - 2)
        else:
            points_2[num] = systematics_dict_total[key][key2]['chi2red_boot']/(len(systematics_dict_total[key][key2]['x']) - 2)
        key2 = '_'+syst_labb[num]
        key2 =  '   '+key+key2
        label.append(key2)
    #print points_2
    
    plt.xticks(x, label, rotation='90')
    plt.xlim([0.,len(syst_labb)+1])
    plt.ylim([0,4.])
    plt.scatter(x,points_2,label = '<e>_0 = <e>')
    plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
    plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
    plt.ylabel('chi2/dof')


    plt.show()
        
        
print ('done')


import skymapper as skm
import matplotlib.cm as cm

def getHealpixCoords(pixels, nside, nest=False):
    # convert healpix cell indices to center ra/dec
    import healpy as hp
    theta, phi = hp.pix2ang(nside, pixels, nest=nest)
    return phi * 180. / np.pi, 90 - theta * 180. / np.pi
    

def skm_plot(mapa, mask, sep=15, ra = None, dec = None, richness = None, title = None, add_rmp = True, small_scale = False, cb_label = 'map value', x_size = 6, y_size = 6, vmin = -0.01, vmax = 0.01):

    cmap = cm.RdYlBu_r
    fig = plt.figure(figsize=(x_size,y_size))
    ax = fig.add_subplot(111, aspect='equal')

    
    
    reticule = sep
    nside = hp.pixelfunc.npix2nside(len(mapa))
    pixels = np.array(range(len(mapa)))
    pixels = pixels[~mask]
    ra_, dec_ = getHealpixCoords(pixels, nside)
    
    #plt.scatter(ra_,dec_)
    ##ig = plt.figure(figsize=(15,10))
    #ax = fig.add_subplot(111, aspect='equal')
    #fig.tight_layout()
    #proj = skm.plotMap(ra_, dec_, mapa[pixels], sep=sep, ax=ax, cb_label= cb_label, cmap='bwr')

    # setup map: define AEA map optimal for given RA/Dec
    proj = skm.createConicMap(ax, ra_, dec_, proj_class=skm.AlbersEqualAreaProjection)

    # add lines and labels for meridians/parallels (separation 5 deg)
    
    parallels = np.arange(0. ,360., reticule)
    meridians = np.arange(-90., 90., reticule)
    skm.setMeridianPatches(ax, proj, meridians, linestyle='-', lw=0.5, alpha=0.3, zorder=10)
    skm.setParallelPatches(ax, proj, parallels, linestyle='-', lw=0.5, alpha=0.3, zorder=10)
    skm.setParallelLabels(ax, proj, parallels, loc="bottom", fmt=skm.pmDegFormatter)
    skm.setMeridianLabels(ax, proj, meridians, loc="left", fmt=skm.pmDegFormatter)
    


    # convert to map coordinates and plot a marker for each point
    x,y = proj(ra_, dec_)
    marker = 's'
    markersize = skm.getMarkerSizeToFill(fig, ax, x, y)

    sc = ax.scatter(x,y, c=mapa[pixels], edgecolors='None', marker=marker, s=markersize, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True, zorder=1)
    
    
    # overplot with another data set
    # here clusters [not implemented]
    if add_rmp:
        x,y = proj(ra, dec)
        ax.scatter(x,y, c='None', edgecolors='k', linewidths=1, s=richness, marker='o', zorder=3)

    if title:
        plt.title(title)

    # add colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes("right", size="5%", pad=0.3)
    cb = fig.colorbar(sc, cax=cax )

    if cb_label != None:
        cb.set_label(cb_label)
    ticks = np.linspace(vmin, vmax, 5)
    cb.set_ticks(ticks)
    if small_scale:
        cb.set_ticklabels([('%.3f' % (t))   for t in ticks])       
    else:
        cb.set_ticklabels([('%.1f' % (t))   for t in ticks])
    cb.solids.set_edgecolor("face")

    # show (and save) ...
    fig.tight_layout()
    fig.show()
    #fig.savefig(imagefile)    
       
        
from scipy.spatial import distance
import os
import numpy as np
import pickle

def update_progress(progress,elapsed_time=0,starting_time=0):

    import time
    import timeit
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))



def dist_cent_2(ra1,dec1,ra2,dec2):

            todeg = np.pi/180.
            ra1 = ra1*todeg
            ra2 = ra2*todeg
            dec1 = dec1*todeg
            dec2 = dec2*todeg

            cos = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
            return np.arccos(cos)/todeg



def save_obj( name,obj ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)






def estimator(w_estimator,DD,DR,RD,RR):

    if w_estimator == 'LS':
        #print (w_estimator)
        results = (DD-DR-RD+RR)/(RR)
    elif w_estimator == 'Natural':
        results = (DD - RR) / (RR)
    elif w_estimator == 'Hamilton':
        results = (DD * RR - RD * DR) / (RD * DR)
    elif w_estimator == 'Natural_noRu':
        results = (DD - DR) / DR
    elif w_estimator == 'Natural_noRr':
        results = (DD - RD) / RD
    return results


import time
import timeit
import treecorr
class Jack(object):
    def __init__(self,w_estimator,conf,pairs_sel,rau,decu,rar,decr,jku,jkr,weight_u,weight_r, dist_path,centers=None, njk=None,verbose=False):
        
        self.jackknife_speedup= True
        self.w_estimator=w_estimator
        self.conf = conf
        
        self.rau=rau
        self.decu=decu

        self.rar=rar
        self.decr=decr

        
        self.jku=jku
        self.jkr=jkr


     
        self.weight_u=weight_u
        self.weight_r=weight_r

        

        self.centers = centers
        self.njk = njk

        self.dist_path = dist_path
        self.pairs_select=pairs_sel

    
    def KKCorrelation(self):
        self.prolog()
        pairs= self.epilog()

        return pairs


    def convert_units(self):
        """Change to 'degrees' the units if necessary.
        """
        if 'sep_units' in self.conf.keys():
            un = self.conf['sep_units']
            if un == 'arcmin':
                todeg = 1./60.
            elif un == 'arcsec':
                todeg = 1./60.**2.
            else:
                todeg = 1.
        else:
            todeg = 1.
        return todeg


    def max_distance_region(self):
        cnt = self.centers
        max_dist_region=load_obj(self.dist_path)
        self.max_dist_region=max_dist_region



    def distance(self):
        """Finds the minimum distance to a center for each center.
           Fixes double of this distance as the criteria for not considering correlations,
           which is a conservative choice. This distance has to be at least 4 times the
           maximum angular separation considered. Centers beyond this distance will not be
           considered in the correlations.
        """

        # Find the minimum distance to a center for each center.

        cnt = self.centers
        dist = np.array([np.sort([self.dist_cent(cnt[i],cnt[j]) for i in range(len(cnt))])[1] for j in range(len(cnt))])
        dist = (dist)*2.

        todeg = self.convert_units()

        if 'max_sep' in self.conf.keys():
            max_sep = self.conf['max_sep'] * todeg
        else:
            raise NotImplementedError("Make use of 'max_sep' in configuration.")

        # Check that the distance is at least 4 times the maximum angular separation.
        self.center_min_dis = np.array( [ 4.*max_sep if x < 4.*max_sep else x for x in dist] )



    def dist_cent(self, a, b):
        """Angular distance between two centers (units: degrees). Makes use of spherical law of cosines.
        """
        todeg = np.pi/180.
        a = a*todeg
        b = b*todeg
        cos = np.sin(a[1])*np.sin(b[1]) + np.cos(a[1])*np.cos(b[1])*np.cos(a[0]-b[0])
        return np.arccos(cos)/todeg


    def cond(self, i, j):

        """Return the maximum conditional distance for a pair of centers to
           determine if the correlation should be computed.
        """

        con = self.center_min_dis
        return max(con[i], con[j])


    def collect(self, pairs):


        shape = (self.njk, self.conf['nbins'])
        
        DD_a,DD_c = np.zeros(shape),np.zeros(shape)
        normm_a,normm_c = np.zeros(shape),np.zeros(shape)
        DD = np.zeros(self.conf['nbins'])
        normm = np.zeros( self.conf['nbins']) 

        if self.jackknife_speedup:
            for n in range(self.conf['nbins']):
                DD[n] = np.sum(pairs[:,0,0,n]) + 0.5 * np.sum(pairs[:,1,0,n])
                normm[n] = np.sum(pairs[:,0,1,n]) + 0.5 * np.sum(pairs[:,1,1,n])    

                # only the jk *****************
                DD_a[:,n] = (pairs[:,0,0,n])
                DD_c[:,n]  =  0.5 * (pairs[:,1,0,n])

                normm_a[:,n] = (pairs[:,0,1,n])
                normm_c[:,n]  =  0.5 * (pairs[:,1,1,n])
        return DD,DD_a, DD_c, normm, normm_a, normm_c


    def parallel(self, i, j):

        [[ra_a, dec_a, jk_a], [ra_b, dec_b, jk_b]] = self.info

        # Create the Catalog object. One for each jackknife region.
        try:
            mask=np.in1d(jk_a,i)
            cat_a = treecorr.Catalog(ra=ra_a[mask], dec=dec_a[mask], ra_units='deg', dec_units='deg',k=self.weight_u[mask])
        except RuntimeError:
            cat_a = None
        try:
            mask=np.in1d(jk_b,j)
            cat_b = treecorr.Catalog(ra=ra_b[mask], dec=dec_b[mask], ra_units='deg', dec_units='deg',k=self.weight_r[mask])
        except RuntimeError:
            cat_b = None

        kk = treecorr.KKCorrelation(self.conf)
        kk.process(cat_a, cat_b)

        return [kk.xi*kk.weight,kk.weight]

    
    def dist_cent_2(self,ra1,dec1,ra2,dec2):

            todeg = np.pi/180.
            ra1 = ra1*todeg
            ra2 = ra2*todeg
            dec1 = dec1*todeg
            dec2 = dec2*todeg

            cos = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
            return np.arccos(cos)/todeg

    def prolog(self):

        #insert njk


        njk = self.njk

        cnt = self.centers

 
        self.distance()
        self.max_distance_region()

        ra_a, dec_a,  jk_a, NA = self.rau, self.decu,self.jku,np.sum(self.weight_u)
        ra_b, dec_b,  jk_b, NB = self.rar, self.decr,self.jkr,np.sum(self.weight_r)
 

        self.info = [[ra_a, dec_a, jk_a],[ra_b, dec_b, jk_b]]

        if 'nbins' in self.conf:
                shape = (2,self.conf['nbins'])

        else:
            raise Exception('Not implemented yet. Please use nbins in config.')
            #shape = (4, self.conf['bin_size'])

        cnt = self.centers

        t1=time.time()

        pairs = [[np.zeros(shape) for i in range(njk)] for j in range(njk)]
   
        pairs_ring = [[np.zeros(shape) for i in range(2)] for j in range(njk)]
        a=np.concatenate((np.array([[(i,j) for i in range(njk)] for j in range(njk)])))

        todeg = self.convert_units()
        max_sep = self.conf['max_sep'] * todeg

        sel = np.array([ max([0.,(self.dist_cent(cnt[i],cnt[j]) - (self.max_dist_region[i]+self.max_dist_region[j]))]) < (3. * max_sep )for (i,j) in a])
        b = a[sel]


        def fun_speedup(othersample,otehrsample1,jackk):
            try:
                pairsCC1=self.parallel(othersample,[jackk])
                pairsCC2=self.parallel([jackk],othersample1)
                pairs_auto=self.parallel([jackk],[jackk])
                for prs in range(2):
                    pairsCC1[prs]+=pairsCC2[prs]
                pairs_ring[jackk][1] = pairsCC1
                pairs_ring[jackk][0] = pairs_auto
                
            except (RuntimeError, e):
                print (e)
                pairs_ring[jackk][0] = np.zeros(shape)
                pairs_ring[jackk][1] = np.zeros(shape)

        start=timeit.default_timer()


        mute=0

        if self.jackknife_speedup:
            startt = timeit.default_timer()
            for counter,jackk in enumerate(np.unique(b[:,1])):
                mask=(b[:,1]==jackk) & (b[:,0]!=jackk)
                othersample=b[mask,0]
                mask=(b[:,0]==jackk) & (b[:,1]!=jackk)
                othersample1=b[mask,1]
                fun_speedup(othersample,othersample1,jackk)
                mute+=len(othersample)+1
                
                #update_progress((float(counter)/len(np.unique(b[:,1]))),timeit.default_timer(),startt)
                if counter == 100:
                    print (counter, timeit.default_timer()-startt)
                if counter == 200:
                    print (counter, timeit.default_timer()-startt)
                if counter == 300:
                    print (counter, timeit.default_timer()-startt)
                if counter == 400:
                    print (counter, timeit.default_timer()-startt)
                if counter == 500:
                    print (counter, timeit.default_timer()-startt)
                if counter == 600:
                    print (counter, timeit.default_timer()-startt)
                if counter == 700:
                    print (counter, timeit.default_timer()-startt)
                if counter == 800:
                    print (counter, timeit.default_timer()-startt)
            self.pairs = np.array(pairs_ring)




    def epilog(self):

        pairs = self.pairs

        DD,DD_a,DD_c, normm, normm_a, normm_c = self.collect(pairs)
        convert=self.convert_units()
        min_sep, max_sep, nedges = self.conf['min_sep']*convert, self.conf['max_sep']*convert, self.conf['nbins']+1
        th = np.linspace(np.log10(min_sep), np.log10(max_sep), nedges)
        theta = 10**np.array([(th[i]+th[i+1])/2 for i in range(len(th)-1)])



        #print (corr.shape,DD.shape,DR_a.shape)#, DD,DR ,RD,RR, DD_a,DR_a ,RD_a,RR_a, DD_c,DR_c ,RD_c,RR_c)

        # ***********************************

        pairs=dict()
        
        pairs.update({'theta':theta})
        pairs.update({'DD':DD})
        pairs.update({'DD_a':DD_a})
        pairs.update({'DD_c':DD_c})
        pairs.update({'normm':normm}) 
        pairs.update({'normm_a':normm_a})
        pairs.update({'normm_c':normm_c})
        
        return pairs

    

# compute maximum distance within jackknife regions *****
def distance_calc(dist_path,ra_m,dec_m,jk_m, njk,centers):

    '''
    based on the maximum angular scale probed by the cross-correlation,
    computes which jackknife regions should be included in the computation
    '''
    
    import timeit
    if not os.path.exists(dist_path+'.pkl'):

        max_dist_region1=np.zeros(njk)

        # convert radec to xyz
        cosdec = np.cos(dec_m)
        aJx_u = cosdec * np.cos(ra_m)
        aJy_u = cosdec * np.sin(ra_m)
        aJz_u = np.sin(dec_m)

        print ('compute maximum distance for each jackknife region')
        start=timeit.default_timer()
        for i in range(njk):
            if len(ra_m[jk_m==i]) ==0 or len(dec_m[jk_m==i])==0:
                max_dist_region1[i,j]=0.
            else:

                ra_c,dec_c=centers[i]

                cosdec = np.cos(dec_c)
                aJx_r = cosdec * np.cos(ra_c)
                aJy_r = cosdec * np.sin(ra_c)
                aJz_r = np.sin(dec_c)
        
                tree_m=spatial.cKDTree(np.c_[aJx_u[jk_m==i], aJy_u[jk_m==i], aJz_u[jk_m==i]])

                max_dist_m,index_dist=tree_m.query([aJx_r,aJy_r,aJz_r],k=len(ra_m[jk_m==i]))

                ra_new=ra_m[jk_m==i]
                dec_new=dec_m[jk_m==i]


                if (len(ra_m[jk_m==i])==1):
                    max_dist_region1[i]=dist_cent_2(ra_c,dec_c,ra_new[index_dist],dec_new[index_dist])
                else:
                    max_dist_region1[i]=dist_cent_2(ra_c,dec_c,ra_new[index_dist[-1]],dec_new[index_dist[-1]])
            update_progress(np.float(i+1)/np.float(njk),timeit.default_timer(),start)
        
        save_obj(dist_path,max_dist_region1)
        


def make_wz_errors(pairs,resampling,number_of_bootstrap,pairs_resampling):
    '''
    take as an input pairs between subregions and creates resmapling and errors.
    '''

    # errors ********************************
    if resampling=='jackknife':
        njk=pairs['DD_a'].shape[0]
    elif resampling=='bootstrap':
        njk=number_of_bootstrap

    DD_j=np.zeros((njk+1,pairs['DD_a'].shape[1]))
    norm_j=np.zeros((njk+1,pairs['DD_a'].shape[1]))
    DD_j[0,:]=pairs['DD']/pairs['normm']
    norm_j[0,:]=pairs['normm']    
         
    for jk in range(njk):
        if resampling=='jackknife':
            
            if pairs_resampling:
                fact=1.
            else:
                fact=2.

            DD_j[jk+1,:]=(pairs['DD'][:]-pairs['DD_a'][jk,:]-fact*pairs['DD_c'][jk,:])
            norm_j[jk+1,:]=(pairs['normm'][:]-pairs['normm_a'][jk,:]-fact*pairs['normm_c'][jk,:])
            

            DD_j[jk+1,:] = DD_j[jk+1,:]/norm_j[jk+1,:]

    return pairs['theta'],DD_j.T,njk




# Initialize ******

def make_plot_corr(title,save_output_folder,kE_label,mute_d):
    plt.figure(figsize=(6,3))
    ax = plt.subplot(111)

    plt.title('systematic map (crosscorr) : {0}'.format(title))
    plt.xlabel('theta (arcmin)'.format(title))        
    plt.ylabel('w(theta)')
    plt.errorbar(mute_d['theta']*60.,mute_d['w'],mute_d['cov']['err'])
    plt.tight_layout()
    plt.text(0.7, 0.1, 'chi2/dof:  {0:.2f}'.format(mute_d['chi2_red']), transform=ax.transAxes, size =12)
    plt.show()     
    plt.savefig(save_output_folder + 'cc_' +kE_label+'_'+ title +'.png')
    

def do_analisys_corr(kE,kE_label,weight,systematic_maps,info,add_label):
    
    systematics_dict_corr = dict()
    for key in systematic_maps.keys():
      if key != kE_label:
        out_file = './output_{0}_{1}/'.format(info['z_min'],info['z_max'])
        if not os.path.exists('./output_{0}_{1}/pairs/'.format(info['z_min'],info['z_max'])):
            os.mkdir('./output_{0}_{1}/pairs/'.format(info['z_min'],info['z_max']))
            
        mute_d = compute_cross_corr(systematic_maps[key]['title'],out_file,kE_label,systematic_maps[key]['map'],kE, info, add_label,  mapa_weight = weight,fractional = systematic_maps[key]['fractional'])
        make_plot_corr(systematic_maps[key]['title'],out_file,kE_label,mute_d)
        systematics_dict_corr.update({systematic_maps[key]['title']  : mute_d}) 
  
    return systematics_dict_corr

def plot_resume_cc(systematics_dict_total,z_min,z_max,kE_label) :
 from scipy import linalg     
 for key in systematics_dict_total.keys():
    plt.figure(1)
    
    plt.subplot(111)
    
    syst_labb = systematics_dict_total[key].keys()


    points = np.zeros(len(syst_labb))
    points_2 = np.zeros(len(syst_labb))
    x = np.array(range(len(syst_labb)))+1
    label = []


  
    for num, key2 in enumerate(syst_labb):
     if key2 != 'mask_tot':
        points[num] = (systematics_dict_total[key][key2]['chi2_red'])
        key2 = '_'+syst_labb[num]
        key2 = '   '+ key+key2
        label.append(key2)
    

    plt.xticks(x, label, rotation='20')
    plt.xlim([0.,len(syst_labb)+1])
    plt.ylim([0.,8])
    plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
    plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
    plt.title('Chi2 null hypothesis (constant shear) : {0}   bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    plt.scatter(x,points)
    plt.ylabel('Chi2/dof')

    #systematics_dict_total[key].update({'mask_tot' : tot_mask})
    plt.tight_layout()
    
    print (kE_label)
    plt.savefig('./output_{1}_{2}/{0}_{3}_total_systematics_cc_{4}.png'.format(key,z_min,z_max, kE_label,fract_limits))

    plt.show()
    
def compute_cross_corr(title,save_output_folder,kE_label,sys_map,mapa, info,  add_label,mapa_weight = None, fractional = False):

    print ('\n****************\n')
    print ('TEST (corr) {1}: -> {0}'.format(title,kE_label))

    mask = info['mask_sims']
    sys_map1 = sys_map[mask]
    mapa1 = mapa[mask]-mapa_weight[mask]

    ave = sys_map1.mean()
    std = sys_map1.std()
    if fractional:
    
        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
        

    if not os.path.exists(save_output_folder+'pairs/pairs_{0}_{1:.3f}_{2:.3f}_{3}_{4}_{5}_{6}_{7}.pkl'.format(info['n_jck'],info['conf']['min_sep'],info['conf']['max_sep'],info['conf']['nbins'],title,kE_label,add_label,fract_limits)):
   
        J = Jack(info['w_estimator'],info['conf'],info['pairs_to_compute'], info['ra'],info['dec'],info['ra'],info['dec'],info['hpix'],info['hpix'],
             mapa1,sys_map1,info['dist_path'],centers=info['centers'],njk=info['n_jck'])
        pairs = J.KKCorrelation()
        save_obj(save_output_folder+'pairs/pairs_{0}_{1:.3f}_{2:.3f}_{3}_{4}_{5}_{6}_{7}'.format(info['n_jck'],info['conf']['min_sep'],info['conf']['max_sep'],info['conf']['nbins'],title,kE_label,add_label,fract_limits),pairs)
    else:
        pairs = load_obj(save_output_folder+'pairs/pairs_{0}_{1:.3f}_{2:.3f}_{3}_{4}_{5}_{6}_{7}'.format(info['n_jck'],info['conf']['min_sep'],info['conf']['max_sep'],info['conf']['nbins'],title,kE_label,add_label,fract_limits))
   
    theta,DD_j,njk=make_wz_errors(pairs,'jackknife',n_jck,True)

    dictu=covariance_jck(DD_j[:,1:],DD_j.shape[1]-1,'jackknife')
    err=dictu['err']
    
    plt.title('systematic map (crosscorr) : {0}'.format(title))
    plt.xlabel('theta (arcmin)'.format(title))        
    plt.ylabel('w(theta)')
    plt.errorbar(theta,DD_j[:,0],err)
    plt.tight_layout()
#plt.yscale('log')
    plt.xscale('log')
    plt.savefig(save_output_folder + 'cc_' + title +'.png')
    plt.show()

    
    inv_cov = linalg.inv(dictu['cov'])
    chi2_val =  (np.matmul((DD_j[:,0]),np.matmul(inv_cov,(DD_j[:,0]))))/len(DD_j[:,0])

    #p_value = 1 - stats.chi2.cdf(np.sum(cv_sol ** 2), len(w.shape(0)))
    
    mute_d = dict()
    mute_d.update({'theta': theta})
    mute_d.update({'w':DD_j[:,0]})
    mute_d.update({'cov' : dictu})
    mute_d.update({'chi2_red' : chi2_val})    
    return mute_d



from multiprocessing import Pool,sharedctypes
from functools import partial
from contextlib import closing
import timeit

def comp(jk,s,m,hp,t):
    if t == 'depth_i' or t == 'depth_g' or t == 'depth_r' or t == 'depth_z':
        
        index = np.random.randint(0,len(s),200000)
        mask1 = (hp[index] != jk)
        mm = np.polyfit(s[index][mask1],m[index][mask1],1)
        
        return jk, mm[0],mm[1]      
    else:
        mask1 = (hp != jk)
        mm = np.polyfit(s[mask1],m[mask1],1)
        return jk, mm[0],mm[1]

def analysis_3(kE_label,sys_map,mask,mapa, hpix, fract_limits=0.,mapa_weight = False, cov_ext = False):
    from scipy import linalg

    print ('\n****************\n')
    title = sys_map['title']
    print ('TEST [method2] {1}: -> {0}'.format(title,kE_label))
    
    n_jck = dict_tot[binx]['info']['n_jck']
    fractional = sys_map['fractional']
    
    # plot typical ranges ***********************
    #area_tot = compute_area(mask, nside)
    y,x = np.histogram(sys_map['map'][mask], bins = 50)
    sys_map1 =  sys_map['map'][mask]
    if mapa_weight:
        mapa1 = mapa[mask]-mapa_weight[mask]
    else:
        mapa1 = mapa[mask]
        
    ave = sys_map1.mean()
    std = sys_map1.std()
    
    min_s, max_s = min(sys_map1),max(sys_map1)


    # plot maps properties *************************

    b_arr = np.zeros(n_jck+1)
    c_arr = np.zeros(n_jck+1)
    if title == 'depth_i' or title == 'depth_g' or title == 'depth_r' or title == 'depth_z':
        print ('fast')
        index = np.random.randint(0,len(sys_map1),200000)

        b_arr[0],c_arr[0] = np.polyfit(sys_map1[index],mapa1[index],1)
        
    else:
        b_arr[0],c_arr[0] = np.polyfit(sys_map1,mapa1,1)
    #print 'jck'
    

    aaa = range(n_jck)
    mas = []
    agents=64
    with closing(Pool(processes=agents)) as pool:
        mas.append(pool.map(partial(comp,s=sys_map1,m=mapa1,hp=hpix,t=title),aaa))
    mas = np.array(mas)

    

    mas =mas[0].reshape(n_jck,3)
    index = np.array(mas[:,0]).astype(int)+1
    b_arr[index]=mas[:,1]
    c_arr[index]=mas[:,2]
    #print mas
    #for i in range(n_jck):
    #    mask1 = (hpix != i)
    #    b_arr[i+1],c_arr[i+1] = np.polyfit(sys_map1[mask1],mapa1[mask1],1)
    
    xx = np.linspace(min_s,max_s,100)
    yy = b_arr[0]*xx+c_arr[0]
    dict_output=dict()
    dict_output.update({'b_arr' : b_arr})
    dict_output.update({'c_arr' : c_arr})
    dict_output.update({'b_err' : covariance_scalar_jck(b_arr[1:],n_jck)['err'] })
    
    dict_output.update({'min_s': min_s})
    dict_output.update({'max_s': max_s})

    dict_output.update({'xx': xx})
    dict_output.update({'yy': yy})

    return dict_output

def do_analysis_3(kE,kE_label,weight,systematic_maps,info,fract_limits=0,len_hist=10,hpix_type = 'normal',cov_ext = False, fit ='linear'):
    systematics_dict = dict() 
    for key in systematic_maps.keys():
      
      #if key =='depth_i'or key =='depth_g'or key =='depth_z'or key =='depth_r':
      # if ((key == 'snr' and kE_label =='kE') or (key == 'size_ratio' and kE_label =='kE') or (key == 'E1' and kE_label =='kE') or (key == 'E2' and kE_label =='kE')):
       # pass
       #else:
        t1 = timeit.default_timer()
        if hpix_type == 'normal':
            dict_output = analysis_3(kE_label,systematic_maps[key],info['mask_sims'],kE,info['hpix'])

        if hpix_type == 'marc':
            dict_output = analysis_3(kE_label,systematic_maps[key],info['mask_sims'],kE,info['hpix_f'])
        systematics_dict.update({systematic_maps[key]['title'] : dict_output})
        t2 = timeit.default_timer()
        print ("slope {0} +- {1}".format(dict_output['b_arr'][0],dict_output['b_err']))
        print (t2-t1)
    return systematics_dict



import config
class field_methods(object):
  """
  Utilities for doing pixel and chip calculations. Information from Mike Jarvis.
  """

  chip_centres = {

  'N7':[16.908,191.670],
  'N6':[16.908,127.780],
  'N5':[16.908,63.890],
  'N4':[16.908,0.],
  'N3':[16.908,-63.890],
  'N2':[16.908,-127.780],
  'N1':[16.908,-191.670],
  'N13':[50.724,159.725],
  'N12':[50.724,95.835],
  'N11':[50.724,31.945],
  'N10':[50.724,-31.945],
  'N9':[50.724,-95.835],
  'N8':[50.724,-159.725],
  'N19':[84.540,159.725],
  'N18':[84.540,95.835],
  'N17':[84.540,31.945],
  'N16':[84.540,-31.945],
  'N15':[84.540,-95.835],
  'N14':[84.540,-159.725],
  'N24':[118.356,127.780],
  'N23':[118.356,63.890],
  'N22':[118.356,0.],
  'N21':[118.356,-63.890],
  'N20':[118.356,-127.780],
  'N28':[152.172,95.835],
  'N27':[152.172,31.945],
  'N26':[152.172,-31.945],
  'N25':[152.172,-95.835],
  'N31':[185.988,63.890],
  'N30':[185.988,0.],
  'N29':[185.988,-63.890],
  'S7':[-16.908,191.670],
  'S6':[-16.908,127.780],
  'S5':[-16.908,63.890],
  'S4':[-16.908,0.],
  'S3':[-16.908,-63.890],
  'S2':[-16.908,-127.780],
  'S1':[-16.908,-191.670],
  'S13':[-50.724,159.725],
  'S12':[-50.724,95.835],
  'S11':[-50.724,31.945],
  'S10':[-50.724,-31.945],
  'S9':[-50.724,-95.835],
  'S8':[-50.724,-159.725],
  'S19':[-84.540,159.725],
  'S18':[-84.540,95.835],
  'S17':[-84.540,31.945],
  'S16':[-84.540,-31.945],
  'S15':[-84.540,-95.835],
  'S14':[-84.540,-159.725],
  'S24':[-118.356,127.780],
  'S23':[-118.356,63.890],
  'S22':[-118.356,0.],
  'S21':[-118.356,-63.890],
  'S20':[-118.356,-127.780],
  'S28':[-152.172,95.835],
  'S27':[-152.172,31.945],
  'S26':[-152.172,-31.945],
  'S25':[-152.172,-95.835],
  'S31':[-185.988,63.890],
  'S30':[-185.988,0.],
  'S29':[-185.988,-63.890]
  }

  bad_ccd_names = ['N30', 'S30', 'S7']
  ccdid = ['S29','S30','S31','S25','S26','S27','S28','S20','S21','S22','S23','S24','S14','S15','S16','S17','S18','S19','S8','S9','S10','S11','S12','S13','S1','S2','S3','S4','S5','S6','S7','N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18','N19','N20','N21','N22','N23','N24','N25','N26','N27','N28','N29','N30','N31']
  
  #bad_ccd_nums = [ ccdid.index(bad_ccd_name) for bad_ccd_name in bad_ccd_names]


  ccdx=2048.*15.e-6*1000. # col
  ccdy=4096.*15.e-6*1000. # row
  
  @staticmethod
  def ccd_centres():

    centrex=[]
    centrey=[]
    for i,x in enumerate(field_methods.ccdid):
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[0])
      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[1])

    return np.vstack((centrex,centrey)).T

  @staticmethod
  def ccd_corners():

    centre=np.zeros((62,4,2))
    for i,x in enumerate(field_methods.ccdid):
      c=field_methods.chip_centres.get(x,None)
      centre[i][0][0]=c[0]-field_methods.ccdx/2. # lower left
      centre[i][1][0]=c[0]-field_methods.ccdx/2. # lower right
      centre[i][2][0]=c[0]+field_methods.ccdx/2. # upper left
      centre[i][3][0]=c[0]+field_methods.ccdx/2. # upper right

      centre[i][0][1]=c[1]-field_methods.ccdy/2.
      centre[i][1][1]=c[1]+field_methods.ccdy/2.
      centre[i][2][1]=c[1]-field_methods.ccdy/2.
      centre[i][3][1]=c[1]+field_methods.ccdy/2.

    return centre

  @staticmethod
  def ccd_to_field(ccd,ccdx,ccdy):

    centre=field_methods.ccd_centres()

    centrex=(centre[:,0])[[ccd]]
    centrey=(centre[:,1])[[ccd]]

    return ccdx*15e-6*1000+centrex-field_methods.ccdx/2.,ccdy*15e-6*1000+centrey-field_methods.ccdy/2.

  @staticmethod
  def get_field_pos(cat):

    x,y=field_methods.ccd_to_field(cat.ccd,cat.col,cat.row)

    return x,y 

  @staticmethod
  def translate_to_wcs(pos,image):

    from esutil import wcsutil
    
    wcs=wcsutil.WCS(image, longpole=180.0, latpole=90.0, theta0=90.0)
    ra,dec=wcs.image2sky(pos[0],pos[1])

    return ra,dec

  @staticmethod
  def get_coadd_tile(ra,dec,tiles=None):

    if tiles is None:
      tiles=fio.FITS(config.coaddtiles)[-1].read()

    tmp=tiles['TILENAME'][(ra<tiles['URAUR'])&(dec<tiles['UDECUR'])&(ra>tiles['URALL'])&(dec>tiles['UDECLL'])]
    if len(tmp)==0:
      tmp=tiles['TILENAME'][((ra+360)<tiles['URAUR'])&(dec<tiles['UDECUR'])&((ra+360)>tiles['URALL'])&(dec>tiles['UDECLL'])]

    return tmp[0].rstrip()

  @staticmethod
  def get_radec_coadd_tiles(tiles=None,tiles0=None,file=config.coaddtiles):

    if tiles is None:
      tiles=fio.FITS(file)[-1].read()

    if tiles0 is None:
      mask=np.ones(len(tiles)).astype(bool)
    else:
      mask=np.in1d(np.core.defchararray.strip(tiles['TILENAME']),tiles0,assume_unique=False)

    return tiles,np.vstack(((tiles['URAUR'][mask]+tiles['URALL'][mask])/2.,(tiles['UDECUR'][mask]+tiles['UDECLL'][mask])/2.)).T



import numpy as np
import timeit
import treecorr # Module for correlation functions, pip install TreeCorr.
from astropy.cosmology import Planck15 as Planck15
import pickle
import sys

cosmol=Planck15

def covariance_jck(TOTAL_PHI,jk_r,type_cov):
  if type_cov=='jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_cov=='bootstrap':
      fact=1./(jk_r)
  #  Covariance estimation

  average=np.zeros(TOTAL_PHI.shape[0])
  cov_jck=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  err_jck=np.zeros(TOTAL_PHI.shape[0])


  for kk in range(jk_r):
    average+=TOTAL_PHI[:,kk]
  average=average/(jk_r)

 # print average
  for ii in range(TOTAL_PHI.shape[0]):
     for jj in range(ii+1):
          for kk in range(jk_r):
            cov_jck[jj,ii]+=TOTAL_PHI[ii,kk]*TOTAL_PHI[jj,kk]

          cov_jck[jj,ii]=(-average[ii]*average[jj]*jk_r+cov_jck[jj,ii])*fact
          cov_jck[ii,jj]=cov_jck[jj,ii]

  for ii in range(TOTAL_PHI.shape[0]):
   err_jck[ii]=np.sqrt(cov_jck[ii,ii])
 # print err_jck

  #compute correlation
  corr=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  for i in range(TOTAL_PHI.shape[0]):
      for j in range(TOTAL_PHI.shape[0]):
        corr[i,j]=cov_jck[i,j]/(np.sqrt(cov_jck[i,i]*cov_jck[j,j]))

  average=average*fact
  return {'cov' : cov_jck,
          'err' : err_jck,
          'corr':corr,
          'mean':average}






def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)#, encoding='latin1')
    
class Jack_shear(object):
    def __init__(self, conf, distance_file, ra_a, dec_a, ra_b, dec_b, hpix_a, hpix_b, w_a,w_b, e1_a, e2_a, e1_b, e2_b, k = None,  corr='GG',  number_of_cores=2, fact_dist=2,
                 centers=None, njk=None,
                 verbose=False, overwrite=False, cov = False,outp_fold='',label =''):

        '''
        Input:
        '''

        self.output_folder =outp_fold
        self.label = label
        self.conf = conf
        self.distance_file = distance_file

        self.number_of_cores = number_of_cores

        self.ra_s_1 = ra_a
        self.dec_s_1 = dec_a
        self.ra_s_2 = ra_b
        self.dec_s_2 = dec_b


        self.e1_a = e1_a#* (-1)
        self.e2_a = e2_a# (-1)
        self.e1_b = e1_b#
        self.e2_b = e2_b#* (-1)



        self.jk_a = hpix_a
        self.jk_b = hpix_b

        self.weight_a = w_a
        self.weight_b = w_b


        self.k = k
        self.corr = corr
        self.FACT_DIST = fact_dist

        self.njk = njk
        self.centers = centers
        self.cov = cov


    def Correlation(self):
        self.start_all = timeit.default_timer()
        self.prolog()
        pairs= self.epilog()
        return pairs



    def max_distance_region(self):

        def load_obj(name):
            with open(name + '.pkl', 'rb') as f:
                return pickle.load(f)#, encoding='latin1')


        max_dist_region = load_obj(self.distance_file)
        self.max_dist_region=max_dist_region


    def convert_units(self):
        """Change to 'degrees' or to Mpc the units if necessary.
        """

        if 'sep_units' in self.conf.keys():
            un = self.conf['sep_units']
            if un == 'arcmin':
                todeg = 1./60.
            elif un == 'arcsec':
                todeg = 1./60.**2.
            elif un == 'kpc':
                todeg= 1./1000.
            else:
                todeg = 1.
        else:
            todeg = 1.
        return todeg

    def dist_cent_2(self,ra1,dec1,ra2,dec2):

            todeg = np.pi/180.
            ra1 = ra1*todeg
            ra2 = ra2*todeg
            dec1 = dec1*todeg
            dec2 = dec2*todeg

            cos = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
            return np.arccos(cos)/todeg


    def distance(self):

        """Finds the minimum distance to a center for each center.
           Fixes double of this distance as the criteria for not considering correlations,
           which is a conservative choice. This distance has to be at least 4 times the
           maximum angular separation considered. Centers beyond this distance will not be
           considered in the correlations.
        """

        # Find the minimum distance to a center for each center.

        cnt = self.centers
        dist = np.array([np.sort([self.dist_cent(cnt[i],cnt[j]) for i in range(len(cnt))])[1] for j in range(len(cnt))])
        dist = (dist)*self.FACT_DIST


        todeg = self.convert_units()

        if 'max_sep' in self.conf.keys():
            max_sep = self.conf['max_sep'] * todeg

        else:
            raise NotImplementedError("Make use of 'max_sep' in configuration.")

        # Check that the distance is at least 4 times the maximum angular separation.

        self.center_min_dis = np.array( [ 4.*max_sep if x < 4.*max_sep else x for x in dist] )


    def dist_cent(self, a, b):
        """Angular distance between two centers (units: degrees). Makes use of spherical law of cosines.
        """
        todeg = np.pi/180.
        a = a*todeg
        b = b*todeg
        cos = np.sin(a[1])*np.sin(b[1]) + np.cos(a[1])*np.cos(b[1])*np.cos(a[0]-b[0])
        return np.arccos(cos)/todeg


    def cond(self, i, j):

        """Return the maximum conditional distance for a pair of centers to
           determine if the correlation should be computed.
        """

        con = self.center_min_dis
        return max(con[i], con[j])



    def collect(self, pairs):

        shape = (self.njk, self.conf['nbins'])
        DD_a, DR_a, RD_a, RR_a = np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(
            shape)
        DD, DR, RD, RR = np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins']), np.zeros(
            self.conf['nbins']), np.zeros(self.conf['nbins'])

        fact = 1.
        FACT_0 = 0.5
        #print pairs
        for n in range(self.conf['nbins']):
            for jk1 in range(len(pairs[:, 0, 0, n])):
                DD[n] += (pairs[jk1, 0, 0, n]) + FACT_0 * (pairs[jk1, 1, 0, n])
                DR[n] += ((pairs[jk1, 0, 1, n]) + FACT_0 * (pairs[jk1, 1, 1, n]))
                #print DR, DD, RD
                if self.corr == 'GG':
                    RD[n] += ((pairs[jk1, 0, 2, n]) + FACT_0 * (pairs[jk1, 1, 2, n]))
                    #print RD
        for n in range(self.conf['nbins']):
            for jk1 in range(len(pairs[:, 0, 0, n])):
                DD_a[jk1, n] = DD[n] - (pairs[jk1, 0, 0, n]) - fact * FACT_0 * (pairs[jk1, 1, 0, n])
                DR_a[jk1, n] = DR[n] - (pairs[jk1, 0, 1, n]) - fact * FACT_0 * (pairs[jk1, 1, 1, n])
                if self.corr == 'GG':
                    RD_a[jk1, n] = RD[n] - (pairs[jk1, 0, 2, n]) - fact * FACT_0 * (pairs[jk1, 1, 2, n])


        if self.corr == 'GG':
            xip = np.zeros(len(RD))
            xim = np.zeros(len(RD))
            masku = RD !=0.
            
            xip[masku] = DD[masku] / RD[masku]
            xim[masku] = DR[masku] / RD[masku]

            masku = RD_a !=0.
            xip_j = np.zeros((RD_a.shape[0],RD_a.shape[1]))
            xim_j = np.zeros((RD_a.shape[0],RD_a.shape[1]))
            xip_j[masku] = DD_a[masku] / RD_a[masku]
            xim_j[masku] = DR_a[masku] / RD_a[masku]


            return xip, xim, xip_j, xim_j
        else:
            xip = DD / DR
            xip_j = DD_a / DR_a

            return xip, xip, xip_j, xip_j





    def parallel(self, i, j):

        def dic(d):
            DD = d['DD']
            DR = d['DR']
            RD = d['RD']
            RR = d['RR']
            return [DD, DR, RD, RR]

        # this should be generalized to 2 different catalogs

        if self.corr == 'GG':
            mask = np.in1d(self.jk_a, i)
            try:
                cat_s1 = treecorr.Catalog(ra=self.ra_s_1[mask], dec=self.dec_s_1[mask], g1=self.e1_a[mask], g2=self.e2_a[mask],
                                      w=self.weight_a[mask],
                                      ra_units='deg', dec_units='deg')
            except:
                cat_s1 = None
            
            #print self.ra_s_1[mask],self.e1_a[mask]
            mask1 = np.in1d(self.jk_b, j)

            try:
                cat_s2 = treecorr.Catalog(ra=self.ra_s_2[mask1], dec=self.dec_s_2[mask1], g1=self.e1_b[mask1], g2=self.e2_b[mask1],
                                      ra_units='deg',
                                      dec_units='deg', w=self.weight_b[mask1])
            except:
                cat_s2 = None
            

            if (cat_s2 == None) or (cat_s1 == None):
                pairs = [np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins'])]
            else:
                gg = treecorr.GGCorrelation(self.conf)
                gg.process(cat_s1, cat_s2)
                ggp = gg.xip
                ggm = gg.xim
                normalization = gg.weight
                pairs = [ggp * normalization, ggm * normalization, normalization]
            return pairs

        if self.corr == 'KK':
            mask = np.in1d(self.jk_s, i)
            try:
                cat_a = treecorr.Catalog(ra=self.ra_s[mask], dec=self.dec_s[mask], ra_units='deg', dec_units='deg',
                                     k=self.k[mask])
            except:
                cat_a = None

            mask = np.in1d(self.jk_s, j)
            try:
                cat_b = treecorr.Catalog(ra=self.ra_s[mask], dec=self.dec_s[mask], ra_units='deg', dec_units='deg',
                                     k=self.k[mask])
            except:
                cat_b = None

            if (cat_a == None) or (cat_b == None):
                return [np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins'])]
            else:
                kk = treecorr.KKCorrelation(self.conf)
                kk.process(cat_a, cat_b)

                return [kk.xi * kk.weight, kk.weight]

    def prolog(self):

        def update_progress(progress, elapsed_time=0, starting_time=0):
            import time

            barLength = 10  # Modify this to change the length of the progress bar
            status = ""
            if isinstance(progress, int):
                progress = float(progress)
            if not isinstance(progress, float):
                progress = 0
                status = "error: progress var must be float\r\n"
            if progress < 0:
                progress = 0
                status = "Halt...\r\n"
            if progress >= 1:
                progress = 1
                status = "Done...\r\n"
            block = int(round(barLength * progress))

            if progress * 100 > 1. and elapsed_time > 0:
                remaining = ((elapsed_time - starting_time) / progress) * (1. - progress)
                text = "\rPercent: [{0}] {1:.2f}% {2}  - elapsed time: {3} - estimated remaining time: {4}".format(
                    "#" * block + "-" * (barLength - block), progress * 100, status,
                    time.strftime('%H:%M:%S', time.gmtime(elapsed_time - starting_time)),
                    time.strftime('%H:%M:%S', time.gmtime(remaining)))
            else:
                text = "\rPercent: [{0}] {1:.2f}% {2}".format("#" * block + "-" * (barLength - block), progress * 100,
                                                              status)
            sys.stdout.write(text)
            sys.stdout.flush()

        njk = self.njk

        self.distance()
        self.distance()
        self.max_distance_region()

        if self.corr == 'GG':
            gm = 3
            shape = (gm, self.conf['nbins'])

        if self.corr == 'KK':
            gm = 2
            shape = (gm, self.conf['nbins'])

        self.pairs_ring = [[np.zeros(shape) for i in range(2)] for j in range(njk)]
        cnt = self.centers
        self.pairs = [[np.zeros(shape) for i in range(njk)] for j in range(njk)]



        a=np.concatenate((np.array([[(i,j) for i in range(njk)] for j in range(njk)])))
        todeg = self.convert_units()
        max_sep = self.conf['max_sep'] * todeg

        sel = np.array([max([0.,(self.dist_cent(cnt[i],cnt[j]) - (self.max_dist_region[i]+self.max_dist_region[j]))]) < (3. * max_sep )for (i,j) in a])
        b = a[sel]
       # b = a

        def fun_speedup(othersample, othersample1, jackk):
            path =self.output_folder + self.label +'_{0}'.format(jackk)
            if os.path.exists(path+'.pkl'):
                dict_m = load_obj(path)
                pairsCC1 = dict_m['c1']
                pairsCC2 = dict_m['c2']
                pairs_auto = dict_m['a']
                for prs in range(gm):
                    pairsCC1[prs] += pairsCC2[prs]
                    

                self.pairs_ring[jackk][1] = pairsCC1
                self.pairs_ring[jackk][0] = pairs_auto                
            else:
                pairsCC1 = self.parallel(othersample, [jackk])

                pairsCC2 = self.parallel([jackk], othersample1)
                pairs_auto = self.parallel([jackk], [jackk])
                dict_m = dict()
                dict_m.update({'c1':pairsCC1})
                dict_m.update({'c2':pairsCC2})
                dict_m.update({'a':pairs_auto})
                save_obj(path,dict_m)
                for prs in range(gm):
                    pairsCC1[prs] += pairsCC2[prs]
                    

                self.pairs_ring[jackk][1] = pairsCC1
                self.pairs_ring[jackk][0] = pairs_auto

            #print 'eee' ,self.pairs_ring[0][1]#, 'ddd'

            #print (self.pairs_ring[jackk])

        startt = timeit.default_timer()
        
        # parallelized?#**********************************
        import os
        import shutil
        from multiprocessing import Pool,Manager,Process,Queue
        import math



        agents =self.number_of_cores

        chunks=int(math.ceil(np.float(len(np.unique(b[:, 1]))))/agents)
        mute_w = 0

        if self.cov:

         def run_s(jackk):
            #print jackk, "iter"
            mask = (b[:, 1] == jackk) & (b[:, 0] != jackk)
            othersample = b[mask, 0]
            mask = (b[:, 0] == jackk) & (b[:, 1] != jackk)
            othersample1 = b[mask, 1]
            fun_speedup(othersample, othersample1, jackk)                   
                   
         for i in range(chunks+1):
           workers=agents
           work_queue = Queue()
           done_queue = Queue()
           processes = []
           
           for w in range(agents ):
              if (mute_w<len(np.unique(b[:, 1]))):
              
               p = Process(target=run_s, args = ((np.unique(b[:, 1])[mute_w]),))
               p.start()
               processes.append(p)
               work_queue.put('STOP')
               mute_w+=1
        

           for p in processes:
               p.join()
            
       #loading
         for mute_w in (np.unique(b[:, 1])):
          path =self.output_folder + self.label +'_{0}'.format(mute_w)

          dict_m = load_obj(path)
          pairsCC1 = dict_m['c1']
          pairsCC2 = dict_m['c2']
          pairs_auto = dict_m['a']
          for prs in range(gm):
            pairsCC1[prs] += pairsCC2[prs]
                    

            self.pairs_ring[mute_w][1] = pairsCC1
            self.pairs_ring[mute_w][0] = pairs_auto       
            
                     
         '''
        for counter, jackk in enumerate(np.unique(b[:, 1])):
            
            if self.cov:
                mask = (b[:, 1] == jackk) & (b[:, 0] != jackk)
                othersample = b[mask, 0]
                mask = (b[:, 0] == jackk) & (b[:, 1] != jackk)
                othersample1 = b[mask, 1]
                fun_speedup(othersample, othersample1, jackk)
            
            update_progress((float(counter)/len(np.unique(b[:,1]))),timeit.default_timer(),startt)
        '''
        self.pairs = np.array(self.pairs_ring)

    def epilog(self):
        
        pairs = self.pairs
        xip1,xim1,xip_j,xim_j = self.collect(pairs)
        convert=self.convert_units()
        min_sep, max_sep, nedges = self.conf['min_sep']*convert, self.conf['max_sep']*convert, self.conf['nbins']+1
        th = np.linspace(np.log10(min_sep), np.log10(max_sep), nedges)
        theta = 10**np.array([(th[i]+th[i+1])/2 for i in range(len(th)-1)])
        end1 = timeit.default_timer()
        ggp_1=0.
        ggm_1=0.

        xip=0.
        xim=0.
        start =  timeit.default_timer()
        end =  timeit.default_timer()
        if not self.cov:
            #crosscheck ***

            if self.corr == 'KK':

                cat_l1 = treecorr.Catalog(ra=self.ra_s, dec=self.dec_s, ra_units='deg', dec_units='deg',k=self.k)
                cat_l2 = treecorr.Catalog(ra=self.ra_s, dec=self.dec_s, ra_units='deg', dec_units='deg',k=self.k)

                dd = treecorr.KKCorrelation(self.conf)

                dd.process(cat_l1, cat_l2)


                xip = dd.xi
                xim = dd.xi
            if self.corr == 'GG':


                start = timeit.default_timer()
                cat_s1 = treecorr.Catalog(ra=self.ra_s_1, dec=self.dec_s_1, g1=self.e1_a, g2=self.e2_a, w=self.weight_a,
                                          ra_units='deg', dec_units='deg')
                cat_s2 = treecorr.Catalog(ra=self.ra_s_2, dec=self.dec_s_2, g1=self.e1_b, g2=self.e2_b, w=self.weight_b,
                                          ra_units='deg', dec_units='deg')
                gg = treecorr.GGCorrelation(self.conf)

                gg.process(cat_s1,cat_s2)
                end = timeit.default_timer()


                xip = gg.xip
                xim = gg.xim
        
        return {'theta': theta, 'xip':xip, 'xim': xim, 'corr_jckp': xip_j,'corr_jckm': xim_j, 'time_fast' : (end1-self.start_all), 'time_slow': (end-start),'xip1':xip1,'xim1':xim1,}

#plot_resume(fold,systematics_dict_e1,z_min,z_max,nside,kE_label='e1', d68=True,keys=keymport numpy as np

def covariance_scalar_jck(TOTAL_PHI,jk_r, type_c = 'jackknife'):

  #  Covariance estimation
  if type_c == 'jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_c=='bootstrap':
      fact=1./(jk_r)
        
  average=0.
  cov_jck=0.
  err_jck=0.


  for kk in range(jk_r):
    average+=TOTAL_PHI[kk]
  average=average/(jk_r)

  for kk in range(jk_r):
    #cov_jck+=TOTAL_PHI[kk]#*TOTAL_PHI[kk]

    cov_jck+=(-average+TOTAL_PHI[kk])*(-average+TOTAL_PHI[kk])


  err_jck=np.sqrt(cov_jck*fact)


  #average=average*(jk_r)/(jk_r-1)
  return {'cov' : cov_jck*fact,
          'err' : err_jck,
          'mean': average}


def analysis(save_output_folder,kE_label,sys_map,mask,mapa, hpix, fract_limits,len_hist=10,mapa_weight = None, cov_ext = False):
    from scipy import linalg
    print ('\n****************\n')
    print ('TEST [method1] {1}: -> {0}'.format(sys_map['title'],kE_label))
    
    if cov_ext:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[cov_ext['indexes']])
                                      
        mask = cov_ext['indexes'] # it can be slightly different
    else:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[mask])     
    fractional = sys_map['fractional']
    
    # plot typical ranges ***********************
    area_tot = compute_area(mask, nside)
    y,x = np.histogram(sys_map['map'][mask], bins = 50)
    sys_map1 =  sys_map['map'][mask]
    mapa1 = mapa[mask]-mapa_weight[mask]

    #print len(mapa1),len(sys_map1)
    #mapa1 = mapa1 - np.mean(mapa1)
    
    ave = sys_map1.mean()
    std = sys_map1.std()
    if fractional:
    
        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
    
    
    min_s, max_s = find_limits(sys_map1, fract_limits)
    if not fractional:
        if 'min_f' in sys_map.keys():
            min_s = max([min_s,sys_map['min_f']-ave])
        if 'max_f' in sys_map.keys():
            max_s = min([max_s,sys_map['max_f']-ave])   
    if fractional:
        min_ss = min_s*ave + ave
        max_ss = max_s*ave + ave    
    else:
        min_ss = min_s + ave
        max_ss = max_s + ave

    
    def histedges_equalN(x, nbin):
        npt = len(x)
        return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

    mh = (sys_map1>min_s) & (sys_map1<max_s)
    bins = histedges_equalN(sys_map1[mh], len_hist)

    
    #print np.histogram(sys_map1, bins= bins)
    
    #define bin centers as the median point.
    x = np.zeros(len(bins)-1)
    for i in range(len(x)):

        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        x[i] = np.median(sys_map1[mask_bin])
        
        
        #print bins[i],x[i],bins[i+1]
    #x = bins[:-1]+0.5*(bins[1]-bins[0])
     
    mask_asd = (sys_map1 > min_s ) & (sys_map1 < max_s )
    mask_pixels = (sys_map1 < min_s) ^ ( sys_map1 > max_s)
    area_rest = compute_area(mask_pixels, nside)
    print ('Area excluded: {0:.1f}/{1:.1f}'.format(area_rest,area_tot))
    

    # plot maps proerties *************************
    value_jck = np.zeros((len_hist,n_jck+1))



    for i in range(len_hist):
        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        value_jck[i,0] = (mapa1[mask_bin]).sum()
        for j in range(n_jck):
            #print ((bins[i],bins[i+1]))
            
            mask_bin1 = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) & (hpix[mask_hpix] == j)
            value_jck[i,j+1] = (value_jck[i,0] - (mapa1[mask_bin1]).sum())/(len(mapa1[mask_bin])-len(mapa1[mask_bin1]))
        value_jck[i,0] = value_jck[i,0]/ len(mapa1[mask_bin])
    dict_jck = covariance_jck(value_jck[:,1:], n_jck, 'jackknife')


    

    params_guess= np.array([0.0,0])
    params_guess2= np.array([.0,0.0,0])
    params_guess3= np.array([.0,.0,0.0,0])
    if not cov_ext:

        ccov = dict_jck#['cov']
        
    else:
        vector_cov = np.zeros((len_hist,cov_ext['cov'].shape[1]))
        for i in range(len_hist):
            mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1])# & cov_ext['indexes']
            vector_cov[i,:] = np.mean(cov_ext['cov'][mask_bin,:],axis = 0)
        mask_bin = (sys_map1 > bins[0]) & (sys_map1 < bins[-1])

        ccov = covariance_jck(vector_cov, cov_ext['cov'].shape[1], 'bootstrap')
        '''
        mute= ccov['cov'].diagonal()
        cov_mute = np.zeros((len(mute),len(mute)))
        for ff in range(len(mute)):
            cov_mute[ff,ff] = mute[ff]
        ccov['cov'] = cov_mute
        '''
        
    dict_output=dict()
    if not cov_ext:

        poptbts, pcovbts = curve_fit(fitting_linear, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess)
        poptbts2, pcovbts2 = curve_fit(fitting_2nd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess2)
        poptbts3, pcovbts3 = curve_fit(fitting_3rd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess3)

        paramsbts=poptbts
        paramsbts2=poptbts2
        paramsbts3=poptbts3
        #print ('slope = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))
        wbt = value_jck[:,0]
        cns = np.zeros(len(wbt))* np.mean(mapa1[mask_asd])
        inv_covbt = linalg.inv(ccov['cov'])
        
        # jackknife values ********************************************************
        cns = fitting_2nd(x,poptbts2[0],poptbts2[1],poptbts2[2]) 
        chi2redbt2 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        cns = fitting_3rd(x,poptbts3[0],poptbts3[1],poptbts3[2],poptbts3[3]) 
        chi2redbt3 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = fitting_linear(x,poptbts[0],poptbts[1]) 
        chi2redbt =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = np.ones(len(wbt))* np.mean(mapa1[mask_asd])
        chi2redbt_null = (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        #print 'pull:'
        ##print (wbt-cns)/ccov['err'],np.sum(((wbt-cns)**2/ccov['err']**2.))
        chi2redbt_diff = chi2redbt_null - chi2redbt 


  
        #********************************************************
        dict_output.update({'chi2red_null': chi2redbt_null})
        dict_output.update({'chi2red_diff': chi2redbt_diff})
        dict_output.update({'chi2red': chi2redbt})
        dict_output.update({'chi2red2': chi2redbt2})
        dict_output.update({'chi2red3': chi2redbt3})

        popt, pcov = curve_fit(fitting_linear, x, value_jck[:,0],sigma = dict_jck['cov'],p0=params_guess)

        params=popt
        print ('slope jck = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))


        
    else:
        
        # all the mocks *****************
        chi2_vect = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect1 = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect2 = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect3 = np.zeros(cov_ext['cov'].shape[1])
        chi2_vect4 = np.zeros(cov_ext['cov'].shape[1]) 
        
        #print 'mean mocks',np.mean(cov_ext['cov'][:,0])
        for jj in range(cov_ext['cov'].shape[1]):

            popt, pcov = curve_fit(fitting_linear, x, vector_cov[:,jj],sigma = ccov['cov'],p0=params_guess)

            params=popt
            #print ('slope jck = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))
            w = vector_cov[:,jj]
            #print np.mean(cov_ext['cov'][:,j])
            cns = np.ones(len(w))*  np.mean(cov_ext['cov'][:,jj])
            inv_cov = linalg.inv(ccov['cov'])
            chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            #plt.plot(x,cns)

            cns = fitting_linear(x,popt[0],popt[1]) 

            chi2_vect[jj] = chi2red - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            chi2_vect1[jj] = cimhi2red# - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            chi2_vect2[jj] = (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
            

            mutw = chi2red/(len(x)-2) - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))//(len(x)-1) 
            df = 1
            p = 1 - stats.t.cdf(mutw,df=df)
            chi2_vect3[jj] = p2s(p)
            
            mutw = chi2red - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))
            df = 1
            p = 1 - stats.t.cdf(mutw,df=df)
            chi2_vect4[jj] = p2s(p)


    
        m1 = np.sort(chi2_vect1)
        m2 = np.sort(chi2_vect2)
        m3 = np.sort(chi2_vect)
        m4 = np.sort(chi2_vect3)
        m5 = np.sort(chi2_vect4)
        chi681 = m1[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        chi682 = m2[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        chi683 = m4[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        chi684 = m5[int(0.68*np.float(cov_ext['cov'].shape[1]))]
        # chi68 out of 1000 sims
        chi68 = m3[int(0.68*np.float(cov_ext['cov'].shape[1]))]

        poptbts, pcovbts = curve_fit(fitting_linear, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess)
        poptbts2, pcovbts2 = curve_fit(fitting_2nd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess2)
        poptbts3, pcovbts3 = curve_fit(fitting_3rd, x, value_jck[:,0],sigma = ccov['cov'],p0=params_guess3)

        paramsbts=poptbts
        paramsbts2=poptbts2
        paramsbts3=poptbts3
        #print ('slope resampling = {0:.4f} +- {1:.4f}').format(paramsbts[0],np.sqrt(pcovbts[0,0]))
        wbt = value_jck[:,0]
        #print 'mean data',np.mean(mapa1[mask_asd])
        
        inv_covbt = linalg.inv(ccov['cov'])
        
        # bootstrap values ********************************************************
        cns = fitting_2nd(x,poptbts2[0],poptbts2[1],poptbts2[2]) 
        chi2redbt2 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        cns = fitting_3rd(x,poptbts3[0],poptbts3[1],poptbts3[2],poptbts3[3]) 
        chi2redbt3 =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = fitting_linear(x,poptbts[0],poptbts[1]) 
        chi2redbt =  (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        
        cns = np.ones(len(wbt))* np.mean(mapa1[mask_asd])
        chi2redbt_null = (np.matmul((wbt-cns),np.matmul(inv_covbt,(wbt-cns))))#/len(wbt)
        #print 'pull:'
        ##print (wbt-cns)/ccov['err'],np.sum(((wbt-cns)**2/ccov['err']**2.))
        chi2redbt_diff = chi2redbt_null - chi2redbt 
        
        #print "                        chi2/dof (d chi2/dof) [d chi2/dof _68]"
        #print ('chi2/dof resampling :  {0:.2f} ({2:.3f}) [{3:.3f}]').format(chi2redbt,sys_map['title'],chi2redbt_diff,chi2redbt_diff/chi68)
         
            
        #********************************************************
        dict_output.update({'chi2red_boot_null': chi2redbt_null})
        dict_output.update({'chi2red_diff_boot': chi2redbt_diff})
        dict_output.update({'chi2red_boot': chi2redbt})
        dict_output.update({'chi2red_boot2': chi2redbt2})
        dict_output.update({'chi2red_boot3': chi2redbt3})
        dict_output.update({'chi2red_diff_boot68': chi2redbt_diff/chi68}) 
        dict_output.update({'chi68': chi68}) 
        popt, pcov = curve_fit(fitting_linear, x, value_jck[:,0],sigma = dict_jck['cov'],p0=params_guess)

        params=popt
        print ('slope jck = {0:.4f} +- {1:.4f}').format(params[0],np.sqrt(pcov[0,0]))
        w = value_jck[:,0]
        cns = np.ones(len(w))* np.mean(mapa1[mask_asd])
        inv_cov = linalg.inv(dict_jck['cov'])
        chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
        
        cns = fitting_linear(x,popt[0],popt[1]) 
        chi2red_diff = chi2red - (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))#/len(w)
        
        print ('chi2/dof jck  <{1}>0 =<{1}>:  {0:.2f} ({2:.3f})  [{3:.3f}]').format(chi2red,sys_map['title'],chi2red_diff,chi2red_diff/chi68)
        
        dict_output.update({'chi2red_diff': chi2red_diff})
        dict_output.update({'chi2red_diff68': chi2red_diff/chi68})
    

  
    
    if not cov_ext:
        dict_output.update({'dict_cov' : ccov})
        dict_output.update({'params_bts' :paramsbts})
        dict_output.update({'params_bts2' :paramsbts2})
        dict_output.update({'params_bts3' :paramsbts3})
        dict_output.update({'params_cov_bts' : pcovbts})  
    else:
        dict_output.update({'dict_cov' : dict_jck})
        dict_output.update({'dict_cov_bts' : ccov})
        dict_output.update({'params_bts' :paramsbts})
        dict_output.update({'params_bts2' :paramsbts2})
        dict_output.update({'params_bts3' :paramsbts3})
        dict_output.update({'params_cov_bts' : pcovbts})        
        
        dict_output.update({'v':vector_cov})
        dict_output.update({'chi68_new':chi683})
        dict_output.update({'chi68_new1':chi684})
        dict_output.update({'chi68_array':m4})
    dict_output.update({'w' :value_jck[:,0]})
    dict_output.update({'x' :x})   
    dict_output.update({'bins' :bins})       

    dict_output.update({'params' :params})
    dict_output.update({'params_cov' : pcov})
    dict_output.update({'min_s': min_s})
    dict_output.update({'max_s': max_s})
    dict_output.update({'min_ss': min_ss})
    dict_output.update({'max_ss': max_ss})
    dict_output.update({'ave_value': ave})
    dict_output.update({'mask_pixel': mask_pixels})
    dict_output.update({'value_jck': value_jck})
    dict_output.update({'Area': area_rest/area_tot})
    dict_output.update({'mean shear': np.mean(mapa1[mask_asd])})


    if not cov_ext:
        label_boot =''
    else:
        label_boot ='_boot'
    if 1==1:
        chi = dict_output['chi2red'+label_boot+'_null'] 
        df1 = len(dict_output['x']) - 1
        chi=chi/df1
        p = 1 - stats.t.cdf(chi,df=df1)
        #print p2s(p)
        dict_output.update({'s_null': p2s(p)})
        chi = dict_output['chi2red'+label_boot] 
        df2 = len(dict_output['x']) - 2
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
    
        dict_output.update({'s_fit': p2s(p)})
        #print p2s(p)
    
        chi = dict_output['chi2red'+label_boot+'2'] 
        df2 = len(dict_output['x']) - 3
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
    
        dict_output.update({'s_fit2': p2s(p)})
        #print p2s(p)
    
        chi = dict_output['chi2red'+label_boot+'3'] 
        df2 = len(dict_output['x']) - 4
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
    
        dict_output.update({'s_fit3': p2s(p)})
        #print p2s(p)
        chi = dict_output['chi2red'+label_boot+'_null']/df1 - dict_output['chi2red'+label_boot+'']/df2 
        df = 1
        p = 1 - stats.t.cdf(chi,df=df)
        

        
        #print p2s(p)
        dict_output.update({'s_H1': p2s(p)})
        
        
        
        
        #******************************************
        chi = dict_output['chi2red'+label_boot] 
        df2 = len(dict_output['x']) - 2
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'_null'] ,dict_output['chi2red'+label_boot] ,len(dict_output['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)
        dict_output.update({'F01': p2s(pf)})
    
        #******************************************
        chi = dict_output['chi2red'+label_boot+'2'] 
        df2 = len(dict_output['x']) - 3
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'_null'] ,dict_output['chi2red'+label_boot+'2'] ,len(dict_output['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot] ,dict_output['chi2red'+label_boot+'2'] ,len(dict_output['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        dict_output.update({'F02': p2s(pf)})
        dict_output.update({'F12': p2s(pf2)})
    
        #******************************************
        chi = dict_output['chi2red'+label_boot+'3'] 
        df2 = len(dict_output['x']) - 4
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'_null'] ,dict_output['chi2red'+label_boot+'3'] ,len(dict_output['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)

        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot] ,dict_output['chi2red'+label_boot+'3'] ,len(dict_output['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
 
        F,v1,v2 = compute_fisher(dict_output['chi2red'+label_boot+'2'] ,dict_output['chi2red'+label_boot+'3'] ,len(dict_output['x']),1,2)  
        pf3 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        print ('ok')
        dict_output.update({'F03': p2s(pf)})
        dict_output.update({'F13': p2s(pf2)})
        dict_output.update({'F23': p2s(pf3)})  
        

    return dict_output

def do_analysis(out_file,add_label,kE,kE_label,weight,systematic_maps,info,fract_limits,len_hist=10,hpix_type = 'normal',cov_ext = False, fit ='linear'):
    systematics_dict = dict() 
    systematics_dict_2 = dict() 
    for key in systematic_maps.keys():
      
      if key != kE_label:
       if ((key == 'snr' and kE_label =='kE') or (key == 'size_ratio' and kE_label =='kE') or (key == 'E1' and kE_label =='kE') or (key == 'E2' and kE_label =='kE')):
        pass
       else:
        if hpix_type == 'normal':
            dict_output = analysis(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix'], fract_limits,len_hist=len_hist,mapa_weight = weight, cov_ext = cov_ext)
            dict_output_2 = analysis_2(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix'], fract_limits,len_hist=len_hist,mapa_weight = weight, cov_ext = cov_ext)
           # print  dict_output_2['b_arr'][0], ' +- ',  dict_output_2['b_err']
      
        if hpix_type == 'marc':
            dict_output = analysis(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix_f'], fract_limits,len_hist=len_hist, mapa_weight = weight, cov_ext = cov_ext)
            dict_output_2 = analysis_2(out_file,kE_label,systematic_maps[key],info['mask_sims'],kE, info['hpix_f'], fract_limits,len_hist=len_hist,mapa_weight = weight, cov_ext = cov_ext)

        make_plot(out_file,add_label,dict_output,systematic_maps[key]['title'],kE_label,systematic_maps[key]['map'],systematic_maps[key]['fractional'],info, fract_limits,len_hist = len_hist,cov_ext = cov_ext,fit=fit,dict_output_2=dict_output_2)
        systematics_dict.update({systematic_maps[key]['title'] : dict_output})
        systematics_dict_2.update({systematic_maps[key]['title'] : dict_output_2})

    return systematics_dict,systematics_dict_2


def derivative(x,f):
            xmi = int(len(x)/2.)
            dx_i = x[xmi+1]-x[xmi]
            dx_mi = x[xmi]-x[xmi-1]
            f_i = f[xmi+1]
            f_0 = f[xmi]
            f_mi = f[xmi-1]
            b = -dx_i*f_mi / (dx_mi*(dx_i+dx_mi)) + (dx_i-dx_mi)*f_0/(dx_i*dx_mi)+dx_mi*f_i/(dx_i*(dx_i+dx_mi))
            return b,xmi

        
def analysis_2(save_output_folder,kE_label,sys_map,mask,mapa, hpix, fract_limits,len_hist=10,mapa_weight = None, cov_ext = False):
    from scipy import linalg
    print ('\n****************\n')
    print ('TEST [method2] {1}: -> {0}'.format(sys_map['title'],kE_label))
    
    if cov_ext:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[cov_ext['indexes']])                        
        mask = cov_ext['indexes'] # it can be slightly different
    else:
        c = np.arange(len(mask))
        mask_hpix = np.in1d(c[mask],c[mask])   
    fractional = sys_map['fractional']
    
    # plot typical ranges ***********************
    area_tot = compute_area(mask, nside)
    y,x = np.histogram(sys_map['map'][mask], bins = 50)
    sys_map1 =  sys_map['map'][mask]
    mapa1 = mapa[mask]-mapa_weight[mask]

    
    ave = sys_map1.mean()
    std = sys_map1.std()
    if fractional:
    
        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
    
    
    min_s, max_s = find_limits(sys_map1, fract_limits)
    if not fractional:
        if 'min_f' in sys_map.keys():
            min_s = max([min_s,sys_map['min_f']-ave])
        if 'max_f' in sys_map.keys():
            max_s = min([max_s,sys_map['max_f']-ave])   
    if fractional:
        min_ss = min_s*ave + ave
        max_ss = max_s*ave + ave    
    else:
        min_ss = min_s + ave
        max_ss = max_s + ave

    
    def histedges_equalN(x, nbin):
        npt = len(x)
        return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

    mh = (sys_map1>min_s) & (sys_map1<max_s)
    bins = histedges_equalN(sys_map1[mh], len_hist)

    
    #print np.histogram(sys_map1, bins= bins)
    
    #define bin centers as the median point.
    x = np.zeros(len(bins)-1)
    for i in range(len(x)):

        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        x[i] = np.median(sys_map1[mask_bin])
        

    mask_asd = (sys_map1 > min_s ) & (sys_map1 < max_s )
    mask_pixels = (sys_map1 < min_s) ^ ( sys_map1 > max_s)
    area_rest = compute_area(mask_pixels, nside)
    print ('Area excluded: {0:.1f}/{1:.1f}'.format(area_rest,area_tot))
    

    # plot maps proerties *************************
    value_jck = np.zeros((len_hist,n_jck+1))



    for i in range(len_hist):
        mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) 
        value_jck[i,0] = (mapa1[mask_bin]).sum()
        for j in range(n_jck):
            #print ((bins[i],bins[i+1]))
            
            mask_bin1 = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1]) & (hpix[mask_hpix] == j)
            value_jck[i,j+1] = (value_jck[i,0] - (mapa1[mask_bin1]).sum())/(len(mapa1[mask_bin])-len(mapa1[mask_bin1]))
        value_jck[i,0] = value_jck[i,0]/ len(mapa1[mask_bin])
    dict_jck = covariance_jck(value_jck[:,1:], n_jck, 'jackknife')



    if not cov_ext:
        ccov = dict_jck#['cov']
    else:
        vector_cov = np.zeros((len_hist,cov_ext['cov'].shape[1]))
        for i in range(len_hist):
            mask_bin = (sys_map1 > bins[i]) & (sys_map1 < bins[i+1])# & cov_ext['indexes']
            vector_cov[i,:] = np.mean(cov_ext['cov'][mask_bin,:],axis = 0)
        mask_bin = (sys_map1 > bins[0]) & (sys_map1 < bins[-1])

        ccov = covariance_jck(vector_cov, cov_ext['cov'].shape[1], 'bootstrap')

        
    
    
    # ***********************************************************************************************
    b_arr = np.zeros(n_jck+1)
    if not cov_ext:

        # compute b with non linear
        b_arr[0],xmi = derivative(x,value_jck[:,0])
        for jk in range(n_jck):
            b_arr[1+jk],xmi = derivative(x,value_jck[:,1+jk])
        b_err = covariance_scalar_jck(b_arr[1:],n_jck)['err']
    
    else:
        # compute b with non linear
        b_arr[0],xmi = derivative(x,value_jck[:,0])
        for jk in range(n_jck):
            b_arr[1+jk],xmi = derivative(x,vector_cov[:,jk])
        b_err = covariance_scalar_jck(b_arr[1:],n_jck, type_c = 'bootstrap')['err']
           
    dict_output=dict()
    dict_output.update({'b_arr' : b_arr})
    dict_output.update({'xmi' : xmi})
    dict_output.update({'x' : x})
    dict_output.update({'b_err' : b_err})   
    dict_output.update({'bins' :bins})       


    dict_output.update({'min_s': min_s})
    dict_output.update({'max_s': max_s})
    dict_output.update({'min_ss': min_ss})
    dict_output.update({'max_ss': max_ss})
    dict_output.update({'ave_value': ave})
    dict_output.update({'mask_pixel': mask_pixels})
    dict_output.update({'value_jck': value_jck})
    if cov_ext:
        dict_output.update({'vector_cov': vector_cov})
    dict_output.update({'Area': area_rest/area_tot})
    dict_output.update({'mean shear': np.mean(mapa1[mask_asd])})
        

    return dict_output





def make_plot_final_method1(dict_systematics_tot_1,label,z_minz,z_maxz,fit ='linear'):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(10,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    if fit =='linear':
        fit_f ='F01'
    elif fit =='quadratic':
        fit_f ='F02'
    elif fit =='cubic':
        fit_f ='F03'
        
    hpix_type = 'normal'

    
    xxx = []
    for kk in dict_systematics_tot_1['{0}_{1}'.format(z_minz[0],z_maxz[0])][label]['_w0'].keys():
        xxx.append(kk)
        
    for i in range(len(z_minz)):

        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        #for count in range(20):
        add_label = '_w0'#.format(count)

        if 1==1:
                mute = dict_systematics_tot_1[binx][label][add_label]
                
                ax[i].set_ylabel('bin : {0} \n F-test'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                yy =[]
                for kk in xxx:
                    yy.append(dict_systematics_tot_1[binx][label][add_label][kk][fit_f])
                    
                ax[i].scatter(np.arange(len(xxx)), yy,label = '{0} fit'.format(fit))
                #ax[i].scatter(mute['x'], mute['original'],label = 'linear_fit (no weight)')
                ll = []
        
                '''
                for xx in (xxx):
                    try:
                        ll.append( '{0} {1}'.format(str(xx.split('_')[2]),str(xx.split('_')[3])))
                    except:
                        ll.append( '{0}'.format(str(xx.split('_')[2])))#,str(xx.split('_')[3])))
                '''
                ax[i].plot([-1.,len(xxx)+1],[2,2], linestyle = 'dashed',color = 'black', alpha = 0.7)
                ax[i].plot([-1.,len(xxx)+1],[1,1], linestyle = 'dashed',color = 'black', alpha = 0.7)
                plt.xticks(np.arange(len(xxx)),xxx, rotation='90')
                ax[0].legend()
                plt.suptitle('mean {0} systematic test'.format(label))
                plt.xlim([-1,len(xxx)])
                #print count, mute.keys()
        #except:
        #        pass
        fig.subplots_adjust(wspace=0.,hspace=0.)
        
def make_plot_final(dict_final1, label,z_minz,z_maxz,type_c = 'jackknife',fact = False,Except=False,only=False):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(10,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    bd = dict()
    for i in range(len(z_minz)):
        #i = 4-i
        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        fold = './output_new1_{0}_{1}_{2}/'.format(z_min,z_max,nside)
        #for count in range(20):
        if 1==1:
                if 1==1:
                    ll = []
                    ll1 = []
                    key = binx
                    m1 = 0
                    count = 0
                    for key1 in dict_final1[key].keys():
                        if only: 
                            new_key = only[key1]
                            m1+=len(only[key1])
                        elif Except:
                            new_key = []
                            for ii,jj in enumerate(~np.in1d(dict_final1[key][key1]['_w0'].keys(),Except[key1])):
                                if jj:
                                    new_key.append(dict_final1[key][key1]['_w0'].keys()[ii])
                            m1 += len(new_key)
                        else:
                            m1 += len(dict_final1[key][key1]['_w0'].keys())
                            new_key = dict_final1[key][key1]['_w0'].keys()
                    b1 = np.zeros((m1,n_jck+1))
                    for key1 in dict_final1[key].keys():
                        
                        if only: 
                            new_key = only[key1]
             
                        elif Except:
                            new_key = []
                            for ii,jj in enumerate(~np.in1d(dict_final1[key][key1]['_w0'].keys(),Except[key1])):
                                if jj:
                                    new_key.append(dict_final1[key][key1]['_w0'].keys()[ii])
                        else:
                            new_key = dict_final1[key][key1]['_w0'].keys()
                            
                        for key2 in new_key:
                            b1[count,:]= dict_final1[key][key1]['_w0'][key2]['b_arr']
                            count+=1
                            ll.append('{0}_{1}'.format(key1,key2))


                if fact:
                    b1 = b1
                    b_dict1 = covariance_jck(b1[:,1:]*fact[binx], n_jck, type_cov = type_c)
                else:
                    b_dict1 = covariance_jck(b1[:,1:], n_jck, type_cov = type_c)
                
                
                ax[i].set_ylabel('bin : {0} \n b/b_err'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                
                ax[i].scatter(np.arange(len( b1[:,0]))+1, b1[:,0]/b_dict1['err'],label = 'b/b_err')
                
                
               
                #ax[i].scatter(np.arange(len( b1[:,0])),  p2s(1 - stats.t.cdf((b1[:,0]/b_dict1['err'])**2.,df=1.)) ,label = 'sigma')

                #print b1[:,0]/b_dict1['err']
                    #print ll
                   # #print '{0} {1}'.format(str(xx.split('_')[2]),str(xx.split('_')[3])) ,mute['label'][ii]
                #ax[i].plot([0.,len( b1[:,0])+1],[2,2], linestyle = 'dashed',color = 'black', alpha = 0.7)
                #ax[i].plot([0.,len( b1[:,0])+1],[1,1], linestyle = 'dashed',color = 'black', alpha = 0.7)
                plt.xticks(np.arange(len( b1[:,0]))+1,ll, rotation='90')
                #ax[0].legend()
                #ax[0].suptitle('{0}'.format(label))
                ax[0].text(0.4,1.05 , '{0}'.format(label),transform=ax[0].transAxes)
                
                w = b1[:,0]
                cns = np.zeros(len(b1[:,0]))
                from scipy import linalg
                inv_cov = linalg.inv(b_dict1['cov'])
                N_p = n_jck
                p_p = len(cns)
                f_hartlap = (N_p-1.)/(N_p-p_p-2.)
                #print f_hartlap
                chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/(len(cns)*f_hartlap)

                ax[i].text(1.1 , 0.5 , 'chi/dof : {0:.2f}; sigma : {1:.2f}'.format(chi2red,p2s(1 - stats.t.cdf(chi2red,df=len(cns)))),transform=ax[i].transAxes)
     #plt.text(0, 5.6, ' Average mean shear = {0:.5f}'.format(np.mean(points_2)))
                plt.xlim([0.,len( b1[:,0])+1])
                #print count, mute.keys()

        fig.subplots_adjust(wspace=0.,hspace=0.)
        bd.update({binx:[b_dict1,b1]})
    return bd, ll


def final_plot_Pearson(dict_final1, label,z_minz,z_maxz,type_c = 'jackknife',fact = False):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(10,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    bd = dict()
    for i in range(len(z_minz)):
        #i = 4-i
        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        fold = './output_new1_{0}_{1}_{2}/'.format(z_min,z_max,nside)
        #for count in range(20):
        if 1==1:
                if 1==1:
                    ll = []
                    
                    m1 = 0
                    for key1 in dict_final1[binx].keys():
                        m1 += len(dict_final1[binx][key1].keys())
                    b1 = np.zeros((m1,n_jck+1))

                    count = 0
                    for key1 in dict_final1[binx].keys():
                        for key2 in dict_final1[binx][key1].keys():
                            b1[count,1:]= dict_final1[binx][key1][key2]['Pj']
                            b1[count,0]= dict_final1[binx][key1][key2]['P']
                            count+=1
                            ll.append('{0}_{1}'.format(key1,key2))
                if fact:
                    b1 = b1
                    b_dict1 = covariance_jck(b1[:,1:]*fact[binx], n_jck, type_cov = type_c)
                else:
                    b_dict1 = covariance_jck(b1[:,1:], n_jck, type_cov = type_c)
                
                
                ax[i].set_ylabel('bin : {0} \n P/P_err'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                
                ax[i].scatter(np.arange(len( b1[:,0])), b1[:,0]/b_dict1['err'],label = 'b/b_err')
                import copy 
                points2 = copy.copy(b1[:,0])
                
                xxt = np.arange(len( b1[:,0]))
                plt.xticks(xxt,ll, rotation='90')
                ax[0].text(0.4,1.05 , '{0}'.format(label),transform=ax[0].transAxes)
                
                w = b1[:,0]
                cns = np.zeros(len(b1[:,0]))
                from scipy import linalg
                inv_cov = linalg.inv(b_dict1['cov'])
                N_p = n_jck
                p_p = len(cns)
                f_hartlap = (N_p-1.)/(N_p-p_p-2.)
                #print f_hartlap
                chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/(len(cns)*f_hartlap)

                ax[i].text(1.1 , 0.5 , 'chi/dof : {0:.2f}; sigma : {1:.2f}'.format(chi2red,p2s(1 - stats.t.cdf(chi2red,df=len(cns)))),transform=ax[i].transAxes)
                plt.xlim([0.,len( b1[:,0])+1])
                plt.ylim([0.,3.5])
                
                # plot above lims
                
                mm = b1[:,0]>3.5

                for hh in range(len(mm[mm])):
                    ax.annotate( '{0:.2f}'.format(ll[mm][hh]), ((xxt[mm][hh]-0.5),3.2),transform=ax.transAxes)
            
            
                if len(mm[mm])>0:
                    points2[mm] = 3.4
                    plt.scatter(xxt[mm],points2[mm], color = 'b', marker = '^')

            
        fig.subplots_adjust(wspace=0.,hspace=0.)
        bd.update({binx:b_dict1})
    return bd, ll

import healpy as hp


def convert_to_pix_coord(ra, dec, nside=1024):
    """
    Converts RA,DEC to hpix coordinates
    """

    theta = (90.0 - dec) * np.pi / 180.
    phi = ra * np.pi / 180.
    pix = hp.ang2pix(nside, theta, phi, nest=False)

    return pix
print ('done') 


# compile also this !


import matplotlib.pyplot as plt
import pyfits as pf
import healpy as hp
import copy
import scipy
from scipy import stats
from scipy.optimize import curve_fit
from scipy import linalg

def covariance_scalar_jck(TOTAL_PHI,jk_r, type_c = 'jackknife'):

  #  Covariance estimation
  if type_c == 'jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_c=='bootstrap':
      fact=1./(jk_r)
        
  average=0.
  cov_jck=0.
  err_jck=0.


  for kk in range(jk_r):
    average+=TOTAL_PHI[kk]
  average=average/(jk_r)

  for kk in range(jk_r):
    #cov_jck+=TOTAL_PHI[kk]#*TOTAL_PHI[kk]

    cov_jck+=(-average+TOTAL_PHI[kk])*(-average+TOTAL_PHI[kk])


  err_jck=np.sqrt(cov_jck*fact)


  #average=average*(jk_r)/(jk_r-1)
  return {'cov' : cov_jck*fact,
          'err' : err_jck,
          'mean': average}



def covariance_jck(TOTAL_PHI,jk_r,type_cov):
  if type_cov=='jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_cov=='bootstrap':
      fact=1./(jk_r)
  #  Covariance estimation

  average=np.zeros(TOTAL_PHI.shape[0])
  cov_jck=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  err_jck=np.zeros(TOTAL_PHI.shape[0])

  for kk in range(jk_r):
    average+=TOTAL_PHI[:,kk]
  average=average/(jk_r)

 # print average
  for ii in range(TOTAL_PHI.shape[0]):
     for jj in range(ii+1):
          for kk in range(jk_r):
            cov_jck[jj,ii]+=TOTAL_PHI[ii,kk]*TOTAL_PHI[jj,kk]

          cov_jck[jj,ii]=(-average[ii]*average[jj]*jk_r+cov_jck[jj,ii])*fact
          cov_jck[ii,jj]=cov_jck[jj,ii]

  for ii in range(TOTAL_PHI.shape[0]):
   err_jck[ii]=np.sqrt(cov_jck[ii,ii])
 # print err_jck

  #compute correlation
  corr=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  for i in range(TOTAL_PHI.shape[0]):
      for j in range(TOTAL_PHI.shape[0]):
        corr[i,j]=cov_jck[i,j]/(np.sqrt(cov_jck[i,i]*cov_jck[j,j]))

  average=average*fact
  return {'cov' : cov_jck,
          'err' : err_jck,
          'corr':corr,
          'mean':average}


def compute_area(map, nside):
    map1 = copy.copy(map)
    #map1[map1 > 0] = 1
    area = np.sum(map1) * 1.0 * hp.pixelfunc.nside2pixarea(nside, degrees=True)
    return area

def find_limits(mapp, fract):
    min_ms = min(mapp)
    max_ms = max(mapp)
    his,arr=np.histogram(mapp, bins = np.linspace(min_ms,max_ms,1000))

    tot = np.sum(his)
    low = tot*(fract)/100.
    mute = 0.
    ll =0
    while ll <(len(his)):
        mute += his[ll]
        
        if mute> low:
            min_s = arr[ll]
            ll = (len(his))
        ll += 1    
    ll = 0
    mute = 0.
    while ll <(len(his)):
        mute += his[len(his)-ll-1]

        if mute> low:
            max_s = arr[len(his)-ll]
            ll = (len(his))
        ll += 1    
        
    return min_s,max_s

# linear fit
def fitting_linear(x, para0,para1):
    return para0*(x)+para1

def fitting_2nd(x, para0,para1,para2):
    return para2*(x*x)+para0*(x)+para1
def fitting_3rd(x, para0,para1,para2,para3):
    return para3*(x*x*x)+para2*(x*x)+para0*(x)+para1

from scipy import stats
from scipy import special


def p2s(p):
    return np.sqrt(2.)*special.erfinv((1-p))
def s2p(sigma):
    return 1-special.erf(sigma/np.sqrt(2.))



def show_only_plots(name_folder,add_label,kE,kE_label,weight,systematic_maps,systematics_dict,info,fract_limits, len_hist = 10, cov_ext = False,fit ='linear'):
    
    for key in systematic_maps.keys():
      if key != kE_label:
        
        out_file = name_folder
        make_plot(out_file,add_label,systematics_dict[key],systematic_maps[key]['title'],kE_label,systematic_maps[key]['map'],systematic_maps[key]['fractional'],info,fract_limits, len_hist = len_hist,cov_ext = cov_ext,fit=fit)
        #make_plot(out_file,add_label,dict_output,systematic_maps[key]['title'],kE_label,systematic_maps[key]['map'],systematic_maps[key]['fractional'],info, fract_limits,len_hist = len_hist,cov_ext = cov_ext,fit=fit,dict_output_2=dict_output_2)
        


from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})



def def_max_val(dictu,lim = 2.0, fit = 'linear'):
    F01 = dictu['F01']
       
    if fit == 'linear':
        h = 0.
        m = 0.
    elif fit == 'quadratic':
        h = 1.
        m = 0.
    elif fit == 'cubic':
        h = 1.
        m = 1.
    
    F02 = dictu['F02']*h
    F03 = dictu['F03']*m
    
    F13 = dictu['F13']*m
    F12 = dictu['F12']*h
    F23 = dictu['F23']*m
    

    
    v = 0
    value = 0
    
    if F01>lim:
        value = F01
        v = 1
        if F02>lim and F12>lim:
            value = F02
            v = 2
            if F03>lim and F23>lim:
                value = F03
                v = 3
        else:
            if F03>lim and F13>lim:
                value = F03
                v = 3
            
    else:
        if F02>lim:   
            value = F02
            v = 2    
            if F03>lim and F23>lim:
                value = F03
                v = 3            
        else:
            if F03>lim:
                value = F03
                v = 3            
    return v,value

def correct_weight(count,systematics_dict_total,systematic_maps, info,weights , lim = 2.0, keys = False, fit = 'linear'):
    #order weights!
    from scipy import linalg     
    if 1==1:
        key = '_w{0}'.format(count)
        if keys:
            syst_labb = keys
        else:
            syst_labb = systematics_dict_total[key].keys()

        points = np.zeros(len(syst_labb))
        vv_arr = np.zeros(len(syst_labb))
        label = []
        string = []
     
        for num, key2 in enumerate(syst_labb):
            if key2 != 'mask_tot':
                points[num] = def_max_val( systematics_dict_total[key][key2],lim = lim,fit=fit)[1]
                vv_arr[num] = def_max_val( systematics_dict_total[key][key2],lim = lim,fit=fit)[0]
                label.append(key2)
                vv=''
                if vv_arr[num] == 1:
                    vv = 'linear'
                if vv_arr[num] == 2:
                    vv = 'quadratic'
                if vv_arr[num] == 3:
                    vv = 'cubic'
                
                string.append('{0} ({1:.2f}, {2}), '.format(key2,points[num],vv))
        
        ind = np.argsort(points)

        mask = points[ind] > lim
        weights_add = np.zeros(len(systematic_maps[key2]['map']))
        if len(points[ind][mask])>0:
            

            syst_l = [str(x) for x in np.array(label)[ind][mask]][::-1][0]
            nn = int(vv_arr[ind][mask][::-1][0])
            
            print (' CORRECTION: {0} {1}'.format(syst_l,nn))
            string = ''.join([str(x) for x in np.array(string)[ind][mask]][::-1])
        
            print ("systematics > lim : ",string)#[mask].flatten()
        
            # define weights:
            import copy
            #weights_add = np.zeros(len(systematic_maps[key2]['map']))
        
            
            
            fractional = systematic_maps[syst_l]['fractional']
            x = systematic_maps[syst_l]['map'][info['mask_sims']]
            ave  = x.mean()                      
            if fractional:
                x -= ave
                x /= ave
            else:
                x -= ave
                  
            #mask_sys = (systematic_maps[syst_l]['map'] > systematics_dict_total[key][key2]['min_ss']) & (systematic_maps[syst_l]['map'] < systematics_dict_total[key][key2]['max_ss'])
            
            mask_nn = (systematic_maps[syst_l]['map'][info['mask_sims']] > systematics_dict_total[key][syst_l]['min_ss']) & (systematic_maps[syst_l]['map'][info['mask_sims']] < systematics_dict_total[key][syst_l]['max_ss'])
            
            mutenn = np.zeros(len(systematic_maps[syst_l]['map'][info['mask_sims']]))
            print ("efficiency: ", (np.float(len(x[mask_nn])))/(np.float(len(x))),len(x[mask_nn]))
            #print len(x[mask_nn]), len(weights_add[info['mask_sims']&mask_sys])
            if nn ==3:
                #cubic correction
                print ('cubic')
                dict_out = systematics_dict_total[key][syst_l]#['params_bts3']
                #the corrections should not correct for the shift
                weights_add[info['mask_sims']][mask_nn] = fitting_3rd(x[mask_nn],dict_out['params_bts3'][0],dict_out['params_bts3'][1],dict_out['params_bts3'][2],dict_out['params_bts3'][3])  

            if nn ==2:
                #cubic correction
                print ('quadratic')
                dict_out = systematics_dict_total[key][syst_l]#['params_bts3']

                nw = fitting_2nd(x[mask_nn],dict_out['params_bts2'][0],dict_out['params_bts2'][1],dict_out['params_bts2'][2])  
                mutenn[mask_nn] = nw
                weights_add[info['mask_sims']] = mutenn
             
            if nn ==1:
                #cubic correction
                print ('linear')
                dict_out = systematics_dict_total[key][syst_l]#['params_bts3']
                #the corrections should not correct for the shift
                #print max(fitting_linear(x[mask_nn],dict_out['params_bts'][0],dict_out['params_bts'][1])  ), min(fitting_linear(x[mask_nn],dict_out['params_bts'][0],dict_out['params_bts'][1])  )
                nw= fitting_linear(x[mask_nn],dict_out['params_bts'][0],dict_out['params_bts'][1])  
                mutenn[mask_nn] = nw
                weights_add[info['mask_sims']] = mutenn
                print (dict_out['params_bts'][0],dict_out['params_bts'][1],len(weights_add[info['mask_sims']]),len(weights_add[info['mask_sims']][mask_nn]), max(weights_add[info['mask_sims']][mask_nn]))
            mute_l = False
        else:
            mute_l = True
        return weights_add,mute_l

    
def compute_fisher(chi_n,chi_fit,N,p_null,p_fit):
    #F = -((N-p_fit)/(p_null-p_fit))*(chi_n/chi_fit-1.),
    F = ((chi_n-chi_fit)/(p_fit-p_null))/((chi_fit)/(N-p_fit))
    v1 = -(p_null-p_fit)
    v2 = (N-p_fit)
    return F, v1,v2


def make_plot(save_output_folder,add_label,dict_out, title, kE_label, sys_map1,fractional,info,fract_limits,len_hist=10,cov_ext = False,fit='linear', dict_output_2 = False):
    
    print ("Jackknife: {0}".format(not cov_ext))
    plt.figure(1,figsize=(6,5))
    ax = plt.subplot(211)
    if fractional:
        plt.xlabel('{0} [{1:.2f},{2:.2f}] (fractional)'.format(title, dict_out['min_s'],dict_out['max_s']))
    else:
        plt.xlabel('{0}[{1:.2f},{2:.2f}] (mean sub)'.format(title, dict_out['min_s'],dict_out['max_s']))
    plt.title('systematic map : {0}'.format(title))
    plt.ylabel('counts')
    
    plt.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        
    sys_map1 = sys_map1[info['mask_sims']]
    ave = sys_map1.mean()
    if fractional:

        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
    plt.hist(sys_map1, bins = dict_out['bins'])
    
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width* 1.4 , box.height])   
    
    plt.figure(1)
    ax = plt.subplot(212)
    plt.title('systematic map : {0}'.format(title))
    if fractional:
        plt.xlabel('{0} (mean subtracted, fractional)'.format(title))
    else:
        plt.xlabel('{0} (mean subtracted)'.format(title))        
    plt.ylabel('<{0}>'.format(kE_label))
    
    bins = dict_out['bins']
    x = dict_out['x']
    
    plt.errorbar(x, dict_out['w'],dict_out['dict_cov']['err'])
    plt.tight_layout()
    
    w = dict_out['w']
    cns = np.ones(len(w))* dict_out['mean shear']
    inv_cov = linalg.inv(dict_out['dict_cov']['cov'])
    chi2red = (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/len(w)
    

    #plt.plot(x, fitting_linear(x,dict_out['params'][0],dict_out['params'][1]),label='fit')
              
    if cov_ext:         
        label_bst ="_boot"
        inv_cov = linalg.inv(dict_out['dict_cov_bts']['cov'])
    else:
        label_bst =""
        inv_cov = linalg.inv(dict_out['dict_cov']['cov'])
    if 1==1:
        
        w = dict_out['w']
        cns = np.ones(len(w))* dict_out['mean shear']
        
        chi2red = (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/len(w)
    
    
            
        scipy.stats.f.cdf
        xx = 1.05
        yy = 1.10

        chi = dict_out['chi2red'+label_bst+'_null'] 
        df1 = len(dict_out['x']) - 1
        chi=chi/df1
        p = 1 - stats.t.cdf(chi,df=df1)
        #print chi,df,p
        #print p2s(p)
        plt.text(xx, yy+1.35, 'summary stat', transform=ax.transAxes)
        plt.text(xx, yy+1.2, '[method 1]', transform=ax.transAxes)

        plt.text(xx, yy+0.95, 'sigma_null: {0:.3f}'.format(p2s(p) ), transform=ax.transAxes)

        chi = dict_out['chi2red'+label_bst] 
        df2 = len(dict_out['x']) - 2
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        #print chi,df,p
        #print p2s(p)
        
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'_null'] ,dict_out['chi2red'+label_bst] ,len(dict_out['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)


        plt.text(xx, yy+0.8, 's_1: {0:.2f}'.format(p2s(p),p2s(pf)), transform=ax.transAxes)
        plt.text(xx, yy+0.65, 'F: [{1:.2f}]'.format(p2s(p),p2s(pf)), transform=ax.transAxes)
  


        chi = dict_out['chi2red'+label_bst+'2'] 
        df2 = len(dict_out['x']) - 3
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        
        
        
        #print chi,df,p
        #print p2s(p)
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'_null'] ,dict_out['chi2red'+label_bst+'2'] ,len(dict_out['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)

        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst] ,dict_out['chi2red'+label_bst+'2'] ,len(dict_out['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        plt.text(xx, yy+0.5, 's_2: {0:.3f}'.format(p2s(p),p2s(pf),p2s(pf2)), transform=ax.transAxes)
        plt.text(xx, yy+0.35, 'F: [{1:.2f},{2:.2f}]'.format(p2s(p),p2s(pf),p2s(pf2)), transform=ax.transAxes)
   
        
        chi = dict_out['chi2red'+label_bst+'3'] 
        df2 = len(dict_out['x']) - 4
        chi=chi/df2
        p = 1 - stats.t.cdf(chi,df=df2)
        #print chi,df,p
        #print p2s(p)
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'_null'] ,dict_out['chi2red'+label_bst+'3'] ,len(dict_out['x']),1,2)  
        pf = 1-scipy.stats.f.cdf(F,v1,v2)


        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst] ,dict_out['chi2red'+label_bst+'3'] ,len(dict_out['x']),1,2)  
        pf2 = 1-scipy.stats.f.cdf(F,v1,v2)
 
        F,v1,v2 = compute_fisher(dict_out['chi2red'+label_bst+'2'] ,dict_out['chi2red'+label_bst+'3'] ,len(dict_out['x']),1,2)  
        pf3 = 1-scipy.stats.f.cdf(F,v1,v2)
        
        plt.text(xx, yy+0.2, 's_3: {0:.3f}'.format(p2s(p),p2s(pf),p2s(pf2),p2s(pf3)), transform=ax.transAxes)
        plt.text(xx, yy+0.05, 'F: [{1:.2f},{2:.2f},{3:.2f}]'.format(p2s(p),p2s(pf),p2s(pf2),p2s(pf3)), transform=ax.transAxes)
                
                 
        chi = dict_out['chi2red'+label_bst+'_null']/df1 - dict_out['chi2red'+label_bst]/df2
        df = 1
        p = 1 - stats.t.cdf(chi,df=df)
        
        '''
        plt.text(xx, yy+0.35, 'sigma_H1 (fit): {0:.3f}'.format(p2s(p)), transform=ax.transAxes)

        plt.text(xx, yy+0.2, 'sigma_68 : {0:.3f}'.format(dict_out['chi68_new']), transform=ax.transAxes)
        plt.text(xx, yy+0.05, 'sigma_68_ratio : {0:.3f}'.format(p2s(p)/dict_out['chi68_new']), color='r',transform=ax.transAxes)
        plt.text(xx, yy-0.1, 'sigma_68 (nd) : {0:.3f}'.format(dict_out['chi68_new1']),transform=ax.transAxes)

        chi = dict_out['chi2red_boot_null'] - dict_out['chi2red_boot']
        df = 1
        p = 1 - stats.t.cdf(chi,df=df)
        plt.text(xx, yy-0.25, 'sigma_68_ratio (nd: {0:.3f}'.format(p2s(p)/dict_out['chi68_new1']), color='r',transform=ax.transAxes)


        '''
        plt.text(xx, 0.65, 'chi2 null: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst+'_null'],len(dict_out['x']) - 1), transform=ax.transAxes)
        plt.text(xx, 0.5, 'chi2 fi1t: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst],len(dict_out['x']) - 2), transform=ax.transAxes)
        plt.text(xx, 0.35, 'chi2 fit2: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst+'2'],len(dict_out['x']) - 2), transform=ax.transAxes)
        plt.text(xx, 0.2, 'chi2 fit3: {0:.2f} (dof: {1})'.format(dict_out['chi2red'+label_bst+'3'],len(dict_out['x']) - 2), transform=ax.transAxes)
      
        '''
        plt.text(xx, 0.35, 'diff chi2: {0:.4f}'.format(dict_out['chi2red_diff_boot'],dict_out['chi68']), transform=ax.transAxes)
        #plt.text(xx, 0.2, 'diff chi2/dof 68:  {0:.4f}'.format(dict_out['chi68']), transform=ax.transAxes)
        #plt.text(xx, 0.05, 'diff chi2/ diff chi2 68: {0:.4f}'.format(dict_out['chi2red_diff_boot68']),color='r', transform=ax.transAxes)
        '''
        plt.text(0.05, 0.05, 'slope [method_1]= {0:.4f} +- {1:.4f}'.format(dict_out['params_bts'][0],np.sqrt(dict_out['params_cov_bts'][0,0])), transform=ax.transAxes)
        
        if dict_output_2:
            plt.text(0.05, 0.2, 'slope [method_2]= {0:.4f} +- {1:.4f}'.format(dict_output_2['b_arr'][0],dict_output_2['b_err']), transform=ax.transAxes)
            
        plt.plot(x, fitting_2nd(x,dict_out['params_bts2'][0],dict_out['params_bts2'][1],dict_out['params_bts3'][2]),label='fit 2nd')
        plt.plot(x, fitting_3rd(x,dict_out['params_bts3'][0],dict_out['params_bts3'][1],dict_out['params_bts3'][2],dict_out['params_bts3'][3]),label='fit 3rd')

        plt.plot(x, fitting_linear(x,dict_out['params_bts'][0],dict_out['params_bts'][1]),label='fit')
        plt.plot(x, cns,label='null')
        



        
        
        
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width* 1.4, box.height])




    plt.ylim([-0.004,0.004])
    plt.savefig(save_output_folder+kE_label+'_'+title+' '+str(fract_limits)+'_'+str(len_hist)+'_'+add_label+'_'+fit+'.png')   

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95),
          ncol=4, fancybox=True)
    #plt.legend()
    plt.show()


def plot_resume(systematics_dict_total,z_min,z_max,nside,kE_label,d68 = False,type_c='jackknife') :
 from scipy import linalg     
 for key in systematics_dict_total.keys():
    fig = plt.figure(figsize=(8, 10)) 
    
    
    
    if type_c == 'jackknife':
        d68 = False
        
    syst_labb = systematics_dict_total[key].keys()
    if d68:
        
        ax = plt.subplot(313)
        
        plt.title('sys map impact wrt null hypotesys - {0}   bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    
        points = np.zeros(len(syst_labb))
        points_2 = np.zeros(len(syst_labb))  
    
        points_err = np.zeros(len(syst_labb))
        x = np.array(range(len(syst_labb)))+1
        label = []

        plt.show
        for num, key2 in enumerate(syst_labb):
          if key2 != 'mask_tot':
        
            points_2[num] = systematics_dict_total[key][key2]['chi2red_diff_boot68']
            key2 = '_'+syst_labb[num]
            key2 =  '   '+key+key2
            label.append(key2)
    
        plt.xticks(x, label, rotation='90')
        plt.xlim([0.,len(syst_labb)+1])
        plt.ylim([0,6.])
        mm = points_2>6

        for hh in range(len(mm[mm])):
            ax.annotate( '{0:.2f}'.format(points_2[mm][hh]), ((x[mm][hh]-0.5),5.4),transform=ax.transAxes)
            
            #    plt.text(x[mm[hh]], points_2[mm[hh]], 'cc')#,transform=ax.transAxes)
        plt.scatter(x,points_2,label = '<e>_0 = <e>', color = 'b')
        if len(mm[mm])>0:
            points_2[mm] = 5.9
            plt.scatter(x[mm],points_2[mm], label = '<e>_0 = 0.', color = 'b', marker = '^')

            
        plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
        plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
        plt.ylabel('delta chi2/delta chi2 (68)')
    

    
    
        
        
        
        
    plt.subplot(312)#,figsize=(10,20))
    #plt.figure(1,figsize=(10,100))
    

    points = np.zeros(len(syst_labb))
    points_2 = np.zeros(len(syst_labb))
    x = np.array(range(len(syst_labb)))+1
    label = []


  
    for num, key2 in enumerate(syst_labb):
     if key2 != 'mask_tot':
        points[num] = np.sqrt((systematics_dict_total[key][key2]['params_bts'][0]/np.sqrt(systematics_dict_total[key][key2]['params_cov'][0,0]))**2.)
        points_2[num] = systematics_dict_total[key][key2]['mean shear'] 
        if num == 0:
            mask_tot = systematics_dict_total[key][key2]['mask_pixel']
        else:
            mask_tot = mask_tot ^ systematics_dict_total[key][key2]['mask_pixel']
        key2 = '_'+syst_labb[num]
        key2 = '   '+ key+key2
        label.append(key2)
    

    plt.xticks(x, label,  rotation='90')
    plt.xlim([0.,len(syst_labb)+1])
    plt.ylim([0.,8])
    plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
    plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
    plt.title('Linear fit slope: {0}   bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    plt.scatter(x,points)
    plt.ylabel('|b|/b_err')
    #fig.tight_layout()
    p_value = 1 - stats.chi2.cdf(np.sum(points**2), len(points))
    #plt.text(0, 7, '  chi2/dof = {0:.2f}'.format(np.sum(points)/len(points)))
    #plt.text(0, 5.6, ' Average mean shear = {0:.5f}'.format(np.mean(points_2)))
    #plt.text(0, 4, '  p-value = {0:.2f}'.format(p_value))
    tot_mask = compute_area(mask_tot,nside)
    #plt.text(10, 7, 'area masked = {0:.1f} deg2'.format(tot_mask))
    #systematics_dict_total[key].update({'mask_tot' : tot_mask})
    #plt.tight_layout()
    
    

    ax = plt.subplot(311)
    #plt.figure(1,figsize=(10,20))
    plt.title('Chi2 null hypothesis {0} - bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    
    points = np.zeros(len(syst_labb))
    points_2 = np.zeros(len(syst_labb))  
    
    points_err = np.zeros(len(syst_labb))
    x = np.array(range(len(syst_labb)))+1
    label = []

    plt.show
    for num, key2 in enumerate(syst_labb):
      if key2 != 'mask_tot':
        
        w = systematics_dict_total[key][key2]['w']
        cns = np.ones(len(w))* systematics_dict_total[key][key2]['mean shear']
        inv_cov = linalg.inv(systematics_dict_total[key][key2]['dict_cov']['cov'])
        points[num] = (np.matmul(w,np.matmul(inv_cov,w)))/len(w)
        if type_c =='jackknife':
            points_2[num] = systematics_dict_total[key][key2]['chi2red']/(len(systematics_dict_total[key][key2]['x']) - 2)
        else:
            points_2[num] = systematics_dict_total[key][key2]['chi2red_boot']/(len(systematics_dict_total[key][key2]['x']) - 2)
        key2 = '_'+syst_labb[num]
        key2 =  '   '+key+key2
        label.append(key2)
    #print points_2
    
    plt.xticks(x, label, rotation='90')
    plt.xlim([0.,len(syst_labb)+1])
    plt.ylim([0,4.])
    plt.scatter(x,points_2,label = '<e>_0 = <e>')
    plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
    plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
    plt.ylabel('chi2/dof')


    plt.show()
        
        
print ('done')


import skymapper as skm
import matplotlib.cm as cm

def getHealpixCoords(pixels, nside, nest=False):
    # convert healpix cell indices to center ra/dec
    import healpy as hp
    theta, phi = hp.pix2ang(nside, pixels, nest=nest)
    return phi * 180. / np.pi, 90 - theta * 180. / np.pi
    

def skm_plot(mapa, mask, sep=15, ra = None, dec = None, richness = None, title = None, add_rmp = True, small_scale = False, cb_label = 'map value', x_size = 6, y_size = 6, vmin = -0.01, vmax = 0.01):

    cmap = cm.RdYlBu_r
    fig = plt.figure(figsize=(x_size,y_size))
    ax = fig.add_subplot(111, aspect='equal')

    
    
    reticule = sep
    nside = hp.pixelfunc.npix2nside(len(mapa))
    pixels = np.array(range(len(mapa)))
    pixels = pixels[~mask]
    ra_, dec_ = getHealpixCoords(pixels, nside)
    
    #plt.scatter(ra_,dec_)
    ##ig = plt.figure(figsize=(15,10))
    #ax = fig.add_subplot(111, aspect='equal')
    #fig.tight_layout()
    #proj = skm.plotMap(ra_, dec_, mapa[pixels], sep=sep, ax=ax, cb_label= cb_label, cmap='bwr')

    # setup map: define AEA map optimal for given RA/Dec
    proj = skm.createConicMap(ax, ra_, dec_, proj_class=skm.AlbersEqualAreaProjection)

    # add lines and labels for meridians/parallels (separation 5 deg)
    
    parallels = np.arange(0. ,360., reticule)
    meridians = np.arange(-90., 90., reticule)
    skm.setMeridianPatches(ax, proj, meridians, linestyle='-', lw=0.5, alpha=0.3, zorder=10)
    skm.setParallelPatches(ax, proj, parallels, linestyle='-', lw=0.5, alpha=0.3, zorder=10)
    skm.setParallelLabels(ax, proj, parallels, loc="bottom", fmt=skm.pmDegFormatter)
    skm.setMeridianLabels(ax, proj, meridians, loc="left", fmt=skm.pmDegFormatter)
    


    # convert to map coordinates and plot a marker for each point
    x,y = proj(ra_, dec_)
    marker = 's'
    markersize = skm.getMarkerSizeToFill(fig, ax, x, y)

    sc = ax.scatter(x,y, c=mapa[pixels], edgecolors='None', marker=marker, s=markersize, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True, zorder=1)
    
    
    # overplot with another data set
    # here clusters [not implemented]
    if add_rmp:
        x,y = proj(ra, dec)
        ax.scatter(x,y, c='None', edgecolors='k', linewidths=1, s=richness, marker='o', zorder=3)

    if title:
        plt.title(title)

    # add colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes("right", size="5%", pad=0.3)
    cb = fig.colorbar(sc, cax=cax )

    if cb_label != None:
        cb.set_label(cb_label)
    ticks = np.linspace(vmin, vmax, 5)
    cb.set_ticks(ticks)
    if small_scale:
        cb.set_ticklabels([('%.3f' % (t))   for t in ticks])       
    else:
        cb.set_ticklabels([('%.1f' % (t))   for t in ticks])
    cb.solids.set_edgecolor("face")

    # show (and save) ...
    fig.tight_layout()
    fig.show()
    #fig.savefig(imagefile)    
       
        
from scipy.spatial import distance
import os
import numpy as np
import pickle

def update_progress(progress,elapsed_time=0,starting_time=0):

    import time
    import timeit
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))



def dist_cent_2(ra1,dec1,ra2,dec2):

            todeg = np.pi/180.
            ra1 = ra1*todeg
            ra2 = ra2*todeg
            dec1 = dec1*todeg
            dec2 = dec2*todeg

            cos = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
            return np.arccos(cos)/todeg



def save_obj( name,obj ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)






def estimator(w_estimator,DD,DR,RD,RR):

    if w_estimator == 'LS':
        #print (w_estimator)
        results = (DD-DR-RD+RR)/(RR)
    elif w_estimator == 'Natural':
        results = (DD - RR) / (RR)
    elif w_estimator == 'Hamilton':
        results = (DD * RR - RD * DR) / (RD * DR)
    elif w_estimator == 'Natural_noRu':
        results = (DD - DR) / DR
    elif w_estimator == 'Natural_noRr':
        results = (DD - RD) / RD
    return results


import time
import timeit
import treecorr
class Jack(object):
    def __init__(self,w_estimator,conf,pairs_sel,rau,decu,rar,decr,jku,jkr,weight_u,weight_r, dist_path,centers=None, njk=None,verbose=False):
        
        self.jackknife_speedup= True
        self.w_estimator=w_estimator
        self.conf = conf
        
        self.rau=rau
        self.decu=decu

        self.rar=rar
        self.decr=decr

        
        self.jku=jku
        self.jkr=jkr


     
        self.weight_u=weight_u
        self.weight_r=weight_r

        

        self.centers = centers
        self.njk = njk

        self.dist_path = dist_path
        self.pairs_select=pairs_sel

    
    def KKCorrelation(self):
        self.prolog()
        pairs= self.epilog()

        return pairs


    def convert_units(self):
        """Change to 'degrees' the units if necessary.
        """
        if 'sep_units' in self.conf.keys():
            un = self.conf['sep_units']
            if un == 'arcmin':
                todeg = 1./60.
            elif un == 'arcsec':
                todeg = 1./60.**2.
            else:
                todeg = 1.
        else:
            todeg = 1.
        return todeg


    def max_distance_region(self):
        cnt = self.centers
        max_dist_region=load_obj(self.dist_path)
        self.max_dist_region=max_dist_region



    def distance(self):
        """Finds the minimum distance to a center for each center.
           Fixes double of this distance as the criteria for not considering correlations,
           which is a conservative choice. This distance has to be at least 4 times the
           maximum angular separation considered. Centers beyond this distance will not be
           considered in the correlations.
        """

        # Find the minimum distance to a center for each center.

        cnt = self.centers
        dist = np.array([np.sort([self.dist_cent(cnt[i],cnt[j]) for i in range(len(cnt))])[1] for j in range(len(cnt))])
        dist = (dist)*2.

        todeg = self.convert_units()

        if 'max_sep' in self.conf.keys():
            max_sep = self.conf['max_sep'] * todeg
        else:
            raise NotImplementedError("Make use of 'max_sep' in configuration.")

        # Check that the distance is at least 4 times the maximum angular separation.
        self.center_min_dis = np.array( [ 4.*max_sep if x < 4.*max_sep else x for x in dist] )



    def dist_cent(self, a, b):
        """Angular distance between two centers (units: degrees). Makes use of spherical law of cosines.
        """
        todeg = np.pi/180.
        a = a*todeg
        b = b*todeg
        cos = np.sin(a[1])*np.sin(b[1]) + np.cos(a[1])*np.cos(b[1])*np.cos(a[0]-b[0])
        return np.arccos(cos)/todeg


    def cond(self, i, j):

        """Return the maximum conditional distance for a pair of centers to
           determine if the correlation should be computed.
        """

        con = self.center_min_dis
        return max(con[i], con[j])


    def collect(self, pairs):


        shape = (self.njk, self.conf['nbins'])
        
        DD_a,DD_c = np.zeros(shape),np.zeros(shape)
        normm_a,normm_c = np.zeros(shape),np.zeros(shape)
        DD = np.zeros(self.conf['nbins'])
        normm = np.zeros( self.conf['nbins']) 

        if self.jackknife_speedup:
            for n in range(self.conf['nbins']):
                DD[n] = np.sum(pairs[:,0,0,n]) + 0.5 * np.sum(pairs[:,1,0,n])
                normm[n] = np.sum(pairs[:,0,1,n]) + 0.5 * np.sum(pairs[:,1,1,n])    

                # only the jk *****************
                DD_a[:,n] = (pairs[:,0,0,n])
                DD_c[:,n]  =  0.5 * (pairs[:,1,0,n])

                normm_a[:,n] = (pairs[:,0,1,n])
                normm_c[:,n]  =  0.5 * (pairs[:,1,1,n])
        return DD,DD_a, DD_c, normm, normm_a, normm_c


    def parallel(self, i, j):

        [[ra_a, dec_a, jk_a], [ra_b, dec_b, jk_b]] = self.info

        # Create the Catalog object. One for each jackknife region.
        try:
            mask=np.in1d(jk_a,i)
            cat_a = treecorr.Catalog(ra=ra_a[mask], dec=dec_a[mask], ra_units='deg', dec_units='deg',k=self.weight_u[mask])
        except RuntimeError:
            cat_a = None
        try:
            mask=np.in1d(jk_b,j)
            cat_b = treecorr.Catalog(ra=ra_b[mask], dec=dec_b[mask], ra_units='deg', dec_units='deg',k=self.weight_r[mask])
        except RuntimeError:
            cat_b = None

        kk = treecorr.KKCorrelation(self.conf)
        kk.process(cat_a, cat_b)

        return [kk.xi*kk.weight,kk.weight]

    
    def dist_cent_2(self,ra1,dec1,ra2,dec2):

            todeg = np.pi/180.
            ra1 = ra1*todeg
            ra2 = ra2*todeg
            dec1 = dec1*todeg
            dec2 = dec2*todeg

            cos = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
            return np.arccos(cos)/todeg

    def prolog(self):

        #insert njk


        njk = self.njk

        cnt = self.centers

 
        self.distance()
        self.max_distance_region()

        ra_a, dec_a,  jk_a, NA = self.rau, self.decu,self.jku,np.sum(self.weight_u)
        ra_b, dec_b,  jk_b, NB = self.rar, self.decr,self.jkr,np.sum(self.weight_r)
 

        self.info = [[ra_a, dec_a, jk_a],[ra_b, dec_b, jk_b]]

        if 'nbins' in self.conf:
                shape = (2,self.conf['nbins'])

        else:
            raise Exception('Not implemented yet. Please use nbins in config.')
            #shape = (4, self.conf['bin_size'])

        cnt = self.centers

        t1=time.time()

        pairs = [[np.zeros(shape) for i in range(njk)] for j in range(njk)]
   
        pairs_ring = [[np.zeros(shape) for i in range(2)] for j in range(njk)]
        a=np.concatenate((np.array([[(i,j) for i in range(njk)] for j in range(njk)])))

        todeg = self.convert_units()
        max_sep = self.conf['max_sep'] * todeg

        sel = np.array([ max([0.,(self.dist_cent(cnt[i],cnt[j]) - (self.max_dist_region[i]+self.max_dist_region[j]))]) < (3. * max_sep )for (i,j) in a])
        b = a[sel]


        def fun_speedup(othersample,otehrsample1,jackk):
            try:
                pairsCC1=self.parallel(othersample,[jackk])
                pairsCC2=self.parallel([jackk],othersample1)
                pairs_auto=self.parallel([jackk],[jackk])
                for prs in range(2):
                    pairsCC1[prs]+=pairsCC2[prs]
                pairs_ring[jackk][1] = pairsCC1
                pairs_ring[jackk][0] = pairs_auto
                
            except (RuntimeError, e):
                print (e)
                pairs_ring[jackk][0] = np.zeros(shape)
                pairs_ring[jackk][1] = np.zeros(shape)

        start=timeit.default_timer()


        mute=0

        if self.jackknife_speedup:
            startt = timeit.default_timer()
            for counter,jackk in enumerate(np.unique(b[:,1])):
                mask=(b[:,1]==jackk) & (b[:,0]!=jackk)
                othersample=b[mask,0]
                mask=(b[:,0]==jackk) & (b[:,1]!=jackk)
                othersample1=b[mask,1]
                fun_speedup(othersample,othersample1,jackk)
                mute+=len(othersample)+1
                
                #update_progress((float(counter)/len(np.unique(b[:,1]))),timeit.default_timer(),startt)
                if counter == 100:
                    print (counter, timeit.default_timer()-startt)
                if counter == 200:
                    print (counter, timeit.default_timer()-startt)
                if counter == 300:
                    print (counter, timeit.default_timer()-startt)
                if counter == 400:
                    print (counter, timeit.default_timer()-startt)
                if counter == 500:
                    print (counter, timeit.default_timer()-startt)
                if counter == 600:
                    print (counter, timeit.default_timer()-startt)
                if counter == 700:
                    print (counter, timeit.default_timer()-startt)
                if counter == 800:
                    print (counter, timeit.default_timer()-startt)
            self.pairs = np.array(pairs_ring)




    def epilog(self):

        pairs = self.pairs

        DD,DD_a,DD_c, normm, normm_a, normm_c = self.collect(pairs)
        convert=self.convert_units()
        min_sep, max_sep, nedges = self.conf['min_sep']*convert, self.conf['max_sep']*convert, self.conf['nbins']+1
        th = np.linspace(np.log10(min_sep), np.log10(max_sep), nedges)
        theta = 10**np.array([(th[i]+th[i+1])/2 for i in range(len(th)-1)])



        #print (corr.shape,DD.shape,DR_a.shape)#, DD,DR ,RD,RR, DD_a,DR_a ,RD_a,RR_a, DD_c,DR_c ,RD_c,RR_c)

        # ***********************************

        pairs=dict()
        
        pairs.update({'theta':theta})
        pairs.update({'DD':DD})
        pairs.update({'DD_a':DD_a})
        pairs.update({'DD_c':DD_c})
        pairs.update({'normm':normm}) 
        pairs.update({'normm_a':normm_a})
        pairs.update({'normm_c':normm_c})
        
        return pairs

    

# compute maximum distance within jackknife regions *****
def distance_calc(dist_path,ra_m,dec_m,jk_m, njk,centers):

    '''
    based on the maximum angular scale probed by the cross-correlation,
    computes which jackknife regions should be included in the computation
    '''
    
    import timeit
    if not os.path.exists(dist_path+'.pkl'):

        max_dist_region1=np.zeros(njk)

        # convert radec to xyz
        cosdec = np.cos(dec_m)
        aJx_u = cosdec * np.cos(ra_m)
        aJy_u = cosdec * np.sin(ra_m)
        aJz_u = np.sin(dec_m)

        print ('compute maximum distance for each jackknife region')
        start=timeit.default_timer()
        for i in range(njk):
            if len(ra_m[jk_m==i]) ==0 or len(dec_m[jk_m==i])==0:
                max_dist_region1[i,j]=0.
            else:

                ra_c,dec_c=centers[i]

                cosdec = np.cos(dec_c)
                aJx_r = cosdec * np.cos(ra_c)
                aJy_r = cosdec * np.sin(ra_c)
                aJz_r = np.sin(dec_c)
        
                tree_m=spatial.cKDTree(np.c_[aJx_u[jk_m==i], aJy_u[jk_m==i], aJz_u[jk_m==i]])

                max_dist_m,index_dist=tree_m.query([aJx_r,aJy_r,aJz_r],k=len(ra_m[jk_m==i]))

                ra_new=ra_m[jk_m==i]
                dec_new=dec_m[jk_m==i]


                if (len(ra_m[jk_m==i])==1):
                    max_dist_region1[i]=dist_cent_2(ra_c,dec_c,ra_new[index_dist],dec_new[index_dist])
                else:
                    max_dist_region1[i]=dist_cent_2(ra_c,dec_c,ra_new[index_dist[-1]],dec_new[index_dist[-1]])
            update_progress(np.float(i+1)/np.float(njk),timeit.default_timer(),start)
        
        save_obj(dist_path,max_dist_region1)
        


def make_wz_errors(pairs,resampling,number_of_bootstrap,pairs_resampling):
    '''
    take as an input pairs between subregions and creates resmapling and errors.
    '''

    # errors ********************************
    if resampling=='jackknife':
        njk=pairs['DD_a'].shape[0]
    elif resampling=='bootstrap':
        njk=number_of_bootstrap

    DD_j=np.zeros((njk+1,pairs['DD_a'].shape[1]))
    norm_j=np.zeros((njk+1,pairs['DD_a'].shape[1]))
    DD_j[0,:]=pairs['DD']/pairs['normm']
    norm_j[0,:]=pairs['normm']    
         
    for jk in range(njk):
        if resampling=='jackknife':
            
            if pairs_resampling:
                fact=1.
            else:
                fact=2.

            DD_j[jk+1,:]=(pairs['DD'][:]-pairs['DD_a'][jk,:]-fact*pairs['DD_c'][jk,:])
            norm_j[jk+1,:]=(pairs['normm'][:]-pairs['normm_a'][jk,:]-fact*pairs['normm_c'][jk,:])
            

            DD_j[jk+1,:] = DD_j[jk+1,:]/norm_j[jk+1,:]

    return pairs['theta'],DD_j.T,njk




# Initialize ******

def make_plot_corr(title,save_output_folder,kE_label,mute_d):
    plt.figure(figsize=(6,3))
    ax = plt.subplot(111)

    plt.title('systematic map (crosscorr) : {0}'.format(title))
    plt.xlabel('theta (arcmin)'.format(title))        
    plt.ylabel('w(theta)')
    plt.errorbar(mute_d['theta']*60.,mute_d['w'],mute_d['cov']['err'])
    plt.tight_layout()
    plt.text(0.7, 0.1, 'chi2/dof:  {0:.2f}'.format(mute_d['chi2_red']), transform=ax.transAxes, size =12)
    plt.show()     
    plt.savefig(save_output_folder + 'cc_' +kE_label+'_'+ title +'.png')
    

def do_analisys_corr(kE,kE_label,weight,systematic_maps,info,add_label):
    
    systematics_dict_corr = dict()
    for key in systematic_maps.keys():
      if key != kE_label:
        out_file = './output_{0}_{1}/'.format(info['z_min'],info['z_max'])
        if not os.path.exists('./output_{0}_{1}/pairs/'.format(info['z_min'],info['z_max'])):
            os.mkdir('./output_{0}_{1}/pairs/'.format(info['z_min'],info['z_max']))
            
        mute_d = compute_cross_corr(systematic_maps[key]['title'],out_file,kE_label,systematic_maps[key]['map'],kE, info, add_label,  mapa_weight = weight,fractional = systematic_maps[key]['fractional'])
        make_plot_corr(systematic_maps[key]['title'],out_file,kE_label,mute_d)
        systematics_dict_corr.update({systematic_maps[key]['title']  : mute_d}) 
  
    return systematics_dict_corr

def plot_resume_cc(systematics_dict_total,z_min,z_max,kE_label) :
 from scipy import linalg     
 for key in systematics_dict_total.keys():
    plt.figure(1)
    
    plt.subplot(111)
    
    syst_labb = systematics_dict_total[key].keys()


    points = np.zeros(len(syst_labb))
    points_2 = np.zeros(len(syst_labb))
    x = np.array(range(len(syst_labb)))+1
    label = []


  
    for num, key2 in enumerate(syst_labb):
     if key2 != 'mask_tot':
        points[num] = (systematics_dict_total[key][key2]['chi2_red'])
        key2 = '_'+syst_labb[num]
        key2 = '   '+ key+key2
        label.append(key2)
    

    plt.xticks(x, label, rotation='20')
    plt.xlim([0.,len(syst_labb)+1])
    plt.ylim([0.,8])
    plt.plot([0.,len(syst_labb)+1],[2,2], linestyle = 'dashed')
    plt.plot([0.,len(syst_labb)+1],[1,1], linestyle = 'dashed')
    plt.title('Chi2 null hypothesis (constant shear) : {0}   bin [{1} , {2}]'.format(kE_label,z_min,z_max))
    plt.scatter(x,points)
    plt.ylabel('Chi2/dof')

    #systematics_dict_total[key].update({'mask_tot' : tot_mask})
    plt.tight_layout()
    
    print (kE_label)
    plt.savefig('./output_{1}_{2}/{0}_{3}_total_systematics_cc_{4}.png'.format(key,z_min,z_max, kE_label,fract_limits))

    plt.show()
    
def compute_cross_corr(title,save_output_folder,kE_label,sys_map,mapa, info,  add_label,mapa_weight = None, fractional = False):

    print ('\n****************\n')
    print ('TEST (corr) {1}: -> {0}'.format(title,kE_label))

    mask = info['mask_sims']
    sys_map1 = sys_map[mask]
    mapa1 = mapa[mask]-mapa_weight[mask]

    ave = sys_map1.mean()
    std = sys_map1.std()
    if fractional:
    
        sys_map1 -= ave
        sys_map1 /= ave
    else:
        sys_map1 -= ave
        

    if not os.path.exists(save_output_folder+'pairs/pairs_{0}_{1:.3f}_{2:.3f}_{3}_{4}_{5}_{6}_{7}.pkl'.format(info['n_jck'],info['conf']['min_sep'],info['conf']['max_sep'],info['conf']['nbins'],title,kE_label,add_label,fract_limits)):
   
        J = Jack(info['w_estimator'],info['conf'],info['pairs_to_compute'], info['ra'],info['dec'],info['ra'],info['dec'],info['hpix'],info['hpix'],
             mapa1,sys_map1,info['dist_path'],centers=info['centers'],njk=info['n_jck'])
        pairs = J.KKCorrelation()
        save_obj(save_output_folder+'pairs/pairs_{0}_{1:.3f}_{2:.3f}_{3}_{4}_{5}_{6}_{7}'.format(info['n_jck'],info['conf']['min_sep'],info['conf']['max_sep'],info['conf']['nbins'],title,kE_label,add_label,fract_limits),pairs)
    else:
        pairs = load_obj(save_output_folder+'pairs/pairs_{0}_{1:.3f}_{2:.3f}_{3}_{4}_{5}_{6}_{7}'.format(info['n_jck'],info['conf']['min_sep'],info['conf']['max_sep'],info['conf']['nbins'],title,kE_label,add_label,fract_limits))
   
    theta,DD_j,njk=make_wz_errors(pairs,'jackknife',n_jck,True)

    dictu=covariance_jck(DD_j[:,1:],DD_j.shape[1]-1,'jackknife')
    err=dictu['err']
    
    plt.title('systematic map (crosscorr) : {0}'.format(title))
    plt.xlabel('theta (arcmin)'.format(title))        
    plt.ylabel('w(theta)')
    plt.errorbar(theta,DD_j[:,0],err)
    plt.tight_layout()
#plt.yscale('log')
    plt.xscale('log')
    plt.savefig(save_output_folder + 'cc_' + title +'.png')
    plt.show()

    
    inv_cov = linalg.inv(dictu['cov'])
    chi2_val =  (np.matmul((DD_j[:,0]),np.matmul(inv_cov,(DD_j[:,0]))))/len(DD_j[:,0])

    #p_value = 1 - stats.chi2.cdf(np.sum(cv_sol ** 2), len(w.shape(0)))
    
    mute_d = dict()
    mute_d.update({'theta': theta})
    mute_d.update({'w':DD_j[:,0]})
    mute_d.update({'cov' : dictu})
    mute_d.update({'chi2_red' : chi2_val})    
    return mute_d



from multiprocessing import Pool,sharedctypes
from functools import partial
from contextlib import closing
import timeit

def comp(jk,s,m,hp,t):
    if t == 'depth_i' or t == 'depth_g' or t == 'depth_r' or t == 'depth_z':
        
        index = np.random.randint(0,len(s),200000)
        mask1 = (hp[index] != jk)
        mm = np.polyfit(s[index][mask1],m[index][mask1],1)
        
        return jk, mm[0],mm[1]      
    else:
        mask1 = (hp != jk)
        mm = np.polyfit(s[mask1],m[mask1],1)
        return jk, mm[0],mm[1]

def analysis_3(kE_label,sys_map,mask,mapa, hpix, fract_limits=0.,mapa_weight = False, cov_ext = False):
    from scipy import linalg

    print ('\n****************\n')
    title = sys_map['title']
    print ('TEST [method2] {1}: -> {0}'.format(title,kE_label))
    
    n_jck = len(np.unique(hpix))
    fractional = sys_map['fractional']
    
    # plot typical ranges ***********************
   # area_tot = compute_area(mask, nside)
    y,x = np.histogram(sys_map['map'][mask], bins = 50)
    sys_map1 =  sys_map['map'][mask]
    if mapa_weight:
        mapa1 = mapa[mask]-mapa_weight[mask]
    else:
        mapa1 = mapa[mask]
        
    ave = sys_map1.mean()
    std = sys_map1.std()
    
    min_s, max_s = min(sys_map1),max(sys_map1)


    # plot maps properties *************************

    b_arr = np.zeros(n_jck+1)
    c_arr = np.zeros(n_jck+1)
    if title == 'depth_i' or title == 'depth_g' or title == 'depth_r' or title == 'depth_z':
        print ('fast')
        index = np.random.randint(0,len(sys_map1),200000)

        b_arr[0],c_arr[0] = np.polyfit(sys_map1[index],mapa1[index],1)
        
    else:
        b_arr[0],c_arr[0] = np.polyfit(sys_map1,mapa1,1)
    #print 'jck'
    

    aaa = range(n_jck)
    mas = []
    agents=64
    with closing(Pool(processes=agents)) as pool:
        mas.append(pool.map(partial(comp,s=sys_map1,m=mapa1,hp=hpix,t=title),aaa))
    mas = np.array(mas)

    

    mas =mas[0].reshape(n_jck,3)
    index = np.array(mas[:,0]).astype(int)+1
    b_arr[index]=mas[:,1]
    c_arr[index]=mas[:,2]
    #print mas
    #for i in range(n_jck):
    #    mask1 = (hpix != i)
    #    b_arr[i+1],c_arr[i+1] = np.polyfit(sys_map1[mask1],mapa1[mask1],1)
    
    xx = np.linspace(min_s,max_s,100)
    yy = b_arr[0]*xx+c_arr[0]
    dict_output=dict()
    dict_output.update({'b_arr' : b_arr})
    dict_output.update({'c_arr' : c_arr})
    dict_output.update({'b_err' : covariance_scalar_jck(b_arr[1:],n_jck)['err'] })
    
    dict_output.update({'min_s': min_s})
    dict_output.update({'max_s': max_s})

    dict_output.update({'xx': xx})
    dict_output.update({'yy': yy})

    return dict_output

def do_analysis_3(kE,kE_label,weight,systematic_maps,info,fract_limits=0,len_hist=10,hpix_type = 'normal',cov_ext = False, fit ='linear'):
    systematics_dict = dict() 
    for key in systematic_maps.keys():
      
      #if key =='depth_i'or key =='depth_g'or key =='depth_z'or key =='depth_r':
      # if ((key == 'snr' and kE_label =='kE') or (key == 'size_ratio' and kE_label =='kE') or (key == 'E1' and kE_label =='kE') or (key == 'E2' and kE_label =='kE')):
       # pass
       #else:
        t1 = timeit.default_timer()
        if hpix_type == 'normal':
            dict_output = analysis_3(kE_label,systematic_maps[key],info['mask_sims'],kE,info['hpix'])

        if hpix_type == 'marc':
            dict_output = analysis_3(kE_label,systematic_maps[key],info['mask_sims'],kE,info['hpix_f'])
        systematics_dict.update({systematic_maps[key]['title'] : dict_output})
        t2 = timeit.default_timer()
        print ("slope {0} +- {1}".format(dict_output['b_arr'][0],dict_output['b_err']))
        print (t2-t1)
    return systematics_dict



import config
class field_methods(object):
  """
  Utilities for doing pixel and chip calculations. Information from Mike Jarvis.
  """

  chip_centres = {

  'N7':[16.908,191.670],
  'N6':[16.908,127.780],
  'N5':[16.908,63.890],
  'N4':[16.908,0.],
  'N3':[16.908,-63.890],
  'N2':[16.908,-127.780],
  'N1':[16.908,-191.670],
  'N13':[50.724,159.725],
  'N12':[50.724,95.835],
  'N11':[50.724,31.945],
  'N10':[50.724,-31.945],
  'N9':[50.724,-95.835],
  'N8':[50.724,-159.725],
  'N19':[84.540,159.725],
  'N18':[84.540,95.835],
  'N17':[84.540,31.945],
  'N16':[84.540,-31.945],
  'N15':[84.540,-95.835],
  'N14':[84.540,-159.725],
  'N24':[118.356,127.780],
  'N23':[118.356,63.890],
  'N22':[118.356,0.],
  'N21':[118.356,-63.890],
  'N20':[118.356,-127.780],
  'N28':[152.172,95.835],
  'N27':[152.172,31.945],
  'N26':[152.172,-31.945],
  'N25':[152.172,-95.835],
  'N31':[185.988,63.890],
  'N30':[185.988,0.],
  'N29':[185.988,-63.890],
  'S7':[-16.908,191.670],
  'S6':[-16.908,127.780],
  'S5':[-16.908,63.890],
  'S4':[-16.908,0.],
  'S3':[-16.908,-63.890],
  'S2':[-16.908,-127.780],
  'S1':[-16.908,-191.670],
  'S13':[-50.724,159.725],
  'S12':[-50.724,95.835],
  'S11':[-50.724,31.945],
  'S10':[-50.724,-31.945],
  'S9':[-50.724,-95.835],
  'S8':[-50.724,-159.725],
  'S19':[-84.540,159.725],
  'S18':[-84.540,95.835],
  'S17':[-84.540,31.945],
  'S16':[-84.540,-31.945],
  'S15':[-84.540,-95.835],
  'S14':[-84.540,-159.725],
  'S24':[-118.356,127.780],
  'S23':[-118.356,63.890],
  'S22':[-118.356,0.],
  'S21':[-118.356,-63.890],
  'S20':[-118.356,-127.780],
  'S28':[-152.172,95.835],
  'S27':[-152.172,31.945],
  'S26':[-152.172,-31.945],
  'S25':[-152.172,-95.835],
  'S31':[-185.988,63.890],
  'S30':[-185.988,0.],
  'S29':[-185.988,-63.890]
  }

  ccdid=['S29','S30','S31','S25','S26','S27','S28','S20','S21','S22','S23','S24','S14','S15','S16','S17','S18','S19','S8','S9','S10','S11','S12','S13','S1','S2','S3','S4','S5','S6','S7','N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18','N19','N20','N21','N22','N23','N24','N25','N26','N27','N28','N29','N30','N31']

  bad_ccd_names = ['N30', 'S30', 'S7']
  #bad_ccd_nums = [ccdid.index(bad_ccd_name) for bad_ccd_name in bad_ccd_names]


  ccdx=2048.*15.e-6*1000. # col
  ccdy=4096.*15.e-6*1000. # row

  @staticmethod
  def ccd_centres():

    centrex=[]
    centrey=[]
    for i,x in enumerate(field_methods.ccdid):
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[0])
      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[1])

    return np.vstack((centrex,centrey)).T

  @staticmethod
  def ccd_corners():

    centre=np.zeros((62,4,2))
    for i,x in enumerate(field_methods.ccdid):
      c=field_methods.chip_centres.get(x,None)
      centre[i][0][0]=c[0]-field_methods.ccdx/2. # lower left
      centre[i][1][0]=c[0]-field_methods.ccdx/2. # lower right
      centre[i][2][0]=c[0]+field_methods.ccdx/2. # upper left
      centre[i][3][0]=c[0]+field_methods.ccdx/2. # upper right

      centre[i][0][1]=c[1]-field_methods.ccdy/2.
      centre[i][1][1]=c[1]+field_methods.ccdy/2.
      centre[i][2][1]=c[1]-field_methods.ccdy/2.
      centre[i][3][1]=c[1]+field_methods.ccdy/2.

    return centre

  @staticmethod
  def ccd_to_field(ccd,ccdx,ccdy):

    centre=field_methods.ccd_centres()

    centrex=(centre[:,0])[[ccd]]
    centrey=(centre[:,1])[[ccd]]

    return ccdx*15e-6*1000+centrex-field_methods.ccdx/2.,ccdy*15e-6*1000+centrey-field_methods.ccdy/2.

  @staticmethod
  def get_field_pos(cat):

    x,y=field_methods.ccd_to_field(cat.ccd,cat.col,cat.row)

    return x,y 

  @staticmethod
  def translate_to_wcs(pos,image):

    from esutil import wcsutil
    
    wcs=wcsutil.WCS(image, longpole=180.0, latpole=90.0, theta0=90.0)
    ra,dec=wcs.image2sky(pos[0],pos[1])

    return ra,dec

  @staticmethod
  def get_coadd_tile(ra,dec,tiles=None):

    if tiles is None:
      tiles=fio.FITS(config.coaddtiles)[-1].read()

    tmp=tiles['TILENAME'][(ra<tiles['URAUR'])&(dec<tiles['UDECUR'])&(ra>tiles['URALL'])&(dec>tiles['UDECLL'])]
    if len(tmp)==0:
      tmp=tiles['TILENAME'][((ra+360)<tiles['URAUR'])&(dec<tiles['UDECUR'])&((ra+360)>tiles['URALL'])&(dec>tiles['UDECLL'])]

    return tmp[0].rstrip()

  @staticmethod
  def get_radec_coadd_tiles(tiles=None,tiles0=None,file=config.coaddtiles):

    if tiles is None:
      tiles=fio.FITS(file)[-1].read()

    if tiles0 is None:
      mask=np.ones(len(tiles)).astype(bool)
    else:
      mask=np.in1d(np.core.defchararray.strip(tiles['TILENAME']),tiles0,assume_unique=False)

    return tiles,np.vstack(((tiles['URAUR'][mask]+tiles['URALL'][mask])/2.,(tiles['UDECUR'][mask]+tiles['UDECLL'][mask])/2.)).T



import timeit
import treecorr # Module for correlation functions, pip install TreeCorr.
from astropy.cosmology import Planck15 as Planck15
import pickle
import sys

cosmol=Planck15

def covariance_jck(TOTAL_PHI,jk_r,type_cov):
  if type_cov=='jackknife':
      fact=(jk_r-1.)/(jk_r)

  elif type_cov=='bootstrap':
      fact=1./(jk_r)
  #  Covariance estimation

  average=np.zeros(TOTAL_PHI.shape[0])
  cov_jck=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  err_jck=np.zeros(TOTAL_PHI.shape[0])


  for kk in range(jk_r):
    average+=TOTAL_PHI[:,kk]
  average=average/(jk_r)

 # print average
  for ii in range(TOTAL_PHI.shape[0]):
     for jj in range(ii+1):
          for kk in range(jk_r):
            cov_jck[jj,ii]+=TOTAL_PHI[ii,kk]*TOTAL_PHI[jj,kk]

          cov_jck[jj,ii]=(-average[ii]*average[jj]*jk_r+cov_jck[jj,ii])*fact
          cov_jck[ii,jj]=cov_jck[jj,ii]

  for ii in range(TOTAL_PHI.shape[0]):
   err_jck[ii]=np.sqrt(cov_jck[ii,ii])
 # print err_jck

  #compute correlation
  corr=np.zeros((TOTAL_PHI.shape[0],TOTAL_PHI.shape[0]))
  for i in range(TOTAL_PHI.shape[0]):
      for j in range(TOTAL_PHI.shape[0]):
        corr[i,j]=cov_jck[i,j]/(np.sqrt(cov_jck[i,i]*cov_jck[j,j]))

  average=average*fact
  return {'cov' : cov_jck,
          'err' : err_jck,
          'corr':corr,
          'mean':average}






def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)#, encoding='latin1')
    
class Jack_shear(object):
    def __init__(self, conf, distance_file, ra_a, dec_a, ra_b, dec_b, hpix_a, hpix_b, w_a,w_b, e1_a, e2_a, e1_b, e2_b, k = None,  corr='GG',  number_of_cores=2, fact_dist=2,
                 centers=None, njk=None,
                 verbose=False, overwrite=False, cov = False,outp_fold='',label =''):

        '''
        Input:
        '''

        self.output_folder =outp_fold
        self.label = label
        self.conf = conf
        self.distance_file = distance_file

        self.number_of_cores = number_of_cores

        self.ra_s_1 = ra_a
        self.dec_s_1 = dec_a
        self.ra_s_2 = ra_b
        self.dec_s_2 = dec_b


        self.e1_a = e1_a#* (-1)
        self.e2_a = e2_a# (-1)
        self.e1_b = e1_b#
        self.e2_b = e2_b#* (-1)



        self.jk_a = hpix_a
        self.jk_b = hpix_b

        self.weight_a = w_a
        self.weight_b = w_b


        self.k = k
        self.corr = corr
        self.FACT_DIST = fact_dist

        self.njk = njk
        self.centers = centers
        self.cov = cov


    def Correlation(self):
        self.start_all = timeit.default_timer()
        self.prolog()
        pairs= self.epilog()
        return pairs



    def max_distance_region(self):

        def load_obj(name):
            with open(name + '.pkl', 'rb') as f:
                return pickle.load(f)#, encoding='latin1')


        max_dist_region = load_obj(self.distance_file)
        self.max_dist_region=max_dist_region


    def convert_units(self):
        """Change to 'degrees' or to Mpc the units if necessary.
        """

        if 'sep_units' in self.conf.keys():
            un = self.conf['sep_units']
            if un == 'arcmin':
                todeg = 1./60.
            elif un == 'arcsec':
                todeg = 1./60.**2.
            elif un == 'kpc':
                todeg= 1./1000.
            else:
                todeg = 1.
        else:
            todeg = 1.
        return todeg

    def dist_cent_2(self,ra1,dec1,ra2,dec2):

            todeg = np.pi/180.
            ra1 = ra1*todeg
            ra2 = ra2*todeg
            dec1 = dec1*todeg
            dec2 = dec2*todeg

            cos = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
            return np.arccos(cos)/todeg


    def distance(self):

        """Finds the minimum distance to a center for each center.
           Fixes double of this distance as the criteria for not considering correlations,
           which is a conservative choice. This distance has to be at least 4 times the
           maximum angular separation considered. Centers beyond this distance will not be
           considered in the correlations.
        """

        # Find the minimum distance to a center for each center.

        cnt = self.centers
        dist = np.array([np.sort([self.dist_cent(cnt[i],cnt[j]) for i in range(len(cnt))])[1] for j in range(len(cnt))])
        dist = (dist)*self.FACT_DIST


        todeg = self.convert_units()

        if 'max_sep' in self.conf.keys():
            max_sep = self.conf['max_sep'] * todeg

        else:
            raise NotImplementedError("Make use of 'max_sep' in configuration.")

        # Check that the distance is at least 4 times the maximum angular separation.

        self.center_min_dis = np.array( [ 4.*max_sep if x < 4.*max_sep else x for x in dist] )


    def dist_cent(self, a, b):
        """Angular distance between two centers (units: degrees). Makes use of spherical law of cosines.
        """
        todeg = np.pi/180.
        a = a*todeg
        b = b*todeg
        cos = np.sin(a[1])*np.sin(b[1]) + np.cos(a[1])*np.cos(b[1])*np.cos(a[0]-b[0])
        return np.arccos(cos)/todeg


    def cond(self, i, j):

        """Return the maximum conditional distance for a pair of centers to
           determine if the correlation should be computed.
        """

        con = self.center_min_dis
        return max(con[i], con[j])



    def collect(self, pairs):

        shape = (self.njk, self.conf['nbins'])
        DD_a, DR_a, RD_a, RR_a = np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(
            shape)
        DD, DR, RD, RR = np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins']), np.zeros(
            self.conf['nbins']), np.zeros(self.conf['nbins'])

        fact = 1.
        FACT_0 = 0.5
        #print pairs
        for n in range(self.conf['nbins']):
            for jk1 in range(len(pairs[:, 0, 0, n])):
                DD[n] += (pairs[jk1, 0, 0, n]) + FACT_0 * (pairs[jk1, 1, 0, n])
                DR[n] += ((pairs[jk1, 0, 1, n]) + FACT_0 * (pairs[jk1, 1, 1, n]))
                #print DR, DD, RD
                if self.corr == 'GG':
                    RD[n] += ((pairs[jk1, 0, 2, n]) + FACT_0 * (pairs[jk1, 1, 2, n]))
                    #print RD
        for n in range(self.conf['nbins']):
            for jk1 in range(len(pairs[:, 0, 0, n])):
                DD_a[jk1, n] = DD[n] - (pairs[jk1, 0, 0, n]) - fact * FACT_0 * (pairs[jk1, 1, 0, n])
                DR_a[jk1, n] = DR[n] - (pairs[jk1, 0, 1, n]) - fact * FACT_0 * (pairs[jk1, 1, 1, n])
                if self.corr == 'GG':
                    RD_a[jk1, n] = RD[n] - (pairs[jk1, 0, 2, n]) - fact * FACT_0 * (pairs[jk1, 1, 2, n])


        if self.corr == 'GG':
            xip = np.zeros(len(RD))
            xim = np.zeros(len(RD))
            masku = RD !=0.
            
            xip[masku] = DD[masku] / RD[masku]
            xim[masku] = DR[masku] / RD[masku]

            masku = RD_a !=0.
            xip_j = np.zeros((RD_a.shape[0],RD_a.shape[1]))
            xim_j = np.zeros((RD_a.shape[0],RD_a.shape[1]))
            xip_j[masku] = DD_a[masku] / RD_a[masku]
            xim_j[masku] = DR_a[masku] / RD_a[masku]


            return xip, xim, xip_j, xim_j
        else:
            xip = DD / DR
            xip_j = DD_a / DR_a

            return xip, xip, xip_j, xip_j





    def parallel(self, i, j):

        def dic(d):
            DD = d['DD']
            DR = d['DR']
            RD = d['RD']
            RR = d['RR']
            return [DD, DR, RD, RR]

        # this should be generalized to 2 different catalogs

        if self.corr == 'GG':
            mask = np.in1d(self.jk_a, i)
            try:
                cat_s1 = treecorr.Catalog(ra=self.ra_s_1[mask], dec=self.dec_s_1[mask], g1=self.e1_a[mask], g2=self.e2_a[mask],
                                      w=self.weight_a[mask],
                                      ra_units='deg', dec_units='deg')
            except:
                cat_s1 = None
            
            #print self.ra_s_1[mask],self.e1_a[mask]
            mask1 = np.in1d(self.jk_b, j)

            try:
                cat_s2 = treecorr.Catalog(ra=self.ra_s_2[mask1], dec=self.dec_s_2[mask1], g1=self.e1_b[mask1], g2=self.e2_b[mask1],
                                      ra_units='deg',
                                      dec_units='deg', w=self.weight_b[mask1])
            except:
                cat_s2 = None
            

            if (cat_s2 == None) or (cat_s1 == None):
                pairs = [np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins'])]
            else:
                gg = treecorr.GGCorrelation(self.conf)
                gg.process(cat_s1, cat_s2)
                ggp = gg.xip
                ggm = gg.xim
                normalization = gg.weight
                pairs = [ggp * normalization, ggm * normalization, normalization]
            return pairs

        if self.corr == 'KK':
            mask = np.in1d(self.jk_s, i)
            try:
                cat_a = treecorr.Catalog(ra=self.ra_s[mask], dec=self.dec_s[mask], ra_units='deg', dec_units='deg',
                                     k=self.k[mask])
            except:
                cat_a = None

            mask = np.in1d(self.jk_s, j)
            try:
                cat_b = treecorr.Catalog(ra=self.ra_s[mask], dec=self.dec_s[mask], ra_units='deg', dec_units='deg',
                                     k=self.k[mask])
            except:
                cat_b = None

            if (cat_a == None) or (cat_b == None):
                return [np.zeros(self.conf['nbins']), np.zeros(self.conf['nbins'])]
            else:
                kk = treecorr.KKCorrelation(self.conf)
                kk.process(cat_a, cat_b)

                return [kk.xi * kk.weight, kk.weight]

    def prolog(self):

        def update_progress(progress, elapsed_time=0, starting_time=0):
            import time

            barLength = 10  # Modify this to change the length of the progress bar
            status = ""
            if isinstance(progress, int):
                progress = float(progress)
            if not isinstance(progress, float):
                progress = 0
                status = "error: progress var must be float\r\n"
            if progress < 0:
                progress = 0
                status = "Halt...\r\n"
            if progress >= 1:
                progress = 1
                status = "Done...\r\n"
            block = int(round(barLength * progress))

            if progress * 100 > 1. and elapsed_time > 0:
                remaining = ((elapsed_time - starting_time) / progress) * (1. - progress)
                text = "\rPercent: [{0}] {1:.2f}% {2}  - elapsed time: {3} - estimated remaining time: {4}".format(
                    "#" * block + "-" * (barLength - block), progress * 100, status,
                    time.strftime('%H:%M:%S', time.gmtime(elapsed_time - starting_time)),
                    time.strftime('%H:%M:%S', time.gmtime(remaining)))
            else:
                text = "\rPercent: [{0}] {1:.2f}% {2}".format("#" * block + "-" * (barLength - block), progress * 100,
                                                              status)
            sys.stdout.write(text)
            sys.stdout.flush()

        njk = self.njk

        self.distance()
        self.distance()
        self.max_distance_region()

        if self.corr == 'GG':
            gm = 3
            shape = (gm, self.conf['nbins'])

        if self.corr == 'KK':
            gm = 2
            shape = (gm, self.conf['nbins'])

        self.pairs_ring = [[np.zeros(shape) for i in range(2)] for j in range(njk)]
        cnt = self.centers
        self.pairs = [[np.zeros(shape) for i in range(njk)] for j in range(njk)]



        a=np.concatenate((np.array([[(i,j) for i in range(njk)] for j in range(njk)])))
        todeg = self.convert_units()
        max_sep = self.conf['max_sep'] * todeg

        sel = np.array([max([0.,(self.dist_cent(cnt[i],cnt[j]) - (self.max_dist_region[i]+self.max_dist_region[j]))]) < (3. * max_sep )for (i,j) in a])
        b = a[sel]
       # b = a

        def fun_speedup(othersample, othersample1, jackk):
            path =self.output_folder + self.label +'_{0}'.format(jackk)
            if os.path.exists(path+'.pkl'):
                dict_m = load_obj(path)
                pairsCC1 = dict_m['c1']
                pairsCC2 = dict_m['c2']
                pairs_auto = dict_m['a']
                for prs in range(gm):
                    pairsCC1[prs] += pairsCC2[prs]
                    

                self.pairs_ring[jackk][1] = pairsCC1
                self.pairs_ring[jackk][0] = pairs_auto                
            else:
                pairsCC1 = self.parallel(othersample, [jackk])

                pairsCC2 = self.parallel([jackk], othersample1)
                pairs_auto = self.parallel([jackk], [jackk])
                dict_m = dict()
                dict_m.update({'c1':pairsCC1})
                dict_m.update({'c2':pairsCC2})
                dict_m.update({'a':pairs_auto})
                save_obj(path,dict_m)
                for prs in range(gm):
                    pairsCC1[prs] += pairsCC2[prs]
                    

                self.pairs_ring[jackk][1] = pairsCC1
                self.pairs_ring[jackk][0] = pairs_auto

            #print 'eee' ,self.pairs_ring[0][1]#, 'ddd'

            #print (self.pairs_ring[jackk])

        startt = timeit.default_timer()
        
        # parallelized?#**********************************
        import os
        import shutil
        from multiprocessing import Pool,Manager,Process,Queue
        import math



        agents =self.number_of_cores

        chunks=int(math.ceil(np.float(len(np.unique(b[:, 1]))))/agents)
        mute_w = 0

        if self.cov:

         def run_s(jackk):
            #print jackk, "iter"
            mask = (b[:, 1] == jackk) & (b[:, 0] != jackk)
            othersample = b[mask, 0]
            mask = (b[:, 0] == jackk) & (b[:, 1] != jackk)
            othersample1 = b[mask, 1]
            fun_speedup(othersample, othersample1, jackk)                   
                   
         for i in range(chunks+1):
           workers=agents
           work_queue = Queue()
           done_queue = Queue()
           processes = []
           
           for w in range(agents ):
              if (mute_w<len(np.unique(b[:, 1]))):
              
               p = Process(target=run_s, args = ((np.unique(b[:, 1])[mute_w]),))
               p.start()
               processes.append(p)
               work_queue.put('STOP')
               mute_w+=1
        

           for p in processes:
               p.join()
            
       #loading
         for mute_w in (np.unique(b[:, 1])):
          path =self.output_folder + self.label +'_{0}'.format(mute_w)

          dict_m = load_obj(path)
          pairsCC1 = dict_m['c1']
          pairsCC2 = dict_m['c2']
          pairs_auto = dict_m['a']
          for prs in range(gm):
            pairsCC1[prs] += pairsCC2[prs]
                    

            self.pairs_ring[mute_w][1] = pairsCC1
            self.pairs_ring[mute_w][0] = pairs_auto       
            
                     
         '''
        for counter, jackk in enumerate(np.unique(b[:, 1])):
            
            if self.cov:
                mask = (b[:, 1] == jackk) & (b[:, 0] != jackk)
                othersample = b[mask, 0]
                mask = (b[:, 0] == jackk) & (b[:, 1] != jackk)
                othersample1 = b[mask, 1]
                fun_speedup(othersample, othersample1, jackk)
            
            update_progress((float(counter)/len(np.unique(b[:,1]))),timeit.default_timer(),startt)
        '''
        self.pairs = np.array(self.pairs_ring)

    def epilog(self):
        
        pairs = self.pairs
        xip1,xim1,xip_j,xim_j = self.collect(pairs)
        convert=self.convert_units()
        min_sep, max_sep, nedges = self.conf['min_sep']*convert, self.conf['max_sep']*convert, self.conf['nbins']+1
        th = np.linspace(np.log10(min_sep), np.log10(max_sep), nedges)
        theta = 10**np.array([(th[i]+th[i+1])/2 for i in range(len(th)-1)])
        end1 = timeit.default_timer()
        ggp_1=0.
        ggm_1=0.

        xip=0.
        xim=0.
        start =  timeit.default_timer()
        end =  timeit.default_timer()
        if not self.cov:
            #crosscheck ***

            if self.corr == 'KK':

                cat_l1 = treecorr.Catalog(ra=self.ra_s, dec=self.dec_s, ra_units='deg', dec_units='deg',k=self.k)
                cat_l2 = treecorr.Catalog(ra=self.ra_s, dec=self.dec_s, ra_units='deg', dec_units='deg',k=self.k)

                dd = treecorr.KKCorrelation(self.conf)

                dd.process(cat_l1, cat_l2)


                xip = dd.xi
                xim = dd.xi
            if self.corr == 'GG':


                start = timeit.default_timer()
                cat_s1 = treecorr.Catalog(ra=self.ra_s_1, dec=self.dec_s_1, g1=self.e1_a, g2=self.e2_a, w=self.weight_a,
                                          ra_units='deg', dec_units='deg')
                cat_s2 = treecorr.Catalog(ra=self.ra_s_2, dec=self.dec_s_2, g1=self.e1_b, g2=self.e2_b, w=self.weight_b,
                                          ra_units='deg', dec_units='deg')
                gg = treecorr.GGCorrelation(self.conf)

                gg.process(cat_s1,cat_s2)
                end = timeit.default_timer()


                xip = gg.xip
                xim = gg.xim
        
        return {'theta': theta, 'xip':xip, 'xim': xim, 'corr_jckp': xip_j,'corr_jckm': xim_j, 'time_fast' : (end1-self.start_all), 'time_slow': (end-start),'xip1':xip1,'xim1':xim1,}

#plot_resume(fold,systematics_dict_e1,z_min,z_max,nside,kE_label='e1', d68=True,keys=keys)s)


def make_plot_final1(dict_final1, label,z_minz,z_maxz,nside,n_jck=1000,type_c = 'jackknife',fact = False,Except=False,only=False,size_x = 10,size_y=10,cchi=True):
    import copy
    fig, ax = plt.subplots(len(z_minz),1,sharex=True, sharey=True, figsize=(size_x,size_y))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    
    bd = dict()
    for i in range(len(z_minz)):
        #i = 4-i
        z_min = z_minz[i]
        z_max = z_maxz[i]
        #print z_min,z_max
        binx = '{0}_{1}'.format(z_min,z_max)
        fold = './output_new1_{0}_{1}_{2}/'.format(z_min,z_max,nside)
        #for count in range(20):
        if 1==1:
                if 1==1:
                    ll = []
                    ll1 = []
                    key = binx
                    m1 = 0
                    count = 0
                    for key1 in dict_final1[key].keys():
                        if only: 
                            new_key = only[key1]
                            m1+=len(only[key1])
                        elif Except:
                            new_key = []
                            for ii,jj in enumerate(~np.in1d(dict_final1[key][key1]['_w0'].keys(),Except[key1])):
                                if jj:
                                    new_key.append(dict_final1[key][key1]['_w0'].keys()[ii])
                            m1 += len(new_key)
                        else:
                            m1 += len(dict_final1[key][key1]['_w0'].keys())
                            new_key = dict_final1[key][key1]['_w0'].keys()
                    b1 = np.zeros((m1,n_jck+1))
                    for key1 in dict_final1[key].keys():
                        
                        if only: 
                            new_key = only[key1]
             
                        elif Except:
                            new_key = []
                            for ii,jj in enumerate(~np.in1d(dict_final1[key][key1]['_w0'].keys(),Except[key1])):
                                if jj:
                                    new_key.append(dict_final1[key][key1]['_w0'].keys()[ii])
                        else:
                            new_key = dict_final1[key][key1]['_w0'].keys()
                            
                        for key2 in new_key:
                            b1[count,:]= dict_final1[key][key1]['_w0'][key2]['b_arr']
                            count+=1
                            ll.append('{0}_{1}'.format(key1,key2))


                if fact:
                    b1 = b1
                    b_dict1 = covariance_jck(b1[:,1:]*fact[binx], n_jck, type_cov = type_c)
                else:
                    b_dict1 = covariance_jck(b1[:,1:], n_jck, type_cov = type_c)
                
                
                ax[i].set_ylabel('bin : {0} \n b/b_err'.format(binx))
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)  
                
                
                ax[i].scatter(np.arange(len( b1[:,0]))+1, b1[:,0]/b_dict1['err'],label = 'b/b_err')

                
                # high correlations!
                mm = ((b1[:,0]/b_dict1['err']>3.8))#^((b1[:,0]/b_dict1['err']<-3.8))
                yc = b1[:,0]/b_dict1['err']
                xc = np.arange(len( b1[:,0]))+1
                for hh in range(len(mm[mm])):
                    ax[i].annotate( '{0:.2f}'.format(yc[mm][hh]), ((xc[mm][hh]-0.5),3.1),transform=ax[i].transAxes)
                    
                
                if len(mm[mm])>0:
                    yc[mm] = 3.8
                    ax[i].scatter(xc[mm],yc[mm],  color = 'b', marker = '^')

            
                mm = ((b1[:,0]/b_dict1['err']<-3.8))#^((b1[:,0]/b_dict1['err']<-3.8))
                yc = b1[:,0]/b_dict1['err']
                xc = np.arange(len( b1[:,0]))+1
                for hh in range(len(mm[mm])):
                    ax[i].annotate( '{0:.2f}'.format(yc[mm][hh]), ((xc[mm][hh]-0.5),-3.1),transform=ax[i].transAxes)
            
            
                if len(mm[mm])>0:
                    yc[mm] = -3.6
                    ax[i].scatter(xc[mm],yc[mm],  color = 'b', marker = 'v')
            
            
                plt.xticks(np.arange(len( b1[:,0]))+1,ll, rotation='90')
                ax[0].text(0.4,1.05 , '{0}'.format(label),transform=ax[0].transAxes)
                
                w = b1[:,0]
                cns = np.zeros(len(b1[:,0]))
                from scipy import linalg
                inv_cov = linalg.inv(b_dict1['cov'])
                N_p = n_jck
                p_p = len(cns)
                f_hartlap = (N_p-1.)/(N_p-p_p-2.)
                #print f_hartlap
                chi2red =  (np.matmul((w-cns),np.matmul(inv_cov,(w-cns))))/(len(cns)*f_hartlap)
                if cchi:
                    ax[i].text(0.05 , 0.85 , '$\chi^2$/dof : {0:.2f}/{1:.0f}'.format(chi2red*len(cns),len(cns)-1),transform=ax[i].transAxes)#,p2s(1 - stats.t.cdf(chi2red,df=len(cns)-1))),transform=ax[i].transAxes)
                
                
                ax[i].plot([0.,len( b1[:,0])+1],[-2,-2],linestyle="dashed",color="r")
                ax[i].plot([0.,len( b1[:,0])+1],[2,2],linestyle="dashed",color="r")
                
                ax[i].xaxis.label.set_size(14)
                ax[i].yaxis.label.set_size(12)
                
                plt.xlim([0.,len( b1[:,0])+1])
                plt.ylim([-4.,4])
                

        fig.subplots_adjust(wspace=0.,hspace=0.)
        bd.update({binx:[b_dict1,b1]})
    return bd, ll
