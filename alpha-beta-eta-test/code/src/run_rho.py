import os
def getVariances( data_stars, data_galaxies, Rs,tau0, tau2, tau5, prefix='piff', mod=True):
    import numpy as np

    e1 = data_stars['obs_e1']
    e2 = data_stars['obs_e2']
    p_e1 = data_stars[prefix+'_e1']
    p_e2 = data_stars[prefix+'_e2']
    T = data_stars['obs_T']
    p_T = data_stars[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T

    w1 = p_e1*dt
    w2 = p_e2*dt

    #Modified ellipticities reserved stars
    if(mod):
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))

    e1gal = data_galaxies['e_1']
    e2gal = data_galaxies['e_2']
    R11 =  data_galaxies['R11']
    R22 =  data_galaxies['R22']
    R11s =  Rs[0]
    R22s =  Rs[1]

    #Modified ellipticities galaxies
    if (mod):
        #e1gal = (e1gal - np.array(np.mean(e1gal)))/(np.mean(R11)) 
        #e2gal = (e2gal - np.array(np.mean(e2gal)))/(np.mean(R22))
        e1gal = (e1gal - np.array(np.mean(e1gal)))/(np.mean(R11) + np.mean(R11s)) 
        e2gal = (e2gal - np.array(np.mean(e2gal)))/(np.mean(R22) + np.mean(R22s))

    varep = np.var(p_e1) + np.var(p_e2)
    varq = np.var(de1) + np.var(de2)
    varw = np.var(w1) + np.var(w2)
    varegal = np.var(e1gal) + np.var(e2gal)
    vartau0 = (varegal*varep)/(2*tau0.npairs);vartau0[vartau0==np.inf] = 0
    vartau2 = (varegal*varq)/(2*tau2.npairs);vartau2[vartau2==np.inf] = 0
    vartau5 = (varegal*varw)/(2*tau5.npairs);vartau5[vartau5==np.inf] = 0
    
    return vartau0, vartau2, vartau5
    
def write_stats( stat_file, rho0, rho1, rho2, rho3, rho4, rho5, shapenoise=False):
    import json

    if(shapenoise):
        stats = [
            rho0.meanlogr.tolist(),
            rho0.xip.tolist(),
            rho0.xip_im.tolist(),
            rho0.xim.tolist(),
            rho0.xim_im.tolist(),
            (2*rho0.varxi).tolist(),
            rho1.xip.tolist(),
            rho1.xip_im.tolist(),
            rho1.xim.tolist(),
            rho1.xim_im.tolist(),
            (2*rho1.varxi).tolist(),
            rho2.xip.tolist(),
            rho2.xip_im.tolist(),
            rho2.xim.tolist(),
            rho2.xim_im.tolist(),
            (2*rho2.varxi).tolist(),
            rho3.xip.tolist(),
            rho3.xip_im.tolist(),
            rho3.xim.tolist(),
            rho3.xim_im.tolist(),
            (2*rho3.varxi).tolist(),
            rho4.xip.tolist(),
            rho4.xip_im.tolist(),
            rho4.xim.tolist(),
            rho4.xim_im.tolist(),
            (2*rho4.varxi).tolist(),
            rho5.xip.tolist(),
            rho5.xip_im.tolist(),
            rho5.xim.tolist(),
            rho5.xim_im.tolist(),
            (2*rho5.varxi).tolist(),
        ]
    else:
        stats = [
            rho0.meanlogr.tolist(),
            rho0.xip.tolist(),
            rho0.xip_im.tolist(),
            rho0.xim.tolist(),
            rho0.xim_im.tolist(),
            rho0.varxi.tolist(),
            rho1.xip.tolist(),
            rho1.xip_im.tolist(),
            rho1.xim.tolist(),
            rho1.xim_im.tolist(),
            rho1.varxi.tolist(),
            rho2.xip.tolist(),
            rho2.xip_im.tolist(),
            rho2.xim.tolist(),
            rho2.xim_im.tolist(),
            rho2.varxi.tolist(),
            rho3.xip.tolist(),
            rho3.xip_im.tolist(),
            rho3.xim.tolist(),
            rho3.xim_im.tolist(),
            rho3.varxi.tolist(),
            rho4.xip.tolist(),
            rho4.xip_im.tolist(),
            rho4.xim.tolist(),
            rho4.xim_im.tolist(),
            rho4.varxi.tolist(),
            rho5.xip.tolist(),
            rho5.xip_im.tolist(),
            rho5.xim.tolist(),
            rho5.xim_im.tolist(),
            rho5.varxi.tolist(),
        ]
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)

def write_tau_stats(vartaus,stat_file, tau0, tau2, tau5, shapenoise=False):
    import json
    tau0var, tau2var, tau5var = vartaus

    if(shapenoise):
        stats = [
            tau0.meanlogr.tolist(),
            tau0.xip.tolist(),
            tau0.xip_im.tolist(),
            tau0.xim.tolist(),
            tau0.xim_im.tolist(),
            tau0var.tolist(),
            tau2.xip.tolist(),
            tau2.xip_im.tolist(),
            tau2.xim.tolist(),
            tau2.xim_im.tolist(),
            tau2var.tolist(),
            tau5.xip.tolist(),
            tau5.xip_im.tolist(),
            tau5.xim.tolist(),
            tau5.xim_im.tolist(),
            tau5var.tolist(),
        ]
    else:
        stats = [
            tau0.meanlogr.tolist(),
            tau0.xip.tolist(),
            tau0.xip_im.tolist(),
            tau0.xim.tolist(),
            tau0.xim_im.tolist(),
            tau0.varxi.tolist(),
            tau2.xip.tolist(),
            tau2.xip_im.tolist(),
            tau2.xim.tolist(),
            tau2.xim_im.tolist(),
            tau2.varxi.tolist(),
            tau5.xip.tolist(),
            tau5.xip_im.tolist(),
            tau5.xim.tolist(),
            tau5.xim_im.tolist(),
            tau5.varxi.tolist(),
        ]
        
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)

def write_xi_stat( stat_file, rho0, shapenoise=False):
    import json

    if(shapenoise):
        stats = [
            rho0.meanlogr.tolist(),
            rho0.xip.tolist(),
            rho0.xip_im.tolist(),
            rho0.xim.tolist(),
            rho0.xim_im.tolist(),
            (2*rho0.varxi).tolist(),
        ]
    else:
        stats = [
            rho0.meanlogr.tolist(),
            rho0.xip.tolist(),
            rho0.xip_im.tolist(),
            rho0.xim.tolist(),
            rho0.xim_im.tolist(),
            rho0.varxi.tolist(),
        ]
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)

def measure_rho(data, max_sep=300, sep_units='arcmin',  tag=None, prefix='piff', mod=True,  obs=False ):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    T = data['obs_T']
    p_T = data[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    w1 = p_e1*dt
    w2 = p_e2*dt
    w1obs = e1*dt 
    w2obs = e2*dt 
        
    
    #Modified ellipticities
    if(mod):
        e1 = e1 - np.array(np.mean(e1))
        e2 = e2 - np.array(np.mean(e2))       
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
        w1obs = w1obs - np.array(np.mean(w1obs))
        w2obs = w2obs - np.array(np.mean(w2obs))
        
    ra = data['ra']
    dec = data['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    if(obs):
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    else:
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        #wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
            cat.name = tag + ":"  + cat.name

    bin_config = dict(
        sep_units = sep_units,
        nbins = 20,
        min_sep = 2.5,
        max_sep = 250,)
        
    #sep_units = 'degrees',
    '''
    bin_slop = 0.1,
    min_sep = 0.5,
    max_sep = max_sep,
    bin_size = 0.2,
    '''   
    

    results = []
    for (cat1, cat2) in [(ecat, ecat), 
                         (decat, decat),
                          (decat, ecat),
                          (wcat, wcat),
                          (decat, wcat),
                          (ecat, wcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    return results

def measure_tau(data_stars, data_galaxies, Rs, max_sep=300, sep_units='arcmin',  tag=None, prefix='piff', mod=True):
    """Compute the tau statistics
    """
    import treecorr
    import numpy as np
    e1 = data_stars['obs_e1']
    e2 = data_stars['obs_e2']
    p_e1 = data_stars[prefix+'_e1']
    p_e2 = data_stars[prefix+'_e2']
    T = data_stars['obs_T']
    p_T = data_stars[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T

    w1 = p_e1*dt
    w2 = p_e2*dt

    #Modified ellipticities reserved stars
    if(mod):
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))

    
    R11 =  data_galaxies['R11']
    R22 =  data_galaxies['R22']
    R11s =  Rs[0]; R22s =  Rs[1]
    e1gal = data_galaxies['e_1']/(R11s + np.mean(R11))
    e2gal = data_galaxies['e_2']/(R22s + np.mean(R22))

    #Modified ellipticities galaxies
    if (mod):
        e1gal = (e1gal - np.array(np.mean(e1gal)))
        e2gal = (e2gal - np.array(np.mean(e2gal)))
        #e1gal = (e1gal - np.array(np.mean(e1gal)))/(np.mean(R11) + np.mean(R11s)) 
        #e2gal = (e2gal - np.array(np.mean(e2gal)))/(np.mean(R22) + np.mean(R22s))
        
    ra = data_stars['ra']
    dec = data_stars['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    ragal = data_galaxies['ra']
    decgal = data_galaxies['dec']
    print('ragal = ',ragal)
    print('decgal = ',decgal)
    
    
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
    wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)

    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    egal_cat.name = 'egal_cat'
    
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
            cat.name = tag + ":"  + cat.name

    bin_config = dict(
        sep_units = sep_units,
        min_sep = 2.5,
        max_sep = 250,
        nbins = 20,)
    #sep_units = 'degrees',
    '''
    bin_slop = 0.1,
    min_sep = 0.5,
    max_sep = max_sep,
    bin_size = 0.2,
    '''
    

    results = []
    for (cat1, cat2) in [(egal_cat, ecat), 
                         (egal_cat, decat),
                          (egal_cat, wcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)
        
    return results

def measure_xi(data_galaxies, Rs, max_sep=300, sep_units='arcmin',  mod=True):
    """Compute the xi statistics
    """
    import treecorr
    import numpy as np
   
    R11 =  data_galaxies['R11']
    R22 =  data_galaxies['R22']
    R11s =  Rs[0]; R22s =  Rs[1]
    e1gal = data_galaxies['e_1']/(R11s + np.mean(R11))
    e2gal = data_galaxies['e_2']/(R22s + np.mean(R22))
    #Modified ellipticities galaxies
    if (mod):
        e1gal = (e1gal - np.array(np.mean(e1gal)))/(np.mean(R11) + np.mean(R11s)) 
        e2gal = (e2gal - np.array(np.mean(e2gal)))/(np.mean(R22) + np.mean(R22s))
        
    ragal = data_galaxies['ra']
    decgal = data_galaxies['dec']
    print('ragal = ',ragal)
    print('decgal = ',decgal)
    
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)
    egal_cat.name = 'egal_cat'

    bin_config = dict(
        sep_units = sep_units,
        min_sep = 2.5,
        max_sep = 250,
        nbins = 20,)
    '''
    bin_slop = 0.1,
    min_sep = 0.5,
    max_sep = max_sep,
    bin_size = 0.2,
    '''
    

    results = []
    for (cat1, cat2) in [(egal_cat, egal_cat)]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    return results

def measure_xi_tomo(datagal_list, Rs_list, i, j,  max_sep=300, sep_units='arcmin',  mod=True):
    """Compute the tau tomo statistics
    """
    import treecorr
    import numpy as np

    data_galaxies_i = datagal_list[i]
    Rs_i = Rs_list[i]
    R11_i =  data_galaxies_i['R11']
    R22_i =  data_galaxies_i['R22']
    R11s_i =  Rs_i[0]; R22s_i =  Rs_i[1]
    e1gal_i = data_galaxies_i['e_1']/(R11s_i + np.mean(R11_i))
    e2gal_i = data_galaxies_i['e_2']/(R22s_i + np.mean(R22_i))
    #Modified ellipticities galaxies
    if (mod):
        e1gal_i = (e1gal_i - np.array(np.mean(e1gal_i)))
        e2gal_i = (e2gal_i - np.array(np.mean(e2gal_i)))
        
    ragal_i = data_galaxies_i['ra']
    decgal_i = data_galaxies_i['dec']
    print('ragal_i = ',ragal_i)
    print('decgal_i = ',decgal_i)
    
    egal_cat_i = treecorr.Catalog(ra=ragal_i, dec=decgal_i, ra_units='deg', dec_units='deg', g1=e1gal_i, g2=e2gal_i)
    egal_cat_i.name = 'egal_cat_i'

    data_galaxies_j = datagal_list[j]
    Rs_j = Rs_list[j]
    R11_j =  data_galaxies_j['R11']
    R22_j =  data_galaxies_j['R22']
    R11s_j =  Rs_j[0]; R22s_j =  Rs_j[1]
    e1gal_j = data_galaxies_j['e_1']/(R11s_j + np.mean(R11_j))
    e2gal_j = data_galaxies_j['e_2']/(R22s_j + np.mean(R22_j))
    #Modified ellipticities galaxies
    if (mod):
        e1gal_j = (e1gal_j - np.array(np.mean(e1gal_j)))
        e2gal_j = (e2gal_j - np.array(np.mean(e2gal_j)))
        
    ragal_j = data_galaxies_j['ra']
    decgal_j = data_galaxies_j['dec']
    print('ragal_j = ',ragal_j)
    print('decgal_j = ',decgal_j)
    
    egal_cat_j = treecorr.Catalog(ra=ragal_j, dec=decgal_j, ra_units='deg', dec_units='deg', g1=e1gal_j, g2=e2gal_j)
    egal_cat_j.name = 'egal_cat_j'

    bin_config = dict(
        sep_units = sep_units,
        min_sep = 2.5,
        max_sep = 250,
        nbins = 20,)
    '''
    bin_slop = 0.1,
    min_sep = 0.5,
    max_sep = max_sep,
    bin_size = 0.2,
    '''
    

    results = []
    for (cat1, cat2) in [(egal_cat_i, egal_cat_j)]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)
        rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    return results

def band_combinations(bands, single=True, combo=True,  allcombo=True):
    if(allcombo):
        if single:
            use_bands = [ [b] for b in bands ]
        else:
            use_bands = []

        if combo:
            if 'r' in bands and 'i' in bands:
                use_bands.append(['r', 'i'])
            if 'r' in bands and 'i' in bands and 'z' in bands:
                use_bands.append(['r', 'i', 'z'])
            if 'g' in bands and 'r' in bands and 'i' in bands and 'z' in bands:
                use_bands.append(['g', 'r', 'i', 'z'])
            if 'g' in bands and 'r' in bands and 'i' in bands and 'z' in bands and 'Y' in bands:
                use_bands.append(['g', 'r', 'i', 'z', 'Y'])
    else:
        letters = [k for k in bands]
        use_bands = []
        use_bands.append(letters)

    print('use_bands = ',use_bands)
    print('tags = ',[ ''.join(band) for band in use_bands ])
    return use_bands
def do_rho_stats(data, bands, tilings, outpath, max_sep=300, sep_units='arcmin', prefix='piff', name='rho_corr', bandcombo=False, mod=True, obs=False,  shapenoise=False):
    import numpy as np 
    print('Start rho CANONICAL: ',prefix,name)
    # Measure the canonical rho stats using all pairs:
    use_bands = band_combinations(bands, allcombo=bandcombo)
    for band in use_bands:
        print('band ',band)
        mask = np.in1d(data['band'],band)
        print('sum(mask) = ',np.sum(mask))
        print('len(data[mask]) = ',len(data[mask]))
        tag = ''.join(band)
        stats = measure_rho(data[mask], max_sep=max_sep, sep_units=sep_units, tag=tag, prefix=prefix,  mod=mod,  obs=obs)
        stat_file = os.path.join(outpath, "rho_%s_%s.json"%(tag, name))
        write_stats(stat_file,*stats, shapenoise=shapenoise)

def do_tau_stats( data_galaxies, Rs, data_stars,  bands, tilings, outpath, max_sep=300, sep_units='arcmin', prefix='piff',  name='tau_corr',   bandcombo=False, mod=True,  shapenoise=True):
    import numpy as np 
    print('Start Tau CANONICAL: ',prefix,name)
    
    # Measure the canonical rho stats using all pairs:
    use_bands = band_combinations(bands, allcombo=bandcombo)
    for band in use_bands:
        print('band ',band)
        mask_stars = np.in1d(data_stars['band'],band)
        #mask_galaxies = np.in1d(data_galaxies['band'],band)
        print('sum(mask) = ',np.sum(mask_stars))
        print('len(data[mask]) = ',len(data_stars[mask_stars]))
        tag = ''.join(band)
        #mod max_sep
        stats = measure_tau(data_stars[mask_stars], data_galaxies, Rs,  max_sep=max_sep, sep_units=sep_units,  tag=tag, prefix=prefix, mod=mod)
        stat_file = os.path.join(outpath, "tau_%s_%s.json"%(tag, name))
        vartau0, vartau2, vartau5 =  getVariances(data_stars, data_galaxies, Rs, *stats, prefix=prefix, mod=mod)
        vartaus = [vartau0, vartau2, vartau5]
        write_tau_stats(vartaus, stat_file,*stats, shapenoise=shapenoise)

def do_xi_stats(data_galaxies, Rs, outpath, max_sep=300, sep_units='arcmin', name='all', bandcombo=True, mod=True,  shapenoise=False):
    import numpy as np 
    print('Start CANONICAL: ',name)
    stats = measure_xi(data_galaxies, Rs,  max_sep=max_sep, sep_units=sep_units,  mod=mod)
    stat_file = os.path.join(outpath, "xi_%s_%s.json"%('riz', name))
    write_xi_stat(stat_file,*stats, shapenoise=shapenoise)

def do_xi_stats_tomo(datagal_list, Rs_list, i, j, outpath, max_sep=300, sep_units='arcmin', name='all', bandcombo=True, mod=True,  shapenoise=False):
    import numpy as np 
    print('Start Xi measurement: ',name)
    stats = measure_xi_tomo(datagal_list, Rs_list, i, j,  max_sep=max_sep, sep_units=sep_units,  mod=mod)
    stat_file = os.path.join(outpath, "xi_%s_%s.json"%('riz', name))
    write_xi_stat(stat_file,*stats, shapenoise=shapenoise)
