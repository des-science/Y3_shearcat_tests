def read_rhos_plots(stat_file, minscale=None, maxscale=None):
    import numpy as np
    import fitsio
    covmat =  fitsio.read(stat_file, ext=1)
    RHO0P =  fitsio.read(stat_file, ext=2); rho0p =  RHO0P['VALUE']; nrows = len(rho0p); covmat0p = covmat[0:nrows,0:nrows] 
    RHO0M =  fitsio.read(stat_file, ext=3); rho0m =  RHO0M['VALUE']; covmat0m = covmat[nrows:2*nrows,nrows:2*nrows ] 
    RHO1P =  fitsio.read(stat_file, ext=4); rho1p =  RHO1P['VALUE']; covmat1p = covmat[2*nrows:3*nrows,2*nrows:3*nrows ] 
    RHO1M =  fitsio.read(stat_file, ext=5); rho1m =  RHO1M['VALUE']; covmat1m = covmat[3*nrows:4*nrows,3*nrows:4*nrows ] 
    RHO2P =  fitsio.read(stat_file, ext=6); rho2p =  RHO2P['VALUE']; covmat2p = covmat[4*nrows:5*nrows,4*nrows:5*nrows ] 
    RHO2M =  fitsio.read(stat_file, ext=7); rho2m =  RHO2M['VALUE']; covmat2m = covmat[5*nrows:6*nrows,5*nrows:6*nrows ] 
    RHO3P =  fitsio.read(stat_file, ext=8); rho3p =  RHO3P['VALUE']; covmat3p = covmat[6*nrows:7*nrows,6*nrows:7*nrows ] 
    RHO3M =  fitsio.read(stat_file, ext=9); rho3m =  RHO3M['VALUE']; covmat3m = covmat[7*nrows:8*nrows,7*nrows:8*nrows ] 
    RHO4P =  fitsio.read(stat_file, ext=10); rho4p =  RHO4P['VALUE']; covmat4p = covmat[8*nrows:9*nrows,8*nrows:9*nrows ] 
    RHO4M =  fitsio.read(stat_file, ext=11); rho4m =  RHO4M['VALUE']; covmat4m = covmat[9*nrows:10*nrows,9*nrows:10*nrows ] 
    RHO5P =  fitsio.read(stat_file, ext=12); rho5p =  RHO5P['VALUE']; covmat5p = covmat[10*nrows:11*nrows,10*nrows:11*nrows ] 
    RHO5M =  fitsio.read(stat_file, ext=13); rho5m =  RHO5M['VALUE']; covmat5m = covmat[11*nrows:12*nrows,11*nrows:12*nrows ] 
    
    meanr = RHO0P['ANG']
    rhos = [rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m,
            rho4p, rho4m, rho5p, rho5m]
    covmats = [covmat0p, covmat0m, covmat1p, covmat1m, covmat2p,
               covmat2m, covmat3p, covmat3m, covmat4p, covmat4m,
               covmat5p, covmat5m]
    if maxscale is not None:
        meanr = meanr[meanr<maxscale]
        idx = len(meanr)    
        for i in range(len(rhos)):
            rhos[i] = rhos[i][:idx]
            covmats[i] = covmats[i][:idx,:idx]
    if minscale is not None:
        idx = len(meanr[meanr<minscale])
        meanr = meanr[idx:]
        for i in range(len(rhos)):
            rhos[i] = rhos[i][idx:]
            covmats[i] = covmats[i][idx:,idx:]
    return  meanr, rhos,  covmats

def read_taus_plots(stat_file, minscale=None, maxscale=None):
    import numpy as np
    import fitsio
    covmat =  fitsio.read(stat_file, ext=1)
    TAU0P =  fitsio.read(stat_file, ext=2); tau0p =  TAU0P['VALUE']; nrows = len(tau0p); covmat0p = covmat[0:nrows,0:nrows] 
    TAU0M =  fitsio.read(stat_file, ext=3); tau0m =  TAU0M['VALUE']; covmat0m = covmat[nrows:2*nrows,nrows:2*nrows ] 
    TAU2P =  fitsio.read(stat_file, ext=4); tau2p =  TAU2P['VALUE']; covmat2p = covmat[2*nrows:3*nrows,2*nrows:3*nrows ] 
    TAU2M =  fitsio.read(stat_file, ext=5); tau2m =  TAU2M['VALUE']; covmat2m = covmat[3*nrows:4*nrows,3*nrows:4*nrows ] 
    TAU5P =  fitsio.read(stat_file, ext=6); tau5p =  TAU5P['VALUE']; covmat5p = covmat[4*nrows:5*nrows,4*nrows:5*nrows ] 
    TAU5M =  fitsio.read(stat_file, ext=7); tau5m =  TAU5M['VALUE']; covmat5m = covmat[5*nrows:6*nrows,5*nrows:6*nrows ] 
    
    meanr = TAU0P['ANG']
    taus = [tau0p, tau0m, tau2p, tau2m, tau5p, tau5m]
    covmats = [covmat0p, covmat0m, covmat2p, covmat2m, covmat5p,
               covmat5m]
    if maxscale is not None:
        meanr = meanr[meanr<maxscale]
        idx = len(meanr)
        for i in range(len(taus)):
            taus[i] = taus[i][:idx]
            covmats[i] = covmats[i][:idx,:idx]
    if minscale is not None:
        idx = len(meanr[meanr<minscale])
        meanr = meanr[idx:]
        for i in range(len(taus)):
            taus[i] = taus[i][idx:]
            covmats[i] = covmats[i][idx:,idx:]
    return  meanr, taus,  covmats

def read_rhos(stat_file, minscale=None, maxscale=None, minbin=None, maxbin=None):
    import numpy as np
    import fitsio
    print("Trying to Read", stat_file)
    covmat=  fitsio.read(stat_file, ext=1)
    RHO0P =  fitsio.read(stat_file, ext=2); rho0p =  RHO0P['VALUE']
    RHO0M =  fitsio.read(stat_file, ext=3); rho0m =  RHO0M['VALUE']
    RHO1P =  fitsio.read(stat_file, ext=4); rho1p =  RHO1P['VALUE']
    RHO1M =  fitsio.read(stat_file, ext=5); rho1m =  RHO1M['VALUE']
    RHO2P =  fitsio.read(stat_file, ext=6); rho2p =  RHO2P['VALUE']
    RHO2M =  fitsio.read(stat_file, ext=7); rho2m =  RHO2M['VALUE']
    RHO3P =  fitsio.read(stat_file, ext=8); rho3p =  RHO3P['VALUE']
    RHO3M =  fitsio.read(stat_file, ext=9); rho3m =  RHO3M['VALUE']
    RHO4P =  fitsio.read(stat_file, ext=10); rho4p =  RHO4P['VALUE']
    RHO4M =  fitsio.read(stat_file, ext=11); rho4m =  RHO4M['VALUE']
    RHO5P =  fitsio.read(stat_file, ext=12); rho5p =  RHO5P['VALUE']
    RHO5M =  fitsio.read(stat_file, ext=13); rho5m =  RHO5M['VALUE']

    meanr = RHO0P['ANG']
    rhos = [rho0p, rho0m, rho1p, rho1m, rho2p, rho2m, rho3p, rho3m,
            rho4p, rho4m, rho5p, rho5m]
    nrhos = len(rhos)
    print(stat_file, "Format seems OK")
    if maxscale is not None and maxbin is None:
        meanr = meanr[meanr<maxscale]
        idx = int(len(meanr))
        ind =  np.arange(idx); size = int(len(covmat)/nrhos)
        indxs = np.concatenate([ind + i*size for i in range(nrhos) ] )
        covmat = covmat[indxs,: ][:,indxs]
        for i in range(len(rhos)):
            rhos[i] = rhos[i][:idx]
    if minscale is not None and minbin is None:
        idx = len(meanr[meanr<minscale])
        meanr = meanr[idx:]
        ind =  np.arange(idx, idx + len(meanr)); size = int(len(covmat)/nrhos)
        indxs = np.concatenate([ind + i*size for i in range(nrhos) ] )
        covmat = covmat[indxs,: ][:,indxs]
        for i in range(len(rhos)):
            rhos[i] = rhos[i][idx:]
    if minbin is not None:
        idx = minbin
        meanr = meanr[idx:]
        ind =  np.arange(idx, idx + len(meanr)); size = int(len(covmat)/nrhos)
        indxs = np.concatenate([ind + i*size for i in range(nrhos) ] )
        covmat = covmat[indxs,: ][:,indxs]
        for i in range(len(rhos)):
            rhos[i] = rhos[i][idx:]
    if maxbin is not None:
        idx = maxbin
        meanr = meanr[:idx]
        ind =  np.arange(idx); size = int(len(covmat)/nrhos)
        indxs = np.concatenate([ind + i*size for i in range(nrhos) ] )
        covmat = covmat[indxs,: ][:,indxs]
        for i in range(len(rhos)):
            rhos[i] = rhos[i][:idx]
        
    return  meanr, rhos,  covmat

def read_taus(stat_file, minscale=None, maxscale=None, minbin=None, maxbin=None):
    import numpy as np
    import fitsio
    #print("Trying to Read", stat_file)
    covmat =  fitsio.read(stat_file, ext=1)
    TAU0P =  fitsio.read(stat_file, ext=2); tau0p =  TAU0P['VALUE']
    TAU0M =  fitsio.read(stat_file, ext=3); tau0m =  TAU0M['VALUE']
    TAU2P =  fitsio.read(stat_file, ext=4); tau2p =  TAU2P['VALUE']
    TAU2M =  fitsio.read(stat_file, ext=5); tau2m =  TAU2M['VALUE']
    TAU5P =  fitsio.read(stat_file, ext=6); tau5p =  TAU5P['VALUE']
    TAU5M =  fitsio.read(stat_file, ext=7); tau5m =  TAU5M['VALUE']
    
    meanr = TAU0P['ANG']
    taus = [tau0p, tau0m, tau2p, tau2m, tau5p, tau5m]
    ntaus = len(taus)
    #print(stat_file, "Format seems OK")
    if maxscale is not None:
        meanr = meanr[meanr<maxscale]
        idx = len(meanr)
        ind =  np.arange(idx); size = int(len(covmat)/ntaus)
        indxs = np.concatenate([ind + i*size for i in range(ntaus) ] )
        covmat = covmat[indxs,: ][:,indxs]   
        for i in range(len(taus)):
            taus[i] = taus[i][:idx]

    if minscale is not None:
        idx = len(meanr[meanr<minscale])
        meanr = meanr[idx:]
        ind =  np.arange(idx, idx + len(meanr)); size = int(len(covmat)/ntaus)
        indxs = np.concatenate([ind + i*size for i in range(ntaus) ] )
        covmat = covmat[indxs,: ][:,indxs]   
        for i in range(len(taus)):
            taus[i] = taus[i][idx:]

    if minbin is not None:
        idx = minbin
        meanr = meanr[idx:]
        ind =  np.arange(idx, idx + len(meanr)); size = int(len(covmat)/ntaus)
        indxs = np.concatenate([ind + i*size for i in range(ntaus) ] )
        covmat = covmat[indxs,: ][:,indxs]
        for i in range(ntaus):
            taus[i] = taus[i][idx:]

    if maxbin is not None:
        idx = maxbin
        meanr = meanr[:idx]
        ind =  np.arange(idx); size = int(len(covmat)/ntaus)
        indxs = np.concatenate([ind + i*size for i in range(ntaus) ] )
        covmat = covmat[indxs,: ][:,indxs]
        for i in range(ntaus):
            taus[i] = taus[i][:idx]
    return  meanr, taus,  covmat

def read_xis(stat_file, minscale=None, maxscale=None):
    import numpy as np
    import fitsio
    covmat =  fitsio.read(stat_file, ext=1)
    XIP =  fitsio.read(stat_file, ext=2); xip =  XIP['VALUE']
    XIM =  fitsio.read(stat_file, ext=3); xim =  XIM['VALUE']
    
    meanr = XIP['ANG']
    xis = [xip, xim]
    nxis = len(xis)
    if maxscale is not None:
        meanr = meanr[meanr<maxscale]
        idx = len(meanr)
        ind =  np.arange(idx); size = int(len(covmat)/nxis)
        indxs = np.concatenate([ind + i*size for i in range(nxis) ] )
        covmat = covmat[indxs,: ][:,indxs]   
        for i in range(nxis):
            xis[i] = xis[i][:idx]

    if minscale is not None:
        idx = len(meanr[meanr<minscale])
        meanr = meanr[idx:]
        ind =  np.arange(idx, idx + len(meanr)); size = int(len(covmat)/nxis)
        indxs = np.concatenate([ind + i*size for i in range(nxis) ] )
        covmat = covmat[indxs,: ][:,indxs]   
        for i in range(len(xis)):
            xis[i] = xis[i][idx:]
    return  meanr, xis,  covmat


