from __future__ import division
from __future__ import print_function
import numpy as np
import timeit
import treecorr # Module for correlation functions, pip install TreeCorr.
from astropy.cosmology import Planck15 as Planck15
import pickle
import sys

#%matplotlib inline
#import matplotlib.pyplot as plt
import sys
import pyfits as pf
import numpy as np
import healpy as hp
import timeit
from astropy.table import Table
import os
import pickle
import os
import timeit
import shutil
from multiprocessing import Pool
import timeit
import numpy as np
from multiprocessing import Pool,sharedctypes
from functools import partial
from contextlib import closing

from mpi4py import MPI


cosmol=Planck15


def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)#, encoding='latin1')
#mute = load_obj("seeing_g_l2_9.pkl")   
class Jack(object):
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


        self.e1_a = e1_a - np.mean(e1_a)
        self.e2_a = e2_a - np.mean(e2_a)
        self.e1_b = e1_b - np.mean(e1_b)#
        self.e2_b = e2_b - np.mean(e2_b)#* (-1)



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

                    #print RD_a
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
                try:
                    dict_m = load_obj(path)
                    pairsCC1 = dict_m['c1']
                    pairsCC2 = dict_m['c2']
                    pairs_auto = dict_m['a']
                    for prs in range(gm):
                        pairsCC1[prs] += pairsCC2[prs]


                    self.pairs_ring[jackk][1] = pairsCC1
                    self.pairs_ring[jackk][0] = pairs_auto  
                except:
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
         #print"pooling dio can"
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
          #print 'loading'
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
        
        return {'theta': theta, 'xip':xip1, 'xim': xim1, 'corr_jckp': xip_j,'corr_jckm': xim_j, 'time_fast' : (end1-self.start_all), 'time_slow': (end-start)}

    
    
    
    


def run_res(path):
    print (path)

    rewrite = False

    J = load_obj(path)
    J.number_of_cores = number_of_cores
    if  ((rewrite) or (not os.path.exists(path[:-2]+".pkl"))):
            print ("running ",path)
            pairs = J.Correlation()
            save_obj(path[:-2],pairs)
   
                
    
if __name__ == '__main__':
    
    number_of_cores=1
    runs_path = load_obj("runs_path")

    runs_to_do = []
    for key in runs_path.keys():
        pathh = key[:-2]+".pkl"
        if not os.path.exists(pathh):
            if not "_0.2_1.3" in key:
                runs_to_do.append(key)


    run_count = 0
    print (len(runs_to_do))
    while run_count<len(runs_to_do):
        comm = MPI.COMM_WORLD
        
	
        print("Hello! I'm rank %d from %d running in total..." % (comm.rank, comm.size))
        try:
            run_res(runs_to_do[(run_count+comm.rank)])#,runs_path)
        except:
            pass
        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 

