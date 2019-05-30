# ALPHA-BETA-ETA test

  Modeling PSF errors and propagate errors in cosmological estimates.

## In the folder code/ you will find all the required code to estimate
   the parameters alpha, beta and eta.  These parameters are obtained
   solving a fitting problem where the model vector is determined by
   rho-stats and the data vector by what we have called tau-stats. See
   https://www.overleaf.com/19082636fhptxgfskhwd for definitions and
   details.

   1) rhos.py measure the rho-stats for reserved stars. You must
   provide:

   a) --piff_cat the stars catalog with PSF measures of shapes
   https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/PSF_Modeling

   b) --exps_file a list of the exposures. Might be you are interested
   in only some few rather the whole survey.

   c) --bands a particular collection of bands

   d) --frac Choose a random fraction of the input stars

   e) --mod By default the mean is substracted to each field
   (e1,e2,T,w1,w2 ...)  before calculating the correlation. If this is
   set False the substraction will no be done

   f) --obs rho stats were defined for observed field rather than PSF
   measurements, by default we use psf measurements, so if this flag
   is set True it will use obs_e1 rather than piff_e1

   g) --outpath the resulting measurements are save in this path. We
   follow the format propoused here https://github.com/joezuntz/2point