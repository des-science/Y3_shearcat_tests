# ALPHA-BETA-ETA test

  Modeling PSF errors and propagate errors in cosmological estimates.
  
## HOW TO RUN
   1) sample_abe.py
   2) create_table.py
   3) produce_2pcfpsfbias.py
   4) plot_samples_and_biases.py
   5) forecast/contaminate.py
   6) forecast/pipeline (COSMOSIS)
   7) plotmarginalised_fidvscont.py

## In the folder code/ 
   you will find all the required code to estimate the parameters alpha, 
   beta and eta and find the total shear PSF bias.  
   These parameters are obtained solving a fitting problem
   where the model vector is determined by the concatenation of
   rho-stats and the data vector by what we have called tau-stats. See
   https://www.overleaf.com/19082636fhptxgfskhwd for definitions and
   details.

   1) abe_test.py main script used for development only (DEPRECATED)
   
   a) --plots & --plotspath: boolean and direction of all the plots the
   correlations rho-stats and tau-staus, covariances matrix both of
   parameters alpha, beta, eta and diagnostics taus, the marginalized
   and countours of the posterior distribution of alpha, beta an eta
   and the final best fits (including a complete error calculation for
   the total bias an each of the bias associated to a particular
   rho-stat)

   b) --outpath the final contamination two point correlation function
   with his diagonal covariance matrix is written

   There are the following modes of analysis:

   a) --eq We have defined a set of three equation and then with the
   concatenation of them we solve the fitting problem. We could use
   less equation to track errors. So to see each equation behaviour we
   can choose a index from 0 to 2. The full right thing is to use the
   whole system of equations then by default is selected --eq=4 which
   represent the whole system. There are also the cases of using two
   equations --eq=10 uses tau0 and tau2, --eq=20 uses tau0 and tau5
   and --eq=21 uses equations with tau2 and tau5

   b) --splitxipxim Calculate independent contamination parameter abe for
   xip and xim tau-stat and rho-stats. The right thing to do is to use
   only one set of parmeter using the whole set of stats
   simultaneously. But again this might be useful to track problems in
   the code.

   c) --overall get the estimated abe values from maximum the overall
   likelihood that correspond to the maximum of the full posterior

   d) --margin get the estimated abe values from the maximum of the
   marginalized posterior.

   e) --abe --ab --ae --be --a --b --e flags to select only a number
   of parameter.  By default it chooses abe. Recalling that a
   represents the leakage, b is proportional to the shape residual and
   e to the size residuals.

   2) sample_abe.py Produce a sample of alphas, betas, etas and chis 
   using emcee and the fitting problem defined above
   
   3) produce_2pcfpsfbias.py Using the sample of abe it produces the bias due
   due to psf errors in the 2pcf. It requires the mean substracted rho-stat which
   will be use to determine the final bias --rhoscosmo. 
   
   Note: --rhoscosmos does not need to be the same --rhos used for the fitting.
   we con go to smaller scales with --rhos and fit better.


   ### In the folder code/essentials 
   There is code to get inputs for the main code. i.e rho-stat file and tau-stat
   file. Particularly getting taus from Flask simulations and Jackknifing 
   require some additional considerations.
   
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

   g) --cosmobin the binning of the rhos used to estimated abe, are
   not necessarily the same to get the PSF bias, since in our
   modelvector in cosmosis might be using a particular binning,
   determine but other physics scale limits. If this flag is present
   then we use the data vector binnig of Y3 analysis

   h) --outpath the resulting measurements are save in this path. We
   follow the format propoused here https://github.com/joezuntz/2point

   2) taus.py Measure the tau-stats i.e, galaxy-stars
   correlations. You must provide:

   a) --metacal_cat A galaxy catalog

   b) --piff_cat the stars catalog with PSF measures of shapes
   https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/PSF_Modeling

   c) --exps_file a list of the exposures. Might be you are interested
   in only some few rather the whole survey.

   d) --bands a particular collection of bands

   e) --frac Choose a random fraction of the input stars

   f) --mod By default the mean is substracted to each field
   (e1,e2,T,w1,w2 ...)  before calculating the correlation. If this is
   set False the substraction will no be done

   g) --outpath the resulting measurements are save in this path. We
   follow the format propoused here https://github.com/joezuntz/2point

   h) --zbin Measure tau in particular redshift bin

   i) --nz_source Indexes catalog to select galaxies in a particular redshiftbin

   3) taus_flask.py Measure tau-stats for each flask realization. Have
   same arguments than taus.py. But additionlally flags associated to
   the flask catalogs

   a) --flask_cat directory containing all the
   flask realization

   b) --seed of the particular flask realization.

   c) --zbin redshift bin of a particular realization

   d) --cookie cookie of a particular realization

   Note: zbin and seed values correspond to lowest limits the code
   will run for a range up to maximum limits defined manually

   4) buildtauscovmat.py Construct the covariance matrix of the taus:
   tau0p tau0m tau2p tau2m tau5p tau5m and write it in a newfile
   together with the taus from metacal. It request:

   a) --tausflask directory with all the measured tau-stats of flask catalogs

   b) --input_tau input file with the taus from Metacal and the
   covariance matrix from this file will be replaced by the covariance
   using flask

   c) --zbin build the covariance matrix for a particular
   tomobin. Warning: use the same tomo bin than the input_file,
   otherwise covariances of different tomo bins will be replaced

  
   ### In the folder code/tests 
   there are some basic test, to validate the catalogs and the code

   a) jk_scaletest.py: run a JK in 4 patches and compare alpha, beta and eta
   obtained solving the fitting problem up to a maximum scale. The scale where
   there is a strong divergence among parameters is where sample variance
   became significant.
   b) overlap_cats.py: stars catalogs areas must be the same that flask catalogs. Plotting the
   areas to compare.
   c) plotallflasktaus.py: dispersion plot of all flask cats and comparison
   with the mean and metacal results.
   d) taus_v1v2.py: compare the mean of correlation functions of  two version of taus flask
   e) covmat_v1v2.py: compare the covariances matrix and inverses of two version of taus flask

   ## In the folder forecast/ 
   there are some scripts to contaminate, test contamination of fiducial cosmology 
   and plot contours plots using the generated samples after running pipeline.

   a) contaminate.py: contaminate tomographycally the fiducial
    cosmology. THere three usage options: by default use the central values of
    the bias, use the --sup of the error bars or use the --inf of the error
    bars

   b) testcontamination.py: plot fiducial, contaminated fiducial and
    contaminant data vector and covariances matrix to test the contamination
    was done properly.

   c) plotmarginalised_fidvscont.py: plot marginalised posteriors and
    contours plots of the samples obtainied after running cosmosis chains with
    datavectors coming from fiducial and contaminated fiducial cosmology

   ### In the folder forecast/pipeline are the used pipelines to run in cosmosis
